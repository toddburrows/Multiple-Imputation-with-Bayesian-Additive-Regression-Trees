library(tidyverse)
library(mice)
library(BART)

missBART <- function(xmis, maxiter = 10, ntree = 50,
                     num_burn_in = 250, num_iterations_after_burn_in = 1000, 
                     alpha_param = 0.95, beta_param = 2, k_param = 2, q_param = 0.9, nu_param = 3,
                     Bverbose = FALSE, verbose = TRUE){
  
  if(verbose){
    start_time <- proc.time()
  }
  
  original_xmis <- xmis
  
  #1. make an initial guess for missing values
  imp_mean <- mice(xmis, method = "mean", m = 1, maxit = 1, printFlag = FALSE)
  xmis <- complete(imp_mean)
  
  #2. get sorted indices of columns in xmis with increasing amount of missingness
  num_variables <- length(original_xmis)
  miss <- c()
  for (i in 1:num_variables){
    miss <- c(miss, sum(is.na(original_xmis[,i])))
  }
  
  k <- data.frame(indices = 1:num_variables,
                  missingness = miss)
  
  k <- k[order(k$missingness),]
  k <- k[,1]
  
  missing_index <- lapply(1:num_variables, function(i)
    which(is.na(original_xmis[,i])))
  
  #3. initialise stopping condition
  stop_criteria <- FALSE
  
  x_imp_new <- xmis
  
  iteration_count <- 0
  
  previous_delta <- NA
  
  max_reached <- FALSE
  
  if(verbose){
    cat("\n\nData matrix initialised, missBART beginning\n\n")
  }
  
  while(!stop_criteria){
    if(verbose){
      cat("\nmissBART iteration ",iteration_count+1," in progress...\n")
      t_start <- proc.time()
    }
    
    #4. store previously imputed matrix
    x_imp_old <- x_imp_new
    
    #5. for each variable, in the order of increasing missingness
    counter <- 1
    for (i in k){
      if(verbose){
        cat("   Column ",counter," in progress...")
      }
      if ( sum(is.na(original_xmis[,i])) > 0){
        #column of variable of current interest from original set with missingness
        y_original <- original_xmis[,i]
        
        #row index of missing values in this variable
        #j <- which(is.na(y_original))
        j <- missing_index[[i]]
        
        #y_obs are the values that are observed for the variable of interest
        y_obs <- na.omit(y_original)
        
        #x_obs are the values (could include imputations) for all the variables 
        #except the one of interest for the rows where the values are observed in
        #the variable of interest
        x_obs <- x_imp_new[-j,-i]
        
        #x_mis are the values (could include imputations) for all the variables
        #except the one of interest for the rows where the values are missing in
        #the variable of interest
        x_mis <- x_imp_new[j,-i]
        
        #6. fit a BART model
        fit <- wbart(x_obs, y_obs,
                     ntree = ntree,
                     nskip = num_burn_in,
                     ndpost = num_iterations_after_burn_in,
                     base = alpha_param,
                     power = beta_param,
                     k = k_param,
                     sigquant = q_param,
                     nu = nu_param,
                     printevery = 0)
        
        if(verbose){
          cat("BART model fitted, predicting...")
        }
        
        x_imp_current <- x_imp_new
        
        #7. impute new values
        x_imp_current[j,i] <- colMeans(predict(fit, newdata = x_mis))
        
        #8. update imputed matrix using imputations
        x_imp_new <- x_imp_current
      } #9. end for loop
      if(verbose){
        cat("done!\n")
        counter <- counter + 1
      }
    }
    
    #10. update stopping condition
    
    #iteration count for component of stopping criterion that is a maximum iteration limit
    iteration_count <- iteration_count + 1
    
    if (iteration_count >= maxiter){
      stop_criteria <- TRUE
      max_reached <- TRUE
    } else{
      #component of stopping criterion that is the relative change of the values
      deltaN <- sum((x_imp_new-x_imp_old)**2) / sum(x_imp_new**2)
      
      #only can do comparison after 2 iterations
      if (iteration_count != 1){
        if (deltaN > previous_delta){
          stop_criteria <- TRUE
          #keep the previous "better" data matrix
          x_imp_new <- x_imp_old
        }
      }
      previous_delta <- deltaN
    }
    if(verbose){
      delta_start <- proc.time() - t_start
      cat("Iteration ",iteration_count," completed\n   Relative change: ",previous_delta,"\n   Time: ",delta_start[3]," seconds\n\n")
      }
  }
  if(verbose){
    if(max_reached){
      cat("Maximum iterations reached")
    }else{
      cat("Converged after ",iteration_count," iterations")
    }
    delta_overall <- proc.time() - start_time
    cat("\nOverall time: ",delta_overall[3]," seconds\n\n")
  }
  return(x_imp_new)
}