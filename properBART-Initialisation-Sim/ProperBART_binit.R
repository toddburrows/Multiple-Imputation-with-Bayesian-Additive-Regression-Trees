library(tidyverse)
library(mice)
library(BART)

library(rJava)
#Sys.setenv(JAVA_TOOL_OPTIONS="--add-modules=jdk.incubator.vector")
options(java.parameters = c("-Xmx20g", "--add-modules=jdk.incubator.vector", "-XX:+UseZGC"))
library(bartMachine)

properBARTbinit <- function(xmis, maxiter = 10, ntree = 50,
                     num_burn_in = 250, num_iterations_after_burn_in = 1000, 
                     alpha_param = 0.95, beta_param = 2, k_param = 2, q_param = 0.9, nu_param = 3,
                     Bverbose = FALSE, verbose = TRUE){
  
  if(verbose){
    start_time <- proc.time()
  }
  
  original_xmis <- xmis
  
  #1. initialise missing values
  num_variables <- ncol(original_xmis)
  
  missing_index <- lapply(1:num_variables, function(i)
    which(is.na(original_xmis[,i])))
  
  #2. get sorted indices of columns in xmis with increasing amount of missingness
  miss <- numeric(num_variables)
  for (i in 1:num_variables){
    miss[i] <- sum(is.na(original_xmis[,i]))
  }
  
  k <- data.frame(indices = 1:num_variables,
                  missingness = miss)
  
  k <- k[order(k$missingness),]
  k <- k[,1]
  
  for(i in k){
    if(sum(is.na(original_xmis[,i]))>0){
      y_original <- original_xmis[,i]
      j <- missing_index[[i]]
      y_obs <- na.omit(y_original)
      x_obs <- original_xmis[-j,-i]
      x_mis <- original_xmis[j,-i]
      
      fit <- bartMachine(X = x_obs, y = y_obs,
                         num_trees = 15, num_burn_in = 25,
                         num_iterations_after_burn_in = 25,
                         alpha = alpha_param, beta = beta_param, k = k_param, q = q_param, nu = nu_param,
                         use_missing_data = TRUE, use_missing_data_dummies_as_covars = TRUE,
                         run_in_sample = FALSE, verbose = FALSE, mem_cache_for_speed = FALSE)
      
      xmis[j,i] <- predict(fit, new_data = x_mis)
      
      rm(fit)
    }
  }
  
  
  #3. initialise stopping condition
  stop_criteria <- FALSE
  
  x_imp_new <- xmis
  
  iteration_count <- 0
  
  max_reached <- FALSE
  
  if(verbose){
    cat("\n\nData matrix initialised, properBART beginning\n\n")
  }
  
  while(!stop_criteria){
    if(verbose){
      cat("\nproperBART iteration ",iteration_count+1," in progress...\n")
      t_start <- proc.time()
    }
    
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
        
        #7. impute new values by posterior predictive draws
        posterior_means <- predict(fit, x_mis)
        
        posterior_errorvars <- fit$sigma[-(1:num_burn_in)]
        
        rm(fit)
        
        S <- nrow(posterior_means)
        
        s <- sample(S, size = 1)
        
        mu_s <- posterior_means[s,]
        #mu_s <- colMeans(posterior_means)
        
        sigsq_s <- posterior_errorvars[s]
        #sigsq_s <- mean(posterior_errorvars)
        
        rm(posterior_means,posterior_errorvars)
        
        #8. update imputed matrix using imputations
        x_imp_new[j,i] <- rnorm(length(mu_s), mu_s, sqrt(sigsq_s))
        
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
    }
    if(verbose){
      delta_start <- proc.time() - t_start
      cat("Iteration ",iteration_count," completed\n   Time: ",delta_start[3]," seconds\n\n")
    }
  }
  if(verbose){
    cat("Maximum iterations reached")
    delta_overall <- proc.time() - start_time
    cat("\nOverall time: ",delta_overall[3]," seconds\n\n")
  }
  return(x_imp_new)
}
