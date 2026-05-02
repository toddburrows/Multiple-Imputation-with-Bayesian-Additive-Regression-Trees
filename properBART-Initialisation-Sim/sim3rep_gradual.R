# Parameters
grad_num <- 5 #simulation number ,i.e. if combining 5 simulations this would be 1 for the first simulation

sim_num <- 2 #Changing Dataset = 2
n_obs <- 1000
missp <- 0.4 #20% Missingness
sims <- 50 #This is number for just this simulation
M <- 12
maxiter <- 10
ntree <- 50
num_burn_in <- 500
num_iterations_after_burn_in <- 200

# Required functions
source('/Users/toddb/Library/Mobile Documents/com~apple~CloudDocs/University/3rd Year/Project/Simulations/SimulationFunctions.R')
source('/Users/toddb/Library/Mobile Documents/com~apple~CloudDocs/University/3rd Year/Project/Simulations/Sim3/ProperBART_binit.R')
source('/Users/toddb/Library/Mobile Documents/com~apple~CloudDocs/University/3rd Year/Project/Simulations/Sim3/ProperBART_minit.R')
source('/Users/toddb/Library/Mobile Documents/com~apple~CloudDocs/University/3rd Year/Project/Simulations/Sim3/ProperBART_rinit.R')
source('/Users/toddb/Library/Mobile Documents/com~apple~CloudDocs/University/3rd Year/Project/Simulations/Sim3/ProperBART_minitmean.R')
library(parallel)

# Reproducibility
RNGkind("L'Ecuyer-CMRG")
set.seed(as.numeric(paste0(193,grad_num)))

# Obtaining True_Beta and Names
dat <- generateSimData(2)

dataset <- dat[[1]]
true_beta <- dat[[2]]

complete_fit <- fitComplete(dataset)
complete_lm_coeff <- complete_fit[[1]]
var_names <- names(complete_lm_coeff)

# Simulation
overall_time <- proc.time()

complete_coeffs <- matrix(NA, ncol = length(true_beta), nrow=sims)
complete_vcovs <- matrix(NA, ncol = length(true_beta), nrow=sims)

# Mean Initialisation
mean_coeffs <- matrix(NA, ncol = length(true_beta), nrow=sims)
mean_withins <- matrix(NA, ncol = length(true_beta), nrow=sims)
mean_betweens <- matrix(NA, ncol = length(true_beta), nrow=sims)
mean_totals <- matrix(NA, ncol = length(true_beta), nrow=sims)
mean_coverage_tot <- matrix(NA, ncol = length(true_beta), nrow=sims)

# Mean Initialisation + Averaged Posterior Draws
meav_coeffs <- matrix(NA, ncol = length(true_beta), nrow=sims)
meav_withins <- matrix(NA, ncol = length(true_beta), nrow=sims)
meav_betweens <- matrix(NA, ncol = length(true_beta), nrow=sims)
meav_totals <- matrix(NA, ncol = length(true_beta), nrow=sims)
meav_coverage_tot <- matrix(NA, ncol = length(true_beta), nrow=sims)

# Random Initialisation
rand_coeffs <- matrix(NA, ncol = length(true_beta), nrow=sims)
rand_withins <- matrix(NA, ncol = length(true_beta), nrow=sims)
rand_betweens <- matrix(NA, ncol = length(true_beta), nrow=sims)
rand_totals <- matrix(NA, ncol = length(true_beta), nrow=sims)
rand_coverage_tot <- matrix(NA, ncol = length(true_beta), nrow=sims)

# BART Initialisation
bart_coeffs <- matrix(NA, ncol = length(true_beta), nrow=sims)
bart_withins <- matrix(NA, ncol = length(true_beta), nrow=sims)
bart_betweens <- matrix(NA, ncol = length(true_beta), nrow=sims)
bart_totals <- matrix(NA, ncol = length(true_beta), nrow=sims)
bart_coverage_tot <- matrix(NA, ncol = length(true_beta), nrow=sims)

# Track Time
sim_time <- numeric(sims)
mean_time <- numeric(sims)
meav_time <- numeric(sims)
rand_time <- numeric(sims)
bart_time <- numeric(sims)

for(j in 1:sims){
  sim_start <- proc.time()
  
  dat <- generateSimData(n_obs)
  
  dataset <- dat[[1]]
  
  complete_fit <- fitComplete(dataset)
  
  complete_lm_coeff <- complete_fit[[1]]
  complete_lm_res <- complete_fit[[2]]
  complete_lm_vcov <- complete_fit[[3]]
  
  complete_coeffs[j,] <- complete_lm_coeff
  complete_vcovs[j,] <- diag(complete_lm_vcov)
  
  dataset_miss <- makeMissing(dataset, missp)
  
  # Mean Initialisation
  mean_start <- proc.time()
  
  mean_results <- mclapply(1:M, function(i) {
    
    mis_imp_mean <- properBARTminit(
      dataset_miss,
      maxiter = maxiter,
      ntree = ntree,
      num_burn_in = num_burn_in,
      num_iterations_after_burn_in = num_iterations_after_burn_in,
      verbose = FALSE
    )
    
    mean_lm <- lm(y ~ 0 + x1 + x2 + I(x3^2) + I(sin(x4)), data = mis_imp_mean)
    
    list(coef = as.numeric(coef(mean_lm)),
         vcov = vcov(mean_lm))
  },
  mc.cores = min(M, detectCores() - 2))
  
  mean_time[j] <- (proc.time() - mean_start)[3]
  
  ## Calculations
  mean_coeffs_mat <- do.call(rbind, lapply(mean_results, `[[`, "coef"))
  
  mean_within <- Reduce(`+`, lapply(mean_results, `[[`, "vcov"))
  
  mean_combined_within <- mean_within / M
  
  mean_combined_coeffs <- colMeans(mean_coeffs_mat)
  
  centre_matrix <- mean_coeffs_mat - matrix(1, ncol = 1, nrow = M) %*% mean_combined_coeffs
  
  mean_combined_between <- 1/(M-1) * t(centre_matrix) %*% centre_matrix
  
  mean_combined_total <- mean_combined_within + (1+1/M) * mean_combined_between
  
  within_diag <- diag(mean_combined_within)
  between_diag <- diag(mean_combined_between)
  total_diag <- diag(mean_combined_total)
  
  lambda <- ((1 + (1/M)) * between_diag) / total_diag
  nu_old <- (M - 1)/(lambda^2)
  nu_com <- complete_lm_res
  nu_obs <- ( (nu_com+1)/(nu_com+3) ) * nu_com * (1-lambda)
  nu <- (nu_old * nu_obs) /(nu_old + nu_obs)
  
  mean_interval_tot <- matrix(c(mean_combined_coeffs - qt(0.975,nu)*sqrt(total_diag),
                                mean_combined_coeffs + qt(0.975,nu)*sqrt(total_diag)),
                              nrow=2,byrow = TRUE)
  mean_coverage_tot[j,] <- ifelse(true_beta > mean_interval_tot[1,] & true_beta < mean_interval_tot[2,]
                                  ,TRUE,FALSE)
  
  mean_coeffs[j,] <- mean_combined_coeffs
  mean_withins[j,] <- within_diag
  mean_betweens[j,] <- between_diag
  mean_totals[j,] <- total_diag
  
  # Mean Initialisation + Averaged Posterior Draws
  meav_start <- proc.time()
  
  meav_results <- mclapply(1:M, function(i) {
    
    mis_imp_meav <- properBARTminitmean(
      dataset_miss,
      maxiter = maxiter,
      ntree = ntree,
      num_burn_in = num_burn_in,
      num_iterations_after_burn_in = num_iterations_after_burn_in,
      verbose = FALSE
    )
    
    meav_lm <- lm(y ~ 0 + x1 + x2 + I(x3^2) + I(sin(x4)), data = mis_imp_meav)
    
    list(coef = as.numeric(coef(meav_lm)),
         vcov = vcov(meav_lm))
  },
  mc.cores = min(M, detectCores() - 2))
  
  meav_time[j] <- (proc.time() - meav_start)[3]
  
  ## Calculations
  meav_coeffs_mat <- do.call(rbind, lapply(meav_results, `[[`, "coef"))
  
  meav_within <- Reduce(`+`, lapply(meav_results, `[[`, "vcov"))
  
  meav_combined_within <- meav_within / M
  
  meav_combined_coeffs <- colMeans(meav_coeffs_mat)
  
  centre_matrix <- meav_coeffs_mat - matrix(1, ncol = 1, nrow = M) %*% meav_combined_coeffs
  
  meav_combined_between <- 1/(M-1) * t(centre_matrix) %*% centre_matrix
  
  meav_combined_total <- meav_combined_within + (1+1/M) * meav_combined_between
  
  within_diag <- diag(meav_combined_within)
  between_diag <- diag(meav_combined_between)
  total_diag <- diag(meav_combined_total)
  
  lambda <- ((1 + (1/M)) * between_diag) / total_diag
  nu_old <- (M - 1)/(lambda^2)
  nu_com <- complete_lm_res
  nu_obs <- ( (nu_com+1)/(nu_com+3) ) * nu_com * (1-lambda)
  nu <- (nu_old * nu_obs) /(nu_old + nu_obs)
  
  meav_interval_tot <- matrix(c(meav_combined_coeffs - qt(0.975,nu)*sqrt(total_diag),
                                meav_combined_coeffs + qt(0.975,nu)*sqrt(total_diag)),
                              nrow=2,byrow = TRUE)
  meav_coverage_tot[j,] <- ifelse(true_beta > meav_interval_tot[1,] & true_beta < meav_interval_tot[2,]
                                  ,TRUE,FALSE)
  
  meav_coeffs[j,] <- meav_combined_coeffs
  meav_withins[j,] <- within_diag
  meav_betweens[j,] <- between_diag
  meav_totals[j,] <- total_diag
  
  # Random Initialisation
  rand_start <- proc.time()
  
  rand_results <- mclapply(1:M, function(i) {
    
    mis_imp_rand <- properBARTrinit(
      dataset_miss,
      maxiter = maxiter,
      ntree = ntree,
      num_burn_in = num_burn_in,
      num_iterations_after_burn_in = num_iterations_after_burn_in,
      verbose = FALSE
    )
    
    rand_lm <- lm(y ~ 0 + x1 + x2 + I(x3^2) + I(sin(x4)), data = mis_imp_rand)
    
    list(coef = as.numeric(coef(rand_lm)),
         vcov = vcov(rand_lm))
  },
  mc.cores = min(M, detectCores() - 2))
  
  rand_time[j] <- (proc.time() - rand_start)[3]
  
  ## Calculations
  rand_coeffs_mat <- do.call(rbind, lapply(rand_results, `[[`, "coef"))
  
  rand_within <- Reduce(`+`, lapply(rand_results, `[[`, "vcov"))
  
  rand_combined_within <- rand_within / M
  
  rand_combined_coeffs <- colMeans(rand_coeffs_mat)
  
  centre_matrix <- rand_coeffs_mat - matrix(1, ncol = 1, nrow = M) %*% rand_combined_coeffs
  
  rand_combined_between <- 1/(M-1) * t(centre_matrix) %*% centre_matrix
  
  rand_combined_total <- rand_combined_within + (1+1/M) * rand_combined_between
  
  within_diag <- diag(rand_combined_within)
  between_diag <- diag(rand_combined_between)
  total_diag <- diag(rand_combined_total)
  
  lambda <- ((1 + (1/M)) * between_diag) / total_diag
  nu_old <- (M - 1)/(lambda^2)
  nu_com <- complete_lm_res
  nu_obs <- ( (nu_com+1)/(nu_com+3) ) * nu_com * (1-lambda)
  nu <- (nu_old * nu_obs) /(nu_old + nu_obs)
  
  rand_interval_tot <- matrix(c(rand_combined_coeffs - qt(0.975,nu)*sqrt(total_diag),
                                rand_combined_coeffs + qt(0.975,nu)*sqrt(total_diag)),
                              nrow=2,byrow = TRUE)
  rand_coverage_tot[j,] <- ifelse(true_beta > rand_interval_tot[1,] & true_beta < rand_interval_tot[2,]
                                  ,TRUE,FALSE)
  
  rand_coeffs[j,] <- rand_combined_coeffs
  rand_withins[j,] <- within_diag
  rand_betweens[j,] <- between_diag
  rand_totals[j,] <- total_diag
  
  # BART Initialisation
  bart_start <- proc.time()
  
  cat("\nBART starting...") #change mclapply to just lapply
  bart_results <- mclapply(1:M, function(i) {
    
    mis_imp_bart <- properBARTbinit(
      dataset_miss,
      maxiter = maxiter,
      ntree = ntree,
      num_burn_in = num_burn_in,
      num_iterations_after_burn_in = num_iterations_after_burn_in,
      verbose = FALSE
    )
    
    bart_lm <- lm(y ~ 0 + x1 + x2 + I(x3^2) + I(sin(x4)), data = mis_imp_bart)
    
    list(coef = as.numeric(coef(bart_lm)),
         vcov = vcov(bart_lm))
  }
  ,mc.cores = min(M, detectCores() - 2)
  )
  #cat("done!")
  
  bart_time[j] <- (proc.time() - bart_start)[3]
  
  ## Calculations
  bart_coeffs_mat <- do.call(rbind, lapply(bart_results, `[[`, "coef"))
  
  bart_within <- Reduce(`+`, lapply(bart_results, `[[`, "vcov"))
  
  bart_combined_within <- bart_within / M
  
  bart_combined_coeffs <- colMeans(bart_coeffs_mat)
  
  centre_matrix <- bart_coeffs_mat - matrix(1, ncol = 1, nrow = M) %*% bart_combined_coeffs
  
  bart_combined_between <- 1/(M-1) * t(centre_matrix) %*% centre_matrix
  
  bart_combined_total <- bart_combined_within + (1+1/M) * bart_combined_between
  
  within_diag <- diag(bart_combined_within)
  between_diag <- diag(bart_combined_between)
  total_diag <- diag(bart_combined_total)
  
  lambda <- ((1 + (1/M)) * between_diag) / total_diag
  nu_old <- (M - 1)/(lambda^2)
  nu_com <- complete_lm_res
  nu_obs <- ( (nu_com+1)/(nu_com+3) ) * nu_com * (1-lambda)
  nu <- (nu_old * nu_obs) /(nu_old + nu_obs)
  
  bart_interval_tot <- matrix(c(bart_combined_coeffs - qt(0.975,nu)*sqrt(total_diag),
                                bart_combined_coeffs + qt(0.975,nu)*sqrt(total_diag)),
                              nrow=2,byrow = TRUE)
  bart_coverage_tot[j,] <- ifelse(true_beta > bart_interval_tot[1,] & true_beta < bart_interval_tot[2,]
                                  ,TRUE,FALSE)
  
  bart_coeffs[j,] <- bart_combined_coeffs
  bart_withins[j,] <- within_diag
  bart_betweens[j,] <- between_diag
  bart_totals[j,] <- total_diag
  
  cat("\nSimulation",j,"complete\n\n")
  sim_time[j] <- (proc.time() - sim_start)[3]
}

ovt <- (proc.time() - overall_time)[3]

# Save Results

sim.output <- list("true_beta" = true_beta,
                   "var_names" = var_names,
                   
                   "complete_coeffs" = complete_coeffs,
                   "complete_vcovs" = complete_vcovs,
                   
                   "mean_coeffs" = mean_coeffs,
                   "mean_coverage_tot" = mean_coverage_tot,
                   "mean_withins" = mean_withins,
                   "mean_betweens" = mean_betweens,
                   "mean_totals" = mean_totals,
                   
                   "meav_coeffs" = meav_coeffs,
                   "meav_coverage_tot" = meav_coverage_tot,
                   "meav_withins" = meav_withins,
                   "meav_betweens" = meav_betweens,
                   "meav_totals" = meav_totals,
                   
                   "rand_coeffs" = rand_coeffs,
                   "rand_coverage_tot" = rand_coverage_tot,
                   "rand_withins" = rand_withins,
                   "rand_betweens" = rand_betweens,
                   "rand_totals" = rand_totals,
                   
                   "bart_coeffs" = bart_coeffs,
                   "bart_coverage_tot" = bart_coverage_tot,
                   "bart_withins" = bart_withins,
                   "bart_betweens" = bart_betweens,
                   "bart_totals" = bart_totals,
                   
                   "overall_time" = ovt,
                   "sim_time" = sim_time,
                   "mean_time" = mean_time,
                   "meav_time" = meav_time,
                   "rand_time" = rand_time,
                   "bart_time" = bart_time)

filename <- paste(grad_num,",",sim_num,",",n_obs,",",missp,",",sims,",",M,",",maxiter,",",ntree,",",num_burn_in,",",num_iterations_after_burn_in,sep = "")
filename <- paste(filename, ".rds", sep = "")
setwd("/Users/toddb/Library/Mobile Documents/com~apple~CloudDocs/University/3rd Year/Project/Simulations/Sim3/SimsGrad")
saveRDS(sim.output, file = filename)