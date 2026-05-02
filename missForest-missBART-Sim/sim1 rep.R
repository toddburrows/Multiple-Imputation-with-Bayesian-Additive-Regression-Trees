# Parameters
sim_num <- 2 #Changing Dataset = 2
n_obs <- 1000
missp <- 0.4 #20% Missingness
sims <- 1000
maxiter <- 10
ntree <- 50
B <- 12
num_burn_in <- 500
num_iterations_after_burn_in <- 200

# Required functions
source('/Users/toddb/Library/Mobile Documents/com~apple~CloudDocs/University/3rd Year/Project/Simulations/SimulationFunctions.R')
source('/Users/toddb/Library/Mobile Documents/com~apple~CloudDocs/University/3rd Year/Project/Simulations/MissBART.R')
library(missForest)
library(parallel)

# Reproducibility
RNGkind("L'Ecuyer-CMRG")
set.seed(191)

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

# For MF
mf_coeffs <- matrix(NA, ncol = length(true_beta), nrow = sims)
mf_dfs <- numeric(sims)
mf_std_errors <- matrix(NA, ncol = length(true_beta), nrow = sims)

# For MB
mb_coeffs <- matrix(NA, ncol = length(true_beta), nrow = sims)
mb_dfs <- numeric(sims)
mb_std_errors <- matrix(NA, ncol = length(true_beta), nrow = sims)

# For MF Bootstrapped
mfBoot_coeffs <- matrix(NA, ncol = length(true_beta), nrow = sims)
mfBoot_withins <- matrix(NA, ncol = length(true_beta), nrow = sims)
mfBoot_betweens <- matrix(NA, ncol = length(true_beta), nrow=sims)
mfBoot_totals <- matrix(NA, ncol = length(true_beta), nrow=sims)
mfBoot_coverage_tot <- matrix(NA, ncol = length(true_beta), nrow=sims)

# For MB Bootstrapped
mbBoot_coeffs <- matrix(NA, ncol = length(true_beta), nrow = sims)
mbBoot_withins <- matrix(NA, ncol = length(true_beta), nrow = sims)
mbBoot_betweens <- matrix(NA, ncol = length(true_beta), nrow=sims)
mbBoot_totals <- matrix(NA, ncol = length(true_beta), nrow=sims)
mbBoot_coverage_tot <- matrix(NA, ncol = length(true_beta), nrow=sims)

# Track Time
sim_time <- numeric(sims)
mf_time <- numeric(sims)
mb_time <- numeric(sims)
mfBoot_time <- numeric(sims)
mbBoot_time <- numeric(sims)

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
  
  # MF
  mf_start <- proc.time()
  
  imp_mf <- missForest(dataset_miss, verbose = FALSE)
  dataset_imp_mf <- imp_mf$ximp
  
  mf_time[j] <- (proc.time() - mf_start)[3]
  
  mf_lm <- lm(y ~ 0 + x1 + x2 + I(x3^2) + I(sin(x4)), data = dataset_imp_mf)
  
  mf_coeffs[j,] <- c(mf_lm$coefficients)
  
  mf_dfs[j] <- mf_lm$df.residual
  
  mf_std_errors[j,] <- c(summary(mf_lm)$coefficients[,2])
  
  # MB
  mb_start <- proc.time()
  
  dataset_imp_mb <- missBART(dataset_miss,
                          maxiter = maxiter,
                          ntree = ntree,
                          num_burn_in = num_burn_in,
                          num_iterations_after_burn_in = num_iterations_after_burn_in,
                          verbose = FALSE)
  
  mb_time[j] <- (proc.time() - mb_start)[3]
  
  mb_lm <- lm(y ~ 0 + x1 + x2 + I(x3^2) + I(sin(x4)), data = dataset_imp_mb)
  
  mb_coeffs[j,] <- c(mb_lm$coefficients)
  
  mb_dfs[j] <- mb_lm$df.residual
  
  mb_std_errors[j,] <- c(summary(mb_lm)$coefficients[,2])
  
  # MF Bootstrapped
  mfBoot_start <- proc.time()
  
  mf_boot_results <- mclapply(1:B, function(i) {
    
    mis_imp_mf <- missForest(
      xmis = dataset_miss[sample(n_obs, replace = TRUE), ],
      verbose = FALSE
    )$ximp
    
    mf_lm <- lm(y ~ 0 + x1 + x2 + I(x3^2) + I(sin(x4)), data = mis_imp_mf)
    
    list(coef = as.numeric(coef(mf_lm)),
         vcov = vcov(mf_lm))
  },
  mc.cores = min(B, detectCores() - 2))
  
  mfBoot_time[j] <- (proc.time() - mfBoot_start)[3]
  
  ## Calculations
  mf_coeffs_mat <- do.call(rbind, lapply(mf_boot_results, `[[`, "coef"))
  
  mf_within <- Reduce(`+`, lapply(mf_boot_results, `[[`, "vcov"))
  
  mf_combined_within <- mf_within / B
  
  mf_combined_coeffs <- colMeans(mf_coeffs_mat)
  
  centre_matrix <- mf_coeffs_mat - matrix(1, ncol = 1, nrow = B) %*% mf_combined_coeffs
  
  mf_combined_between <- 1/(B-1) * t(centre_matrix) %*% centre_matrix
  
  mf_combined_total <- mf_combined_within + (1+1/B) * mf_combined_between
  
  within_diag <- diag(mf_combined_within)
  between_diag <- diag(mf_combined_between)
  total_diag <- diag(mf_combined_total)
  
  lambda <- ((1 + (1/B)) * between_diag) / total_diag
  nu_old <- (B - 1)/(lambda^2)
  nu_com <- complete_lm_res
  nu_obs <- ( (nu_com+1)/(nu_com+3) ) * nu_com * (1-lambda)
  nu <- (nu_old * nu_obs) /(nu_old + nu_obs)
  
  mf_interval_tot <- matrix(c(mf_combined_coeffs - qt(0.975,nu)*sqrt(total_diag),
                              mf_combined_coeffs + qt(0.975,nu)*sqrt(total_diag)),
                            nrow=2,byrow = TRUE)
  mfBoot_coverage_tot[j,] <- ifelse(true_beta > mf_interval_tot[1,] & true_beta < mf_interval_tot[2,]
                                ,TRUE,FALSE)
  
  mfBoot_coeffs[j,] <- mf_combined_coeffs
  mfBoot_withins[j,] <- within_diag
  mfBoot_betweens[j,] <- between_diag
  mfBoot_totals[j,] <- total_diag
  ##
  
  # MB Bootstrapped
  mbBoot_start <- proc.time()
  
  mb_boot_results <- mclapply(1:B, function(i) {
    
    mis_imp_mb <- missBART(
      dataset_miss[sample(n_obs, replace = TRUE), ],
      maxiter = maxiter,
      ntree = ntree,
      num_burn_in = num_burn_in,
      num_iterations_after_burn_in = num_iterations_after_burn_in,
      verbose = FALSE)
    
    mb_lm <- lm(y ~ 0 + x1 + x2 + I(x3^2) + I(sin(x4)), data = mis_imp_mb)
    
    list(coef = as.numeric(coef(mb_lm)),
         vcov = vcov(mb_lm))
  },
  mc.cores = min(B, detectCores() - 2))
  
  mbBoot_time[j] <- (proc.time() - mbBoot_start)[3]
  
  ## Calculations
  mb_coeffs_mat <- do.call(rbind, lapply(mb_boot_results, `[[`, "coef"))
  
  mb_within <- Reduce(`+`, lapply(mb_boot_results, `[[`, "vcov"))
  
  mb_combined_within <- mb_within / B
  
  mb_combined_coeffs <- colMeans(mb_coeffs_mat)
  
  centre_matrix <- mb_coeffs_mat - matrix(1, ncol = 1, nrow = B) %*% mb_combined_coeffs
  
  mb_combined_between <- 1/(B-1) * t(centre_matrix) %*% centre_matrix
  
  mb_combined_total <- mb_combined_within + (1+1/B) * mb_combined_between
  
  within_diag <- diag(mb_combined_within)
  between_diag <- diag(mb_combined_between)
  total_diag <- diag(mb_combined_total)
  
  lambda <- ((1 + (1/B)) * between_diag) / total_diag
  nu_old <- (B - 1)/(lambda^2)
  nu_com <- complete_lm_res
  nu_obs <- ( (nu_com+1)/(nu_com+3) ) * nu_com * (1-lambda)
  nu <- (nu_old * nu_obs) /(nu_old + nu_obs)
  
  mb_interval_tot <- matrix(c(mb_combined_coeffs - qt(0.975,nu)*sqrt(total_diag),
                              mb_combined_coeffs + qt(0.975,nu)*sqrt(total_diag)),
                            nrow=2,byrow = TRUE)
  mbBoot_coverage_tot[j,] <- ifelse(true_beta > mb_interval_tot[1,] & true_beta < mb_interval_tot[2,]
                                    ,TRUE,FALSE)
  
  mbBoot_coeffs[j,] <- mb_combined_coeffs
  mbBoot_withins[j,] <- within_diag
  mbBoot_betweens[j,] <- between_diag
  mbBoot_totals[j,] <- total_diag
  ##
  
  cat("\nSimulation",j,"complete\n\n")
  sim_time[j] <- (proc.time() - sim_start)[3]
}

complete_lm_coeff <- colMeans(complete_coeffs)
complete_lm_vcov <- colMeans(complete_vcovs)

names(complete_lm_coeff) <- var_names

# MF Calculations

mf_lower <- mf_coeffs - qt(0.975,mf_dfs)*mf_std_errors
mf_upper <- mf_coeffs + qt(0.975,mf_dfs)*mf_std_errors
mf_coverage_tot <- ifelse(true_beta>=mf_lower & true_beta<=mf_upper,TRUE,FALSE)

mf_bias_tru <- colMeans(mf_coeffs)-true_beta
mf_pb_tru <- 100*abs(mf_bias_tru/true_beta)
mf_coverage_tot <- colMeans(mf_coverage_tot)

# MB Calculations

mb_lower <- mb_coeffs - qt(0.975,mb_dfs)*mb_std_errors
mb_upper <- mb_coeffs + qt(0.975,mb_dfs)*mb_std_errors
mb_coverage_tot <- ifelse(true_beta>=mb_lower & true_beta<=mb_upper,TRUE,FALSE)

mb_bias_tru <- colMeans(mb_coeffs)-true_beta
mb_pb_tru <- 100*abs(mb_bias_tru/true_beta)
mb_coverage_tot <- colMeans(mb_coverage_tot)

# MF Bootstrapped Calculations

mfBoot_bias_tru <- colMeans(mfBoot_coeffs)-true_beta
mfBoot_pb_tru <- 100*abs(mfBoot_bias_tru/true_beta)

mfBoot_within_final <- colMeans(mfBoot_withins)
mfBoot_between_final <- (1+1/B)*colMeans(mfBoot_betweens)
mfBoot_total_final <- colMeans(mfBoot_totals)

mfBoot_coverage_tot <- colMeans(mfBoot_coverage_tot)

# MB Bootstrapped Calculations

mbBoot_bias_tru <- colMeans(mbBoot_coeffs)-true_beta
mbBoot_pb_tru <- 100*abs(mbBoot_bias_tru/true_beta)

mbBoot_within_final <- colMeans(mbBoot_withins)
mbBoot_between_final <- (1+1/B)*colMeans(mbBoot_betweens)
mbBoot_total_final <- colMeans(mbBoot_totals)

mbBoot_coverage_tot <- colMeans(mbBoot_coverage_tot)


ovt <- (proc.time() - overall_time)[3]

# Save Results
sim.output <- list("true_beta" = true_beta,
                   
                   "complete_lm_coeff" = complete_lm_coeff,
                   "complete_lm_vcov" = complete_lm_vcov,
                   
                   "mf_coeffs" = mf_coeffs,
                   "mf_bias_tru" = mf_bias_tru,
                   "mf_pb_tru" = mf_pb_tru,
                   "mf_coverage_tot" = mf_coverage_tot,
                   "mf_std_errors" = mf_std_errors,
                   
                   "mb_coeffs" = mb_coeffs,
                   "mb_bias_tru" = mb_bias_tru,
                   "mb_pb_tru" = mb_pb_tru,
                   "mb_coverage_tot" = mb_coverage_tot,
                   "mb_std_errors" = mb_std_errors,
                   
                   "mfBoot_coeffs" = mfBoot_coeffs,
                   "mfBoot_bias_tru" = mfBoot_bias_tru,
                   "mfBoot_pb_tru" = mfBoot_pb_tru,
                   "mfBoot_coverage_tot" = mfBoot_coverage_tot,
                   "mfBoot_within_final" = mfBoot_within_final,
                   "mfBoot_between_final" = mfBoot_between_final,
                   "mfBoot_total_final" = mfBoot_total_final,
                   
                   "mbBoot_coeffs" = mbBoot_coeffs,
                   "mbBoot_bias_tru" = mbBoot_bias_tru,
                   "mbBoot_pb_tru" = mbBoot_pb_tru,
                   "mbBoot_coverage_tot" = mbBoot_coverage_tot,
                   "mbBoot_within_final" = mbBoot_within_final,
                   "mbBoot_between_final" = mbBoot_between_final,
                   "mbBoot_total_final" = mbBoot_total_final,
                   
                   "overall_time" = ovt,
                   "sim_time" = sim_time,
                   "mf_time" = mf_time,
                   "mb_time" = mb_time,
                   "mfBoot_time" = mfBoot_time,
                   "mbBoot_time" = mbBoot_time)

filename <- paste(sim_num,",",n_obs,",",missp,",",sims,",",maxiter,",",ntree,",",B,",",num_burn_in,",",num_iterations_after_burn_in,sep = "")
filename <- paste(filename, ".rds", sep = "")
setwd("/Users/toddb/Library/Mobile Documents/com~apple~CloudDocs/University/3rd Year/Project/Simulations/Sim1/Sims")
saveRDS(sim.output, file = filename)