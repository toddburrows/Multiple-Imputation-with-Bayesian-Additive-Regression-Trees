# Parameters
grad_num <- 5 #simulation number ,i.e. if combining 5 simulations this would be 1 for the first simulation

sim_num <- 2 #Changing Dataset = 2
n_obs <- 1000
missp <- 0.4 #20% Missingness
sims <- 100 #100 #This is number for just this simulation
M <- 12
maxiter <- 10
ntree <- 50
num_burn_in <- 500
num_iterations_after_burn_in <- 200

# Required functions
source('/Users/toddb/Library/Mobile Documents/com~apple~CloudDocs/University/3rd Year/Project/Simulations/SimulationFunctions.R')
source('/Users/toddb/Library/Mobile Documents/com~apple~CloudDocs/University/3rd Year/Project/Simulations/ProperBART.R')
library(parallel)

# Reproducibility
RNGkind("L'Ecuyer-CMRG")
set.seed(as.numeric(paste0(192,grad_num)))

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

# For MICE PMM
pmm_coeffs <- matrix(NA, ncol = length(true_beta), nrow=sims)
pmm_withins <- matrix(NA, ncol = length(true_beta), nrow=sims)
pmm_betweens <- matrix(NA, ncol = length(true_beta), nrow=sims)
pmm_totals <- matrix(NA, ncol = length(true_beta), nrow=sims)
pmm_coverage_tot <- matrix(NA, ncol = length(true_beta), nrow=sims)

# For MICE CART
cart_coeffs <- matrix(NA, ncol = length(true_beta), nrow=sims)
cart_withins <- matrix(NA, ncol = length(true_beta), nrow=sims)
cart_betweens <- matrix(NA, ncol = length(true_beta), nrow=sims)
cart_totals <- matrix(NA, ncol = length(true_beta), nrow=sims)
cart_coverage_tot <- matrix(NA, ncol = length(true_beta), nrow=sims)

# For MICE RF
rf_coeffs <- matrix(NA, ncol = length(true_beta), nrow=sims)
rf_withins <- matrix(NA, ncol = length(true_beta), nrow=sims)
rf_betweens <- matrix(NA, ncol = length(true_beta), nrow=sims)
rf_totals <- matrix(NA, ncol = length(true_beta), nrow=sims)
rf_coverage_tot <- matrix(NA, ncol = length(true_beta), nrow=sims)

# For Proper BART
bart_coeffs <- matrix(NA, ncol = length(true_beta), nrow=sims)
bart_withins <- matrix(NA, ncol = length(true_beta), nrow=sims)
bart_betweens <- matrix(NA, ncol = length(true_beta), nrow=sims)
bart_totals <- matrix(NA, ncol = length(true_beta), nrow=sims)
bart_coverage_tot <- matrix(NA, ncol = length(true_beta), nrow=sims)

# Track Time
sim_time <- numeric(sims)
pmm_time <- numeric(sims)
cart_time <- numeric(sims)
rf_time <- numeric(sims)
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
  
  # MICE PMM
  pmm_start <- proc.time()
  
  pmm_imp <- futuremice(dataset_miss, n.core = min(M, detectCores() - 2),
                        m = M, method = 'pmm', maxit = maxiter)
  
  pmm_time[j] <- (proc.time() - pmm_start)[3]
  
  ## Calculations
  pmm_fit <- with(pmm_imp, lm(y ~ 0 + x1 + x2 + I(x3^2) + I(sin(x4))))
  
  pmm_pool <- pool(pmm_fit)$pooled
  
  pmm_combined_coeffs <- pmm_pool$estimate
  pmm_within <- pmm_pool$ubar
  pmm_between <- pmm_pool$b
  pmm_total <- pmm_pool$t
  pmm_nu <- pmm_pool$df
  
  pmm_interval_tot <- matrix(c(pmm_combined_coeffs - qt(0.975,pmm_nu)*sqrt(pmm_total),
                               pmm_combined_coeffs + qt(0.975,pmm_nu)*sqrt(pmm_total)),
                             nrow=2,byrow = TRUE)
  pmm_coverage_tot[j,] <- ifelse(true_beta > pmm_interval_tot[1,] & true_beta < pmm_interval_tot[2,]
                                  ,TRUE,FALSE)
  
  pmm_coeffs[j,] <- pmm_combined_coeffs
  pmm_withins[j,] <- pmm_within
  pmm_betweens[j,] <- pmm_between
  pmm_totals[j,] <- pmm_total
  
  # MICE CART
  cart_start <- proc.time()
  
  cart_imp <- futuremice(dataset_miss, n.core = min(M, detectCores() - 2),
                        m = M, method = 'cart', maxit = maxiter)
  
  cart_time[j] <- (proc.time() - cart_start)[3]
  
  ## Calculations
  cart_fit <- with(cart_imp, lm(y ~ 0 + x1 + x2 + I(x3^2) + I(sin(x4))))
  
  cart_pool <- pool(cart_fit)$pooled
  
  cart_combined_coeffs <- cart_pool$estimate
  cart_within <- cart_pool$ubar
  cart_between <- cart_pool$b
  cart_total <- cart_pool$t
  cart_nu <- cart_pool$df

  cart_interval_tot <- matrix(c(cart_combined_coeffs - qt(0.975,cart_nu)*sqrt(cart_total),
                               cart_combined_coeffs + qt(0.975,cart_nu)*sqrt(cart_total)),
                             nrow=2,byrow = TRUE)
  cart_coverage_tot[j,] <- ifelse(true_beta > cart_interval_tot[1,] & true_beta < cart_interval_tot[2,]
                                 ,TRUE,FALSE)
  
  cart_coeffs[j,] <- cart_combined_coeffs
  cart_withins[j,] <- cart_within
  cart_betweens[j,] <- cart_between
  cart_totals[j,] <- cart_total
  
  # MICE RF
  rf_start <- proc.time()
  
  rf_imp <- futuremice(dataset_miss, n.core = min(M, detectCores() - 2),
                        m = M, method = 'rf', maxit = maxiter)
  
  rf_time[j] <- (proc.time() - rf_start)[3]
  
  ## Calculations
  rf_fit <- with(rf_imp, lm(y ~ 0 + x1 + x2 + I(x3^2) + I(sin(x4))))
  
  rf_pool <- pool(rf_fit)$pooled
  
  rf_combined_coeffs <- rf_pool$estimate
  rf_within <- rf_pool$ubar
  rf_between <- rf_pool$b
  rf_total <- rf_pool$t
  rf_nu <- rf_pool$df

  rf_interval_tot <- matrix(c(rf_combined_coeffs - qt(0.975,rf_nu)*sqrt(rf_total),
                               rf_combined_coeffs + qt(0.975,rf_nu)*sqrt(rf_total)),
                             nrow=2,byrow = TRUE)
  rf_coverage_tot[j,] <- ifelse(true_beta > rf_interval_tot[1,] & true_beta < rf_interval_tot[2,]
                                 ,TRUE,FALSE)
  
  rf_coeffs[j,] <- rf_combined_coeffs
  rf_withins[j,] <- rf_within
  rf_betweens[j,] <- rf_between
  rf_totals[j,] <- rf_total
  
  # Proper BART
  bart_start <- proc.time()
  
  bart_results <- mclapply(1:M, function(i) {
    
    mis_imp_bart <- properBART(
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
  },
  mc.cores = min(M, detectCores() - 2))
  
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
                   
                   "pmm_coeffs" = pmm_coeffs,
                   "pmm_coverage_tot" = pmm_coverage_tot,
                   "pmm_withins" = pmm_withins,
                   "pmm_betweens" = pmm_betweens,
                   "pmm_totals" = pmm_totals,
                   
                   "cart_coeffs" = cart_coeffs,
                   "cart_coverage_tot" = cart_coverage_tot,
                   "cart_withins" = cart_withins,
                   "cart_betweens" = cart_betweens,
                   "cart_totals" = cart_totals,
                   
                   "rf_coeffs" = rf_coeffs,
                   "rf_coverage_tot" = rf_coverage_tot,
                   "rf_withins" = rf_withins,
                   "rf_betweens" = rf_betweens,
                   "rf_totals" = rf_totals,
                   
                   "bart_coeffs" = bart_coeffs,
                   "bart_coverage_tot" = bart_coverage_tot,
                   "bart_withins" = bart_withins,
                   "bart_betweens" = bart_betweens,
                   "bart_totals" = bart_totals,
                   
                   "overall_time" = ovt,
                   "sim_time" = sim_time,
                   "pmm_time" = pmm_time,
                   "cart_time" = cart_time,
                   "rf_time" = rf_time,
                   "bart_time" = bart_time)

filename <- paste(grad_num,",",sim_num,",",n_obs,",",missp,",",sims,",",M,",",maxiter,",",ntree,",",num_burn_in,",",num_iterations_after_burn_in,sep = "")
filename <- paste(filename, ".rds", sep = "")
setwd("/Users/toddb/Library/Mobile Documents/com~apple~CloudDocs/University/3rd Year/Project/Simulations/Sim2/SimsGrad")
saveRDS(sim.output, file = filename)