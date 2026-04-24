generateSimData <- function(n_obs){
  
  library(MASS)
  
  Sigma <- matrix(c(
    1, 0.5, 0.3, 0.2, 0.3,
    0.5, 1, 0.4, 0.3, 0.4,
    0.3, 0.4, 1, 0.4, 0.3,
    0.2, 0.3, 0.4, 1, 0.4,
    0.3, 0.4, 0.3, 0.4, 1
  ),5,5)
  
  X <- MASS::mvrnorm(n_obs, mu = rep(0,5), Sigma = Sigma)
  
  x1 <- X[,1]
  x2 <- X[,2]
  x3 <- X[,3]
  x4 <- X[,4]
  x5 <- X[,5]
  
  epsilon <- rnorm(n_obs)
  
  true_beta <- c(1,1.5,0.5,1)
  
  y <- true_beta[1]*x1 + true_beta[2]*x2 + true_beta[3]*(x3^2) + true_beta[4]*sin(x4) + epsilon
  
  dataset <- data.frame(y,x1,x2,x3,x4,x5)
  
  return(list(dataset, true_beta))
}

makeMissing <- function(dataset,missp){
  dataset_miss <- dataset
  n <- nrow(dataset)
  
  dataset_miss$x1_mar <- ifelse(runif(n)<missp & dataset_miss$x5>0, NA, dataset_miss$x1)
  
  dataset_miss$x2_mcar <- ifelse(runif(n)<(missp/2), NA, dataset_miss$x2)
  
  dataset_miss$x3_mar <- ifelse(runif(n)<missp & dataset_miss$x5>0, NA, dataset_miss$x3)
  
  dataset_miss$x4_mar <- ifelse(runif(n)<missp & dataset_miss$x5<0, NA, dataset_miss$x4)
  
  dataset_miss <- dataset_miss[,c(1,7,8,9,10,6)]
  
  names(dataset_miss) <- c("y","x1","x2","x3","x4","x5")
  
  return(dataset_miss)
}

fitComplete <- function(dataset){
  
  complete_lm <- lm(y ~ 0 + x1 + x2 + I(x3^2) + I(sin(x4)), data = dataset)
  
  complete_lm_coeff <- c(complete_lm$coefficients)
  
  complete_lm_res <- complete_lm$df.residual
  
  complete_lm_vcov <- vcov(complete_lm)
  
  return(list(complete_lm_coeff, complete_lm_res, complete_lm_vcov))
}