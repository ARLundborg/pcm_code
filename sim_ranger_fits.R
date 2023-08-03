suppressPackageStartupMessages({
  library(ranger)
})

eval_f <- function(X, Z, setting) {
  d <- 7
  f <- function(z) sin(2 * pi * z)
  
  if (setting == 0) {
    true_f <- f(Z[, 1]) * (Z[, 3] + 1)
  } else if (setting == 1) {
    true_f <- f(Z[, 1]) * (Z[, 3] + 1) + 0.04 * X^2
  } else if (setting == 2) {
    true_f <- f(Z[, 1]) * (1 + Z[, 3]) + 0.04 * X^2 
  } else if (setting == 3) {
    true_f <- f(Z[, 1]) * (1 + Z[, 3]) + 0.04 * X^2 * Z[, 2]
  } else {
    stop("Setting must be either 0, 1, 2 or 3.")
  }
  return(true_f)
}

sim_data <- function(n, setting) {
  d <- 7
  f <- function(z) sin(2 * pi * z)
  
  if (setting == 0) {
    Z <- MASS::mvrnorm(n, rep(0, d), diag(1, d))
    X <- f(Z[, 1]) * (Z[, 3] + 1) + rnorm(n, sd = 1)
    eps <- sqrt(0.5 + (X > 0)) * rnorm(n)
    Y <- f(Z[, 1]) * (Z[, 3] + 1) + eps
  } else if (setting == 1) {
    Z <- MASS::mvrnorm(n, rep(0, d), diag(1, d))
    X <- f(Z[, 1]) * (1 + Z[, 3]) + rnorm(n, sd = 1)
    eps <- sqrt(0.5 + (X > 0)) * rnorm(n)
    Y <- f(Z[, 1]) * (1 + Z[, 3]) + 0.04 * X^2 + eps
  } else if (setting == 2) {
    alpha <- beta <- 1
    Z <- MASS::mvrnorm(n, rep(0, d), diag(1, d))
    sigma_X <- -beta * f(Z[, 1]) * (1 + Z[, 3])
    xi <- rgamma(n, alpha, beta) - alpha / beta
    X <- f(Z[, 1]) * (1 + Z[, 3]) + sigma_X * xi
    eps <- sqrt(0.5 + (X > 0)) * rnorm(n)
    Y <- f(Z[, 1]) * (1 + Z[, 3]) + 0.04 * X^2 + eps
  } else if (setting == 3) {
    Z <- MASS::mvrnorm(n, rep(0, d), diag(1, d))
    X <- f(Z[, 1]) * (1 + Z[, 3]) + rnorm(n, sd = 1)
    eps <- sqrt(0.5 + (X > 0)) * rnorm(n)
    Y <- f(Z[, 1]) * (1 + Z[, 3]) + 0.04 * X^2 * Z[, 2] + eps
  } else {
    stop("Setting must be either 0, 1, 2 or 3.")
  }
  return(list(Y=Y, X=X, Z=Z))
}


ranger_reg_method <- function(X, y, mtry) {
  n <- length(y)
  X <- as.matrix(X, nrow = n)
  d <- dim(X)[2]
  W <- cbind(y, X)
  colnames(W) <- c("y", 1:d)
  m <- ranger::ranger(data = W, dependent.variable.name = "y",
                      num.tree = 500, mtry = mtry)
  pred_func <- function(X_new) {
    X_new <- as.matrix(X_new)
    d <- dim(X_new)[2]
    colnames(X_new) <- as.character(1:d)
    as.numeric(predict(m, X_new)$predictions)
  }
  return(pred_func)
}

sim_ranger_performance <- function(n, setting) {
  data_list <- sim_data(n, setting)
  Z <- as.matrix(data_list$Z)
  X <- as.matrix(data_list$X)
  Y <- as.matrix(data_list$Y)
  n <- length(Y)
  
  fitted_reg_XZs <- list()
  fitted_reg_Zs <- list()
  for (i in seq_len(d+1)) {
    fitted_reg_XZs[[i]] <- ranger_reg_method(cbind(X, Z), Y, mtry=i)
  }
  for (i in seq_len(d)) {
    fitted_reg_Zs[[i]] <- ranger_reg_method(Z, Y, mtry=i)
  }
  
  new_data_list <- sim_data(n, setting)
  Z_new <- as.matrix(new_data_list$Z)
  X_new <- as.matrix(new_data_list$X)
  true_f <- eval_f(X_new, Z_new, setting)
  
  MSE_XZ <- lapply(fitted_reg_XZs, function(reg) mean((true_f-reg(cbind(X, Z)))^2))
  MSE_Z <- lapply(fitted_reg_Zs, function(reg) mean((true_f-reg(Z))^2))
  return(list(MSE_XZ=MSE_XZ, MSE_Z=MSE_Z))
}


# no tests, run rangers, check perofmrance, different mtrys