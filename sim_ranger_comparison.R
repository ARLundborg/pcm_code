source("test_functions.R")
suppressPackageStartupMessages({
  library(ranger)
})


sim_ranger_comparison <- function(
    n, setting, pcm_reg_params = list(),
    gcm_reg_params = list(), wgsc_reg_params = list(),
    wgcm_est_reg_params = list()) {
  d <- 7
  f <- function(z) sin(pi * z)

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

  p_pcms <- sapply(1:6, function(i) {
    pcm_test(Y, X, Z, ranger_reg_method, reg_params = pcm_reg_params)
  })

  p_gcm <- gcm_test(Y, X, Z, ranger_reg_method, reg_params = gcm_reg_params)

  p_wgsc <- wgsc(Y, X, Z, ranger_reg_method, reg_params = wgsc_reg_params)
  p_wgcm_est <- wGCM_est(Y, X, Z, ranger_reg_method,
    reg_params = wgcm_est_reg_params
  )
  p_wgcm_fix <- wGCM_fix(Y, X, Z, ranger_reg_method,
    reg_params = gcm_reg_params
  )
  return(c(
    p_pcm1 = p_pcms[1], p_pcm2 = p_pcms[2], p_pcm3 = p_pcms[3],
    p_pcm4 = p_pcms[4], p_pcm5 = p_pcms[5], p_pcm6 = p_pcms[6], p_gcm = p_gcm,
    p_wgsc = p_wgsc, p_wgcm_est = p_wgcm_est, p_wgcm_fix = p_wgcm_fix
  ))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  print(sim_ranger_comparison(as.integer(args[1]), as.integer(args[2])))
}