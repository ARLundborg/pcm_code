source("test_functions.R")
suppressPackageStartupMessages({
  library(ranger)
})


sim_ranger_comparison <- function(n, setting) {
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

  pcms <- sapply(1:6, function(i) pcm_test(Y, X, Z, ranger_reg_method))

  p_pcm <- 1 - pnorm(pcms[1])
  p_pcm_avg <- 1 - pnorm(mean(pcms))


  eps <- X - ranger_reg_method(Z, X)(Z)
  xi <- Y - ranger_reg_method(Z, Y)(Z)
  R <- eps * xi

  p_gcm <- 2 * pnorm(-abs(mean(R) / sd(R) * sqrt(n)))


  p_wgsc_sep <- wgsc(Y, X, Z, ranger_reg_method)
  p_wgsc_seq <- wgsc(Y, X, Z, ranger_reg_method, sequential = TRUE)
  p_wgcm_est <- wGCM_est(Y, X, Z, ranger_reg_method)
  p_wgcm_fix <- wGCM_fix(eps, xi, Z, 7) ### same weight.num as wgcm paper sim
  return(c(
    p_pcm = p_pcm, p_pcm_avg = p_pcm_avg, p_gcm = p_gcm,
    p_wgsc_sep = p_wgsc_sep, p_wgsc_seq = p_wgsc_seq,
    p_wgcm_est = p_wgcm_est, p_wgcm_fix = p_wgcm_fix
  ))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  print(sim_ranger_comparison(as.integer(args[1]), as.integer(args[2])))
}