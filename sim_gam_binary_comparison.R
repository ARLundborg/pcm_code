source("test_functions.R")
suppressPackageStartupMessages({
  library(mgcv)
})

sim_gam_binary_comparison <- function(n, setting) {
  d <- 7
  f <- function(z) sin(2 * pi * z)
  expit <- function(z) 1 / (1 + exp(-z))
  if (setting == 0) {
    Z <- MASS::mvrnorm(n, rep(0, d), diag(1, d))
    X <- f(Z[, 1]) + 0.1 * rnorm(n, sd = 1)
    Y <- rbinom(n, 1, expit(f(Z[, 1])))
  } else if (setting == 1) {
    Z <- MASS::mvrnorm(n, rep(0, d), diag(1, d))
    X <- f(Z[, 1]) + rnorm(n, sd = 1)
    Y <- rbinom(n, 1, expit(f(Z[, 1]) + 0.25 * X^2))
  } else if (setting == 2) {
    alpha <- beta <- 1
    Z <- MASS::mvrnorm(n, rep(0, d), diag(1, d))
    xi <- rgamma(n, alpha, beta) - alpha / beta
    X <- f(Z[, 1]) - beta * f(Z[, 1]) * xi
    Y <- rbinom(n, 1, expit(f(Z[, 1]) + 0.5 * X^2))
  } else if (setting == 3) {
    Z <- MASS::mvrnorm(n, rep(0, d), diag(1, d))
    X <- f(Z[, 1]) + rnorm(n, sd = 1)
    Y <- rbinom(n, 1, expit(f(Z[, 1]) + 0.5 * X^2 * Z[, 2]))
  } else {
    stop("Setting must be either 0, 1, 2 or 3.")
  }

  k <- floor((n - 1) / 2 / (d + 1))

  m_formula <- paste0(
    "Y ~ 1+", "s(X, k=k) +",
    paste0(sapply(1:d, function(x) {
      paste0("s(Z[,", x, "], k=k)")
    }), collapse = "+")
  )

  m <- gam(as.formula(m_formula), family = binomial)
  p_gam <- summary(m)$s.table[1, 4]


  p_pcms <- sapply(1:6, function(i) {
    pcm_test_binary(
      Y, X, Z, gam_reg_method,
      gam_reg_method_binary
    )
  })


  p_gcm <- gcm_test_binary(Y, X, Z, gam_reg_method, gam_reg_method_binary)


  p_wgsc <- wgsc_binary(Y, X, Z, gam_reg_method, gam_reg_method_binary)
  p_wgcm_est <- wGCM_est_binary(
    Y, X, Z, gam_reg_method,
    gam_reg_method_binary
  )
  p_wgcm_fix <- wGCM_fix_binary(Y, X, Z, gam_reg_method, gam_reg_method_binary)
  return(c(
    p_gam = p_gam, p_pcm1 = p_pcms[1], p_pcm2 = p_pcms[2], p_pcm3 = p_pcms[3],
    p_pcm4 = p_pcms[4], p_pcm5 = p_pcms[5], p_pcm6 = p_pcms[6], p_gcm = p_gcm,
    p_wgsc = p_wgsc, p_wgcm_est = p_wgcm_est, p_wgcm_fix = p_wgcm_fix
  ))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  print(sim_gam_binary_comparison(as.integer(args[1]), as.integer(args[2])))
}