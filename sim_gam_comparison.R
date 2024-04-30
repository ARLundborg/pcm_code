source("test_functions.R")
suppressPackageStartupMessages({
  library(mgcv)
})

sim_gam_comparison <- function(n, setting) {
  d <- 7
  f <- function(z) sin(2 * pi * z)
  if (setting == 0) {
    Z <- MASS::mvrnorm(n, rep(0, d), diag(1, d))
    X <- f(Z[, 1]) + 0.1 * rnorm(n, sd = 1)
    Y <- f(Z[, 1]) + rnorm(n)
  } else if (setting == 1) {
    Z <- MASS::mvrnorm(n, rep(0, d), diag(1, d))
    X <- f(Z[, 1]) + rnorm(n, sd = 1)
    Y <- f(Z[, 1]) + 0.2 * X^2 + rnorm(n)
  } else if (setting == 2) {
    alpha <- beta <- 1
    Z <- MASS::mvrnorm(n, rep(0, d), diag(1, d))
    X <- f(Z[, 1]) - beta * f(Z[, 1]) * (rgamma(n, alpha, beta) - alpha / beta)
    Y <- f(Z[, 1]) + 0.4 * X^2 + rnorm(n)
  } else if (setting == 3) {
    Z <- MASS::mvrnorm(n, rep(0, d), diag(1, d))
    X <- f(Z[, 1]) + rnorm(n, sd = 1)
    Y <- f(Z[, 1]) + 0.4 * X^2 * Z[, 2] + rnorm(n)
  } else {
    stop("Setting must be either 0, 1, 2 or 3.")
  }
  k <- max(floor((n - 1) / (d + 1)), d + 1)
  m_formula <- paste0(
    "Y ~ 1+", "s(X, k=k) +",
    paste0(sapply(1:d, function(x) {
      paste0("s(Z[,", x, "], k=k)")
    }), collapse = "+")
  )

  m <- gam(as.formula(m_formula))
  p_gam <- summary(m)$s.table[1, 4]

  p_pcms <- sapply(1:6, function(i) {
    pcm_test(Y, X, Z, gam_reg_method, gam_gtilde, gam_vhat_reg_method)
  })

  p_gcm <- gcm_test(Y, X, Z, gam_reg_method)

  suppressWarnings(
    if (n > 250) {
      p_kci <- kci_test(Y, X, Z, GP = FALSE)
    } else {
      p_kci <- kci_test(Y, X, Z, GP = TRUE)
    }
  )

  p_wgsc <- wgsc(Y, X, Z, gam_reg_method)
  p_wgcm_est <- wGCM_est(Y, X, Z, gam_reg_method)
  p_wgcm_fix <- wGCM_fix(Y, X, Z, gam_reg_method)
  return(c(
    p_gam = p_gam, p_pcm1 = p_pcms[1], p_pcm2 = p_pcms[2], p_pcm3 = p_pcms[3],
    p_pcm4 = p_pcms[4], p_pcm5 = p_pcms[5], p_pcm6 = p_pcms[6], p_gcm = p_gcm,
    p_kci = p_kci, p_wgsc = p_wgsc, p_wgcm_est = p_wgcm_est,
    p_wgcm_fix = p_wgcm_fix
  ))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  print(sim_gam_comparison(as.integer(args[1]), as.integer(args[2])))
}