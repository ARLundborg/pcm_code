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

  m_formula <- paste0("Y ~ 1+", "s(X, k=k) +",
    paste0(sapply(1:d, function(x) {
      paste0("s(Z[,", x, "], k=k)")
    }), collapse = "+"))

  m <- gam(as.formula(m_formula), family = binomial)
  p_gam <- summary(m)$s.table[1, 4]


  pcms <- sapply(1:6, function(i) {
    pcm_test_binary(
      Y, X, Z, gam_reg_method,
      gam_reg_method_binary
    )
  })

  p_pcm <- 1 - pnorm(pcms[1])
  p_pcm_avg <- 1 - pnorm(mean(pcms))


  eps <- X - gam_reg_method(Z, X)(Z)
  xi <- Y - gam_reg_method_binary(Z, Y)(Z)
  R <- eps * xi

  p_gcm <- 2 * pnorm(-abs(mean(R) / sd(R) * sqrt(n)))


  p_wgsc_sep <- wgsc_binary(Y, X, Z, gam_reg_method, gam_reg_method_binary)
  p_wgsc_seq <- wgsc_binary(Y, X, Z, gam_reg_method, gam_reg_method_binary,
                             sequential = TRUE)
  p_wgcm_est <- wGCM_est_binary(Y, X, Z, gam_reg_method,
    gam_reg_method_binary)
  p_wgcm_fix <- wGCM_fix(eps, xi, Z, 7) ### same weight.num as wgcm paper sims
  return(c(
    p_gam = p_gam, p_pcm = p_pcm, p_pcm_avg = p_pcm_avg, p_gcm = p_gcm,
    p_wgsc_sep = p_wgsc_sep, p_wgsc_seq = p_wgsc_seq,
    p_wgcm_est = p_wgcm_est, p_wgcm_fix = p_wgcm_fix
  ))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  print(sim_gam_binary_comparison(as.integer(args[1]), as.integer(args[2])))
}