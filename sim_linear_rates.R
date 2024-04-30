source("test_functions.R")
suppressPackageStartupMessages({
  library(vimp)
  library(lmtest)
  library(sandwich)
  library(ranger)
})


sim_linear_rates <- function(n) {
  my_expit <- function(x) 2 / (1 + exp(-3 * x))

  d_Z <- 5
  d_X <- 5
  A <- diag(1, d_X, d_Z)

  Z <- MASS::mvrnorm(n, rep(0, d_Z), diag(1, d_Z))
  X <- tcrossprod(Z, A) + MASS::mvrnorm(n, rep(0, d_X), diag(1, d_X))
  Y <- X %*% 1:5 / (sqrt(n)) + (my_expit(X[, 1]) + my_expit(Z[, 1])) * rnorm(n)
  Y <- as.numeric(Y)

  p_pcm <- pcm_test(Y, X, Z, lm_reg_method, lm_gtilde, ranger_reg_method)
  p_wgsc <- wgsc(Y, X, Z, lm_reg_method)
  p_wgsc_no_x <- wgsc(Y, X, Z, lm_reg_method, no_crossfit = TRUE)

  lm_Z <- lm(Y ~ Z)
  lm_XZ <- lm(Y ~ X + Z)
  p_lmtest <- waldtest(lm_Z, lm_XZ,
    vcov = function(x) vcovHC(x, type = "HC0")
  )$"Pr(>F"[2]
  return(c(p_pcm = p_pcm, p_wgsc = p_wgsc, p_wgsc_no_x = p_wgsc_no_x,
   p_lmtest = p_lmtest))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  print(sim_linear_rates(as.integer(args[1])))
}
