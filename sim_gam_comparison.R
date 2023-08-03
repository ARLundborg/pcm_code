source("test_functions.R")
suppressPackageStartupMessages({
    library(mgcv)
})

sim_gam_comparison <- function(n, setting) {
    d <- 7
    f <- function(z) sin(2 * pi * z)
    expit <- function(z) 1 / (1 + exp(-z))
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
    k <- max(floor((n - 1) / (d + 1)), d+1)
    m_formula <- paste0("Y ~ 1+", "s(X, k=k) +",
        paste0(sapply(1:d, function(x) {
            paste0("s(Z[,", x, "], k=k)")
        }), collapse = "+"))

    m <- gam(as.formula(m_formula))
    p_gam <- summary(m)$s.table[1, 4]

    pcms <- sapply(1:6, function(i) {
        pcm_test(Y, X, Z, gam_reg_method, gam_g_hat)
    })


    p_pcm <- 1 - pnorm(pcms[1])
    p_pcm_avg <- 1 - pnorm(mean(pcms))


    p_gcm <- gcm_test(Y, X, Z, gam_reg_method)

    suppressWarnings(
        if (n > 250) {
            p_kci <- kci_test(Y, X, Z, GP = FALSE)
        } else {
            p_kci <- kci_test(Y, X, Z, GP = TRUE)
        }
    )

    p_wgsc_sep <- wgsc(Y, X, Z, gam_reg_method)
    p_wgsc_seq <- wgsc(Y, X, Z, gam_reg_method, sequential = TRUE)
    p_wgcm_est <- wGCM_est(Y, X, Z, gam_reg_method)
    p_wgcm_fix <- wGCM_fix(eps, xi, Z, 7) # same weight.num as wgcm paper sims
    return(c(
        p_gam = p_gam, p_pcm = p_pcm, p_pcm_avg = p_pcm_avg, p_gcm = p_gcm,
        p_kci = p_kci, p_wgsc_sep = p_wgsc_sep,
        p_wgsc_seq = p_wgsc_seq, p_wgcm_est = p_wgcm_est,
        p_wgcm_fix = p_wgcm_fix
    ))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
    print(sim_gam_comparison(as.integer(args[1]), as.integer(args[2])))
}