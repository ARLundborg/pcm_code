suppressPackageStartupMessages({
  library(mgcv)
  library(ranger)
})

sim_regression_errors <- function(n, method) {
  f <- function(z) sin(2 * pi * z)
  d <- 7

  Z <- MASS::mvrnorm(n, rep(0, d), diag(1, d))
  Y <- f(Z[, 1]) + rnorm(n)

  Z_ <- MASS::mvrnorm(n, rep(0, d), diag(1, d))
  Y_ <- f(Z_[, 1]) + rnorm(n)

  if (method == "gam") {
    k <- floor((n - 1) / d)
    m_formula <- paste0("Y ~ 1 +", paste0(sapply(1:d,
      function(x) paste0("s(Z[,", x, "], k=k)")), collapse = "+"))
    m <- gam(as.formula(m_formula))
    in_mse <- mean((f(Z[, 1]) - fitted(m))^2)
    out_mse <- mean((f(Z_[, 1]) - predict(m, list(Z = Z_)))^2)
  } else if (method == "ranger") {
    W <- cbind(Y, Z)
    colnames(W) <- c("Y", 1:d)
    m <- ranger(data = W, dependent.variable.name = "Y",
                num.tree = 500, mtry = d)
    in_mse <- mean((f(Z[, 1]) - m$predictions)^2)
    W_ <- Z_
    colnames(W_) <- as.character(1:d)
    out_mse <- mean((f(Z_[, 1]) - predict(m, W_)$predictions)^2)
  } else{
    stop("Method must be either gam or ranger.")
  }
  return(c(in_mse = in_mse, out_mse = out_mse))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  print(sim_regression_error(as.integer(args[1]), args[2]))
}
