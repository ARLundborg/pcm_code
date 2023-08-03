pcm_test <- function(Y, X, Z, reg_method, ghat_method = NULL,
                     vhat_reg_method = NULL) {
  Y <- as.numeric(Y)
  n <- length(Y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  direction_indices <- sample(1:n, n / 2)
  main_indices <- sample(setdiff(1:n, direction_indices))

  X_dir <- X[direction_indices, ]
  Z_dir <- Z[direction_indices, ]
  Y_dir <- Y[direction_indices]

  m_Y_XZ <- reg_method(cbind(X_dir, Z_dir), Y_dir)
  if (is.null(ghat_method)) {
    ghat <- function(Z, X) m_Y_XZ(cbind(X, Z))
  } else {
    ghat <- ghat_method(Z_dir, X_dir, Y_dir)
  }
  m_Y_Z <- reg_method(Z_dir, ghat(Z_dir, X_dir))

  rho <- mean((Y_dir - m_Y_XZ(cbind(Z_dir, X_dir)) +
    ghat(Z_dir, X_dir) -  m_Y_Z(Z_dir)) * (ghat(Z_dir, X_dir) -  m_Y_Z(Z_dir)))
  sgn <- ifelse((rho < 0), -1, 1)

  if (is.null(vhat_reg_method)) {
    vhat <- reg_method(cbind(X_dir, Z_dir),
                       (Y_dir - m_Y_XZ(cbind(Z_dir, X_dir)))^2)
  } else {
    vhat <- vhat_reg_method(cbind(X_dir, Z_dir),
                            (Y_dir - m_Y_XZ(cbind(Z_dir, X_dir)))^2)
  }

  Z_main <- Z[main_indices, ]
  X_main <- X[main_indices, ]
  Y_main <- Y[main_indices]

  hhat_vals <- sgn * (ghat(Z_main, X_main) - m_Y_Z(Z_main)) /
    vhat(cbind(X_main, Z_main))

  xi <- Y_main - reg_method(Z_main, Y_main)(Z_main)
  eps <- hhat_vals - reg_method(Z_main, hhat_vals)(Z_main)


  R <- xi * eps
  test_statistic <- sqrt(length(R)) * mean(R) / stats::sd(R)
  return(test_statistic)
}

pcm_test_binary <- function(Y, X, Z, reg_method, binary_reg_method,
                            vhat_reg_method = NULL) {
  Y <- as.numeric(Y)
  n <- length(Y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  direction_indices <- sample(1:n, n / 2)
  main_indices <- sample(setdiff(1:n, direction_indices))

  X_dir <- X[direction_indices, ]
  Z_dir <- Z[direction_indices, ]
  Y_dir <- Y[direction_indices]

  m_Y_XZ <- binary_reg_method(cbind(X_dir, Z_dir), Y_dir)
  ghat <- function(Z, X) m_Y_XZ(cbind(X, Z))

  m_Y_Z <- reg_method(Z_dir, ghat(Z_dir, X_dir))

  rho <- mean((Y_dir - m_Y_XZ(cbind(Z_dir, X_dir)) +
    ghat(Z_dir, X_dir) - m_Y_Z(Z_dir)) *
    (ghat(Z_dir, X_dir) - m_Y_Z(Z_dir)))
  sgn <- ifelse((rho < 0), -1, 1)

  if (is.null(vhat_reg_method)) {
    vhat <- reg_method(
      cbind(X_dir, Z_dir),
      (Y_dir - m_Y_XZ(cbind(Z_dir, X_dir)))^2
    )
  } else {
    vhat <- vhat_reg_method(
      cbind(X_dir, Z_dir),
      (Y_dir - m_Y_XZ(cbind(Z_dir, X_dir)))^2
    )
  }

  Z_main <- Z[main_indices, ]
  X_main <- X[main_indices, ]
  Y_main <- Y[main_indices]

  hhat_vals <- sgn * (ghat(Z_main, X_main) - m_Y_Z(Z_main)) /
    vhat(cbind(X_main, Z_main))

  xi <- Y_main - binary_reg_method(Z[main_indices, ], Y_main)(Z_main)
  eps <- hhat_vals - reg_method(Z[main_indices, ], hhat_vals)(Z_main)


  R <- xi * eps
  test_statistic <- sqrt(length(R)) * mean(R) / stats::sd(R)
  return(test_statistic)
}


wgsc <- function(Y, X, Z, reg_method, sequential = FALSE, no_crossfit = FALSE) {
  n <- length(Y)
  Z <- as.matrix(Z)
  X <- as.matrix(X)

  full_fitted <- numeric(n)
  reduced_fitted <- numeric(n)
  data_folds <- sample(c(rep(1, floor(n / 4)), rep(2, floor(n / 4)),
                         rep(3, floor(n / 4)), rep(4, n - 3*floor(n / 4))))
  sample_splitting_folds <- vimp::make_folds(unique(data_folds), V = 2)
  if (no_crossfit) {
    full_fitted <- reg_method(cbind(X, Z), Y)(cbind(X, Z))
    reduced_fitted <- reg_method(Z, Y)(Z)
  } else {
    for (j in 1:4) {
      full_fit <- reg_method(cbind(X, Z)[data_folds != j, ], Y[data_folds != j])
      full_fitted[data_folds == j] <- full_fit(cbind(X, Z)[data_folds == j, ])
      if (sequential) {
        reduced_fitted[data_folds == j] <- reg_method(Z[data_folds != j, ],
      full_fit(cbind(X, Z)[data_folds != j, ]))(Z[data_folds == j, ])
      } else {
        reduced_fitted[data_folds == j] <- reg_method(Z[data_folds != j, ],
      Y[data_folds != j])(Z[data_folds == j, ])
      }
    }
  }


  suppressWarnings(
  est <- vimp::cv_vim(Y = Y, cross_fitted_f1 = full_fitted,
                cross_fitted_f2 = reduced_fitted, V = 2, type = "r_squared",
                cross_fitting_folds = data_folds,
                sample_splitting_folds = sample_splitting_folds,
                run_regression = FALSE, alpha = 0.05)
  )
  return(est$p_value)
}

wgsc_binary <- function(Y, X, Z, reg_method, binary_reg_method,
                        sequential = FALSE) {
  n <- length(Y)
  Z <- as.matrix(Z)
  X <- as.matrix(X)

  full_fitted <- numeric(n)
  reduced_fitted <- numeric(n)
  data_folds <- sample(c(rep(1, floor(n / 4)), rep(2, floor(n / 4)),
                         rep(3, floor(n / 4)), rep(4, n - 3*floor(n / 4))))
  sample_splitting_folds <- vimp::make_folds(unique(data_folds), V = 2)

  for (j in 1:4) {
    full_fit <- binary_reg_method(cbind(X, Z)[data_folds != j, ],
                                  Y[data_folds != j])
    full_fitted[data_folds == j] <- full_fit(cbind(X, Z)[data_folds == j, ])
    if (sequential) {
      reduced_fitted[data_folds == j] <- reg_method(Z[data_folds != j, ],
     full_fit(cbind(X, Z)[data_folds != j, ]))(Z[data_folds == j, ])
    } else {
      reduced_fitted[data_folds == j] <- binary_reg_method(Z[data_folds != j, ],
     Y[data_folds != j])(Z[data_folds == j, ])
    }

  }

  suppressWarnings(
  est <- vimp::cv_vim(Y = Y, cross_fitted_f1 = full_fitted,
                cross_fitted_f2 = reduced_fitted, V = 2, type = "r_squared",
                cross_fitting_folds = data_folds,
                sample_splitting_folds = sample_splitting_folds,
                run_regression = FALSE, alpha = 0.05)
  )
  return(est$p_value)
}

wGCM_est <- function(Y, X, Z, reg_method) {
  n <- length(Y)

  data_folds <- sample(c(rep(1, floor(0.3 * n)), rep(2, n - floor(0.3 * n))))

  Z1 <- Z[data_folds == 1, ]
  Z2 <- Z[data_folds == 2, ]
  X1 <- X[data_folds == 1]
  X2 <- X[data_folds == 2]
  Y1 <- Y[data_folds == 1]
  Y2 <- Y[data_folds == 2]

  eps1 <- X1 - reg_method(Z1, X1)(Z1)
  xi1 <- Y1 - reg_method(Z1, Y1)(Z1)
  W <- sign(reg_method(Z1, eps1 * xi1)(Z2))

  eps2 <- X2 - reg_method(Z2, X2)(Z2)
  xi2 <- Y2 - reg_method(Z2, Y2)(Z2)
  R <- eps2 * xi2 * W

  return(1 - pnorm(mean(R) / sd(R) * sqrt(n)))
}

wGCM_est_binary <- function(Y, X, Z, reg_method, binary_reg_method) {
  n <- length(Y)

  data_folds <- sample(c(rep(1, floor(0.3 * n)), rep(2, n - floor(0.3 * n))))

  Z1 <- Z[data_folds == 1, ]
  Z2 <- Z[data_folds == 2, ]
  X1 <- X[data_folds == 1]
  X2 <- X[data_folds == 2]
  Y1 <- Y[data_folds == 1]
  Y2 <- Y[data_folds == 2]

  eps1 <- X1 - reg_method(Z1, X1)(Z1)
  xi1 <- Y1 - binary_reg_method(Z1, Y1)(Z1)
  W <- sign(reg_method(Z1, eps1 * xi1)(Z2))

  eps2 <- X2 - reg_method(Z2, X2)(Z2)
  xi2 <- Y2 - binary_reg_method(Z2, Y2)(Z2)
  R <- eps2 * xi2 * W

  return(1 - pnorm(mean(R) / sd(R) * sqrt(n)))
}


wGCM_fix <- function(eps, xi, Z, weight.num) {
  ## Copied from wgcm.fix
  n <- length(xi)
  nsim <- 499
  W <- weightedGCM:::weight_matrix(Z, weight.num, "sign")
  R <- eps * xi * W
  R <- t(R)
  R.norm <- R / sqrt(rowMeans(R^2) - rowMeans(R)^2)
  T.stat <- sqrt(n) * max(abs(rowMeans(R.norm)))
  T.stat.sim <- apply(abs(R.norm %*% matrix(rnorm(n *
    nsim), n, nsim)), 2, max) / sqrt(n)
  p.value <- (sum(T.stat.sim >= T.stat) + 1) / (nsim +
    1)
  return(p.value)
}

gam_g_hat <- function(Z, X, Y) {
  n <- length(Y)
  Z <- as.matrix(Z, nrow = n)
  X <- as.matrix(X, nrow = n)
  d_Z <- dim(Z)[2]
  d_X <- dim(X)[2]
  k <- floor((n - 1) / (d_Z + d_X))
  m_formula <- paste0(
    c(
      "Y ~ 1",
      paste0(sapply(1:d_Z,
       function(x) paste0("s(Z[,", x, "], k=", k, ")")), collapse = "+"),
      paste0(sapply(1:d_X,
       function(x) paste0("s(X[,", x, "], k=", k, ")")), collapse = "+")
    ),
    collapse = "+"
  )
  m <- mgcv::gam(as.formula(m_formula))
  return(
    function(Z_new, X_new) {
    rowSums(as.matrix(predict(m, list(X = as.matrix(X_new),
    Z = as.matrix(Z_new)), type = "terms")[, -(1:d_Z)]))
    }
  )
}

gam_reg_method <- function(X, y) {
  n <- length(y)
  X <- as.matrix(X, nrow = n)
  d <- dim(X)[2]
  k <- floor((n - 1) / d)
  m_formula <- paste0("y ~ 1 +", paste0(sapply(1:d,
   function(x) paste0("s(X[,", x, "], k=", k, ")")), collapse = "+"))
  m <- mgcv::gam(as.formula(m_formula))
  return(
    function(X_new) as.numeric(predict(m, list(X = X_new)))
  )
}


lm_reg_method <- function(X, y) {
  n <- length(y)
  X <- as.matrix(X, nrow = n)
  m <- lm(y ~ X)
  return(
    function(X_new) {
      X_new <- as.matrix(X_new)
      as.numeric(predict(m, list(X = X_new)))
    }
  )
}

lm_ghat <- function(Z, X, Y) {
  n <- length(Y)
  Z <- as.matrix(Z, nrow = n)
  X <- as.matrix(X, nrow = n)
  m <- lm(Y ~ Z + X)
  return(
    function(Z_new, X_new) {
    predict(m, list(X = as.matrix(X_new), Z = as.matrix(Z_new)),
    type = "terms")[, 2]
    }
  )
}

ranger_reg_method <- function(X, y) {
  n <- length(y)
  X <- as.matrix(X, nrow = n)
  d <- dim(X)[2]
  W <- cbind(y, X)
  colnames(W) <- c("y", 1:d)
  m <- ranger::ranger(data = W, dependent.variable.name = "y",
    num.tree = 500, mtry = d)
  pred_func <- function(X_new) {
    X_new <- as.matrix(X_new)
    d <- dim(X_new)[2]
    colnames(X_new) <- as.character(1:d)
    as.numeric(predict(m, X_new)$predictions)
  }
  return(pred_func)
}

gam_reg_method_binary <- function(X, y) {
  n <- length(y)
  X <- as.matrix(X, nrow = n)
  d <- dim(X)[2]
  k <- floor((n - 1) / 2 / d)

  m_formula <- paste0("y ~ 1 +", paste0(sapply(1:d,
   function(x) paste0("s(X[,", x, "], k=", k, ")")), collapse = "+"))
  m <- mgcv::gam(as.formula(m_formula), family = binomial)
  pred_func <- function(X_new) {
    as.numeric(predict(m, list(X = X_new), type = "response"))
  }
  attr(pred_func, "model") <- m
  return(
    pred_func
  )
}