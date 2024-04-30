suppressPackageStartupMessages({
  library(vimp)
  library(CondIndTests)
})

pcm_test <- function(Y, X, Z, reg_method, gtilde_method = NULL,
                     vhat_reg_method = NULL, var_min = 0.01,
                     reg_params = list()) {
  #' reg_params is a list containing optional regression parameters for "mhat",
  #' "mtilde", "ghat", "gtilde", "vhat" and "mhat_fhat"


  Y <- as.numeric(Y)
  n <- length(Y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  direction_indices <- sample(1:n, n / 2)
  main_indices <- sample(setdiff(1:n, direction_indices))

  X_dir <- X[direction_indices, ]
  Z_dir <- Z[direction_indices, ]
  Y_dir <- Y[direction_indices]

  ghat <- do.call(reg_method, c(
    list(X = cbind(X_dir, Z_dir), y = Y_dir),
    reg_params[["ghat"]]
  ))

  if (is.null(gtilde_method)) {
    gtilde <- function(X, Z) ghat(cbind(X, Z))
  } else {
    gtilde <- do.call(gtilde_method, c(
      list(X = X_dir, Z = Z_dir, Y = Y_dir),
      reg_params[["gtilde"]]
    ))
  }
  ghat_dir <- ghat(cbind(X_dir, Z_dir))
  gtilde_dir <- gtilde(X_dir, Z_dir)

  mtilde <- do.call(reg_method, c(
    list(X = Z_dir, y = ghat_dir),
    reg_params[["mtilde"]]
  ))
  mtilde_dir <- mtilde(Z_dir)

  htilde_dir <- gtilde_dir - mtilde_dir

  rho <- mean((Y_dir - ghat_dir + gtilde_dir - mtilde_dir) * htilde_dir)
  sgn <- ifelse((rho < 0), -1, 1)

  sqr_resid_dir <- (Y_dir - ghat_dir)^2
  if (is.null(vhat_reg_method)) {
    vtilde <- do.call(
      reg_method,
      c(
        list(X = cbind(X_dir, Z_dir), y = sqr_resid_dir),
        reg_params[["vhat"]]
      )
    )
  } else {
    vtilde <- do.call(
      vhat_reg_method,
      c(
        list(X = cbind(X_dir, Z_dir), y = sqr_resid_dir),
        reg_params[["vhat"]]
      )
    )
  }
  vtilde_dir <- vtilde(cbind(X_dir, Z_dir))
  a <- function(c) mean(sqr_resid_dir / (pmax(vtilde_dir, 0) + c))

  if (a(0) <= 1) {
    chat <- 0
  } else {
    chat <- uniroot(function(c) a(c) - 1, c(0, 10), extendInt = "yes")$root
  }
  vhat <- function(X, Z) pmax(vtilde(cbind(X, Z)) + chat, var_min)

  Z_main <- Z[main_indices, ]
  X_main <- X[main_indices, ]
  Y_main <- Y[main_indices]

  fhat_main <- sgn * (gtilde(X_main, Z_main) - mtilde(Z_main)) /
    vhat(X_main, Z_main)

  mhat <- do.call(reg_method, c(
    list(X = Z_main, y = Y_main),
    reg_params[["mhat"]]
  ))
  eps <- Y_main - mhat(Z_main)

  mhat_fhat <- do.call(reg_method, c(
    list(X = Z_main, y = fhat_main),
    reg_params[["mhat_fhat"]]
  ))
  xi <- fhat_main - mhat_fhat(Z_main)


  R <- xi * eps
  test_statistic <- sqrt(length(R)) * mean(R) / stats::sd(R)
  return(1 - pnorm(test_statistic))
}

pcm_test_binary <- function(Y, X, Z, reg_method, binary_reg_method,
                            var_min = 0.01, reg_params = list()) {
  #' reg_params is a list containing optional regression parameters for "mhat",
  #' "mtilde", "ghat" and "mhat_fhat"

  Y <- as.numeric(Y)
  n <- length(Y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  direction_indices <- sample(1:n, n / 2)
  main_indices <- sample(setdiff(1:n, direction_indices))

  X_dir <- X[direction_indices, ]
  Z_dir <- Z[direction_indices, ]
  Y_dir <- Y[direction_indices]

  ghat <- do.call(binary_reg_method, c(
    list(X = cbind(X_dir, Z_dir), y = Y_dir),
    reg_params[["ghat"]]
  ))
  gtilde <- function(X, Z) ghat(cbind(X, Z))

  ghat_dir <- ghat(cbind(X_dir, Z_dir))
  gtilde_dir <- gtilde(X_dir, Z_dir)

  mtilde <- do.call(reg_method, c(
    list(X = Z_dir, y = ghat_dir),
    reg_params[["mtilde"]]
  ))
  mtilde_dir <- mtilde(Z_dir)

  htilde_dir <- gtilde_dir - mtilde_dir

  rho <- mean((Y_dir - ghat_dir + gtilde_dir - mtilde_dir) * htilde_dir)
  sgn <- ifelse((rho < 0), -1, 1)

  sqr_resid_dir <- (Y_dir - ghat_dir)^2
  vtilde_dir <- ghat_dir * (1 - ghat_dir)
  a <- function(c) mean(sqr_resid_dir / (pmax(vtilde_dir, 0) + c))
  if (a(0) <= 1) {
    chat <- 0
  } else {
    chat <- uniroot(function(c) a(c) - 1, c(0, 10), extendInt = "yes")$root
  }
  vhat <- function(X, Z) {
    ghat_eval <- ghat(cbind(X, Z))
    return(pmax(ghat_eval * (1 - ghat_eval) + chat, var_min))
  }

  Z_main <- Z[main_indices, ]
  X_main <- X[main_indices, ]
  Y_main <- Y[main_indices]

  fhat_main <- sgn * (gtilde(X_main, Z_main) - mtilde(Z_main)) /
    vhat(X_main, Z_main)

  mhat <- do.call(binary_reg_method, c(
    list(X = Z_main, y = Y_main),
    reg_params[["mhat"]]
  ))
  eps <- Y_main - mhat(Z_main)

  mhat_fhat <- do.call(reg_method, c(
    list(X = Z_main, y = fhat_main),
    reg_params[["mhat_fhat"]]
  ))
  xi <- fhat_main - mhat_fhat(Z_main)

  R <- xi * eps
  test_statistic <- sqrt(length(R)) * mean(R) / stats::sd(R)
  return(1 - pnorm(test_statistic))
}

wgsc <- function(
    Y, X, Z, reg_method, no_crossfit = FALSE,
    reg_params = list()) {
  #' reg params is a list containing optional regression parameters for "ghat",
  #' "mtilde"
  n <- length(Y)
  Z <- as.matrix(Z)
  X <- as.matrix(X)

  full_fitted <- numeric(n)
  reduced_fitted <- numeric(n)
  data_folds <- sample(c(
    rep(1, floor(n / 4)), rep(2, floor(n / 4)),
    rep(3, floor(n / 4)), rep(4, n - 3 * floor(n / 4))
  ))
  sample_splitting_folds <- vimp::make_folds(unique(data_folds), V = 2)
  if (no_crossfit) {
    full_fitted <- do.call(reg_method, c(
      list(X = cbind(X, Z), y = Y),
      reg_params[["ghat"]]
    ))(cbind(X, Z))
    reduced_fitted <- do.call(reg_method, c(
      list(X = Z, y = full_fitted),
      reg_params[["mtilde"]]
    ))(Z)
  } else {
    for (j in 1:4) {
      full_fit <- do.call(
        reg_method,
        c(
          list(
            X = cbind(X, Z)[data_folds != j, ],
            y = Y[data_folds != j]
          ),
          reg_params[["ghat"]]
        )
      )
      full_fitted[data_folds == j] <- full_fit(cbind(X, Z)[data_folds == j, ])
      reduced_fit <- do.call(
        reg_method,
        c(
          list(
            X = Z[data_folds != j, ],
            y = full_fit(cbind(X, Z)[data_folds != j, ])
          ),
          reg_params[["mtilde"]]
        )
      )
      reduced_fitted[data_folds == j] <- reduced_fit(Z[data_folds == j, ])
    }
  }


  suppressWarnings(
    est <- vimp::cv_vim(
      Y = Y, cross_fitted_f1 = full_fitted,
      cross_fitted_f2 = reduced_fitted, V = 2, type = "r_squared",
      cross_fitting_folds = data_folds,
      sample_splitting_folds = sample_splitting_folds,
      run_regression = FALSE, alpha = 0.05
    )
  )
  return(est$p_value)
}


wgsc_binary <- function(
    Y, X, Z, reg_method, binary_reg_method,
    reg_params = list()) {
  #' reg params is a list containing optional regression parameters for "ghat",
  #' "mtilde"
  n <- length(Y)
  Z <- as.matrix(Z)
  X <- as.matrix(X)

  full_fitted <- numeric(n)
  reduced_fitted <- numeric(n)
  data_folds <- sample(c(
    rep(1, floor(n / 4)), rep(2, floor(n / 4)),
    rep(3, floor(n / 4)), rep(4, n - 3 * floor(n / 4))
  ))
  sample_splitting_folds <- vimp::make_folds(unique(data_folds), V = 2)

  for (j in 1:4) {
    full_fit <- do.call(
      binary_reg_method,
      c(
        list(
          X = cbind(X, Z)[data_folds != j, ],
          y = Y[data_folds != j]
        ),
        reg_params[["ghat"]]
      )
    )
    full_fitted[data_folds == j] <- full_fit(cbind(X, Z)[data_folds == j, ])
    reduced_fit <- do.call(
      reg_method,
      c(
        list(
          X = Z[data_folds != j, ],
          y = full_fit(cbind(X, Z)[data_folds != j, ])
        ),
        reg_params[["mtilde"]]
      )
    )
    reduced_fitted[data_folds == j] <- reduced_fit(Z[data_folds == j, ])
  }

  suppressWarnings(
    est <- vimp::cv_vim(
      Y = Y, cross_fitted_f1 = full_fitted,
      cross_fitted_f2 = reduced_fitted, V = 2, type = "r_squared",
      cross_fitting_folds = data_folds,
      sample_splitting_folds = sample_splitting_folds,
      run_regression = FALSE, alpha = 0.05
    )
  )
  return(est$p_value)
}

wGCM_est <- function(Y, X, Z, reg_method, reg_params = list()) {
  #' reg params is a list containing optional regression parameters for "mhat",
  #' "X_on_Z", "covYX_Z"
  n <- length(Y)

  data_folds <- sample(c(rep(1, floor(0.3 * n)), rep(2, n - floor(0.3 * n))))

  Z1 <- Z[data_folds == 1, ]
  Z2 <- Z[data_folds == 2, ]
  X1 <- X[data_folds == 1]
  X2 <- X[data_folds == 2]
  Y1 <- Y[data_folds == 1]
  Y2 <- Y[data_folds == 2]

  X_on_Z1 <- do.call(reg_method, c(
    list(X = Z1, y = X1),
    reg_params[["X_on_Z"]]
  ))
  eps1 <- X1 - X_on_Z1(Z1)
  mhat1 <- do.call(reg_method, c(list(X = Z1, y = Y1), reg_params[["mhat"]]))
  xi1 <- Y1 - mhat1(Z1)
  covYX_Z <- do.call(reg_method, c(
    list(X = Z1, y = eps1 * xi1),
    reg_params[["covYX_Z"]]
  ))
  W <- sign(covYX_Z(Z2))


  X_on_Z2 <- do.call(reg_method, c(
    list(X = Z2, y = X2),
    reg_params[["X_on_Z"]]
  ))
  eps2 <- X2 - X_on_Z2(Z2)
  mhat2 <- do.call(reg_method, c(list(X = Z2, y = Y2), reg_params[["mhat"]]))
  xi2 <- Y2 - mhat2(Z2)
  R <- eps2 * xi2 * W

  return(1 - pnorm(mean(R) / sd(R) * sqrt(n)))
}

wGCM_est_binary <- function(
    Y, X, Z, reg_method, binary_reg_method,
    reg_params = list()) {
  #' reg params is a list containing optional regression parameters for "mhat",
  #' "X_on_Z", "covYX_Z"
  n <- length(Y)

  data_folds <- sample(c(rep(1, floor(0.3 * n)), rep(2, n - floor(0.3 * n))))

  Z1 <- Z[data_folds == 1, ]
  Z2 <- Z[data_folds == 2, ]
  X1 <- X[data_folds == 1]
  X2 <- X[data_folds == 2]
  Y1 <- Y[data_folds == 1]
  Y2 <- Y[data_folds == 2]

  X_on_Z1 <- do.call(reg_method, c(
    list(X = Z1, y = X1),
    reg_params[["X_on_Z"]]
  ))
  eps1 <- X1 - X_on_Z1(Z1)
  mhat1 <- do.call(binary_reg_method, c(
    list(X = Z1, y = Y1),
    reg_params[["mhat"]]
  ))
  xi1 <- Y1 - mhat1(Z1)
  covYX_Z <- do.call(reg_method, c(
    list(X = Z1, y = eps1 * xi1),
    reg_params[["covYX_Z"]]
  ))
  W <- sign(covYX_Z(Z2))


  X_on_Z2 <- do.call(reg_method, c(
    list(X = Z2, y = X2),
    reg_params[["X_on_Z"]]
  ))
  eps2 <- X2 - X_on_Z2(Z2)
  mhat2 <- do.call(binary_reg_method, c(
    list(X = Z2, y = Y2),
    reg_params[["mhat"]]
  ))
  xi2 <- Y2 - mhat2(Z2)
  R <- eps2 * xi2 * W

  return(1 - pnorm(mean(R) / sd(R) * sqrt(n)))
}


wGCM_fix <- function(Y, X, Z, reg_method, weight.num = 7, reg_params = list()) {
  #' Copied from wgcm.fix
  #' reg params is a list containing optional regression parameters for "mhat",
  #' "X_on_Z"
  n <- length(Y)
  X_on_Z <- do.call(reg_method, c(list(X = Z, y = X), reg_params[["X_on_Z"]]))
  eps <- X - X_on_Z(Z)
  mhat <- do.call(reg_method, c(list(X = Z, y = Y), reg_params[["mhat"]]))
  xi <- Y - mhat(Z)
  nsim <- 2499
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

wGCM_fix_binary <- function(Y, X, Z, reg_method, binary_reg_method,
                            weight.num = 7, reg_params = list()) {
  #' Copied from wgcm.fix
  #' reg params is a list containing optional regression parameters for "mhat",
  #' "X_on_Z"
  n <- length(Y)
  X_on_Z <- do.call(reg_method, c(list(X = Z, y = X), reg_params[["X_on_Z"]]))
  eps <- X - X_on_Z(Z)
  mhat <- do.call(binary_reg_method, c(
    list(X = Z, y = Y),
    reg_params[["mhat"]]
  ))
  xi <- Y - mhat(Z)

  nsim <- 2499
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

kci_test <- function(Y, X, Z, GP = TRUE) {
  return(CondIndTests::KCI(Y, X, Z, GP = GP)$pvalue)
}

gcm_test <- function(Y, X, Z, reg_method, reg_params = list()) {
  #' reg params is a list containing optional regression parameters for "mhat",
  #' "X_on_Z"
  n <- length(Y)
  X_on_Z <- do.call(reg_method, c(list(X = Z, y = X), reg_params[["X_on_Z"]]))
  eps <- X - X_on_Z(Z)
  mhat <- do.call(reg_method, c(list(X = Z, y = Y), reg_params[["mhat"]]))
  xi <- Y - mhat(Z)
  R <- eps * xi

  p_gcm <- 2 * pnorm(-abs(mean(R) / sd(R) * sqrt(n)))
  return(p_gcm)
}

gcm_test_binary <- function(
    Y, X, Z, reg_method, binary_reg_method,
    reg_params = list()) {
  #' Copied from wgcm.fix
  #' reg params is a list containing optional regression parameters for "mhat",
  #' "X_on_Z"
  n <- length(Y)
  X_on_Z <- do.call(reg_method, c(list(X = Z, y = X), reg_params[["X_on_Z"]]))
  eps <- X - X_on_Z(Z)
  mhat <- do.call(binary_reg_method, c(
    list(X = Z, y = Y),
    reg_params[["mhat"]]
  ))
  xi <- Y - mhat(Z)
  R <- eps * xi

  p_gcm <- 2 * pnorm(-abs(mean(R) / sd(R) * sqrt(n)))
  return(p_gcm)
}

gam_gtilde <- function(X, Z, Y, ...) {
  n <- length(Y)
  Z <- as.matrix(Z, nrow = n)
  X <- as.matrix(X, nrow = n)
  d_Z <- dim(Z)[2]
  d_X <- dim(X)[2]
  k <- floor((n - 1) / (d_Z + d_X))
  m_formula <- paste0(
    c(
      "Y ~ 1",
      paste0(sapply(
        1:d_Z,
        function(x) paste0("s(Z[,", x, "], k=", k, ")")
      ), collapse = "+"),
      paste0(sapply(
        1:d_X,
        function(x) paste0("s(X[,", x, "], k=", k, ")")
      ), collapse = "+")
    ),
    collapse = "+"
  )
  m <- mgcv::gam(as.formula(m_formula))
  return(
    function(X_new, Z_new) {
      rowSums(as.matrix(predict(m, list(
        X = as.matrix(X_new),
        Z = as.matrix(Z_new)
      ), type = "terms")[, -(1:d_Z)]))
    }
  )
}

gam_reg_method <- function(X, y, ...) {
  n <- length(y)
  X <- as.matrix(X, nrow = n)
  d <- dim(X)[2]
  k <- floor((n - 1) / d)
  m_formula <- paste0("y ~ 1 +", paste0(sapply(
    1:d,
    function(x) paste0("s(X[,", x, "], k=", k, ")")
  ), collapse = "+"))
  m <- mgcv::gam(as.formula(m_formula))
  return(
    function(X_new) {
      as.numeric(
        predict(m, list(X = as.matrix(X_new, nrow = NROW(X_new))))
      )
    }
  )
}

gam_vhat_reg_method <- function(X, y, ...) {
  n <- length(y)
  X <- as.matrix(X, nrow = n)
  d <- dim(X)[2]
  k <- floor((n - 1) / 2 / d)
  m_formula <- paste0("y ~ 1 +", paste0(sapply(
    1:d,
    function(x) paste0("s(X[,", x, "], k=", k, ")")
  ), collapse = "+"))
  m <- mgcv::gam(as.formula(m_formula), family = Gamma(link = "log"))
  return(
    function(X_new) as.numeric(predict(m, list(X = X_new), type = "response"))
  )
}


lm_reg_method <- function(X, y, ...) {
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

lm_gtilde <- function(X, Z, Y) {
  n <- length(Y)
  Z <- as.matrix(Z, nrow = n)
  X <- as.matrix(X, nrow = n)
  m <- lm(Y ~ Z + X)
  return(
    function(X_new, Z_new) {
      predict(m, list(X = as.matrix(X_new), Z = as.matrix(Z_new)),
        type = "terms"
      )[, 2]
    }
  )
}

ranger_reg_method <- function(X, y, mtry = NULL, max.depth = NULL, ...) {
  n <- length(y)
  X <- as.matrix(X, nrow = n)
  d <- dim(X)[2]
  if (is.null(mtry)) {
    mtry <- d
  }
  W <- cbind(y, X)
  colnames(W) <- c("y", 1:d)
  m <- ranger::ranger(
    data = W, dependent.variable.name = "y",
    num.tree = 500, mtry = mtry, max.depth = max.depth
  )
  pred_func <- function(X_new) {
    X_new <- as.matrix(X_new)
    d <- dim(X_new)[2]
    colnames(X_new) <- as.character(1:d)
    as.numeric(predict(m, X_new)$predictions)
  }
  return(pred_func)
}

gam_reg_method_binary <- function(X, y, ...) {
  n <- length(y)
  X <- as.matrix(X, nrow = n)
  d <- dim(X)[2]
  k <- floor((n - 1) / 2 / d)

  m_formula <- paste0("y ~ 1 +", paste0(sapply(
    1:d,
    function(x) paste0("s(X[,", x, "], k=", k, ")")
  ), collapse = "+"))
  m <- mgcv::gam(as.formula(m_formula), family = binomial)
  pred_func <- function(X_new) {
    as.numeric(predict(m, list(X = X_new), type = "response"))
  }
  attr(pred_func, "model") <- m
  return(
    pred_func
  )
}