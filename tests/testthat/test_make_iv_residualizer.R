manual_iv_residualize <- function(X, Z, Y, fixef_id = NULL) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  Y <- as.matrix(Y)

  if (!is.null(fixef_id)) {
    X <- as.matrix(fixest::demean(X, f = fixef_id, nthreads = 1L))
    Z <- as.matrix(fixest::demean(Z, f = fixef_id, nthreads = 1L))
    Y <- as.matrix(fixest::demean(Y, f = fixef_id, nthreads = 1L))
  }

  qrZ <- qr(Z)
  PzX <- qr.fitted(qrZ, X)
  PzY <- qr.fitted(qrZ, Y)
  Y - X %*% solve(crossprod(X, PzX), crossprod(X, PzY))
}

test_that("make_iv_residualizer matches the direct QR-based 2SLS formula", {
  X <- cbind(
    "(Intercept)" = 1,
    w = c(-1, 0, 1, 2, 0.5, -0.5),
    x = c(0.5, 1.0, 1.4, 2.1, -0.1, 0.3)
  )
  Z <- cbind(
    "(Intercept)" = 1,
    w = c(-1, 0, 1, 2, 0.5, -0.5),
    z = c(1.5, 0.2, -0.4, 1.1, 0.3, -0.8)
  )
  Y <- c(2.0, 1.3, 3.1, 4.2, 0.4, 1.0)

  residualize <- .make_iv_residualizer(X, Z)
  expected <- manual_iv_residualize(X, Z, Y)

  expect_equal(residualize(Y), as.numeric(expected), tolerance = 1e-10)
})

test_that("make_iv_residualizer handles multiple Y columns", {
  X <- cbind(
    "(Intercept)" = 1,
    w = c(-1, 0, 1, 2, 0.5, -0.5),
    x = c(0.5, 1.0, 1.4, 2.1, -0.1, 0.3)
  )
  Z <- cbind(
    "(Intercept)" = 1,
    w = c(-1, 0, 1, 2, 0.5, -0.5),
    z = c(1.5, 0.2, -0.4, 1.1, 0.3, -0.8)
  )
  Y <- cbind(
    y1 = c(2.0, 1.3, 3.1, 4.2, 0.4, 1.0),
    y2 = c(-0.5, 0.1, 0.7, 1.4, 0.3, -1.2),
    y3 = c(1.0, 1.0, 2.0, 3.0, 5.0, 8.0)
  )

  residualize <- .make_iv_residualizer(X, Z)
  expected <- manual_iv_residualize(X, Z, Y)

  expect_equal(residualize(Y), expected, tolerance = 1e-10)
})

test_that("make_iv_residualizer supports fixed-effect demeaning", {
  skip_if_not_installed("fixest")

  X <- cbind(
    w = c(-1, 0, 1, 2, 0.5, -0.5),
    x = c(0.5, 1.0, 1.4, 2.1, -0.1, 0.3)
  )
  Z <- cbind(
    w = c(-1, 0, 1, 2, 0.5, -0.5),
    z = c(1.5, 0.2, -0.4, 1.1, 0.3, -0.8)
  )
  Y <- cbind(
    y1 = c(2.0, 1.3, 3.1, 4.2, 0.4, 1.0),
    y2 = c(-0.5, 0.1, 0.7, 1.4, 0.3, -1.2)
  )
  fixef_id <- list(fe = c(1L, 1L, 2L, 2L, 3L, 3L))

  residualize <- .make_iv_residualizer(X, Z, fixef_id = fixef_id)
  expected <- manual_iv_residualize(X, Z, Y, fixef_id = fixef_id)

  expect_equal(residualize(Y), expected, tolerance = 1e-10)
})

test_that("make_iv_residualizer errors when X'P_ZX is rank deficient", {
  X <- cbind(
    "(Intercept)" = 1,
    x1 = c(0, 1, 2, 3),
    x2 = c(0, 1, 2, 3)
  )
  Z <- cbind(
    "(Intercept)" = 1,
    z1 = c(1, 0, 1, 0),
    z2 = c(0, 1, 0, 1)
  )

  expect_error(
    .make_iv_residualizer(X, Z),
    "rank deficient"
  )
})

test_that("IV projection uses the same identifying span after rescaling or duplicate instruments", {
  X <- cbind(1, c(-3, 0, 2, 2, 6))
  Z <- cbind(1, c(0, 4, 1, 2, 7), c(1, 0, 5, 4, 0))
  Q <- qr.Q(qr(Z))
  X_hat <- Q %*% crossprod(Q, X)
  expected <- diag(5) - t(X %*% solve(crossprod(X_hat), t(X_hat)))
  U <- sweep(X, 2L, c(1e130, 1e-130), "*")
  V <- sweep(cbind(Z, Z[, 3]), 2L, c(1e-70, 1e95, 1e-110, 1e95), "*")
  expect_equal(.make_iv_residualizer(U, V, transpose = TRUE)(diag(5)), expected)
})

test_that("FE projection preserves small within-group variation and drops absorbed instruments", {
  fe <- list(group = rep(1:4, each = 3))
  x <- cbind(rep(c(-2, -1, 3), 4))
  z <- cbind(rep(c(-1, -1, 2), 4))
  y <- sin(seq_len(nrow(x)))
  absorbed <- cbind(2^40 * fe$group)

  ordinary <- .make_iv_residualizer(x, z, fixef_id = fe, transpose = TRUE)
  offset_units <- .make_iv_residualizer(
    (x + absorbed) * 2^100, cbind((z + absorbed) * 2^-100, absorbed),
    fixef_id = fe, transpose = TRUE
  )
  expect_equal(offset_units(y), ordinary(y), tolerance = 1e-10)
  expect_error(.make_iv_residualizer(x, absorbed, fixef_id = fe), "rank deficient")

  ols <- .make_ols_residualizer(x, fixef_id = fe)
  ols_scaled <- .make_ols_residualizer((x + absorbed) * 2^100, fixef_id = fe)
  expect_equal(ols_scaled(y), ols(y), tolerance = 1e-10)
  expect_error(.make_ols_residualizer(absorbed, fixef_id = fe), "rank deficient")
})
