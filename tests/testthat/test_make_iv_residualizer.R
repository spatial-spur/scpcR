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
