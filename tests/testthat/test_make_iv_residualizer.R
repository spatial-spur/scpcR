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

test_that("IV rank rejects roundoff when no projected column identifies the design", {
  x <- c(-1, -1, 1, 1, -1, -1, 1, 1)
  z <- c(-1, 1, -1, 1, -1, 1, -1, 1)
  x2 <- rep(c(-1, 1), each = 4) * x
  z2 <- rep(c(-1, 1), each = 4) * z
  expect_error(.build_iv_projector(x + 1e-14 * z, z), "rank deficient")
  expect_error(.build_iv_projector(cbind(x + 1e-14 * z, x2 + 1e-14 * z2), cbind(z, z2)),
               "rank deficient")
  expect_no_error(.build_iv_projector(x + 1e-8 * z, z))
})

test_that("FE projection retains within variation at the input spacing", {
  fe <- list(group = rep(0:2, each = 3))
  x <- c(0, 1, 2, 0, 2, 1, 1, 0, 2)
  y <- seq_along(x)^2
  expected <- .make_ols_residualizer(x, fixef_id = fe)(y)
  expect_equal(.make_ols_residualizer(x + 2^52, fixef_id = fe)(y), expected,
               tolerance = 1e-12)
})

test_that("two-way projection keeps translated variation and rejects absorbed columns", {
  l <- 80L
  a <- c(seq_len(l), seq_len(l - 1L), seq_len(l))
  b <- c(seq_len(l), seq.int(2L, l), seq_len(l))
  x <- c(rep(c(1, 2), length.out = l), numeric(l - 1L), -rep(c(1, 2), length.out = l))
  z <- c(rep(c(2, 1), length.out = l), numeric(l - 1L), -rep(c(2, 1), length.out = l))
  y <- sin(seq_along(x))
  fe <- list(a, b)
  translated <- x + 2^40 * b

  expected <- lm.fit(cbind(model.matrix(~ factor(a) + factor(b)), x), y)$residuals
  expect_equal(.make_ols_residualizer(translated, fixef_id = fe)(y), expected, tolerance = 1e-9)
  expect_equal(.make_ols_residualizer(translated, fixef_id = rev(fe))(y), expected, tolerance = 1e-9)
  expect_equal(.make_ols_residualizer(x + 2^42 * (a + b), fixef_id = fe)(y),
               expected, tolerance = 1e-9)
  expect_error(.make_ols_residualizer(a + 0.75 * b, fixef_id = fe), "rank deficient")
  clean <- .make_iv_residualizer(x, z, fixef_id = fe, transpose = TRUE)(y)
  dirty <- .make_iv_residualizer(
    translated, cbind(z, a + 0.75 * b), fixef_id = fe, transpose = TRUE
  )(y)
  expect_equal(dirty, clean, tolerance = 1e-9)
})

test_that("selected OLS rank is invariant to regressor units", {
  x <- c(-1, 1, -1, 1)
  z <- c(-1, -1, 1, 1)
  for (units in c(1e-200, 1, 1e200)) {
    fit <- .make_ols_residualizer(cbind(x, units * (x + 1e-9*z)), FALSE)
    expect_lt(max(abs(fit(z))), 1e-6)
    for (gap in c(1e-11, 1e-13)) {
      expect_error(.make_ols_residualizer(cbind(x, units * (x + gap*z)), FALSE),
                   "rank deficient")
    }
  }
})

test_that("within contrasts use the spacing and scale of their own group", {
  top <- .Machine$double.xmax
  ulp <- 2^971
  tiny <- 1e-300
  x <- c(top, tiny, top - ulp, 2*tiny, top - 2*ulp, 3*tiny)
  residuals <- .demean_for_scpc(cbind(x), list(c(1,2,1,2,1,2)), "x", drop_absorbed=TRUE)
  expect_equal(residuals[c(1,3,5)] / ulp, c(1,0,-1))
  expect_equal(residuals[c(2,4,6)] / tiny, c(-1,0,1))
  group <- factor(c(3,3,1), levels=1:3)
  expect_equal(.normalize_cluster_direction(c(4,3,0), group), c(.8,.6,1))
  expect_equal(.normalize_cluster_direction(c(0,0,0), group), c(0,0,1))
  expect_equal(.cluster_l2_norms(c(4,3,0), group), c(0,0,5))
})
