test_that("IV weights use the transpose of an oblique residual maker", {
  x <- c(-2, 1, 3, -1, 4, 0, 2, -3)
  z <- c(1, 2, -1, 0, 3, 4, -2, 1)
  X <- cbind(1, x)
  Z <- cbind(1, z, c(2, 0, 1, 3, -2, 4, 1, -1))
  pzx <- Z %*% solve(crossprod(Z), crossprod(Z, X))
  m <- diag(8) - X %*% solve(crossprod(X, pzx), t(pzx))
  expect_gt(max(abs(m - t(m))), 0.1)
  w <- cbind(rep(1 / sqrt(8), 8), sin(1:8), cos(1:8))
  direction <- c(-1, -4, 2, 2, 1, 3, -2, 5)
  residualize <- .make_iv_residualizer(X, Z, transpose = TRUE)
  expected <- w
  expected[, 1] <- w[, 1] * abs(direction)
  expected[, -1] <- sign(direction) * (t(m) %*% (direction * w[, -1]))
  expect_equal(.conditional_weights(w, direction, sign(direction), residualize),
               expected, tolerance = 1e-12)

  cl <- factor(rep(1:4, each = 2))
  wc <- w[c(1, 3, 5, 7), ]
  wi <- wc[as.integer(cl), ]
  normalized <- direction / sqrt(as.numeric(rowsum(direction^2, cl))[as.integer(cl)])
  expected_cluster <- cbind(
    as.numeric(rowsum(wi[, 1] * normalized * direction, cl)),
    rowsum(normalized * (t(m) %*% (direction * wi[, -1])), cl)
  )
  expect_equal(.conditional_cluster_weights(wc, cl, direction, residualize),
               unname(expected_cluster), ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("identity-instrument weights reproduce the independent five-site cutoff", {
  d <- data.frame(s = c(0, 1, 2, 4, 6), x = c(-1, -4, 2, 2, 1))
  d$y <- d$x + as.numeric(d$s == 0)
  ols <- lm(y ~ x, d)
  X <- model.matrix(ols)
  direction <- as.numeric(nrow(d) * (sandwich::bread(ols) / nrow(d))["x", ] %*% t(X))
  spc <- .setOmsWfin(cbind(d$s, 0), avc0 = 0.5, latlong = FALSE, method = "exact")
  residualize <- .make_iv_residualizer(X, X, transpose = TRUE)
  weights <- .conditional_weights(spc$Wfin, direction, sign(direction), residualize)
  q <- ncol(weights) - 1L
  oms <- .getOms(spc$distmat, spc$c0, spc$cmax, weights, 1.2)
  cv <- max(spc$cvfin, .getcv(oms, q, 0.05))
  expect_equal(as.numeric(cv), 17.62678224, tolerance = 1e-7)
  backends <- c("fallback", if (is.loaded("wrap__rust_run_model_scpc", PACKAGE = "scpcR")) "rust")
  for (backend in backends) {
    a <- scpc(ols, d, coords_euclidean = "s", avc = 0.5, cvs = TRUE,
              backend = backend, nthreads = 1L)
    expect_equal(unname(a$scpccvs["x", 3]), cv, tolerance = 1e-8)
  }
})
test_that("clustered absorbed FE weights match explicit dummies", {
  skip_if_not_installed("fixest")
  set.seed(167)
  n <- 60
  d <- data.frame(g = rep(1:6, times = 10), cl = rep(1:30, each = 2),
                  x = rnorm(n), y = rnorm(n))
  d$s <- rep(runif(30), each = 2)
  absorbed <- fixest::feols(y ~ x | g, d, nthreads = 1)
  dummy <- lm(y ~ x + factor(g), d)
  backends <- c("fallback", if (is.loaded("wrap__rust_run_model_scpc", PACKAGE = "scpcR")) "rust")
  for (backend in backends) {
    a <- scpc(absorbed, d, coords_euclidean = "s", cluster = "cl", cvs = TRUE,
              backend = backend, nthreads = 1L)
    b <- scpc(dummy, d, coords_euclidean = "s", cluster = "cl", cvs = TRUE,
              backend = backend, nthreads = 1L)
    expect_equal(unname(a$scpcstats["x", ]), unname(b$scpcstats["x", ]), tolerance = 1e-7)
    expect_equal(unname(a$scpccvs["x", ]), unname(b$scpccvs["x", ]), tolerance = 1e-7)
  }
})

test_that("clustered oblique IV with crossed FE matches explicit FE dummies", {
  skip_if_not_installed("fixest")
  set.seed(168)
  n <- 60L
  d <- data.frame(g = rep(1:6, times = 10), cl = rep(1:30, each = 2),
                  z = rnorm(n), u = rnorm(n))
  d$x <- 0.8 * d$z + d$u
  d$y <- 0.3 * d$x + sin(d$g) + 0.5 * d$u + rnorm(n)
  d$s <- rep(runif(30), each = 2)
  a_fit <- fixest::feols(y ~ 1 | g | x ~ z, d, nthreads = 1L)
  b_fit <- fixest::feols(y ~ factor(g) | x ~ z, d, nthreads = 1L)
  backends <- c("fallback", if (is.loaded("wrap__rust_run_model_scpc", PACKAGE = "scpcR")) "rust")
  for (backend in backends) {
    a <- scpc(a_fit, d, coords_euclidean = "s", cluster = "cl", cvs = TRUE,
              backend = backend, nthreads = 1L)
    b <- scpc(b_fit, d, coords_euclidean = "s", cluster = "cl", cvs = TRUE,
              backend = backend, nthreads = 1L)
    expect_equal(unname(a$scpcstats["fit_x", ]), unname(b$scpcstats["fit_x", ]), tolerance = 1e-7)
    expect_equal(unname(a$scpccvs["fit_x", ]), unname(b$scpccvs["fit_x", ]), tolerance = 1e-7)
  }
})
test_that("a large coefficient level preserves the SCPC variance", {
  set.seed(15)
  noise <- round(16 * rnorm(20)) / 4
  d <- data.frame(s = runif(40), y = c(noise, -noise))
  for (backend in c("fallback", if (is.loaded("wrap__rust_run_model_scpc", PACKAGE = "scpcR")) "rust")) {
    fit <- function(data) scpc(lm(y ~ 1, data), data, coords_euclidean = "s",
                              avc = 0.1, backend = backend, nthreads = 1L)$scpcstats
    base <- fit(d)
    d_shifted <- transform(d, y = y + 1e12)
    shifted <- fit(d_shifted)
    expect_equal(shifted[1, 2], base[1, 2], tolerance = 1e-4)
    expect_equal(shifted[1, 3], shifted[1, 1] / shifted[1, 2])
  }
})

test_that("singleton clusters preserve intervals at an exact-zero direction", {
  x <- c(-2, -1, 0, 1, 2)
  d <- data.frame(x = x, y = x + c(1, -2, 0, 2, -1),
                  s = c(0.774, 0.439, 0.859, 0.697, 0.094), id = seq_along(x))
  fit <- fixest::feols(y ~ x, d)
  for (backend in c("fallback", if (is.loaded("wrap__rust_run_model_scpc", PACKAGE = "scpcR")) "rust")) {
    compute <- function(cluster) scpc(fit, d, coords_euclidean = "s", avc = 0.3,
                                     cluster = cluster, backend = backend, nthreads = 1L)
    expect_equal(compute(NULL)$scpcstats, compute("id")$scpcstats, tolerance = 1e-10)
  }
})
