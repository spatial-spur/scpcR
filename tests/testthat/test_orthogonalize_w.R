test_that("orthogonalize_W leaves an already orthogonal basis unchanged", {
  W <- matrix(
    c(1 / sqrt(2), 1 / sqrt(2),
      1 / sqrt(2), -1 / sqrt(2)),
    nrow = 2,
    byrow = TRUE
  )
  xj <- c(1, 1)
  xjs <- c(1, 1)
  model_mat <- matrix(1, nrow = 2, ncol = 1)

  expect_equal(
    .orthogonalize_W(W, xj, xjs, model_mat, include_intercept = FALSE),
    W,
    tolerance = 1e-12
  )
})

test_that("demean_for_scpc preserves binary 1-ULP, single-group offset, and balanced 2-way FE", {
  skip_if_not_installed("fixest")

  # 1. Binary 0/1 contrast shifted by 2^52
  fe_bin <- c(1L, 1L, 2L, 2L, 3L, 3L)
  x_bin <- matrix(c(0, 1, 0, 1, 1, 0), ncol = 1L)
  vals_bin <- seq_len(6L)^2
  exp_bin <- .make_ols_residualizer(x_bin, include_intercept = FALSE, fixef_id = fe_bin)(vals_bin)
  act_bin <- .make_ols_residualizer(x_bin + 2^52, include_intercept = FALSE, fixef_id = fe_bin)(vals_bin)
  expect_equal(act_bin, exp_bin, tolerance = 1e-12)

  # 2. Single-group huge offset (1e15 in group 1, 0.001*z in groups 2..10)
  set.seed(42)
  fe_sg <- rep(1:10, each = 10L)
  z <- stats::rnorm(100L)
  x_clean <- matrix(0.001 * z * (fe_sg > 1L), ncol = 1L)
  x_sg <- matrix(0.001 * z * (fe_sg > 1L) + 1e15 * (fe_sg == 1L), ncol = 1L)
  vals_sg <- stats::rnorm(100L)
  exp_sg <- .make_ols_residualizer(x_clean, include_intercept = FALSE, fixef_id = fe_sg)(vals_sg)
  act_sg <- .make_ols_residualizer(x_sg, include_intercept = FALSE, fixef_id = fe_sg)(vals_sg)
  expect_equal(act_sg, exp_sg, tolerance = 1e-12)

  # 3. Balanced 2-way FE with large offsets in both factors (1e12 * g + 1e6 * h)
  g <- rep(0:4, each = 5L)
  h <- rep(0:4, times = 5L)
  c_vec <- as.numeric((g - 2) * (h - 2))
  x_2w <- matrix(c_vec + 1e12 * g + 1e6 * h, ncol = 1L)
  res_2w <- .demean_for_scpc(x_2w, list(g = g, h = h), "test", drop_absorbed = TRUE)
  expect_equal(as.numeric(res_2w[, 1L]), c_vec, tolerance = 1e-9)
})

test_that("cluster direction normalization handles 1e200 and 1e-200 without overflow or underflow", {
  set.seed(199)
  nclust <- 6L
  cl_vec <- factor(rep(seq_len(nclust), each = 3L))
  n <- length(cl_vec)
  W <- qr.Q(qr(matrix(stats::rnorm(nclust * 3L), nrow = nclust, ncol = 3L)))
  X <- matrix(stats::rnorm(n * 2L), nrow = n, ncol = 2L)
  Z <- matrix(stats::rnorm(n * 2L), nrow = n, ncol = 2L) + 0.5 * X
  xj_base <- stats::rnorm(n)
  iv_res <- .make_iv_residualizer(X, Z, transpose = TRUE)

  base_ols <- .orthogonalize_W_cluster(W, cl_vec, xj_base, X)
  base_iv <- .conditional_cluster_weights(W, cl_vec, xj_base, iv_res)
  for (sc in c(1e200, 1e-200)) {
    scaled_ols <- .orthogonalize_W_cluster(W, cl_vec, xj_base * sc, X)
    scaled_iv <- .conditional_cluster_weights(W, cl_vec, xj_base * sc, iv_res)
    expect_equal(scaled_ols / sc, base_ols, tolerance = 1e-12)
    expect_equal(scaled_iv / sc, base_iv, tolerance = 1e-12)
  }
})

test_that("glm with na.exclude aligns estfun rows with model.matrix", {
  set.seed(303)
  d <- data.frame(
    y = stats::rnorm(25L),
    x = stats::rnorm(25L),
    s1 = stats::runif(25L),
    s2 = stats::runif(25L)
  )
  d$x[c(3L, 11L)] <- NA_real_
  fit_lm <- stats::lm(y ~ x, data = d, na.action = stats::na.exclude)
  fit_glm <- stats::glm(y ~ x, data = d, family = stats::gaussian(), na.action = stats::na.exclude)
  res_lm <- scpc(fit_lm, d, coords_euclidean = c("s1", "s2"), backend = "fallback")
  res_glm <- scpc(fit_glm, d, coords_euclidean = c("s1", "s2"), backend = "fallback")
  expect_equal(res_glm$scpcstats, res_lm$scpcstats, tolerance = 1e-10)
})
