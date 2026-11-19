test_that("getOms builds the expected omega grid", {
  distmat <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  W <- diag(2)
  c0 <- log(2)
  cmax <- 2 * c0
  cgridfac <- 2

  Oms <- .getOms(distmat, c0, cmax, W, cgridfac)

  expect_length(Oms, 2)
  expect_equal(Oms[[1]], matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = TRUE), tolerance = 1e-12)
  expect_equal(Oms[[2]], matrix(c(1, 0.25, 0.25, 1), nrow = 2, byrow = TRUE), tolerance = 1e-12)
})
test_that("sampled covariance validation is scale invariant and accepts rank deficiency", {
  for (scale in c(1e-110, 1, 1e110)) {
    expect_no_error(.check_sampled_covariance(scale * matrix(1, 2, 2)))
    expect_error(
      .check_sampled_covariance(scale * matrix(c(0, 1, 1, 0), 2, 2)),
      class = "spatial_invalid_sampled_covariance"
    )
  }
  expect_identical(.scpc_score_sign(c(-1e-100, -0.0, 0, 1e-100)), c(-1, 1, 1, 1))
})
