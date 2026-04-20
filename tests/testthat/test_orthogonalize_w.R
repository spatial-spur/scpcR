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
