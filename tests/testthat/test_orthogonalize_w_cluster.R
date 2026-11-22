test_that("clustered conditional projection leaves the simple cluster basis unchanged", {
  W <- matrix(
    c(1 / sqrt(2), 1 / sqrt(2),
      1 / sqrt(2), -1 / sqrt(2)),
    nrow = 2,
    byrow = TRUE
  )
  cl_vec <- factor(c("a", "b"))
  xj_indiv <- c(1, 1)
  model_mat_indiv <- matrix(1, nrow = 2, ncol = 1)
  residualize <- .make_ols_residualizer(model_mat_indiv, include_intercept = FALSE)

  expect_equal(
    .conditional_cluster_weights(W, cl_vec, xj_indiv, residualize),
    W,
    tolerance = 1e-12
  )
})
