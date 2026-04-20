test_that("orthogonalize_W_cluster leaves the simple cluster basis unchanged", {
  W <- matrix(
    c(1 / sqrt(2), 1 / sqrt(2),
      1 / sqrt(2), -1 / sqrt(2)),
    nrow = 2,
    byrow = TRUE
  )
  cl_vec <- factor(c("a", "b"))
  xj_indiv <- c(1, 1)
  model_mat_indiv <- matrix(1, nrow = 2, ncol = 1)

  expect_equal(
    .orthogonalize_W_cluster(
      W,
      cl_vec,
      xj_indiv,
      model_mat_indiv,
      include_intercept = FALSE
    ),
    W,
    tolerance = 1e-12
  )
})
