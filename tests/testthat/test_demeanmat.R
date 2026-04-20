test_that("demeanmat double-demeans the matrix", {
  mat <- matrix(c(1, 2, 4, 8), nrow = 2, byrow = TRUE)
  expected <- matrix(c(0.75, -0.75, -0.75, 0.75), nrow = 2, byrow = TRUE)

  expect_equal(.demeanmat(mat), expected, tolerance = 1e-12)
})
