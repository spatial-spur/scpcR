test_that("lvech returns the strict lower triangle", {
  mat <- matrix(1:9, nrow = 3, byrow = TRUE)

  expect_equal(.lvech(mat), c(4, 7, 8), tolerance = 1e-12)
})
