test_that("gettau returns the expected statistic", {
  y <- c(2, 1)
  W <- matrix(
    c(1 / sqrt(2), 1 / sqrt(2),
      1 / sqrt(2), -1 / sqrt(2)),
    nrow = 2,
    byrow = TRUE
  )

  expect_equal(as.numeric(.gettau(y, W)), 3, tolerance = 1e-12)
})
