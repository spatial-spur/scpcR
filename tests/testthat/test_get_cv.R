test_that("getcv matches the t critical value for zero omega matrices", {
  Oms <- list(matrix(0, 2, 2), matrix(0, 2, 2))
  q <- 1
  level <- 0.05

  expect_equal(.getcv(Oms, q, level), stats::qt(1 - level / 2, df = q), tolerance = 1e-10)
})
