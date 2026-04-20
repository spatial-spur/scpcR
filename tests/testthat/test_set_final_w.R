test_that("setfinalW returns the only admissible projection when qmax is one", {
  Oms <- list(matrix(0, 2, 2), matrix(0, 2, 2))
  W <- diag(2)

  out <- .setfinalW(Oms, W, 1)

  expect_equal(out$W, W, tolerance = 1e-12)
  expect_equal(out$cv, stats::qt(0.975, df = 1), tolerance = 1e-10)
  expect_equal(out$q, 1L)
})
