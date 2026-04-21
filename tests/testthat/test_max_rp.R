test_that("maxrp returns the first zero case when all cases match", {
  Oms <- list(matrix(0, 2, 2), matrix(0, 2, 2))

  out <- .maxrp(Oms, 1, 2)

  expect_equal(out$max, 0, tolerance = 1e-12)
  expect_equal(out$i, 1L)
})
