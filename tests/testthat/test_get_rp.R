test_that("getrp is zero for the zero omega matrix", {
  Om <- matrix(0, 2, 2)

  expect_equal(.getrp(Om, 2), 0, tolerance = 1e-12)
})
