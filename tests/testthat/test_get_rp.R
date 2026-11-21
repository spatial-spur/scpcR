test_that("getrp is zero for the zero omega matrix", {
  Om <- matrix(0, 2, 2)

  expect_equal(.getrp(Om, 2), 0, tolerance = 1e-12)
})

test_that("strict Gaussian rejection handles ties and tiny positive directions", {
  for (r in c(1 - 5e-15, -1 + .Machine$double.eps)) {
    om <- matrix(c(1, r, r, 1), 2)
    expect_identical(.getrp(om, 1), 0.5)
    expect_identical(.getrp(rbind(cbind(om, 0), 0), 1), 0.5)
  }
  for (v in list(c(1, 1), c(1, -1), c(5, 3, 4))) {
    expect_identical(.getrp(tcrossprod(v), 1), 0)
    expect_identical(.getrp(tcrossprod(v), 1 - 1e-15), 1)
  }
  for (scale in c(2^-900, 1, 2^900)) {
    a <- 2401 * scale
    b <- 49 * scale
    expect_identical(.getrp(matrix(c(a, b, b, scale), 2), 49), 0)
    b <- b * (1 - .Machine$double.eps)
    expect_identical(.getrp(matrix(c(a, b, b, scale), 2), 49), 0.5)
  }
  near <- tcrossprod(c(5, 3, 4))
  near[1, 1] <- 25 + 2^-48
  expect_equal(.getrp(near, 1), 0.5, tolerance = 2e-7)
  expect_error(.getrp(matrix(c(0, 1, 1, 0), 2), 1), "positive semidefinite")
})
