test_that("getOms builds the expected omega grid", {
  distmat <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  W <- diag(2)
  c0 <- log(2)
  cmax <- 2 * c0
  cgridfac <- 2

  Oms <- .getOms(distmat, c0, cmax, W, cgridfac)

  expect_length(Oms, 2)
  expect_equal(Oms[[1]], matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = TRUE), tolerance = 1e-12)
  expect_equal(Oms[[2]], matrix(c(1, 0.25, 0.25, 1), nrow = 2, byrow = TRUE), tolerance = 1e-12)
})
