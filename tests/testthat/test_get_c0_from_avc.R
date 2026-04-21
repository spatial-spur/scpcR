test_that("getc0fromavc inverts the target average correlation", {
  dist <- c(1, 2)
  avc0 <- 0.375

  c0 <- .getc0fromavc(dist, avc0)

  expect_gt(c0, 0)
  expect_equal(c0, log(2), tolerance = 1e-3)
  expect_lt(abs(.getavc(c0, dist) - avc0), 4e-4)
})
