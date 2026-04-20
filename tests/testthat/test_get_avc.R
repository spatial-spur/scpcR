test_that("getavc matches the exponential average", {
  c <- log(2)
  dist <- c(1, 2)

  expect_equal(.getavc(c, dist), 0.375, tolerance = 1e-12)
})
