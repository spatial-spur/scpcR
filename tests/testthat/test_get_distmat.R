test_that("getdistmat returns Euclidean distances", {
  coords <- matrix(c(0, 0, 3, 4, 6, 8), ncol = 2, byrow = TRUE)
  expected <- matrix(c(0, 5, 10, 5, 0, 5, 10, 5, 0), ncol = 3, byrow = TRUE)

  expect_equal(unname(.getdistmat(coords, FALSE)), expected, tolerance = 1e-12)
})

test_that("getdistmat returns a symmetric geodesic distance matrix", {
  skip_if_not_installed("geodist")
  coords <- matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE)

  distmat <- .getdistmat(coords, TRUE)

  expect_equal(dim(distmat), c(2, 2))
  expect_equal(unname(diag(distmat)), c(0, 0), tolerance = 1e-12)
  expect_true(all(distmat >= 0))
  expect_equal(distmat[1, 2], distmat[2, 1], tolerance = 1e-12)
})
