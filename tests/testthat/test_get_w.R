test_that("getW builds the expected two-point basis up to sign", {
  distmat <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  expected <- matrix(
    c(1 / sqrt(2), 1 / sqrt(2),
      1 / sqrt(2), -1 / sqrt(2)),
    nrow = 2,
    byrow = TRUE
  )

  assert_columns_allclose_up_to_sign(.getW(distmat, log(2), 1), expected)
})
