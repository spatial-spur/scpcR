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

test_that("getW matches the leading algebraic eigenspace when RSpectra is used", {
  skip_if_not_installed("RSpectra")

  pts <- structure(
    c(
      0.91480604, 0.93707541, 0.28613953, 0.83044763, 0.64174552,
      0.51909595, 0.73658831, 0.13466660, 0.65699229, 0.70506478
    ),
    dim = c(5L, 2L)
  )
  c0 <- 9.209061347343958
  distmat <- as.matrix(stats::dist(pts, diag = TRUE, upper = TRUE))

  sig_d <- .demeanmat(exp(-c0 * distmat))
  eig <- eigen(sig_d, symmetric = TRUE)
  ord <- order(Re(eig$values), decreasing = TRUE)
  expected <- cbind(
    rep(1, nrow(distmat)) / sqrt(nrow(distmat)),
    eig$vectors[, ord[1:2], drop = FALSE]
  )

  assert_columns_allclose_up_to_sign(.getW(distmat, c0, 2), expected)
})
