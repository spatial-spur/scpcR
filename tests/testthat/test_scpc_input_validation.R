test_that("scpc validates scalar option inputs", {
  dat <- make_scpc_data()
  fit <- stats::lm(y ~ x, data = dat)

  expect_error(
    scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = c(0.1, 0.2), uncond = TRUE),
    "`avc` must be a single finite numeric value"
  )
  expect_error(
    scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = NA_real_, uncond = TRUE),
    "`avc` must be a single finite numeric value"
  )
  expect_error(
    scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = 0.001, uncond = TRUE),
    "`avc` must lie in \\(0.001, 0.99\\)"
  )
  expect_error(
    scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = 0.1, ncoef = 1.5, uncond = TRUE),
    "`ncoef` must be NULL or a single positive integer"
  )
  expect_error(
    scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = 0.1, ncoef = 0, uncond = TRUE),
    "`ncoef` must be NULL or a single positive integer"
  )
  expect_error(
    scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = 0.1, uncond = NA),
    "`uncond` must be a single TRUE/FALSE value"
  )
  expect_error(
    scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = 0.1, cvs = NA, uncond = TRUE),
    "`cvs` must be a single TRUE/FALSE value"
  )
  expect_error(
    scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = 0.1, cluster = "", uncond = TRUE),
    "`cluster` must be a single non-empty column name"
  )
})
