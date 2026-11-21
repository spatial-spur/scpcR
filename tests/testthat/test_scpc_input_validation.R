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
test_that("GLM prior and working weights are distinguished", {
  d <- data.frame(y = c(1, 3, 2, 6, 8), x = 1:5, w = c(1, 2, 1, 1, 2))

  expect_false(.has_regression_weights(stats::glm(y ~ x, data = d, weights = NULL)))
  expect_false(.has_regression_weights(stats::glm(y ~ x, data = d, family = stats::poisson())))
  expect_true(.has_regression_weights(stats::glm(y ~ x, data = d, weights = w)))
  expect_true(.has_regression_weights(stats::glm(y ~ x, data = d, weights = rep(1, 5))))
  grouped <- stats::glm(cbind(y, 9 - y) ~ x, data = d, family = stats::binomial())
  expect_true(.has_regression_weights(grouped))
  expect_error(scpc(grouped, data = d, coords_euclidean = "x"), "regression weights")
})

test_that("GLM weight detection uses only fitted prior weights and accepts a named NULL", {
  d <- data.frame(y = c(1, 3, 2, NA, 4), x = 1:5, w = 1)
  w <- NULL
  w_null <- NULL
  excluded <- stats::glm(y ~ x, data = d, na.action = stats::na.exclude)
  expect_false(.has_regression_weights(excluded))
  no_frame <- stats::glm(y ~ x, data = d, weights = w_null, model = FALSE)
  expect_false(.has_regression_weights(no_frame, d))
  with_column <- stats::glm(y ~ x, data = d, weights = w, model = FALSE)
  expect_true(.has_regression_weights(with_column, d))
})

test_that("cluster coordinates are checked observation by observation", {
  data <- data.frame(y = sin(1:12), x = cos(1:12), site = rep(1:4, each = 3),
                     loc = c(0, 1, -1, rep(2:4, each = 3)))
  fit <- lm(y ~ x, data)
  for (backend in c("fallback", "rust")) {
    expect_error(scpc(fit, data, cluster = "site", coords_euclidean = "loc", backend = backend),
                 "Coordinates vary within clusters")
    missing_id <- data
    missing_id$site[1] <- NA_integer_
    expect_error(scpc(fit, missing_id, cluster = "site", coords_euclidean = "loc", backend = backend),
                 "Cluster IDs must be present")
  }
})
