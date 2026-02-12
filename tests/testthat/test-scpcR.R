make_scpc_data <- function(n = 10, with_na = FALSE) {
  set.seed(42)
  x <- stats::rnorm(n)
  y <- 1 + 0.5 * x + stats::rnorm(n, sd = 0.2)
  lon <- seq(-125, -66, length.out = n) + stats::rnorm(n, sd = 0.01)
  lat <- seq(35, 45, length.out = n) + stats::rnorm(n, sd = 0.01)
  if (with_na) {
    x[3] <- NA_real_
  }
  data.frame(
    y = y,
    x = x,
    lon = lon,
    lat = lat,
    cluster = rep(seq_len(ceiling(n / 2)), each = 2)[seq_len(n)]
  )
}

test_that("scpc returns expected structure for lm in geodesic mode", {
  skip_if_not_installed("geodist")
  dat <- make_scpc_data()
  fit <- stats::lm(y ~ x, data = dat)
  out <- scpc(
    model = fit,
    data = dat,
    lon = "lon",
    lat = "lat",
    k = 2,
    avc = 0.1,
    uncond = TRUE
  )

  expect_s3_class(out, "scpc")
  expect_equal(dim(out$scpcstats), c(2, 6))
  expect_true(all(is.finite(out$scpcstats)))
  expect_null(out$scpccvs)
})

test_that("scpc returns expected structure in euclidean mode", {
  dat <- make_scpc_data()
  fit <- stats::lm(y ~ x, data = dat)

  out_2d <- scpc(
    model = fit,
    data = dat,
    coord_euclidean = c("lon", "lat"),
    k = 2,
    avc = 0.1,
    uncond = TRUE
  )
  out_1d <- scpc(
    model = fit,
    data = dat,
    coord_euclidean = "lon",
    k = 1,
    avc = 0.1,
    uncond = TRUE
  )

  expect_equal(dim(out_2d$scpcstats), c(2, 6))
  expect_equal(dim(out_1d$scpcstats), c(1, 6))
  expect_true(all(is.finite(out_2d$scpcstats)))
  expect_true(all(is.finite(out_1d$scpcstats)))
})

test_that("scpc computes critical values when requested", {
  skip_if_not_installed("geodist")
  dat <- make_scpc_data()
  fit <- stats::lm(y ~ x, data = dat)
  out_geo <- scpc(
    model = fit,
    data = dat,
    lon = "lon",
    lat = "lat",
    k = 1,
    avc = 0.1,
    uncond = TRUE,
    cvs = TRUE
  )
  out_euc <- scpc(
    model = fit,
    data = dat,
    coord_euclidean = c("lon", "lat"),
    k = 1,
    avc = 0.1,
    uncond = TRUE,
    cvs = TRUE
  )

  expect_equal(dim(out_geo$scpccvs), c(1, 4))
  expect_equal(dim(out_euc$scpccvs), c(1, 4))
  expect_true(all(is.finite(out_geo$scpccvs)))
  expect_true(all(is.finite(out_euc$scpccvs)))
})

test_that("scpc supports full-length cluster vectors with missing rows", {
  dat <- make_scpc_data(with_na = TRUE)
  fit <- stats::lm(y ~ x, data = dat)
  out <- scpc(
    model = fit,
    data = dat,
    coord_euclidean = c("lon", "lat"),
    cluster = dat$cluster,
    k = 1,
    avc = 0.1,
    uncond = TRUE
  )

  expect_s3_class(out, "scpc")
  expect_equal(nrow(out$scpcstats), 1)
})

test_that("cvs option does not change scpcstats", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)

  out_no <- scpc(
    model = fit,
    data = dat,
    coord_euclidean = c("lon", "lat"),
    k = 2,
    avc = 0.1,
    uncond = FALSE,
    cvs = FALSE
  )
  out_yes <- scpc(
    model = fit,
    data = dat,
    coord_euclidean = c("lon", "lat"),
    k = 2,
    avc = 0.1,
    uncond = FALSE,
    cvs = TRUE
  )

  expect_equal(out_no$scpcstats, out_yes$scpcstats, tolerance = 1e-10)
})

test_that("geodesic mode works with explicit custom lon/lat names", {
  skip_if_not_installed("geodist")
  dat <- make_scpc_data(n = 24)

  out_ref <- scpc(
    model = stats::lm(y ~ x, data = dat),
    data = dat,
    lon = "lon",
    lat = "lat",
    k = 2,
    avc = 0.1,
    uncond = TRUE,
    cvs = TRUE
  )

  dat_named <- dat
  names(dat_named)[names(dat_named) == "lon"] <- "my_longitude"
  names(dat_named)[names(dat_named) == "lat"] <- "my_latitude"
  out_named <- scpc(
    model = stats::lm(y ~ x, data = dat_named),
    data = dat_named,
    lon = "my_longitude",
    lat = "my_latitude",
    k = 2,
    avc = 0.1,
    uncond = TRUE,
    cvs = TRUE
  )

  expect_equal(out_ref$scpcstats, out_named$scpcstats, tolerance = 1e-8)
  expect_equal(out_ref$scpccvs, out_named$scpccvs, tolerance = 1e-8)
})

test_that("scpc validates coordinate mode selection and avc bounds", {
  dat <- make_scpc_data()
  fit <- stats::lm(y ~ x, data = dat)

  expect_error(
    scpc(fit, data = dat, avc = 0.1, uncond = TRUE),
    "Specify coordinates via `lon`/`lat` or `coord_euclidean`"
  )
  expect_error(
    scpc(fit, data = dat, lon = "lon", lat = "lat", coord_euclidean = c("lon", "lat"), avc = 0.1, uncond = TRUE),
    "Specify either `lon`/`lat` or `coord_euclidean`, not both"
  )
  expect_error(
    scpc(fit, data = dat, lon = "lon", avc = 0.1, uncond = TRUE),
    "For geodesic coordinates, provide both `lon` and `lat`"
  )
  expect_error(
    scpc(fit, data = dat, coord_euclidean = "lon", avc = 1, uncond = TRUE),
    "Option avc\\(\\) must be in \\(0.001, 0.99\\)"
  )
})

test_that("scpc validates coordinate columns and ranges", {
  skip_if_not_installed("geodist")
  dat <- make_scpc_data()
  fit <- stats::lm(y ~ x, data = dat)

  dat_bad_lon <- dat
  dat_bad_lon$lon[1] <- 400
  expect_error(
    scpc(fit, data = dat_bad_lon, lon = "lon", lat = "lat", avc = 0.1, uncond = TRUE),
    "Longitude values must be in \\[-180, 180\\]"
  )

  dat_bad_lat <- dat
  dat_bad_lat$lat[1] <- 120
  expect_error(
    scpc(fit, data = dat_bad_lat, lon = "lon", lat = "lat", avc = 0.1, uncond = TRUE),
    "Latitude values must be in \\[-90, 90\\]"
  )

  dat_bad_type <- dat
  dat_bad_type$lon <- as.character(dat_bad_type$lon)
  expect_error(
    scpc(fit, data = dat_bad_type, lon = "lon", lat = "lat", avc = 0.1, uncond = TRUE),
    "`lon` and `lat` must reference numeric columns"
  )

  expect_error(
    scpc(fit, data = dat, lon = "missing_lon", lat = "lat", avc = 0.1, uncond = TRUE),
    "Coordinate variables not found in data: missing_lon"
  )
  expect_error(
    scpc(fit, data = dat, coord_euclidean = c("lon", "missing"), avc = 0.1, uncond = TRUE),
    "Coordinate variables not found in data: missing"
  )
  expect_error(
    scpc(fit, data = dat, coord_euclidean = 1, avc = 0.1, uncond = TRUE),
    "`coord_euclidean` must be a character vector with at least one column name"
  )
})
