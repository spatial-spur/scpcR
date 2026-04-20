test_that("resolve_coords_input selects the euclidean coordinate columns", {
  dat <- data.frame(
    coord_x = c(0, 1, 2, 3),
    coord_y = c(0, 1, 1.5, 2)
  )
  obs_index <- c(1L, 2L, 4L)

  out <- .resolve_coords_input(dat, obs_index, NULL, NULL, c("coord_x", "coord_y"))

  expect_equal(
    unname(out$coords),
    matrix(c(0, 0, 1, 1, 3, 2), ncol = 2, byrow = TRUE),
    tolerance = 1e-12
  )
  expect_false(out$latlong)
})

test_that("resolve_coords_input validates coordinate mode selection", {
  dat <- data.frame(
    lon = c(0, 1),
    lat = c(0, 1),
    coord_x = c(0, 1),
    coord_y = c(0, 1)
  )
  obs_index <- c(1L, 2L)

  expect_error(
    .resolve_coords_input(dat, obs_index, NULL, NULL, NULL),
    "Specify coordinates via `lon`/`lat` or `coords_euclidean`"
  )
  expect_error(
    .resolve_coords_input(dat, obs_index, "lon", "lat", c("coord_x", "coord_y")),
    "Specify either `lon`/`lat` or `coords_euclidean`, not both"
  )
  expect_error(
    .resolve_coords_input(dat, obs_index, "lon", NULL, NULL),
    "For geodesic coordinates, provide both `lon` and `lat`"
  )
})

test_that("resolve_coords_input validates coordinate columns and ranges", {
  dat <- data.frame(
    lon = c(0, 400),
    lat = c(0, 1),
    coord_x = c(0, 1),
    coord_y = c(0, 1)
  )
  obs_index <- c(1L, 2L)

  expect_error(
    .resolve_coords_input(dat, obs_index, "lon", "lat", NULL),
    "Longitude values must be in \\[-180, 180\\]"
  )

  dat_bad_lat <- dat
  dat_bad_lat$lon[2] <- 1
  dat_bad_lat$lat[2] <- 120
  expect_error(
    .resolve_coords_input(dat_bad_lat, obs_index, "lon", "lat", NULL),
    "Latitude values must be in \\[-90, 90\\]"
  )

  dat_bad_type <- dat
  dat_bad_type$lon <- as.character(dat_bad_type$lon)
  expect_error(
    .resolve_coords_input(dat_bad_type, obs_index, "lon", "lat", NULL),
    "`lon` and `lat` must reference numeric columns"
  )

  expect_error(
    .resolve_coords_input(dat, obs_index, "missing_lon", "lat", NULL),
    "Coordinate variables not found in data: missing_lon"
  )
  expect_error(
    .resolve_coords_input(dat, obs_index, NULL, NULL, c("coord_x", "missing")),
    "Coordinate variables not found in data: missing"
  )
  expect_error(
    .resolve_coords_input(dat, obs_index, NULL, NULL, 1),
    "`coords_euclidean` must be a character vector with at least one column name"
  )
})
