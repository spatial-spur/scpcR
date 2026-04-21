test_that("get_obs_index returns the surviving data index labels", {
  dat <- data.frame(
    y = c(1, 2, 3, 4),
    x = c(0, 1, NA, 3)
  )
  fit <- stats::lm(y ~ x, data = dat)

  expect_equal(.get_obs_index(fit, dat), c(1L, 2L, 4L))
})

test_that("get_obs_index maps custom lm row names back to data positions", {
  dat <- data.frame(
    y = c(1, 2, 3, 4),
    x = c(0, 1, NA, 3),
    row.names = c("a", "b", "c", "d")
  )
  fit <- stats::lm(y ~ x, data = dat)

  expect_equal(.get_obs_index(fit, dat), c(1L, 2L, 4L))
})

test_that("get_obs_index matches numeric-looking row names before integer coercion", {
  dat <- data.frame(
    y = c(1, 2, 3, 4),
    x = c(0, 1, NA, 3),
    row.names = c("10", "20", "30", "40")
  )
  fit <- stats::lm(y ~ x, data = dat)

  expect_equal(.get_obs_index(fit, dat), c(1L, 2L, 4L))
})

test_that("get_obs_index errors when lm rows cannot be aligned back to data", {
  dat <- data.frame(
    y = c(1, 2, 3),
    x = c(0, 1, 2),
    row.names = c("a", "b", "c")
  )
  fit <- stats::lm(y ~ x, data = dat)
  dat_bad <- dat
  rownames(dat_bad) <- c("x", "y", "z")

  expect_error(
    .get_obs_index(fit, dat_bad),
    "Could not align model estimation sample to data rows"
  )
})

test_that("get_obs_index supports fixest model objects", {
  skip_if_not_installed("fixest")

  dat <- data.frame(
    y = c(1, 2, 3, 4),
    x = c(0, 1, NA, 3)
  )
  fit <- fixest::feols(y ~ x, data = dat)

  expect_equal(.get_obs_index(fit, dat), fixest::obs(fit))
})
