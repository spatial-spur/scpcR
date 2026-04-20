test_that("has_fixest_fe is false for lm", {
  dat <- data.frame(y = c(1, 2, 3), x = c(0, 1, 2))
  fit <- stats::lm(y ~ x, data = dat)

  expect_false(.has_fixest_fe(fit))
})

test_that("has_fixest_fe distinguishes fixest models with and without absorbed FE", {
  skip_if_not_installed("fixest")

  dat <- data.frame(
    y = c(1, 2, 3, 4),
    x = c(0, 1, 2, 3),
    fe = c(1, 1, 2, 2)
  )

  fit_plain <- fixest::feols(y ~ x, data = dat)
  fit_fe <- fixest::feols(y ~ x | fe, data = dat)

  expect_false(.has_fixest_fe(fit_plain))
  expect_true(.has_fixest_fe(fit_fe))
})
