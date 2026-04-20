test_that("get_scpc_model_matrix returns intercept and regressor columns", {
  dat <- data.frame(
    y = c(1, 2, 3, 4),
    x = c(0, 1, NA, 3)
  )
  fit <- stats::lm(y ~ x, data = dat)
  expected <- matrix(c(1, 0, 1, 1, 1, 3), nrow = 3, byrow = TRUE)
  mm <- .get_scpc_model_matrix(fit)
  attr(mm, "assign") <- NULL
  dimnames(mm) <- NULL

  expect_equal(mm, expected, tolerance = 1e-12)
})

test_that("get_scpc_model_matrix aligns fixest IV regressors with coefficients", {
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    set.seed(102)
    n <- 120
    z <- stats::rnorm(n)
    w <- stats::rnorm(n)
    u <- stats::rnorm(n)
    x <- 0.9 * z + 0.4 * w + 0.7 * u + stats::rnorm(n, sd = 0.2)
    y <- 1 + 1.2 * x + 0.5 * w + u
    dat <- data.frame(y = y, x = x, w = w, z = z)

    fit <- fixest::feols(y ~ w | x ~ z, data = dat)
    mm <- .get_scpc_model_matrix(fit)

    expect_equal(nrow(mm), nobs(fit))
    expect_equal(ncol(mm), length(stats::coef(fit)))
  })
})
