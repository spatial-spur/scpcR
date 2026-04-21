test_that("get_fixest_iv_design extracts aligned design matrices with exogenous regressors", {
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    set.seed(1001)
    n <- 80
    z <- stats::rnorm(n)
    w <- stats::rnorm(n)
    u <- stats::rnorm(n)
    x <- 0.8 * z + 0.4 * w + u
    y <- 1 + 1.2 * x + 0.5 * w + u
    dat <- data.frame(y = y, x = x, w = w, z = z)

    fit <- fixest::feols(y ~ w | x ~ z, data = dat)
    design <- .get_fixest_iv_design(fit)

    expect_equal(nrow(design$X), nobs(fit))
    expect_equal(nrow(design$Z), nobs(fit))
    expect_equal(nrow(design$model_mat), nobs(fit))
    expect_equal(colnames(design$model_mat), names(stats::coef(fit)))
    expect_equal(colnames(design$X), c("(Intercept)", "w", "x"))
    expect_equal(colnames(design$Z), c("(Intercept)", "w", "z"))
    expect_identical(design$coef_names, names(stats::coef(fit)))
    expect_false(design$has_fixef)
    expect_null(design$fixef_id)
  })
})

test_that("get_fixest_iv_design handles intercept-plus-endogenous specifications", {
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    set.seed(1002)
    n <- 60
    z <- stats::rnorm(n)
    u <- stats::rnorm(n)
    x <- 0.9 * z + u
    y <- 2 + 1.1 * x + u
    dat <- data.frame(y = y, x = x, z = z)

    fit <- fixest::feols(y ~ 1 | x ~ z, data = dat)
    design <- .get_fixest_iv_design(fit)

    expect_equal(colnames(design$model_mat), names(stats::coef(fit)))
    expect_equal(colnames(design$X), c("(Intercept)", "x"))
    expect_equal(colnames(design$Z), c("(Intercept)", "z"))
  })
})

test_that("get_fixest_iv_design retains fixed-effect identifiers", {
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    set.seed(1003)
    n_fe <- 10
    t_per_fe <- 4
    n <- n_fe * t_per_fe
    fe <- rep(seq_len(n_fe), each = t_per_fe)
    z <- stats::rnorm(n)
    w <- stats::rnorm(n)
    u <- stats::rnorm(n)
    x <- 0.7 * z + 0.3 * w + u
    y <- 1 + 1.1 * x + 0.4 * w + stats::rnorm(n_fe)[fe] + u
    dat <- data.frame(y = y, x = x, w = w, z = z, fe = fe)

    fit <- fixest::feols(y ~ w | fe | x ~ z, data = dat)
    design <- .get_fixest_iv_design(fit)

    expect_true(design$has_fixef)
    expect_false(is.null(design$fixef_id))
    expect_equal(colnames(design$model_mat), names(stats::coef(fit)))
    expect_equal(colnames(design$X), c("w", "x"))
    expect_equal(colnames(design$Z), c("w", "z"))
  })
})

test_that("get_fixest_iv_design supports multi-way fixed effects", {
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    set.seed(1004)
    n1 <- 4
    n2 <- 5
    reps <- 3
    n <- n1 * n2 * reps
    fe1 <- rep(rep(seq_len(n1), each = n2), each = reps)
    fe2 <- rep(rep(seq_len(n2), times = n1), each = reps)
    z <- stats::rnorm(n)
    w <- stats::rnorm(n)
    u <- stats::rnorm(n)
    x <- 0.7 * z + 0.2 * w + u
    y <- 1 + 1.2 * x + 0.3 * w + stats::rnorm(n1)[fe1] + stats::rnorm(n2)[fe2] + u
    dat <- data.frame(y = y, x = x, w = w, z = z, fe1 = fe1, fe2 = fe2)

    fit <- fixest::feols(y ~ w | fe1 + fe2 | x ~ z, data = dat)
    design <- .get_fixest_iv_design(fit)

    expect_true(design$has_fixef)
    expect_length(design$fixef_id, 2L)
    expect_equal(colnames(design$model_mat), names(stats::coef(fit)))
  })
})
