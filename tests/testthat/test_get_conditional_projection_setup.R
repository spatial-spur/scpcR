test_that("get_conditional_projection_setup keeps the baseline matrix for lm", {
  dat <- data.frame(y = c(1, 2, 4), x = c(0, 1, 3))
  fit <- stats::lm(y ~ x, data = dat)
  model_mat <- matrix(c(1, 0, 1, 1, 1, 3), nrow = 3, byrow = TRUE)

  setup <- .get_conditional_projection_setup(fit, model_mat, n = 3, uncond = FALSE)

  expect_equal(setup$model_mat, model_mat, tolerance = 1e-12)
  expect_true(setup$include_intercept)
  expect_null(setup$fixef_id)
})

test_that("get_conditional_projection_setup prepares demeaned regressors for fixest FE models", {
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    set.seed(2026)
    n_fe <- 10
    t_per_fe <- 4
    n <- n_fe * t_per_fe
    fe <- rep(seq_len(n_fe), each = t_per_fe)
    x1 <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    y <- 0.8 * x1 - 0.4 * x2 + stats::rnorm(n_fe)[fe] + stats::rnorm(n, sd = 0.5)
    dat <- data.frame(y = y, x1 = x1, x2 = x2, fe = fe)

    fit <- fixest::feols(y ~ x1 + x2 | fe, data = dat)
    model_mat <- .get_scpc_model_matrix(fit)
    setup <- .get_conditional_projection_setup(fit, model_mat, n = nrow(model_mat), uncond = FALSE)

    expect_equal(nrow(setup$model_mat), nrow(model_mat))
    expect_equal(ncol(setup$model_mat), length(stats::coef(fit)))
    expect_false(setup$include_intercept)
    expect_false(is.null(setup$fixef_id))
  })
})
