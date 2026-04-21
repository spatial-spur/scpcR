test_that("is_fixest_iv_second_stage recognizes second-stage fixest iv fits", {
  dat <- data.frame(
    y = c(1.0, 2.1, 2.8, 4.2, 5.1, 6.0),
    x = c(0.2, 0.7, 1.0, 1.5, 1.8, 2.2),
    w = c(-0.4, 0.1, 0.6, -0.2, 0.5, 0.9),
    z = c(0.1, 0.5, 0.8, 1.0, 1.3, 1.7),
    fe = c(1, 1, 2, 2, 3, 3)
  )

  with_fixest_single_thread({
    fit_ols <- fixest::feols(y ~ w, data = dat)
    fit_iv <- fixest::feols(y ~ w | x ~ z, data = dat)
    fit_fe_iv <- fixest::feols(y ~ w | fe | x ~ z, data = dat)

    # fix: first-stage fits are iv-related, but they are not the main stage
    # that scpc should treat as the final coefficient problem.
    first_stage <- fit_iv$iv_first_stage[[1L]]

    expect_false(.is_fixest_iv_second_stage(fit_ols))
    expect_true(.is_fixest_iv_second_stage(fit_iv))
    expect_true(.is_fixest_iv_second_stage(fit_fe_iv))
    expect_false(.is_fixest_iv_second_stage(first_stage))
  })
})
