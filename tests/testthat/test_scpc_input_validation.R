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

test_that("print.scpc handles single-coefficient models without dropping dimensions", {
  dat <- make_scpc_data(12L)
  fit <- stats::lm(y ~ x - 1, data = dat)
  out <- scpc(fit, dat, coords_euclidean = c("lon", "lat"),
              backend = "fallback", uncond = TRUE, cvs = FALSE)
  expect_identical(nrow(out$scpcstats), 1L)
  expect_no_error(printed <- utils::capture.output(print(out)))
  expect_true(any(grepl("\\bx\\b", printed)))
})

test_that("GLM rejects non-Gaussian families and non-identity links", {
  set.seed(11)
  d <- data.frame(
    y_bin = c(0, 1, 0, 1, 1, 0, 1, 0, 1, 0),
    y_cnt = c(1, 2, 0, 3, 2, 1, 4, 2, 3, 1),
    y_pos = seq(1.5, 6.0, length.out = 10L),
    x = seq_len(10L),
    s1 = seq(0, 1, length.out = 10L),
    s2 = sin(seq_len(10L))
  )
  fit_logit <- stats::glm(y_bin ~ x, data = d, family = stats::binomial())
  fit_pois <- stats::glm(y_cnt ~ x, data = d, family = stats::poisson())
  fit_loggauss <- stats::glm(y_pos ~ x, data = d, family = stats::gaussian(link = "log"))

  for (fit in list(fit_logit, fit_pois, fit_loggauss)) {
    expect_error(
      scpc(fit, d, coords_euclidean = c("s1", "s2"), backend = "fallback", uncond = TRUE),
      "Gaussian identity-link"
    )
  }
})

test_that("fixest varying-slope FEs match explicit lm interactions and pure-FE models error cleanly", {
  skip_if_not_installed("fixest")
  set.seed(101)
  n <- 40L
  d <- data.frame(
    fe1 = factor(rep(1:2, each = 20L)),
    fe2 = factor(rep(1:8, times = 5L)),
    z = stats::rnorm(n),
    z2 = stats::rnorm(n),
    x = stats::rnorm(n),
    s1 = stats::runif(n),
    s2 = stats::runif(n)
  )
  d$y <- 1 + 0.7 * d$x +
    as.numeric(d$fe1) * (0.4 + 0.5 * d$z - 0.3 * d$z2) +
    0.2 * as.numeric(d$fe2) * (1 + d$z) +
    stats::rnorm(n, sd = 0.3)

  fit_pure_fe <- fixest::feols(y ~ 1 | fe1, data = d)
  expect_error(
    scpc(fit_pure_fe, d, coords_euclidean = c("s1", "s2"), backend = "fallback", uncond = TRUE),
    "at least one estimated model coefficient"
  )

  spec_pairs <- list(
    list(
      fe = fixest::feols(y ~ x | fe1[z], data = d, fixef.rm = "none"),
      lm = stats::lm(y ~ x + fe1 + fe1:z, data = d),
      target = "x"
    ),
    list(
      fe = fixest::feols(y ~ x | fe1[[z]], data = d, fixef.rm = "none"),
      lm = stats::lm(y ~ x + fe1:z, data = d),
      target = c("(Intercept)", "x")
    ),
    list(
      fe = fixest::feols(y ~ x | fe1[z, z2] + fe2[z], data = d, fixef.rm = "none",
                       fixef.tol = 1e-11),
      lm = stats::lm(y ~ x + fe1 + fe1:z + fe1:z2 + fe2 + fe2:z, data = d),
      target = "x"
    )
  )

  for (pair in spec_pairs) {
    for (uncond_flag in c(FALSE, TRUE)) {
      res_fe <- scpc(pair$fe, d, coords_euclidean = c("s1", "s2"),
                     avc = 0.05, uncond = uncond_flag, cvs = FALSE)
      res_lm <- scpc(pair$lm, d, coords_euclidean = c("s1", "s2"),
                     avc = 0.05, ncoef = length(stats::coef(pair$lm)),
                     uncond = uncond_flag, cvs = FALSE)
      expect_equal(res_fe$c0, res_lm$c0, tolerance = 1e-10)
      expect_equal(res_fe$cv, res_lm$cv, tolerance = 1e-10)
      expect_equal(
        unname(res_fe$scpcstats[pair$target, , drop = FALSE]),
        unname(res_lm$scpcstats[pair$target, , drop = FALSE]),
        tolerance = 1e-10
      )
    }
  }
})
