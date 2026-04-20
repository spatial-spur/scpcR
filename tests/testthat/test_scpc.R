test_that("scpc returns a result with the expected structure", {
  dat <- make_python_scpc_data()
  fit <- stats::lm(y ~ x, data = dat)

  out <- scpc(
    fit,
    dat,
    coords_euclidean = c("coord_x", "coord_y"),
    avc = 0.1,
    uncond = TRUE,
    cvs = FALSE
  )

  expect_s3_class(out, "scpc")
  expect_equal(dim(out$scpcstats), c(2, 6))
  expect_null(out$scpccvs)
  expect_equal(out$avc, 0.1, tolerance = 1e-12)
  expect_true(out$q >= 1)
  expect_equal(nrow(out$W), 5)
})

test_that("scpc supports geodesic and euclidean modes", {
  skip_if_not_installed("geodist")
  dat <- make_scpc_data()
  fit <- stats::lm(y ~ x, data = dat)

  out_geo <- scpc(
    fit,
    dat,
    lon = "lon",
    lat = "lat",
    ncoef = 2,
    avc = 0.1,
    uncond = TRUE
  )
  out_euc <- scpc(
    fit,
    dat,
    coords_euclidean = c("lon", "lat"),
    ncoef = 2,
    avc = 0.1,
    uncond = TRUE
  )

  expect_equal(dim(out_geo$scpcstats), c(2, 6))
  expect_equal(dim(out_euc$scpcstats), c(2, 6))
  expect_true(all(is.finite(out_geo$scpcstats)))
  expect_true(all(is.finite(out_euc$scpcstats)))
})

test_that("scpc computes critical values when requested", {
  skip_if_not_installed("geodist")
  dat <- make_scpc_data()
  fit <- stats::lm(y ~ x, data = dat)

  out_geo <- scpc(
    fit,
    dat,
    lon = "lon",
    lat = "lat",
    ncoef = 1,
    avc = 0.1,
    uncond = TRUE,
    cvs = TRUE
  )
  out_euc <- scpc(
    fit,
    dat,
    coords_euclidean = c("lon", "lat"),
    ncoef = 1,
    avc = 0.1,
    uncond = TRUE,
    cvs = TRUE
  )

  expect_equal(dim(out_geo$scpccvs), c(1, 4))
  expect_equal(dim(out_euc$scpccvs), c(1, 4))
  expect_true(all(is.finite(out_geo$scpccvs)))
  expect_true(all(is.finite(out_euc$scpccvs)))
})

test_that("scpc supports clustering with NA rows in the model", {
  dat <- make_scpc_data(with_na = TRUE)
  fit <- stats::lm(y ~ x, data = dat)

  expect_warning(
    out <- scpc(
      fit,
      dat,
      coords_euclidean = c("lon", "lat"),
      cluster = "cluster",
      ncoef = 1,
      avc = 0.1,
      uncond = TRUE
    ),
    "Coordinates vary within clusters"
  )

  expect_s3_class(out, "scpc")
  expect_equal(nrow(out$scpcstats), 1)
})

test_that("scpc validates cluster input and coordinate mode selection", {
  dat <- make_scpc_data()
  fit <- stats::lm(y ~ x, data = dat)

  expect_error(
    scpc(fit, data = dat, coords_euclidean = c("lon", "lat"), cluster = dat$cluster, avc = 0.1, uncond = TRUE),
    "`cluster` must be a single column name"
  )
  expect_error(
    scpc(fit, data = dat, coords_euclidean = c("lon", "lat"), cluster = "nonexistent", avc = 0.1, uncond = TRUE),
    "Cluster variable not found in data"
  )
  expect_error(
    scpc(fit, data = dat, avc = 0.1, uncond = TRUE),
    "Specify coordinates via `lon`/`lat` or `coords_euclidean`"
  )
  expect_error(
    scpc(fit, data = dat, lon = "lon", lat = "lat", coords_euclidean = c("lon", "lat"), avc = 0.1, uncond = TRUE),
    "Specify either `lon`/`lat` or `coords_euclidean`, not both"
  )
  expect_error(
    scpc(fit, data = dat, lon = "lon", avc = 0.1, uncond = TRUE),
    "For geodesic coordinates, provide both `lon` and `lat`"
  )
  expect_error(
    scpc(fit, data = dat, coords_euclidean = "lon", avc = 1, uncond = TRUE),
    "Option avc\\(\\) must be in \\(0.001, 0.99\\)"
  )
})

test_that("scpc does not accept the old k argument", {
  dat <- make_python_scpc_data()
  fit <- stats::lm(y ~ x, data = dat)

  expect_error(
    scpc(fit, dat, coords_euclidean = c("coord_x", "coord_y"), k = 1),
    "unused argument"
  )
})

test_that("cvs does not change scpcstats", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)

  out_no <- scpc(
    fit,
    dat,
    coords_euclidean = c("lon", "lat"),
    ncoef = 2,
    avc = 0.1,
    uncond = FALSE,
    cvs = FALSE
  )
  out_yes <- scpc(
    fit,
    dat,
    coords_euclidean = c("lon", "lat"),
    ncoef = 2,
    avc = 0.1,
    uncond = FALSE,
    cvs = TRUE
  )

  expect_equal(out_no$scpcstats, out_yes$scpcstats, tolerance = 1e-10)
})

test_that("geodesic mode works with custom lon and lat names", {
  skip_if_not_installed("geodist")
  dat <- make_scpc_data(n = 24)

  out_ref <- scpc(
    stats::lm(y ~ x, data = dat),
    dat,
    lon = "lon",
    lat = "lat",
    ncoef = 2,
    avc = 0.1,
    uncond = TRUE,
    cvs = TRUE
  )

  dat_named <- dat
  names(dat_named)[names(dat_named) == "lon"] <- "my_longitude"
  names(dat_named)[names(dat_named) == "lat"] <- "my_latitude"
  out_named <- scpc(
    stats::lm(y ~ x, data = dat_named),
    dat_named,
    lon = "my_longitude",
    lat = "my_latitude",
    ncoef = 2,
    avc = 0.1,
    uncond = TRUE,
    cvs = TRUE
  )

  expect_equal(out_ref$scpcstats, out_named$scpcstats, tolerance = 1e-8)
  expect_equal(out_ref$scpccvs, out_named$scpccvs, tolerance = 1e-8)
})

test_that("scpc supports fixest IV models in unconditional and conditional modes", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("geodist")

  with_fixest_single_thread({
    set.seed(102)
    n <- 120
    z <- stats::rnorm(n)
    w <- stats::rnorm(n)
    u <- stats::rnorm(n)
    x <- 0.9 * z + 0.4 * w + 0.7 * u + stats::rnorm(n, sd = 0.2)
    y <- 1 + 1.2 * x + 0.5 * w + u
    dat <- data.frame(
      y = y,
      x = x,
      w = w,
      z = z,
      lon = runif(n, -125, -66),
      lat = runif(n, 25, 49)
    )

    fit <- fixest::feols(y ~ w | x ~ z, data = dat)

    out_u <- scpc(
      fit,
      dat,
      coords_euclidean = c("lon", "lat"),
      ncoef = 2,
      avc = 0.1,
      uncond = TRUE,
      cvs = TRUE
    )
    out_c <- scpc(
      fit,
      dat,
      lon = "lon",
      lat = "lat",
      ncoef = 2,
      avc = 0.1,
      uncond = FALSE,
      cvs = TRUE
    )

    expect_equal(dim(out_u$scpcstats), c(2, 6))
    expect_equal(dim(out_c$scpcstats), c(2, 6))
    expect_equal(dim(out_u$scpccvs), c(2, 4))
    expect_equal(dim(out_c$scpccvs), c(2, 4))
    expect_true(all(is.finite(out_u$scpcstats)))
    expect_true(all(is.finite(out_c$scpcstats)))
  })
})

test_that("conditional SCPC for fixest FE matches lm with explicit FE dummies", {
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    set.seed(2026)
    n_fe <- 40
    t_per_fe <- 6
    n <- n_fe * t_per_fe
    fe <- rep(seq_len(n_fe), each = t_per_fe)
    x1 <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    alpha <- stats::rnorm(n_fe)[fe]
    y <- 0.8 * x1 - 0.4 * x2 + alpha + stats::rnorm(n, sd = 0.7)
    dat <- data.frame(
      y = y,
      x1 = x1,
      x2 = x2,
      fe = fe,
      lon = runif(n, -125, -66),
      lat = runif(n, 25, 49)
    )

    fit_fixest <- fixest::feols(y ~ x1 + x2 | fe, data = dat)
    fit_lm <- stats::lm(y ~ x1 + x2 + factor(fe), data = dat)

    out_fixest <- scpc(
      fit_fixest,
      dat,
      coords_euclidean = c("lon", "lat"),
      ncoef = 2,
      avc = 0.05,
      uncond = FALSE,
      cvs = TRUE
    )
    out_lm <- scpc(
      fit_lm,
      dat,
      coords_euclidean = c("lon", "lat"),
      ncoef = 3,
      avc = 0.05,
      uncond = FALSE,
      cvs = TRUE
    )

    stats_fixest <- as.data.frame(out_fixest$scpcstats)
    stats_fixest$term <- rownames(out_fixest$scpcstats)
    stats_lm <- as.data.frame(out_lm$scpcstats)
    stats_lm$term <- rownames(out_lm$scpcstats)

    merged_stats <- merge(stats_fixest, stats_lm, by = "term", suffixes = c("_fixest", "_lm"), sort = FALSE)
    expect_equal(sort(merged_stats$term), c("x1", "x2"))
    for (v in c("Coef", "Std_Err", "t", "P>|t|", "2.5 %", "97.5 %")) {
      d <- max(abs(merged_stats[[paste0(v, "_fixest")]] - merged_stats[[paste0(v, "_lm")]]))
      expect_true(d < 2e-4, info = paste("max abs diff for", v, "=", format(d, scientific = TRUE)))
    }

    cvs_fixest <- as.data.frame(out_fixest$scpccvs)
    cvs_fixest$term <- rownames(out_fixest$scpcstats)
    cvs_lm <- as.data.frame(out_lm$scpccvs)
    cvs_lm$term <- rownames(out_lm$scpcstats)[seq_len(nrow(cvs_lm))]

    merged_cvs <- merge(cvs_fixest, cvs_lm, by = "term", suffixes = c("_fixest", "_lm"), sort = FALSE)
    expect_equal(sort(merged_cvs$term), c("x1", "x2"))
    for (v in c("32%", "10%", "5%", "1%")) {
      d <- max(abs(merged_cvs[[paste0(v, "_fixest")]] - merged_cvs[[paste0(v, "_lm")]]))
      expect_true(d < 2e-4, info = paste("max abs diff for cvs", v, "=", format(d, scientific = TRUE)))
    }
  })
})

test_that("conditional SCPC for fixest absorbed FE IV matches explicit FE dummies", {
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    set.seed(2027)
    n <- 200
    fe <- rep(seq_len(40), each = 5)
    z <- stats::rnorm(n)
    u <- stats::rnorm(n)
    w <- stats::rnorm(n)
    x <- 0.8 * z + 0.6 * u + stats::rnorm(n, sd = 0.2)
    y <- 1 + 1.2 * x + 0.4 * w + stats::rnorm(n) + stats::rnorm(40)[fe]
    dat <- data.frame(
      y = y,
      x = x,
      z = z,
      w = w,
      fe = fe,
      lon = runif(n, -125, -66),
      lat = runif(n, 25, 49)
    )

    fit_absorbed <- fixest::feols(y ~ w | fe | x ~ z, data = dat)
    fit_explicit <- fixest::feols(y ~ w + i(fe) | x ~ z, data = dat)

    out_absorbed <- scpc(
      fit_absorbed,
      dat,
      coords_euclidean = c("lon", "lat"),
      ncoef = 2,
      avc = 0.1,
      uncond = FALSE,
      cvs = TRUE
    )
    out_explicit <- scpc(
      fit_explicit,
      dat,
      coords_euclidean = c("lon", "lat"),
      ncoef = 30,
      avc = 0.1,
      uncond = FALSE,
      cvs = TRUE
    )

    st_abs <- as.data.frame(out_absorbed$scpcstats)
    st_abs$term <- rownames(out_absorbed$scpcstats)
    st_exp <- as.data.frame(out_explicit$scpcstats)
    st_exp$term <- rownames(out_explicit$scpcstats)
    merged_stats <- merge(st_abs, st_exp, by = "term", suffixes = c("_abs", "_exp"), sort = FALSE)
    expect_equal(sort(merged_stats$term), sort(c("fit_x", "w")))
    for (v in c("Coef", "Std_Err", "t", "P>|t|", "2.5 %", "97.5 %")) {
      d <- max(abs(merged_stats[[paste0(v, "_abs")]] - merged_stats[[paste0(v, "_exp")]]))
      expect_true(d < 3e-4, info = paste("max abs diff for", v, "=", format(d, scientific = TRUE)))
    }

    cv_abs <- as.data.frame(out_absorbed$scpccvs)
    cv_abs$term <- rownames(out_absorbed$scpcstats)
    cv_exp <- as.data.frame(out_explicit$scpccvs)
    cv_exp$term <- rownames(out_explicit$scpcstats)[seq_len(nrow(cv_exp))]
    merged_cvs <- merge(cv_abs, cv_exp, by = "term", suffixes = c("_abs", "_exp"), sort = FALSE)
    expect_equal(sort(merged_cvs$term), sort(c("fit_x", "w")))
    for (v in c("32%", "10%", "5%", "1%")) {
      d <- max(abs(merged_cvs[[paste0(v, "_abs")]] - merged_cvs[[paste0(v, "_exp")]]))
      expect_true(d < 3e-4, info = paste("max abs diff for cvs", v, "=", format(d, scientific = TRUE)))
    }
  })
})

test_that("coef.scpc returns a named vector of coefficients", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)
  out <- scpc(fit, data = dat, coords_euclidean = c("lon", "lat"), ncoef = 2, avc = 0.1, uncond = TRUE)
  co <- coef(out)

  expect_true(is.numeric(co))
  expect_equal(length(co), 2)
  expect_equal(names(co), c("(Intercept)", "x"))
  expect_equal(co, stats::coef(fit), tolerance = 1e-10)
})

test_that("confint.scpc returns the default 95 percent interval and supports subsetting", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)
  out <- scpc(fit, data = dat, coords_euclidean = c("lon", "lat"), ncoef = 2, avc = 0.1, uncond = TRUE)

  ci <- confint(out)
  ci_x <- confint(out, parm = "x")

  expect_equal(nrow(ci), 2)
  expect_equal(ncol(ci), 2)
  expect_equal(rownames(ci), c("(Intercept)", "x"))
  expect_equal(unname(ci[, 1]), unname(out$scpcstats[, "2.5 %"]))
  expect_equal(unname(ci[, 2]), unname(out$scpcstats[, "97.5 %"]))
  expect_equal(nrow(ci_x), 1)
  expect_equal(rownames(ci_x), "x")
})

test_that("confint.scpc supports additional levels only when cvs are stored", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)
  out_no <- scpc(fit, data = dat, coords_euclidean = c("lon", "lat"), ncoef = 2, avc = 0.1, uncond = TRUE, cvs = FALSE)
  out_yes <- scpc(fit, data = dat, coords_euclidean = c("lon", "lat"), ncoef = 2, avc = 0.1, uncond = TRUE, cvs = TRUE)

  expect_error(confint(out_no, level = 0.90), "not available")

  ci90 <- confint(out_yes, level = 0.90)
  ci95 <- confint(out_yes, level = 0.95)
  expect_true(all(ci90[, 2] - ci90[, 1] < ci95[, 2] - ci95[, 1]))
})

test_that("print.scpc, summary.scpc, and stored cvs all behave as expected", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)
  out <- scpc(fit, data = dat, coords_euclidean = c("lon", "lat"), ncoef = 2, avc = 0.1, uncond = TRUE, cvs = TRUE)

  expect_output(print(out), "SCPC Inference")
  expect_output(summary(out), "SCPC Inference")
  expect_equal(rownames(out$scpccvs), c("(Intercept)", "x"))
  expect_equal(colnames(out$scpccvs), c("32%", "10%", "5%", "1%"))
})

test_that("ncoef = NULL reports all coefficients", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)
  out <- scpc(fit, data = dat, coords_euclidean = c("lon", "lat"), avc = 0.1, uncond = TRUE)

  expect_equal(nrow(out$scpcstats), 2)
})
