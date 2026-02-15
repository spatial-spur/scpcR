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
    ncoef =2,
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
    ncoef =2,
    avc = 0.1,
    uncond = TRUE
  )
  out_1d <- scpc(
    model = fit,
    data = dat,
    coord_euclidean = "lon",
    ncoef =1,
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
    ncoef =1,
    avc = 0.1,
    uncond = TRUE,
    cvs = TRUE
  )
  out_euc <- scpc(
    model = fit,
    data = dat,
    coord_euclidean = c("lon", "lat"),
    ncoef =1,
    avc = 0.1,
    uncond = TRUE,
    cvs = TRUE
  )

  expect_equal(dim(out_geo$scpccvs), c(1, 4))
  expect_equal(dim(out_euc$scpccvs), c(1, 4))
  expect_true(all(is.finite(out_geo$scpccvs)))
  expect_true(all(is.finite(out_euc$scpccvs)))
})

test_that("scpc supports clustering with NA rows in model", {
  dat <- make_scpc_data(with_na = TRUE)
  fit <- stats::lm(y ~ x, data = dat)
  expect_warning(
    out <- scpc(
      model = fit,
      data = dat,
      coord_euclidean = c("lon", "lat"),
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

test_that("cluster must be a column name string", {
  dat <- make_scpc_data()
  fit <- stats::lm(y ~ x, data = dat)
  expect_error(
    scpc(fit, data = dat, coord_euclidean = c("lon", "lat"),
         cluster = dat$cluster, avc = 0.1, uncond = TRUE),
    "`cluster` must be a single column name"
  )
  expect_error(
    scpc(fit, data = dat, coord_euclidean = c("lon", "lat"),
         cluster = "nonexistent", avc = 0.1, uncond = TRUE),
    "Cluster variable not found in data"
  )
})

test_that("cvs option does not change scpcstats", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)

  out_no <- scpc(
    model = fit,
    data = dat,
    coord_euclidean = c("lon", "lat"),
    ncoef =2,
    avc = 0.1,
    uncond = FALSE,
    cvs = FALSE
  )
  out_yes <- scpc(
    model = fit,
    data = dat,
    coord_euclidean = c("lon", "lat"),
    ncoef =2,
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
    ncoef =2,
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
    ncoef =2,
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

test_that("scpc supports fixest IV models (unconditional and conditional)", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("geodist")
  old_threads <- fixest::getFixest_nthreads()
  fixest::setFixest_nthreads(1)
  on.exit(fixest::setFixest_nthreads(old_threads), add = TRUE)

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
    model = fit,
    data = dat,
    coord_euclidean = c("lon", "lat"),
    ncoef =2,
    avc = 0.1,
    uncond = TRUE,
    cvs = TRUE
  )
  out_c <- scpc(
    model = fit,
    data = dat,
    lon = "lon",
    lat = "lat",
    ncoef =2,
    avc = 0.1,
    uncond = FALSE,
    cvs = TRUE
  )

  expect_s3_class(out_u, "scpc")
  expect_s3_class(out_c, "scpc")
  expect_equal(dim(out_u$scpcstats), c(2, 6))
  expect_equal(dim(out_c$scpcstats), c(2, 6))
  expect_equal(dim(out_u$scpccvs), c(2, 4))
  expect_equal(dim(out_c$scpccvs), c(2, 4))
  expect_true(all(is.finite(out_u$scpcstats)))
  expect_true(all(is.finite(out_c$scpcstats)))
})

test_that("conditional SCPC for fixest FE matches lm with explicit FE dummies", {
  skip_if_not_installed("fixest")
  old_threads <- fixest::getFixest_nthreads()
  fixest::setFixest_nthreads(1)
  on.exit(fixest::setFixest_nthreads(old_threads), add = TRUE)

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
    model = fit_fixest,
    data = dat,
    coord_euclidean = c("lon", "lat"),
    ncoef =2,
    avc = 0.05,
    uncond = FALSE,
    cvs = TRUE
  )
  out_lm <- scpc(
    model = fit_lm,
    data = dat,
    coord_euclidean = c("lon", "lat"),
    ncoef =3,
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

test_that("conditional SCPC for fixest absorbed FE IV matches explicit FE dummies", {
  skip_if_not_installed("fixest")
  old_threads <- fixest::getFixest_nthreads()
  fixest::setFixest_nthreads(1)
  on.exit(fixest::setFixest_nthreads(old_threads), add = TRUE)

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
    model = fit_absorbed,
    data = dat,
    coord_euclidean = c("lon", "lat"),
    ncoef =2,
    avc = 0.1,
    uncond = FALSE,
    cvs = TRUE
  )
  out_explicit <- scpc(
    model = fit_explicit,
    data = dat,
    coord_euclidean = c("lon", "lat"),
    ncoef =30,
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

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

test_that("coef.scpc returns named vector of coefficients", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)
  out <- scpc(fit, data = dat, coord_euclidean = c("lon", "lat"),
              ncoef =2, avc = 0.1, uncond = TRUE)
  co <- coef(out)
  expect_true(is.numeric(co))
  expect_equal(length(co), 2)
  expect_equal(names(co), c("(Intercept)", "x"))
  expect_equal(co, stats::coef(fit), tolerance = 1e-10)
})

test_that("confint.scpc returns 95% CI by default", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)
  out <- scpc(fit, data = dat, coord_euclidean = c("lon", "lat"),
              ncoef =2, avc = 0.1, uncond = TRUE)
  ci <- confint(out)
  expect_equal(nrow(ci), 2)
  expect_equal(ncol(ci), 2)
  expect_equal(rownames(ci), c("(Intercept)", "x"))
  expect_equal(unname(ci[, 1]), unname(out$scpcstats[, "2.5 %"]))
  expect_equal(unname(ci[, 2]), unname(out$scpcstats[, "97.5 %"]))
})

test_that("confint.scpc supports parm subsetting", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)
  out <- scpc(fit, data = dat, coord_euclidean = c("lon", "lat"),
              ncoef =2, avc = 0.1, uncond = TRUE)
  ci <- confint(out, parm = "x")
  expect_equal(nrow(ci), 1)
  expect_equal(rownames(ci), "x")
})

test_that("confint.scpc errors for unsupported level without cvs", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)
  out <- scpc(fit, data = dat, coord_euclidean = c("lon", "lat"),
              ncoef =2, avc = 0.1, uncond = TRUE, cvs = FALSE)
  expect_error(confint(out, level = 0.90), "not available")
})

test_that("confint.scpc supports other levels when cvs = TRUE", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)
  out <- scpc(fit, data = dat, coord_euclidean = c("lon", "lat"),
              ncoef =2, avc = 0.1, uncond = TRUE, cvs = TRUE)
  ci90 <- confint(out, level = 0.90)
  ci95 <- confint(out, level = 0.95)
  # 90% CI should be narrower than 95% CI
  expect_true(all(ci90[, 2] - ci90[, 1] < ci95[, 2] - ci95[, 1]))
})

test_that("summary.scpc runs without error", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)
  out <- scpc(fit, data = dat, coord_euclidean = c("lon", "lat"),
              ncoef =2, avc = 0.1, uncond = TRUE)
  expect_output(summary(out), "SCPC Inference")
})

test_that("scpccvs has proper row and column names", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)
  out <- scpc(fit, data = dat, coord_euclidean = c("lon", "lat"),
              ncoef =2, avc = 0.1, uncond = TRUE, cvs = TRUE)
  expect_equal(rownames(out$scpccvs), c("(Intercept)", "x"))
  expect_equal(colnames(out$scpccvs), c("32%", "10%", "5%", "1%"))
})

test_that("ncoef = NULL reports all coefficients", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)
  out <- scpc(fit, data = dat, coord_euclidean = c("lon", "lat"),
              avc = 0.1, uncond = TRUE)
  expect_equal(nrow(out$scpcstats), 2)
})

test_that("k alias works for backwards compatibility", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)
  out <- scpc(fit, data = dat, coord_euclidean = c("lon", "lat"),
              k = 1, avc = 0.1, uncond = TRUE)
  expect_equal(nrow(out$scpcstats), 1)
})

test_that("using both ncoef and k errors", {
  dat <- make_scpc_data(n = 30)
  fit <- stats::lm(y ~ x, data = dat)
  expect_error(
    scpc(fit, data = dat, coord_euclidean = c("lon", "lat"),
         ncoef = 2, k = 1, avc = 0.1, uncond = TRUE),
    "Specify either `ncoef` or `k`, not both"
  )
})
