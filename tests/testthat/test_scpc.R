manual_fixest_iv_conditional_reference <- function(model,
                                                   data,
                                                   coord_cols,
                                                   avc,
                                                   coef_names,
                                                   method = "exact",
                                                   large_n_seed = 1) {
  design <- .get_fixest_iv_design(model)

  n <- nrow(data)
  coef_vec <- stats::coef(model)
  S <- sandwich::estfun(model)
  bread_inv <- sandwich::bread(model) / n
  model_mat <- design$model_mat
  if (!is.null(design$fixef_id)) {
    model_mat <- as.matrix(fixest::demean(model_mat, f = design$fixef_id, nthreads = 1L))
  }
  spc <- .setOmsWfin(
    as.matrix(data[, coord_cols, drop = FALSE]),
    avc0 = avc,
    latlong = FALSE,
    method = method,
    large_n_seed = large_n_seed
  )
  Wfin <- spc$Wfin
  Omsfin <- spc$Omsfin
  perm <- spc$perm
  q <- ncol(Wfin) - 1L
  levs <- c(0.32, 0.10, 0.05, 0.01)
  cvs_uncond <- vapply(levs, function(lv) .getcv(Omsfin, q, lv), 0.0)
  large_n_random_state <- spc$random_state

  stats_out <- matrix(
    NA_real_,
    nrow = length(coef_names),
    ncol = 6,
    dimnames = list(coef_names, c("Coef", "Std_Err", "t", "P>|t|", "2.5 %", "97.5 %"))
  )
  cvs_out <- matrix(
    NA_real_,
    nrow = length(coef_names),
    ncol = 4,
    dimnames = list(coef_names, c("32%", "10%", "5%", "1%"))
  )

  for (i in seq_along(coef_names)) {
    coef_name <- coef_names[[i]]
    pos <- match(coef_name, colnames(bread_inv))
    if (is.na(pos)) {
      stop("Missing coefficient in manual IV reference: ", coef_name, call. = FALSE)
    }

    wj <- as.numeric(n * bread_inv[pos, , drop = TRUE] %*% t(S)) + coef_vec[[coef_name]]
    wj_perm <- wj[perm]
    tau_u <- as.numeric(
      sqrt(q) * crossprod(Wfin[, 1], wj_perm) /
        sqrt(sum((t(Wfin[, -1]) %*% wj_perm)^2))
    )
    se <- as.numeric(
      sqrt(sum((t(Wfin[, -1]) %*% wj_perm)^2)) / (sqrt(q) * sqrt(n))
    )
    p_u <- .maxrp(Omsfin, q, abs(tau_u) / sqrt(q))$max

    xj_raw <- as.numeric(n * bread_inv[pos, , drop = TRUE] %*% t(model_mat))
    if (spc$large_n) {
      residualize <- .make_iv_residualizer(
        design$X[perm, , drop = FALSE],
        design$Z[perm, , drop = FALSE],
        fixef_id = .permute_fixef_id(design$fixef_id, perm)
      )
      xj <- xj_raw[perm]
    } else {
      residualize <- .make_iv_residualizer(
        design$X,
        design$Z,
        fixef_id = design$fixef_id
      )
      xj <- xj_raw
    }
    xjs <- sign(xj)
    Wx <- .orthogonalize_W_iv(Wfin, xj, xjs, residualize = residualize)

    if (spc$large_n) {
      omsx_res <- .lnget_Oms(
        spc$coords,
        spc$c0,
        spc$cmax,
        Wx,
        1.2,
        capM = 1000000L,
        random_t = large_n_random_state,
        latlong = FALSE
      )
      Omsx <- omsx_res$Oms
      large_n_random_state <- omsx_res$state
    } else {
      Omsx <- .getOms(spc$distmat, spc$c0, spc$cmax, Wx, 1.2)
    }
    p_c <- .maxrp(Omsx, q, abs(tau_u) / sqrt(q))$max
    cvs_cond <- vapply(levs, function(lv) .getcv(Omsx, q, lv), 0.0)
    cv <- max(spc$cvfin, cvs_cond[[3]])

    stats_out[i, ] <- c(
      coef_vec[[coef_name]],
      se,
      tau_u,
      max(p_u, p_c),
      coef_vec[[coef_name]] - cv * se,
      coef_vec[[coef_name]] + cv * se
    )
    cvs_out[i, ] <- pmax(cvs_uncond, cvs_cond)
  }

  list(stats = stats_out, cvs = cvs_out)
}

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
  dat$lon <- ave(dat$lon, dat$cluster, FUN = mean)
  dat$lat <- ave(dat$lat, dat$cluster, FUN = mean)
  fit <- stats::lm(y ~ x, data = dat)

  out <- scpc(
    fit,
    dat,
    coords_euclidean = c("lon", "lat"),
    cluster = "cluster",
    ncoef = 1,
    avc = 0.1,
    uncond = TRUE
  )

  expect_s3_class(out, "scpc")
  expect_equal(nrow(out$scpcstats), 1)
})

test_that("scpc errors when coordinates vary within clusters", {
  dat <- make_scpc_data()
  fit <- stats::lm(y ~ x, data = dat)

  expect_error(
    scpc(
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
})

test_that("scpc validates cluster input and coordinate mode selection", {
  dat <- make_scpc_data()
  fit <- stats::lm(y ~ x, data = dat)

  expect_error(
    scpc(fit, data = dat, coords_euclidean = c("lon", "lat"), cluster = dat$cluster, avc = 0.1, uncond = TRUE),
    "`cluster` must be a single non-empty column name"
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
    "`avc` must lie in \\(0.001, 0.99\\)"
  )
})

test_that("scpc validates method and large_n_seed", {
  dat <- make_scpc_data()
  fit <- stats::lm(y ~ x, data = dat)

  expect_error(
    scpc(
      fit,
      data = dat,
      coords_euclidean = c("lon", "lat"),
      avc = 0.1,
      method = "bad",
      uncond = TRUE
    ),
    "`method` must be one of \"auto\", \"exact\", or \"approx\""
  )
  expect_error(
    scpc(
      fit,
      data = dat,
      coords_euclidean = c("lon", "lat"),
      avc = 0.1,
      large_n_seed = 1.5,
      uncond = TRUE
    ),
    "`large_n_seed` must be a single integer-valued number in \\[0, 2\\^32\\)"
  )
})

test_that("scpc method override can force exact or approximation branches", {
  dat_small <- make_scpc_data(n = 40)
  fit_small <- stats::lm(y ~ x, data = dat_small)

  out_auto_small <- scpc(
    fit_small,
    data = dat_small,
    coords_euclidean = c("lon", "lat"),
    ncoef = 1,
    avc = 0.1,
    method = "auto",
    large_n_seed = 1,
    uncond = TRUE
  )
  out_exact_small <- scpc(
    fit_small,
    data = dat_small,
    coords_euclidean = c("lon", "lat"),
    ncoef = 1,
    avc = 0.1,
    method = "exact",
    large_n_seed = 1,
    uncond = TRUE
  )
  out_approx_small_a <- scpc(
    fit_small,
    data = dat_small,
    coords_euclidean = c("lon", "lat"),
    ncoef = 1,
    avc = 0.1,
    method = "approx",
    large_n_seed = 1,
    uncond = TRUE
  )
  out_approx_small_b <- scpc(
    fit_small,
    data = dat_small,
    coords_euclidean = c("lon", "lat"),
    ncoef = 1,
    avc = 0.1,
    method = "approx",
    large_n_seed = 1,
    uncond = TRUE
  )

  expect_identical(out_auto_small$method, "exact")
  expect_identical(out_exact_small$method, "exact")
  expect_identical(out_approx_small_a$method, "approx")
  expect_equal(out_auto_small$scpcstats, out_exact_small$scpcstats, tolerance = 1e-10)
  expect_equal(out_auto_small$c0, out_exact_small$c0, tolerance = 1e-10)
  expect_equal(out_auto_small$cv, out_exact_small$cv, tolerance = 1e-10)
  expect_equal(out_auto_small$q, out_exact_small$q)
  expect_equal(out_approx_small_a$scpcstats, out_approx_small_b$scpcstats, tolerance = 1e-12)
  expect_equal(out_approx_small_a$c0, out_approx_small_b$c0, tolerance = 1e-12)
  expect_equal(out_approx_small_a$cv, out_approx_small_b$cv, tolerance = 1e-12)
  expect_equal(out_approx_small_a$q, out_approx_small_b$q)
  expect_true(all(is.finite(out_approx_small_a$scpcstats)))
})

test_that("scpc routes auto to the large-n path in a full run", {
  skip_if(
    !identical(Sys.getenv("SCPCR_RUN_HEAVY_TESTS", "false"), "true"),
    "set SCPCR_RUN_HEAVY_TESTS=true to run the full large-n routing integration test"
  )

  dat_large <- make_scpc_data(n = 4500)
  fit_large <- stats::lm(y ~ x, data = dat_large)

  out_auto_large <- scpc(
    fit_large,
    data = dat_large,
    coords_euclidean = c("lon", "lat"),
    ncoef = 1,
    avc = 0.1,
    method = "auto",
    large_n_seed = 1,
    uncond = TRUE
  )

  expect_identical(out_auto_large$method, "approx")
  expect_true(all(is.finite(out_auto_large$scpcstats)))
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

test_that("conditional SCPC for fixest IV matches a direct 2SLS residualization reference", {
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    set.seed(2028)
    n_fe <- 30
    t_per_fe <- 6
    n <- n_fe * t_per_fe
    fe <- rep(seq_len(n_fe), each = t_per_fe)
    z <- stats::rnorm(n)
    w <- stats::rnorm(n)
    u <- stats::rnorm(n)
    x <- 0.8 * z + 0.4 * w + 0.6 * u + stats::rnorm(n, sd = 0.2)
    y <- 1 + 1.1 * x + 0.5 * w + stats::rnorm(n_fe)[fe] + u
    dat <- data.frame(
      y = y,
      x = x,
      z = z,
      w = w,
      fe = fe,
      coord_x = runif(n),
      coord_y = runif(n)
    )

    fit <- fixest::feols(y ~ w + i(fe) | x ~ z, data = dat)
    out <- scpc(
      fit,
      dat,
      coords_euclidean = c("coord_x", "coord_y"),
      ncoef = 3,
      avc = 0.05,
      uncond = FALSE,
      cvs = TRUE
    )
    ref <- manual_fixest_iv_conditional_reference(
      fit,
      dat,
      coord_cols = c("coord_x", "coord_y"),
      avc = 0.05,
      coef_names = c("fit_x", "w")
    )

    stats_out <- as.data.frame(out$scpcstats)
    stats_out$term <- rownames(out$scpcstats)
    stats_out <- stats_out[stats_out$term %in% rownames(ref$stats), , drop = FALSE]
    stats_mat <- as.matrix(stats_out[, c("Coef", "Std_Err", "t", "P>|t|", "2.5 %", "97.5 %")])
    rownames(stats_mat) <- stats_out$term

    expect_equal(stats_mat[rownames(ref$stats), , drop = FALSE], ref$stats, tolerance = 1e-8)
    expect_equal(out$scpccvs[rownames(ref$cvs), , drop = FALSE], ref$cvs, tolerance = 1e-8)
  })
})

test_that("conditional SCPC for fixest IV matches a direct large-n reference", {
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    set.seed(2031)
    n_fe <- 30
    t_per_fe <- 6
    n <- n_fe * t_per_fe
    fe <- rep(seq_len(n_fe), each = t_per_fe)
    z <- stats::rnorm(n)
    w <- stats::rnorm(n)
    u <- stats::rnorm(n)
    x <- 0.8 * z + 0.4 * w + 0.6 * u + stats::rnorm(n, sd = 0.2)
    y <- 1 + 1.1 * x + 0.5 * w + stats::rnorm(n_fe)[fe] + u
    dat <- data.frame(
      y = y,
      x = x,
      z = z,
      w = w,
      fe = fe,
      coord_x = runif(n),
      coord_y = runif(n)
    )

    fit <- fixest::feols(y ~ w | fe | x ~ z, data = dat)
    out <- scpc(
      fit,
      dat,
      coords_euclidean = c("coord_x", "coord_y"),
      ncoef = 2,
      avc = 0.05,
      method = "approx",
      large_n_seed = 1,
      uncond = FALSE,
      cvs = TRUE
    )
    ref <- manual_fixest_iv_conditional_reference(
      fit,
      dat,
      coord_cols = c("coord_x", "coord_y"),
      avc = 0.05,
      coef_names = c("fit_x", "w"),
      method = "approx",
      large_n_seed = 1
    )

    stats_out <- as.data.frame(out$scpcstats)
    stats_out$term <- rownames(out$scpcstats)
    stats_out <- stats_out[stats_out$term %in% rownames(ref$stats), , drop = FALSE]
    stats_mat <- as.matrix(stats_out[, c("Coef", "Std_Err", "t", "P>|t|", "2.5 %", "97.5 %")])
    rownames(stats_mat) <- stats_out$term

    expect_identical(out$method, "approx")
    expect_equal(out$large_n_seed, 1)
    expect_equal(stats_mat[rownames(ref$stats), , drop = FALSE], ref$stats, tolerance = 1e-8)
    expect_equal(out$scpccvs[rownames(ref$cvs), , drop = FALSE], ref$cvs, tolerance = 1e-8)
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
