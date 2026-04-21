# Parity against Stata should only tolerate floating-point-scale noise.
# Larger gaps should fail loudly and be treated as implementation drift.
STATA_STATS_TOL <- c(
  coef = 1e-4,
  std_err = 1e-4,
  t = 1e-4,
  p = 1e-4,
  ci_low = 1e-3,
  ci_high = 1e-3
)

STATA_CVS_TOL <- c(
  cv32 = 1e-4,
  cv10 = 1e-4,
  cv5 = 1e-4,
  cv1 = 1e-4
)

write_stata_parity_csv <- function(data) {
  path <- tempfile("scpc_parity_data_", fileext = ".csv")
  utils::write.csv(data, path, row.names = FALSE)
  path
}

build_custom_stata_scpc_lines <- function(data_path,
                                          reg_line,
                                          scpc_opts,
                                          extra_lines = character(),
                                          coord_cols = c("lat", "lon"),
                                          latlong = TRUE) {
  coord_x <- coord_cols[[1]]
  coord_y <- coord_cols[[2]]
  scpc_line <- if (length(scpc_opts) > 0L) {
    sprintf("scpc, %s", paste(scpc_opts, collapse = " "))
  } else {
    "scpc"
  }

  c(
    "clear all",
    "set more off",
    sprintf("import delimited \"%s\", varnames(1) clear", data_path),
    "rename *, lower",
    extra_lines,
    sprintf("gen s_1 = %s", coord_x),
    sprintf("gen s_2 = %s", coord_y),
    reg_line,
    scpc_line
  )
}

normalize_iv_r_terms <- function(tab) {
  if (is.null(tab)) {
    return(NULL)
  }
  tab$term <- sub("^fit_", "", tab$term)
  tab
}

expect_custom_stata_scpc_match <- function(stata_lines,
                                           out_r,
                                           terms_target,
                                           label,
                                           stata_bin,
                                           stats_tol = STATA_STATS_TOL,
                                           cvs_tol = STATA_CVS_TOL) {
  stata_out <- run_stata_scpc_lines(
    stata_lines,
    write_cvs = !is.null(out_r$scpccvs),
    stata_bin = stata_bin
  )

  r_stats <- normalize_iv_r_terms(normalize_r_scpcstats(out_r))
  r_cvs <- normalize_iv_r_terms(normalize_r_scpccvs(out_r))

  stata_stats <- stata_out$stats[stata_out$stats$term %in% terms_target, , drop = FALSE]
  r_stats <- r_stats[r_stats$term %in% terms_target, , drop = FALSE]
  expect_equal(sort(stata_stats$term), sort(terms_target))
  expect_equal(sort(r_stats$term), sort(terms_target))

  expect_table_close(
    stata_stats,
    r_stats,
    vars = c("coef", "std_err", "t", "p", "ci_low", "ci_high"),
    tolerance = stats_tol,
    label = paste0(label, " stats")
  )

  if (!is.null(r_cvs)) {
    stata_cvs <- stata_out$cvs[stata_out$cvs$term %in% terms_target, , drop = FALSE]
    r_cvs <- r_cvs[r_cvs$term %in% terms_target, , drop = FALSE]
    expect_equal(sort(stata_cvs$term), sort(terms_target))
    expect_equal(sort(r_cvs$term), sort(terms_target))

    expect_table_close(
      stata_cvs,
      r_cvs,
      vars = c("cv32", "cv10", "cv5", "cv1"),
      tolerance = cvs_tol,
      label = paste0(label, " cvs")
    )
  }

  invisible(list(stata = stata_out, r = list(stats = r_stats, cvs = r_cvs)))
}

make_iv_parity_data <- function(seed,
                                n = NULL,
                                n_cluster = NULL,
                                n_fe = NULL,
                                t_per_fe = NULL,
                                n_fe2 = NULL,
                                reps = NULL,
                                include_w = TRUE) {
  set.seed(seed)

  if (!is.null(n_fe) && !is.null(n_fe2)) {
    if (is.null(reps)) {
      stop("`reps` must be supplied for two-way FE parity data.", call. = FALSE)
    }
    n <- n_fe * n_fe2 * reps
    fe1 <- rep(rep(seq_len(n_fe), each = n_fe2), each = reps)
    fe2 <- rep(rep(seq_len(n_fe2), times = n_fe), each = reps)
    fe_effect <- stats::rnorm(n_fe, sd = 0.6)[fe1] + stats::rnorm(n_fe2, sd = 0.5)[fe2]
  } else if (!is.null(n_fe)) {
    if (is.null(t_per_fe)) {
      stop("`t_per_fe` must be supplied for one-way FE parity data.", call. = FALSE)
    }
    n <- n_fe * t_per_fe
    fe <- rep(seq_len(n_fe), each = t_per_fe)
    fe_effect <- stats::rnorm(n_fe, sd = 0.7)[fe]
  } else {
    if (is.null(n)) {
      n <- 180L
    }
    fe_effect <- rep(0, n)
  }

  idx <- seq_len(n)
  z_base <- stats::rnorm(n)
  w <- if (isTRUE(include_w)) stats::rnorm(n) else NULL
  u_det <- sin(idx / 7) + 0.5 * cos(idx / 11)
  v_det <- cos(idx / 5) - 0.25 * sin(idx / 13)
  z_iv <- 0.9 * z_base + 0.15 * v_det
  x_endog <- 0.8 * z_iv + 0.25 * v_det + 0.2 * u_det
  if (isTRUE(include_w)) {
    x_endog <- x_endog + 0.35 * w
  }
  x_endog <- x_endog + stats::rnorm(n, sd = 0.25)

  y <- 1.1 * x_endog + fe_effect + 0.45 * u_det + stats::rnorm(n, sd = 0.35)
  if (isTRUE(include_w)) {
    y <- y + 0.55 * w
  } else {
    y <- y + 0.8
  }

  if (!is.null(n_cluster)) {
    clust_id <- rep(seq_len(n_cluster), length.out = n)
    lon_cl <- runif(n_cluster, -125, -66)
    lat_cl <- runif(n_cluster, 25, 49)
    lon <- lon_cl[clust_id]
    lat <- lat_cl[clust_id]
  } else {
    clust_id <- NULL
    lon <- runif(n, -125, -66)
    lat <- runif(n, 25, 49)
  }

  out <- data.frame(
    y = y,
    x_endog = x_endog,
    z_iv = z_iv,
    lat = lat,
    lon = lon
  )

  if (isTRUE(include_w)) {
    out$w <- w
  }
  if (!is.null(n_cluster)) {
    out$clust_id <- clust_id
  }
  if (exists("fe", inherits = FALSE)) {
    out$fe <- fe
  }
  if (exists("fe1", inherits = FALSE)) {
    out$fe1 <- fe1
    out$fe2 <- fe2
  }

  out
}

STANDARD_STATA_SCENARIOS <- list(
  list(
    id = "default_cond_latlong",
    dep = "am",
    rhs = c("tlfpr"),
    latlong = TRUE,
    uncond = FALSE,
    avc = NULL,
    k = NULL,
    cvs = FALSE,
    cluster_mode = "none"
  ),
  list(
    id = "uncond_latlong",
    dep = "am",
    rhs = c("tlfpr"),
    latlong = TRUE,
    uncond = TRUE,
    avc = 0.03,
    k = 2L,
    cvs = FALSE,
    cluster_mode = "none"
  ),
  list(
    id = "multivar_cond_cvs",
    dep = "am",
    rhs = c("tlfpr", "fracblack", "gini"),
    latlong = TRUE,
    uncond = FALSE,
    avc = 0.03,
    k = 3L,
    cvs = TRUE,
    cluster_mode = "none"
  ),
  list(
    id = "fe_cond_latlong",
    dep = "am",
    rhs = c("tlfpr", "fracblack", "gini"),
    include_fe = TRUE,
    latlong = TRUE,
    uncond = FALSE,
    avc = 0.03,
    k = 3L,
    cvs = TRUE,
    cluster_mode = "none"
  ),
  list(
    id = "multivar_uncond_avc005",
    dep = "am",
    rhs = c("tlfpr", "fracblack", "gini", "colpc"),
    latlong = TRUE,
    uncond = TRUE,
    avc = 0.05,
    k = 4L,
    cvs = FALSE,
    cluster_mode = "none"
  ),
  list(
    id = "euclidean_cond_avc001",
    dep = "am",
    rhs = c("tlfpr", "fracdiv", "fracmar"),
    latlong = FALSE,
    uncond = FALSE,
    avc = 0.01,
    k = 3L,
    cvs = FALSE,
    cluster_mode = "none"
  ),
  list(
    id = "euclidean_uncond_klarge",
    dep = "am",
    rhs = c("tlfpr", "fracdiv", "fracmar", "fracblack"),
    latlong = FALSE,
    uncond = TRUE,
    avc = 0.01,
    k = 50L,
    cvs = FALSE,
    cluster_mode = "none"
  ),
  list(
    id = "k1_uncond",
    dep = "am",
    rhs = c("tlfpr", "fracblack"),
    latlong = TRUE,
    uncond = TRUE,
    avc = 0.03,
    k = 1L,
    cvs = FALSE,
    cluster_mode = "none"
  ),
  list(
    id = "avc_minish_cond",
    dep = "am",
    rhs = c("tlfpr", "fracblack", "gini"),
    latlong = TRUE,
    uncond = FALSE,
    avc = 0.005,
    k = 3L,
    cvs = FALSE,
    cluster_mode = "none"
  )
)

CLUSTER_STATA_SCENARIOS <- list(
  list(
    id = "cluster_pairmean_cvs",
    dep = "am",
    rhs = c("tlfpr", "fracblack", "gini"),
    latlong = TRUE,
    uncond = TRUE,
    avc = 0.03,
    k = 3L,
    cvs = TRUE,
    cluster_mode = "pair_mean"
  ),
  list(
    id = "cluster_pairmean_cond",
    dep = "am",
    rhs = c("tlfpr", "fracblack", "gini"),
    latlong = TRUE,
    uncond = FALSE,
    avc = 0.03,
    k = 3L,
    cvs = TRUE,
    cluster_mode = "pair_mean"
  )
)

test_that("scpcR matches Stata scpc across the standard scenario grid", {
  ctx <- skip_if_stata_scpc_unavailable()

  for (scenario in STANDARD_STATA_SCENARIOS) {
    expect_stata_scpc_match(
      scenario,
      stats_tol = STATA_STATS_TOL,
      cvs_tol = if (isTRUE(scenario$cvs)) STATA_CVS_TOL else NULL,
      stata_bin = ctx$stata_bin,
      data_path = ctx$data_path
    )
  }
})

test_that("scpcR matches Stata scpc for clustered pair-mean scenarios", {
  ctx <- skip_if_stata_scpc_unavailable()

  for (scenario in CLUSTER_STATA_SCENARIOS) {
    expect_stata_scpc_match(
      scenario,
      stats_tol = STATA_STATS_TOL,
      cvs_tol = STATA_CVS_TOL,
      stata_bin = ctx$stata_bin,
      data_path = ctx$data_path
    )
  }
})

test_that("absorbed FE IV conditional SCPC matches Stata within tight tolerances", {
  ctx <- skip_if_stata_scpc_unavailable()
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    stata_out <- run_stata_scpc_lines(
      c(
        "clear all",
        "set more off",
        sprintf("import delimited \"%s\", varnames(1) clear", ctx$data_path),
        "rename *, lower",
        "drop if missing(am) | missing(fracblack) | missing(tlfpr) | missing(gini) | missing(lat) | missing(lon)",
        "gen long fe = mod(_n - 1, 20) + 1",
        "gen double u_det = sin(_n / 7)",
        "gen double v_det = cos(_n / 11)",
        "gen double z_iv = tlfpr + 0.25*fracblack + 0.1*v_det",
        "gen double x_endog = 0.6*z_iv + 0.1*gini + 0.2*u_det",
        "gen s_1 = lat",
        "gen s_2 = lon",
        "ivregress 2sls am fracblack i.fe (x_endog = z_iv), robust",
        "scpc, latlong avc(0.03) k(2) cvs"
      ),
      write_cvs = TRUE,
      stata_bin = ctx$stata_bin
    )

    d <- read.csv(ctx$data_path, stringsAsFactors = FALSE, check.names = FALSE)
    names(d) <- tolower(names(d))
    keep <- stats::complete.cases(d[, c("am", "fracblack", "tlfpr", "gini", "lat", "lon"), drop = FALSE])
    d <- d[keep, , drop = FALSE]
    rownames(d) <- seq_len(nrow(d))
    d$fe <- (seq_len(nrow(d)) - 1L) %% 20L + 1L
    d$u_det <- sin(seq_len(nrow(d)) / 7)
    d$v_det <- cos(seq_len(nrow(d)) / 11)
    d$z_iv <- d$tlfpr + 0.25 * d$fracblack + 0.1 * d$v_det
    d$x_endog <- 0.6 * d$z_iv + 0.1 * d$gini + 0.2 * d$u_det

    fit_r <- fixest::feols(am ~ fracblack | fe | x_endog ~ z_iv, data = d)
    out_r <- scpc(
      model = fit_r,
      data = d,
      lon = "lon",
      lat = "lat",
      avc = 0.03,
      ncoef = 2,
      uncond = FALSE,
      cvs = TRUE
    )

    r_stats <- normalize_r_scpcstats(out_r)
    r_stats$term <- sub("^fit_", "", r_stats$term)
    r_cvs <- normalize_r_scpccvs(out_r)
    r_cvs$term <- sub("^fit_", "", r_cvs$term)

    terms_target <- c("x_endog", "fracblack")
    stata_stats <- stata_out$stats[stata_out$stats$term %in% terms_target, , drop = FALSE]
    r_stats <- r_stats[r_stats$term %in% terms_target, , drop = FALSE]
    expect_equal(sort(stata_stats$term), sort(terms_target))
    expect_equal(sort(r_stats$term), sort(terms_target))

    stata_cvs <- stata_out$cvs[stata_out$cvs$term %in% terms_target, , drop = FALSE]
    r_cvs <- r_cvs[r_cvs$term %in% terms_target, , drop = FALSE]
    expect_equal(sort(stata_cvs$term), sort(terms_target))
    expect_equal(sort(r_cvs$term), sort(terms_target))

    expect_table_close(
      stata_stats,
      r_stats,
      vars = c("coef", "std_err", "t", "p", "ci_low", "ci_high"),
      tolerance = STATA_STATS_TOL,
      label = "iv_fe stats"
    )
    expect_table_close(
      stata_cvs,
      r_cvs,
      vars = c("cv32", "cv10", "cv5", "cv1"),
      tolerance = STATA_CVS_TOL,
      label = "iv_fe cvs"
    )
  })
})

test_that("additional IV parity cases match Stata within tight tolerances", {
  ctx <- skip_if_stata_scpc_unavailable()
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    cases <- list(
      list(
        id = "iv_no_fe_conditional_cvs",
        data = make_iv_parity_data(seed = 3101, n = 180L, include_w = TRUE),
        fit_fun = function(d) fixest::feols(y ~ w | x_endog ~ z_iv, data = d),
        reg_line = "ivregress 2sls y w (x_endog = z_iv), robust",
        terms = c("(Intercept)", "x_endog", "w"),
        avc = 0.03,
        uncond = FALSE,
        cluster = NULL
      ),
      list(
        id = "iv_intercept_only_unconditional_cvs",
        data = make_iv_parity_data(seed = 3102, n = 160L, include_w = FALSE),
        fit_fun = function(d) fixest::feols(y ~ 1 | x_endog ~ z_iv, data = d),
        reg_line = "ivregress 2sls y (x_endog = z_iv), robust",
        terms = c("(Intercept)", "x_endog"),
        avc = 0.05,
        uncond = TRUE,
        cluster = NULL
      ),
      list(
        id = "iv_cluster_conditional_cvs",
        data = make_iv_parity_data(seed = 3103, n = 180L, n_cluster = 60L, include_w = TRUE),
        fit_fun = function(d) fixest::feols(y ~ w | x_endog ~ z_iv, data = d),
        reg_line = "ivregress 2sls y w (x_endog = z_iv), cluster(clust_id)",
        terms = c("(Intercept)", "x_endog", "w"),
        avc = 0.03,
        uncond = FALSE,
        cluster = "clust_id"
      ),
      list(
        id = "iv_oneway_fe_unconditional_cvs",
        data = make_iv_parity_data(seed = 3104, n_fe = 18L, t_per_fe = 6L, include_w = TRUE),
        fit_fun = function(d) fixest::feols(y ~ w | fe | x_endog ~ z_iv, data = d),
        reg_line = "ivregress 2sls y w i.fe (x_endog = z_iv), robust",
        terms = c("x_endog", "w"),
        avc = 0.03,
        uncond = TRUE,
        cluster = NULL
      ),
      list(
        id = "iv_oneway_fe_cluster_conditional_cvs",
        data = make_iv_parity_data(
          seed = 3205,
          n_fe = 24L,
          t_per_fe = 8L,
          n_cluster = 96L,
          include_w = TRUE
        ),
        fit_fun = function(d) fixest::feols(y ~ w | fe | x_endog ~ z_iv, data = d),
        reg_line = "ivregress 2sls y w i.fe (x_endog = z_iv), cluster(clust_id)",
        terms = c("x_endog", "w"),
        avc = 0.03,
        uncond = FALSE,
        cluster = "clust_id"
      ),
      list(
        id = "iv_twoway_fe_conditional_cvs",
        data = make_iv_parity_data(
          seed = 3106,
          n_fe = 5L,
          n_fe2 = 4L,
          reps = 5L,
          include_w = TRUE
        ),
        fit_fun = function(d) fixest::feols(y ~ w | fe1 + fe2 | x_endog ~ z_iv, data = d),
        reg_line = "ivregress 2sls y w i.fe1 i.fe2 (x_endog = z_iv), robust",
        terms = c("x_endog", "w"),
        avc = 0.03,
        uncond = FALSE,
        cluster = NULL
      )
    )

    for (case in cases) {
      data_path <- write_stata_parity_csv(case$data)
      on.exit(unlink(data_path), add = TRUE)

      fit_r <- case$fit_fun(case$data)
      out_r <- scpc(
        model = fit_r,
        data = case$data,
        lon = "lon",
        lat = "lat",
        cluster = case$cluster,
        avc = case$avc,
        ncoef = length(stats::coef(fit_r)),
        uncond = case$uncond,
        cvs = TRUE
      )

      stata_opts <- c(
        "latlong",
        sprintf("avc(%s)", format(case$avc, scientific = FALSE, trim = TRUE)),
        sprintf("k(%d)", length(case$terms)),
        if (isTRUE(case$uncond)) "uncond",
        "cvs"
      )
      stata_lines <- build_custom_stata_scpc_lines(
        data_path = data_path,
        reg_line = case$reg_line,
        scpc_opts = stata_opts
      )

      expect_custom_stata_scpc_match(
        stata_lines = stata_lines,
        out_r = out_r,
        terms_target = case$terms,
        label = case$id,
        stata_bin = ctx$stata_bin,
        stats_tol = STATA_STATS_TOL,
        cvs_tol = STATA_CVS_TOL
      )
    }
  })
})
