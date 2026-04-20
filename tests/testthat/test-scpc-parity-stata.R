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
