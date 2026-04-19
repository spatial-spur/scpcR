find_stata_binary <- function() {
  candidates <- c(
    Sys.which("stata-se"),
    Sys.which("stata-mp"),
    Sys.which("stata"),
    "/Applications/Stata/StataSE.app/Contents/MacOS/stata-se"
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0) {
    return(NA_character_)
  }
  hits[[1]]
}

run_stata_do <- function(stata_bin, lines) {
  do_file <- tempfile("scpc_parity_", fileext = ".do")
  base <- tools::file_path_sans_ext(basename(do_file))
  log_candidates <- c(
    sub("\\.do$", ".log", do_file),
    file.path(getwd(), paste0(base, ".log"))
  )
  for (p in log_candidates) {
    if (file.exists(p)) {
      unlink(p)
    }
  }

  writeLines(lines, do_file)
  output <- suppressWarnings(system2(stata_bin, c("-b", "do", do_file), stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  found_logs <- unique(log_candidates[file.exists(log_candidates)])
  list(
    status = status,
    output = output,
    do_file = do_file,
    log_file = if (length(found_logs)) found_logs[[1]] else NA_character_,
    log_files = found_logs
  )
}

stata_has_scpc <- function(stata_bin) {
  res <- run_stata_do(
    stata_bin,
    c(
      "clear all",
      "set more off",
      "which scpc"
    )
  )
  ok <- isTRUE(res$status == 0L)
  unlink(c(res$do_file, res$log_files))
  ok
}

resolve_data_path <- function() {
  p <- normalizePath(file.path(getwd(), "chetty_data_1.csv"), mustWork = FALSE)
  if (file.exists(p)) {
    return(p)
  }
  p <- normalizePath(file.path("..", "..", "chetty_data_1.csv"), mustWork = FALSE)
  if (file.exists(p)) {
    return(p)
  }
  NA_character_
}

get_scpc_fun <- local({
  cache <- NULL
  function(data_path) {
    if (!is.null(cache)) {
      return(cache)
    }
    if (exists("scpc", mode = "function")) {
      cache <<- get("scpc", mode = "function")
      return(cache)
    }
    cache <<- tryCatch(getExportedValue("scpcR", "scpc"), error = function(e) NULL)
    if (is.null(cache)) {
      source(file.path(dirname(data_path), "R", "scpcR.R"))
      cache <<- get("scpc", mode = "function")
    }
    cache
  }
})

prepare_data <- function(data_path, scenario) {
  d <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)
  names(d) <- tolower(names(d))
  vars_needed <- unique(c(scenario$dep, scenario$rhs, "lat", "lon"))
  miss <- setdiff(vars_needed, names(d))
  if (length(miss) > 0) {
    stop("Missing variables in chetty_data_1.csv: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  keep <- stats::complete.cases(d[, vars_needed, drop = FALSE])
  d <- d[keep, , drop = FALSE]
  rownames(d) <- seq_len(nrow(d))

  if (isTRUE(scenario$include_fe)) {
    d$fe <- (seq_len(nrow(d)) - 1L) %% 20L + 1L
  }

  if (identical(scenario$cluster_mode, "pair_mean")) {
    d$clust_id <- ceiling(seq_len(nrow(d)) / 2)
    d$lat <- ave(d$lat, d$clust_id, FUN = mean)
    d$lon <- ave(d$lon, d$clust_id, FUN = mean)
  }
  d
}

stata_export_helpers <- function(out_stats_path, out_cvs_path, write_cvs) {
  lines <- c(
    "capture program drop export_scpcstats",
    "program define export_scpcstats",
    "    syntax , outfile(string)",
    "    tempname M",
    "    matrix `M' = e(scpcstats)",
    "    local rn : rownames `M'",
    "    preserve",
    "    clear",
    "    set obs `=rowsof(`M')'",
    "    gen str64 term = \"\"",
    "    gen double coef = .",
    "    gen double std_err = .",
    "    gen double t = .",
    "    gen double p = .",
    "    gen double ci_low = .",
    "    gen double ci_high = .",
    "    forvalues i = 1/`=rowsof(`M')' {",
    "        replace coef = `M'[`i',1] in `i'",
    "        replace std_err = `M'[`i',2] in `i'",
    "        replace t = `M'[`i',3] in `i'",
    "        replace p = `M'[`i',4] in `i'",
    "        replace ci_low = `M'[`i',5] in `i'",
    "        replace ci_high = `M'[`i',6] in `i'",
    "    }",
    "    local i = 1",
    "    foreach nm of local rn {",
    "        replace term = \"`nm'\" in `i'",
    "        local ++i",
    "    }",
    "    export delimited using \"`outfile'\", replace",
    "    restore",
    "end",
    sprintf("export_scpcstats, outfile(\"%s\")", out_stats_path)
  )

  if (write_cvs) {
    lines <- c(
      lines,
      "capture program drop export_scpccvs",
      "program define export_scpccvs",
      "    syntax , outfile(string)",
      "    tempname M C",
      "    matrix `M' = e(scpcstats)",
      "    local rn : rownames `M'",
      "    matrix `C' = e(scpccvs)",
      "    preserve",
      "    clear",
      "    set obs `=rowsof(`C')'",
      "    gen str64 term = \"\"",
      "    gen double cv32 = .",
      "    gen double cv10 = .",
      "    gen double cv5 = .",
      "    gen double cv1 = .",
      "    forvalues i = 1/`=rowsof(`C')' {",
      "        replace cv32 = `C'[`i',1] in `i'",
      "        replace cv10 = `C'[`i',2] in `i'",
      "        replace cv5 = `C'[`i',3] in `i'",
      "        replace cv1 = `C'[`i',4] in `i'",
      "    }",
      "    local i = 1",
      "    foreach nm of local rn {",
      "        replace term = \"`nm'\" in `i'",
      "        local ++i",
      "    }",
      "    export delimited using \"`outfile'\", replace",
      "    restore",
      "end",
      sprintf("export_scpccvs, outfile(\"%s\")", out_cvs_path)
    )
  }
  lines
}

run_stata_scenario <- function(stata_bin, data_path, scenario, out_stats_path, out_cvs_path) {
  drop_expr <- paste(sprintf("missing(%s)", c(scenario$dep, scenario$rhs, "lat", "lon")), collapse = " | ")
  rhs_terms <- scenario$rhs
  if (isTRUE(scenario$include_fe)) {
    rhs_terms <- c(rhs_terms, "i.fe")
  }
  rhs <- paste(rhs_terms, collapse = " ")

  reg_line <- if (identical(scenario$cluster_mode, "pair_mean")) {
    sprintf("regress %s %s, cluster(clust_id)", scenario$dep, rhs)
  } else {
    sprintf("regress %s %s, robust", scenario$dep, rhs)
  }

  opts <- character()
  if (!is.null(scenario$avc)) {
    opts <- c(opts, sprintf("avc(%s)", format(scenario$avc, scientific = FALSE, trim = TRUE)))
  }
  if (!is.null(scenario$k)) {
    opts <- c(opts, sprintf("k(%d)", scenario$k))
  }
  if (isTRUE(scenario$latlong)) {
    opts <- c(opts, "latlong")
  }
  if (isTRUE(scenario$uncond)) {
    opts <- c(opts, "uncond")
  }
  if (isTRUE(scenario$cvs)) {
    opts <- c(opts, "cvs")
  }
  scpc_line <- if (length(opts) > 0) {
    sprintf("scpc, %s", paste(opts, collapse = " "))
  } else {
    "scpc"
  }

  lines <- c(
    "clear all",
    "set more off",
    sprintf("import delimited \"%s\", varnames(1) clear", data_path),
    "rename *, lower",
    sprintf("drop if %s", drop_expr),
    if (isTRUE(scenario$include_fe)) "gen long fe = mod(_n - 1, 20) + 1",
    if (identical(scenario$cluster_mode, "pair_mean")) "gen long clust_id = ceil(_n/2)",
    if (identical(scenario$cluster_mode, "pair_mean")) "bysort clust_id: egen s_1 = mean(lat)",
    if (identical(scenario$cluster_mode, "pair_mean")) "bysort clust_id: egen s_2 = mean(lon)",
    if (!identical(scenario$cluster_mode, "pair_mean")) "gen s_1 = lat",
    if (!identical(scenario$cluster_mode, "pair_mean")) "gen s_2 = lon",
    reg_line,
    scpc_line,
    stata_export_helpers(out_stats_path, out_cvs_path, isTRUE(scenario$cvs))
  )

  res <- run_stata_do(stata_bin, lines)
  failed <- res$status != 0L || !file.exists(out_stats_path)
  if (!failed && isTRUE(scenario$cvs) && !file.exists(out_cvs_path)) {
    failed <- TRUE
  }
  if (failed) {
    log_tail <- ""
    if (!is.na(res$log_file) && file.exists(res$log_file)) {
      lg <- readLines(res$log_file, warn = FALSE)
      log_tail <- paste(utils::tail(lg, 60), collapse = "\n")
    }
    unlink(c(res$do_file, res$log_files))
    stop(paste("Stata scenario failed:", scenario$id, log_tail), call. = FALSE)
  }
  unlink(c(res$do_file, res$log_files))
}

run_r_scenario <- function(data_path, scenario, out_stats_path, out_cvs_path) {
  d <- prepare_data(data_path, scenario)
  rhs <- paste(scenario$rhs, collapse = " + ")
  if (isTRUE(scenario$include_fe)) {
    rhs <- paste0(rhs, " + factor(fe)")
  }
  form <- stats::as.formula(sprintf("%s ~ %s", scenario$dep, rhs))
  fit <- stats::lm(form, data = d)
  scpc_fun <- get_scpc_fun(data_path)

  args <- list(
    model = fit,
    data = d,
    cluster = if (identical(scenario$cluster_mode, "pair_mean")) "clust_id" else NULL,
    uncond = scenario$uncond,
    cvs = scenario$cvs
  )
  if (isTRUE(scenario$latlong)) {
    args$lon <- "lon"
    args$lat <- "lat"
  } else {
    args$coord_euclidean <- c("lon", "lat")
  }
  if (!is.null(scenario$avc)) {
    args$avc <- scenario$avc
  }
  if (!is.null(scenario$k)) {
    # Stata keeps _cons last in e(b), while R keeps intercept first.
    # Increase k when intercept is present so Stata's first k terms are
    # represented in the R output for term-wise parity checks.
    k_adj <- as.integer(scenario$k)
    if ("(Intercept)" %in% names(stats::coef(fit))) {
      k_adj <- k_adj + 1L
    }
    args$ncoef <- k_adj
  }

  out <- do.call(scpc_fun, args)
  tab <- as.data.frame(out$scpcstats)
  tab$term <- rownames(out$scpcstats)
  rownames(tab) <- NULL
  names(tab) <- c("coef", "std_err", "t", "p", "ci_low", "ci_high", "term")
  tab <- tab[, c("term", "coef", "std_err", "t", "p", "ci_low", "ci_high")]
  write.csv(tab, out_stats_path, row.names = FALSE)

  if (isTRUE(scenario$cvs)) {
    cv <- as.data.frame(out$scpccvs)
    names(cv) <- c("cv32", "cv10", "cv5", "cv1")
    cv$term <- tab$term[seq_len(nrow(cv))]
    cv <- cv[, c("term", "cv32", "cv10", "cv5", "cv1")]
    write.csv(cv, out_cvs_path, row.names = FALSE)
  }
}

read_stata_scpcstats <- function(path) {
  tab <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  tab$term[tab$term == "_cons"] <- "(Intercept)"
  tab
}

assert_stats_close <- function(stata_tab, r_tab, tol = 2e-4) {
  merged <- merge(stata_tab, r_tab, by = "term", all.x = TRUE, suffixes = c("_stata", "_r"), sort = FALSE)
  expect_gt(nrow(merged), 0L)
  expect_true(all(!is.na(merged$coef_r)), info = "R output is missing one or more Stata terms.")
  for (v in c("coef", "std_err", "t", "p", "ci_low", "ci_high")) {
    d <- max(abs(merged[[paste0(v, "_stata")]] - merged[[paste0(v, "_r")]]))
    expect_true(d < tol, info = paste("max abs diff for", v, "=", format(d, scientific = TRUE)))
  }
}

assert_cvs_close <- function(stata_tab, r_tab, tol = 2e-4) {
  merged <- merge(stata_tab, r_tab, by = "term", all.x = TRUE, suffixes = c("_stata", "_r"), sort = FALSE)
  expect_gt(nrow(merged), 0L)
  expect_true(all(!is.na(merged$cv32_r)), info = "R cvs output is missing one or more Stata terms.")
  for (v in c("cv32", "cv10", "cv5", "cv1")) {
    d <- max(abs(merged[[paste0(v, "_stata")]] - merged[[paste0(v, "_r")]]))
    expect_true(d < tol, info = paste("max abs diff for", v, "=", format(d, scientific = TRUE)))
  }
}

test_that("scpcR matches Stata scpc across models and options (non-clustered)", {
  data_path <- resolve_data_path()
  skip_if(is.na(data_path), "Missing chetty_data_1.csv in repository base.")

  stata_bin <- find_stata_binary()
  skip_if(is.na(stata_bin), "Stata binary not found.")
  skip_if_not(stata_has_scpc(stata_bin), "Stata scpc command not installed.")

  scenarios <- list(
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

  tmpdir <- tempfile("scpc_stata_parity_grid_")
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE), add = TRUE)

  for (sc in scenarios) {
    stata_stats <- file.path(tmpdir, paste0("stata_stats_", sc$id, ".csv"))
    r_stats <- file.path(tmpdir, paste0("r_stats_", sc$id, ".csv"))
    stata_cvs <- file.path(tmpdir, paste0("stata_cvs_", sc$id, ".csv"))
    r_cvs <- file.path(tmpdir, paste0("r_cvs_", sc$id, ".csv"))

    run_stata_scenario(stata_bin, data_path, sc, stata_stats, stata_cvs)
    run_r_scenario(data_path, sc, r_stats, r_cvs)

    s_stats <- read_stata_scpcstats(stata_stats)
    r_stats_tab <- read.csv(r_stats, stringsAsFactors = FALSE, check.names = FALSE)
    assert_stats_close(s_stats, r_stats_tab, tol = 2e-4)

    if (isTRUE(sc$cvs)) {
      s_cvs <- read.csv(stata_cvs, stringsAsFactors = FALSE, check.names = FALSE)
      s_cvs$term[s_cvs$term == "_cons"] <- "(Intercept)"
      r_cvs_tab <- read.csv(r_cvs, stringsAsFactors = FALSE, check.names = FALSE)
      assert_cvs_close(s_cvs, r_cvs_tab, tol = 2e-4)
    }
  }
})

test_that("clustered pair-mean scenario matches Stata", {
  data_path <- resolve_data_path()
  skip_if(is.na(data_path), "Missing chetty_data_1.csv in repository base.")

  stata_bin <- find_stata_binary()
  skip_if(is.na(stata_bin), "Stata binary not found.")
  skip_if_not(stata_has_scpc(stata_bin), "Stata scpc command not installed.")

  sc <- list(
    id = "cluster_pairmean_cvs",
    dep = "am",
    rhs = c("tlfpr", "fracblack", "gini"),
    latlong = TRUE,
    uncond = TRUE,
    avc = 0.03,
    k = 3L,
    cvs = TRUE,
    cluster_mode = "pair_mean"
  )

  tmpdir <- tempfile("scpc_stata_cluster_gap_")
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE), add = TRUE)

  stata_stats <- file.path(tmpdir, "stata_stats_cluster.csv")
  r_stats <- file.path(tmpdir, "r_stats_cluster.csv")
  stata_cvs <- file.path(tmpdir, "stata_cvs_cluster.csv")
  r_cvs <- file.path(tmpdir, "r_cvs_cluster.csv")

  run_stata_scenario(stata_bin, data_path, sc, stata_stats, stata_cvs)
  run_r_scenario(data_path, sc, r_stats, r_cvs)

  s_stats <- read_stata_scpcstats(stata_stats)
  r_stats_tab <- read.csv(r_stats, stringsAsFactors = FALSE, check.names = FALSE)
  assert_stats_close(s_stats, r_stats_tab, tol = 2e-4)

  s_cvs <- read.csv(stata_cvs, stringsAsFactors = FALSE, check.names = FALSE)
  s_cvs$term[s_cvs$term == "_cons"] <- "(Intercept)"
  r_cvs_tab <- read.csv(r_cvs, stringsAsFactors = FALSE, check.names = FALSE)
  assert_cvs_close(s_cvs, r_cvs_tab, tol = 2e-4)
})

test_that("clustered pair-mean conditional scenario matches Stata", {
  data_path <- resolve_data_path()
  skip_if(is.na(data_path), "Missing chetty_data_1.csv in repository base.")

  stata_bin <- find_stata_binary()
  skip_if(is.na(stata_bin), "Stata binary not found.")
  skip_if_not(stata_has_scpc(stata_bin), "Stata scpc command not installed.")

  sc <- list(
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

  tmpdir <- tempfile("scpc_stata_cluster_cond_")
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE), add = TRUE)

  stata_stats <- file.path(tmpdir, "stata_stats_cluster_cond.csv")
  r_stats <- file.path(tmpdir, "r_stats_cluster_cond.csv")
  stata_cvs <- file.path(tmpdir, "stata_cvs_cluster_cond.csv")
  r_cvs <- file.path(tmpdir, "r_cvs_cluster_cond.csv")

  run_stata_scenario(stata_bin, data_path, sc, stata_stats, stata_cvs)
  run_r_scenario(data_path, sc, r_stats, r_cvs)

  s_stats <- read_stata_scpcstats(stata_stats)
  r_stats_tab <- read.csv(r_stats, stringsAsFactors = FALSE, check.names = FALSE)
  assert_stats_close(s_stats, r_stats_tab, tol = 2e-4)

  s_cvs <- read.csv(stata_cvs, stringsAsFactors = FALSE, check.names = FALSE)
  s_cvs$term[s_cvs$term == "_cons"] <- "(Intercept)"
  r_cvs_tab <- read.csv(r_cvs, stringsAsFactors = FALSE, check.names = FALSE)
  assert_cvs_close(s_cvs, r_cvs_tab, tol = 2e-4)
})

test_that("absorbed FE IV conditional SCPC matches Stata ivregress with i.fe", {
  skip_if_not_installed("fixest")
  old_threads <- fixest::getFixest_nthreads()
  fixest::setFixest_nthreads(1)
  on.exit(fixest::setFixest_nthreads(old_threads), add = TRUE)

  data_path <- resolve_data_path()
  skip_if(is.na(data_path), "Missing chetty_data_1.csv in repository base.")

  stata_bin <- find_stata_binary()
  skip_if(is.na(stata_bin), "Stata binary not found.")
  skip_if_not(stata_has_scpc(stata_bin), "Stata scpc command not installed.")

  tmpdir <- tempfile("scpc_stata_iv_fe_")
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE), add = TRUE)

  stata_stats <- file.path(tmpdir, "stata_stats_iv_fe.csv")
  stata_cvs <- file.path(tmpdir, "stata_cvs_iv_fe.csv")

  lines <- c(
    "clear all",
    "set more off",
    sprintf("import delimited \"%s\", varnames(1) clear", data_path),
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
    "scpc, latlong avc(0.03) k(2) cvs",
    stata_export_helpers(stata_stats, stata_cvs, write_cvs = TRUE)
  )

  stata_res <- run_stata_do(stata_bin, lines)
  failed <- stata_res$status != 0L || !file.exists(stata_stats) || !file.exists(stata_cvs)
  if (failed) {
    log_tail <- ""
    if (!is.na(stata_res$log_file) && file.exists(stata_res$log_file)) {
      lg <- readLines(stata_res$log_file, warn = FALSE)
      log_tail <- paste(utils::tail(lg, 80), collapse = "\n")
    }
    unlink(c(stata_res$do_file, stata_res$log_files))
    stop(paste("Stata IV+FE scenario failed:", log_tail), call. = FALSE)
  }
  unlink(c(stata_res$do_file, stata_res$log_files))

  d <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)
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

  r_stats <- as.data.frame(out_r$scpcstats)
  r_stats$term <- rownames(out_r$scpcstats)
  rownames(r_stats) <- NULL
  names(r_stats) <- c("coef", "std_err", "t", "p", "ci_low", "ci_high", "term")
  r_stats$term <- sub("^fit_", "", r_stats$term)
  r_stats <- r_stats[, c("term", "coef", "std_err", "t", "p", "ci_low", "ci_high")]

  r_cvs <- as.data.frame(out_r$scpccvs)
  names(r_cvs) <- c("cv32", "cv10", "cv5", "cv1")
  r_cvs$term <- sub("^fit_", "", rownames(out_r$scpcstats)[seq_len(nrow(r_cvs))])
  r_cvs <- r_cvs[, c("term", "cv32", "cv10", "cv5", "cv1")]

  s_stats <- read_stata_scpcstats(stata_stats)
  s_cvs <- read.csv(stata_cvs, stringsAsFactors = FALSE, check.names = FALSE)

  terms_target <- c("x_endog", "fracblack")
  s_stats <- s_stats[s_stats$term %in% terms_target, , drop = FALSE]
  r_stats <- r_stats[r_stats$term %in% terms_target, , drop = FALSE]
  expect_equal(sort(s_stats$term), sort(terms_target))
  expect_equal(sort(r_stats$term), sort(terms_target))
  merged_stats <- merge(s_stats, r_stats, by = "term", suffixes = c("_stata", "_r"), sort = FALSE)
  tol_by_var <- c(
    coef = 1e-3,
    std_err = 1e-3,
    t = 1e-4,
    p = 2e-3,
    ci_low = 0.5,
    ci_high = 0.5
  )
  for (v in names(tol_by_var)) {
    d <- max(abs(merged_stats[[paste0(v, "_stata")]] - merged_stats[[paste0(v, "_r")]]))
    expect_true(d < tol_by_var[[v]], info = paste("max abs diff for", v, "=", format(d, scientific = TRUE)))
  }

  s_cvs <- s_cvs[s_cvs$term %in% terms_target, , drop = FALSE]
  r_cvs <- r_cvs[r_cvs$term %in% terms_target, , drop = FALSE]
  expect_equal(sort(s_cvs$term), sort(terms_target))
  expect_equal(sort(r_cvs$term), sort(terms_target))
  merged_cvs <- merge(s_cvs, r_cvs, by = "term", suffixes = c("_stata", "_r"), sort = FALSE)
  for (v in c("cv32", "cv10", "cv5", "cv1")) {
    d <- max(abs(merged_cvs[[paste0(v, "_stata")]] - merged_cvs[[paste0(v, "_r")]]))
    expect_true(d < 0.12, info = paste("max abs diff for", v, "=", format(d, scientific = TRUE)))
  }
})

test_that("large-n synthetic scenarios match Stata in both unconditional and conditional modes", {
  stata_bin <- find_stata_binary()
  skip_if(is.na(stata_bin), "Stata binary not found.")
  skip_if_not(stata_has_scpc(stata_bin), "Stata scpc command not installed.")

  set.seed(1)
  n <- 4600
  dat <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    lon = runif(n, min = -125, max = -66),
    lat = runif(n, min = 25, max = 49)
  )
  data_path <- tempfile("scpc_large_n_", fileext = ".csv")
  write.csv(dat, data_path, row.names = FALSE)

  tmpdir <- tempfile("scpc_stata_large_n_")
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE), add = TRUE)

  run_large_n_case <- function(id, uncond) {
    stata_stats <- file.path(tmpdir, paste0("stata_stats_", id, ".csv"))
    r_stats <- file.path(tmpdir, paste0("r_stats_", id, ".csv"))
    opts <- c("latlong", "avc(0.1)", "k(1)")
    if (isTRUE(uncond)) {
      opts <- c(opts, "uncond")
    }

    lines <- c(
      "clear all",
      "set more off",
      sprintf("import delimited \"%s\", varnames(1) clear", data_path),
      "regress y x, robust",
      "gen s_1 = lat",
      "gen s_2 = lon",
      sprintf("scpc, %s", paste(opts, collapse = " ")),
      stata_export_helpers(stata_stats, tempfile(fileext = ".csv"), write_cvs = FALSE)
    )

    stata_res <- run_stata_do(stata_bin, lines)
    failed <- stata_res$status != 0L || !file.exists(stata_stats)
    if (failed) {
      log_tail <- ""
      if (!is.na(stata_res$log_file) && file.exists(stata_res$log_file)) {
        lg <- readLines(stata_res$log_file, warn = FALSE)
        log_tail <- paste(utils::tail(lg, 80), collapse = "\n")
      }
      unlink(c(stata_res$do_file, stata_res$log_files))
      stop(paste("Stata large-n scenario failed:", id, log_tail), call. = FALSE)
    }
    unlink(c(stata_res$do_file, stata_res$log_files))

    fit <- stats::lm(y ~ x, data = dat)
    out <- scpc(
      model = fit,
      data = dat,
      lon = "lon",
      lat = "lat",
      ncoef = 2,
      avc = 0.1,
      method = "approx",
      large_n_seed = 1,
      uncond = uncond,
      cvs = FALSE
    )

    r_stats_tab <- as.data.frame(out$scpcstats)
    r_stats_tab$term <- rownames(out$scpcstats)
    rownames(r_stats_tab) <- NULL
    names(r_stats_tab) <- c("coef", "std_err", "t", "p", "ci_low", "ci_high", "term")
    r_stats_tab <- r_stats_tab[, c("term", "coef", "std_err", "t", "p", "ci_low", "ci_high")]
    write.csv(r_stats_tab, r_stats, row.names = FALSE)

    s_stats <- read_stata_scpcstats(stata_stats)
    r_stats_tab <- read.csv(r_stats, stringsAsFactors = FALSE, check.names = FALSE)
    assert_stats_close(s_stats, r_stats_tab, tol = 2e-4)
  }

  run_large_n_case("large_n_uncond", uncond = TRUE)
  run_large_n_case("large_n_cond", uncond = FALSE)
})
