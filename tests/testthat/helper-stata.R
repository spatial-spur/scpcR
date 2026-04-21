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
  for (path in log_candidates) {
    if (file.exists(path)) {
      unlink(path)
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

stata_log_tail <- function(result, n = 60L) {
  if (!is.na(result$log_file) && file.exists(result$log_file)) {
    log_lines <- readLines(result$log_file, warn = FALSE)
    return(paste(utils::tail(log_lines, n), collapse = "\n"))
  }
  paste(result$output, collapse = "\n")
}

stata_has_scpc <- function(stata_bin = find_stata_binary()) {
  if (is.na(stata_bin)) {
    return(FALSE)
  }
  result <- run_stata_do(
    stata_bin,
    c(
      "clear all",
      "set more off",
      "which scpc"
    )
  )
  ok <- isTRUE(result$status == 0L)
  unlink(c(result$do_file, result$log_files))
  ok
}

skip_if_stata_scpc_unavailable <- function() {
  skip_on_cran()

  data_path <- scpc_example_data_path()
  skip_if(is.na(data_path), "Missing chetty_data_1.csv in package extdata.")

  stata_bin <- find_stata_binary()
  skip_if(is.na(stata_bin), "Stata binary not found.")
  skip_if_not(stata_has_scpc(stata_bin), "Stata scpc command not installed.")

  list(stata_bin = stata_bin, data_path = data_path)
}

stata_export_helpers <- function(out_stats_path, out_cvs_path = NULL, write_cvs = FALSE) {
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

  if (isTRUE(write_cvs)) {
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

read_stata_scpcstats <- function(path) {
  tab <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  tab$term[tab$term == "_cons"] <- "(Intercept)"
  tab
}

read_stata_scpccvs <- function(path) {
  tab <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  tab$term[tab$term == "_cons"] <- "(Intercept)"
  tab
}

normalize_r_scpcstats <- function(out) {
  tab <- as.data.frame(out$scpcstats)
  tab$term <- rownames(out$scpcstats)
  rownames(tab) <- NULL
  names(tab) <- c("coef", "std_err", "t", "p", "ci_low", "ci_high", "term")
  tab[, c("term", "coef", "std_err", "t", "p", "ci_low", "ci_high")]
}

normalize_r_scpccvs <- function(out) {
  if (is.null(out$scpccvs)) {
    return(NULL)
  }
  tab <- as.data.frame(out$scpccvs)
  names(tab) <- c("cv32", "cv10", "cv5", "cv1")
  tab$term <- rownames(out$scpccvs)
  tab[, c("term", "cv32", "cv10", "cv5", "cv1")]
}

prepare_parity_data <- function(data_path, scenario) {
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

stata_scpc_option_tokens <- function(scenario) {
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
  opts
}

build_stata_scpc_lines <- function(data_path, scenario) {
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

  opts <- stata_scpc_option_tokens(scenario)
  scpc_line <- if (length(opts) > 0) {
    sprintf("scpc, %s", paste(opts, collapse = " "))
  } else {
    "scpc"
  }

  c(
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
    scpc_line
  )
}

run_stata_scpc_lines <- function(lines, write_cvs = FALSE, stata_bin = find_stata_binary()) {
  tmpdir <- tempfile("scpc_stata_run_")
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE), add = TRUE)

  out_stats_path <- file.path(tmpdir, "stata_stats.csv")
  out_cvs_path <- file.path(tmpdir, "stata_cvs.csv")
  result <- run_stata_do(
    stata_bin,
    c(lines, stata_export_helpers(out_stats_path, out_cvs_path, write_cvs = write_cvs))
  )
  on.exit(unlink(c(result$do_file, result$log_files)), add = TRUE)

  failed <- result$status != 0L || !file.exists(out_stats_path)
  if (!failed && isTRUE(write_cvs) && !file.exists(out_cvs_path)) {
    failed <- TRUE
  }
  if (failed) {
    stop("Stata scpc execution failed:\n", stata_log_tail(result, n = 80L), call. = FALSE)
  }

  list(
    stats = read_stata_scpcstats(out_stats_path),
    cvs = if (isTRUE(write_cvs)) read_stata_scpccvs(out_cvs_path) else NULL
  )
}

run_stata_scpc_scenario <- function(scenario, data_path = scpc_example_data_path(), stata_bin = find_stata_binary()) {
  run_stata_scpc_lines(
    build_stata_scpc_lines(data_path, scenario),
    write_cvs = isTRUE(scenario$cvs),
    stata_bin = stata_bin
  )
}

run_r_scpc_scenario <- function(scenario, data_path = scpc_example_data_path()) {
  d <- prepare_parity_data(data_path, scenario)
  rhs <- paste(scenario$rhs, collapse = " + ")
  if (isTRUE(scenario$include_fe)) {
    rhs <- paste0(rhs, " + factor(fe)")
  }
  form <- stats::as.formula(sprintf("%s ~ %s", scenario$dep, rhs))
  fit <- stats::lm(form, data = d)

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
    args$coords_euclidean <- c("lon", "lat")
  }
  if (!is.null(scenario$avc)) {
    args$avc <- scenario$avc
  }
  if (!is.null(scenario$k)) {
    k_adj <- as.integer(scenario$k)
    if ("(Intercept)" %in% names(stats::coef(fit))) {
      k_adj <- k_adj + 1L
    }
    args$ncoef <- k_adj
  }

  out <- do.call(scpc, args)
  list(
    stats = normalize_r_scpcstats(out),
    cvs = normalize_r_scpccvs(out)
  )
}

normalize_tolerance <- function(vars, tolerance) {
  if (length(tolerance) == 1L && is.null(names(tolerance))) {
    return(stats::setNames(rep(as.numeric(tolerance), length(vars)), vars))
  }

  if (is.null(names(tolerance))) {
    stop("Tolerance vectors must be named when they are not scalar.", call. = FALSE)
  }

  missing_vars <- setdiff(vars, names(tolerance))
  if (length(missing_vars) > 0L) {
    stop("Missing tolerance values for: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  as.numeric(tolerance[vars]) |> stats::setNames(vars)
}

table_max_abs_diffs <- function(stata_tab, r_tab, vars) {
  merged <- merge(stata_tab, r_tab, by = "term", all.x = TRUE, suffixes = c("_stata", "_r"), sort = FALSE)
  if (nrow(merged) == 0L) {
    stop("No overlapping terms to compare.", call. = FALSE)
  }
  if (anyNA(merged[[paste0(vars[[1]], "_r")]])) {
    stop("R output is missing one or more Stata terms.", call. = FALSE)
  }

  stats::setNames(vapply(
    vars,
    function(v) max(abs(merged[[paste0(v, "_stata")]] - merged[[paste0(v, "_r")]])),
    numeric(1)
  ), vars)
}

expect_table_close <- function(stata_tab, r_tab, vars, tolerance, label) {
  diffs <- table_max_abs_diffs(stata_tab, r_tab, vars)
  tolerance <- normalize_tolerance(vars, tolerance)
  for (v in vars) {
    expect_true(
      diffs[[v]] <= tolerance[[v]],
      info = paste0(
        label, ": max abs diff for ", v, " = ",
        format(diffs[[v]], scientific = TRUE),
        " > ", format(tolerance[[v]], scientific = TRUE)
      )
    )
  }
  invisible(diffs)
}

compare_stata_scpc_scenario <- function(scenario, stata_bin = find_stata_binary(), data_path = scpc_example_data_path()) {
  stata_out <- run_stata_scpc_scenario(scenario, data_path = data_path, stata_bin = stata_bin)
  r_out <- run_r_scpc_scenario(scenario, data_path = data_path)
  list(stata = stata_out, r = r_out)
}

expect_stata_scpc_match <- function(scenario, stats_tol, cvs_tol = NULL,
                                    stata_bin = find_stata_binary(),
                                    data_path = scpc_example_data_path()) {
  out <- compare_stata_scpc_scenario(scenario, stata_bin = stata_bin, data_path = data_path)
  stats_diffs <- expect_table_close(
    out$stata$stats,
    out$r$stats,
    vars = c("coef", "std_err", "t", "p", "ci_low", "ci_high"),
    tolerance = stats_tol,
    label = scenario$id
  )

  cvs_diffs <- NULL
  if (isTRUE(scenario$cvs)) {
    cvs_diffs <- expect_table_close(
      out$stata$cvs,
      out$r$cvs,
      vars = c("cv32", "cv10", "cv5", "cv1"),
      tolerance = cvs_tol,
      label = paste0(scenario$id, " cvs")
    )
  }

  invisible(list(stats = stats_diffs, cvs = cvs_diffs))
}
