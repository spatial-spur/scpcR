parse_json_minimal <- function(text) {
  chars <- strsplit(text, "", fixed = TRUE)[[1L]]
  n <- length(chars)
  pos <- 1L

  skip_ws <- function() {
    while (pos <= n && chars[pos] %in% c(" ", "\t", "\n", "\r")) {
      pos <<- pos + 1L
    }
  }

  parse_string <- function() {
    pos <<- pos + 1L
    out <- character()
    while (pos <= n) {
      ch <- chars[pos]
      if (ch == "\"") {
        pos <<- pos + 1L
        return(paste0(out, collapse = ""))
      }
      if (ch == "\\") {
        pos <<- pos + 1L
        esc <- chars[pos]
        out <- c(out, switch(esc,
          "\"" = "\"",
          "\\" = "\\",
          "/" = "/",
          "b" = "\b",
          "f" = "\f",
          "n" = "\n",
          "r" = "\r",
          "t" = "\t",
          esc
        ))
        pos <<- pos + 1L
      } else {
        out <- c(out, ch)
        pos <<- pos + 1L
      }
    }
    stop("Unterminated JSON string.")
  }

  parse_number <- function() {
    start <- pos
    while (pos <= n && grepl("[-0-9eE.+]", chars[pos])) {
      pos <<- pos + 1L
    }
    as.numeric(paste0(chars[start:(pos - 1L)], collapse = ""))
  }

  parse_value <- function() {
    skip_ws()
    ch <- chars[pos]
    if (ch == "{") {
      pos <<- pos + 1L
      skip_ws()
      res <- list()
      if (chars[pos] == "}") {
        pos <<- pos + 1L
        return(res)
      }
      repeat {
        skip_ws()
        key <- parse_string()
        skip_ws()
        pos <<- pos + 1L # :
        res[[key]] <- parse_value()
        skip_ws()
        if (chars[pos] == ",") {
          pos <<- pos + 1L
        } else if (chars[pos] == "}") {
          pos <<- pos + 1L
          break
        } else {
          stop("Malformed JSON object.")
        }
      }
      return(res)
    }
    if (ch == "[") {
      pos <<- pos + 1L
      skip_ws()
      items <- list()
      if (chars[pos] == "]") {
        pos <<- pos + 1L
        return(items)
      }
      repeat {
        items[[length(items) + 1L]] <- parse_value()
        skip_ws()
        if (chars[pos] == ",") {
          pos <<- pos + 1L
        } else if (chars[pos] == "]") {
          pos <<- pos + 1L
          break
        } else {
          stop("Malformed JSON array.")
        }
      }
      return(items)
    }
    if (ch == "\"") {
      return(parse_string())
    }
    if (ch == "t") {
      pos <<- pos + 4L
      return(TRUE)
    }
    if (ch == "f") {
      pos <<- pos + 5L
      return(FALSE)
    }
    if (ch == "n") {
      pos <<- pos + 4L
      return(NULL)
    }
    parse_number()
  }

  parse_value()
}

find_spur_file <- function(rel_path) {
  candidates <- c(
    normalizePath(file.path("..", "..", "..", rel_path), mustWork = FALSE),
    normalizePath(file.path("..", "..", rel_path), mustWork = FALSE),
    normalizePath(file.path(getwd(), rel_path), mustWork = FALSE),
    normalizePath(file.path("/Users/pauldavidboll/Github/spatial-spur", rel_path), mustWork = FALSE)
  )
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0L) {
    return(NA_character_)
  }
  hits[[1L]]
}

json_rows_to_matrix <- function(rows) {
  do.call(rbind, lapply(rows, as.numeric))
}

test_that("conditional OLS cutoffs match the fitted design with and without an intercept", {
  sites <- c(0, 1, 2, 4, 6)
  x <- c(-1, -4, 2, 2, 1)
  d <- data.frame(s = sites, x = x, y = x + (sites == 0))
  backends <- "fallback"
  if (is.loaded("wrap__rust_run_exact_scpc", PACKAGE = "scpcR")) {
    backends <- c(backends, "rust")
  }
  for (backend in backends) {
    with_intercept <- scpc(lm(y ~ x, d), d, coords_euclidean = "s", avc = .5, cvs = TRUE, backend = backend)
    without_intercept <- scpc(lm(y ~ x - 1, d), d, coords_euclidean = "s", avc = .5, cvs = TRUE, backend = backend)
    expect_lt(abs(with_intercept$scpccvs["x", "5%"] - 17.625169), .002)
    expect_lt(abs(without_intercept$scpccvs["x", "5%"] - 9.973326), .002)
  }

  d <- d[rep(seq_len(nrow(d)), each = 2), ]
  d$x <- d$x + c(-.25, .25)
  d$y <- d$x + (d$s == 0)
  with_intercept <- scpc(lm(y ~ x, d), d, coords_euclidean = "s", cluster = "s", avc = .5, cvs = TRUE, backend = "fallback")
  without_intercept <- scpc(lm(y ~ x - 1, d), d, coords_euclidean = "s", cluster = "s", avc = .5, cvs = TRUE, backend = "fallback")
  expect_lt(abs(with_intercept$scpccvs["x", "5%"] - 17.380977), .002)
  expect_lt(abs(without_intercept$scpccvs["x", "5%"] - 9.902023), .002)
})

assert_relative_close <- function(actual, expected, rel_tol) {
  expect_equal(dim(actual), dim(expected))
  rel_err <- abs(actual - expected) / pmax(abs(expected), 1e-15)
  expect_true(all(rel_err <= rel_tol), info = paste("max rel_err =", max(rel_err)))
}

assert_scpc_matches_fixture <- function(out, fixture, uncond, tol = 1e-10, cal_rel_tol = 2e-3) {
  expected_stats <- json_rows_to_matrix(fixture$scpcstats)
  rownames(expected_stats) <- unlist(fixture$coef_names)
  colnames(expected_stats) <- c("Coef", "Std_Err", "t", "P>|t|", "2.5 %", "97.5 %")

  expect_equal(
    unname(out$scpcstats[, 1L]),
    unname(expected_stats[, 1L]),
    tolerance = tol
  )
  # Frozen fixtures reflect the historical 0.001 bisection bracket for c0/cv/SE/t.
  assert_relative_close(
    unname(out$scpcstats[, 2:3, drop = FALSE]),
    unname(expected_stats[, 2:3, drop = FALSE]),
    cal_rel_tol
  )
  if (uncond) {
    assert_relative_close(
      unname(out$scpcstats[, 5:6, drop = FALSE]),
      unname(expected_stats[, 5:6, drop = FALSE]),
      cal_rel_tol
    )
    expect_lte(
      max(abs(unname(out$scpcstats[, 4L]) - unname(expected_stats[, 4L]))),
      1e-6
    )
  }
  expect_equal(rownames(out$scpcstats), rownames(expected_stats))
  expect_lte(abs(out$c0 - fixture$c0) / abs(fixture$c0), cal_rel_tol)
  expect_lte(abs(out$cv - fixture$cv) / abs(fixture$cv), cal_rel_tol)
  expect_equal(out$q, fixture$q)

  if (is.null(fixture$scpccvs)) {
    expect_null(out$scpccvs)
  } else if (uncond) {
    expected_cvs <- json_rows_to_matrix(fixture$scpccvs)
    assert_relative_close(unname(out$scpccvs), unname(expected_cvs), cal_rel_tol)
  }
}

test_that("Rust and fallback agree and preserve the stable Milestone 1 results", {
  fixture_path <- find_spur_file("spatial-rust/fixtures/m1_exact_scpc.json")
  bench_200_path <- find_spur_file("verification/data/bench_200.csv")
  bench_1000_path <- find_spur_file("verification/data/bench_1000.csv")
  skip_if(
    is.na(fixture_path) || is.na(bench_200_path) || is.na(bench_1000_path),
    "spatial-rust fixture or verification CSV files not available"
  )

  fixtures <- parse_json_minimal(paste(readLines(fixture_path, warn = FALSE), collapse = "\n"))$r
  d200 <- utils::read.csv(bench_200_path)
  fit200 <- stats::lm(y ~ x1 + x2, data = d200)

  d1000 <- utils::read.csv(bench_1000_path)
  fit1000 <- stats::lm(y ~ x1 + x2, data = d1000)

  rust_available <- is.loaded("wrap__rust_run_exact_scpc", PACKAGE = "scpcR")
  backends <- if (rust_available) c("rust", "fallback") else "fallback"

  cases <- list(
    list(name = "bench_200_euclid_cond", fit = fit200, data = d200, latlong = FALSE, uncond = FALSE, cvs = TRUE),
    list(name = "bench_200_euclid_uncond", fit = fit200, data = d200, latlong = FALSE, uncond = TRUE, cvs = FALSE),
    list(name = "bench_200_latlong_cond", fit = fit200, data = d200, latlong = TRUE, uncond = FALSE, cvs = TRUE),
    list(name = "bench_1000_euclid_cond", fit = fit1000, data = d1000, latlong = FALSE, uncond = FALSE, cvs = TRUE),
    list(name = "bench_1000_euclid_uncond", fit = fit1000, data = d1000, latlong = FALSE, uncond = TRUE, cvs = FALSE),
    list(name = "bench_1000_latlong_cond", fit = fit1000, data = d1000, latlong = TRUE, uncond = FALSE, cvs = TRUE)
  )

  for (cs in cases) {
    outputs <- list()
    for (bk in backends) {
      if (cs$latlong && bk == "fallback" && !requireNamespace("geodist", quietly = TRUE)) {
        next
      }
      out <- if (cs$latlong) {
        scpc(
          cs$fit, cs$data,
          lon = "s_1", lat = "s_2",
          avc = 0.03, uncond = cs$uncond, cvs = cs$cvs,
          backend = bk
        )
      } else {
        scpc(
          cs$fit, cs$data,
          coords_euclidean = c("s_1", "s_2"),
          avc = 0.03, uncond = cs$uncond, cvs = cs$cvs,
          backend = bk
        )
      }
      expect_identical(out$backend, bk)
      expect_false(isTRUE(out$tied))
      if (!cs$latlong || requireNamespace("geodist", quietly = TRUE)) {
        distv <- .lvech(.getdistmat(cs$data[, c("s_1", "s_2"), drop = FALSE], cs$latlong))
        expect_lte(abs(mean(exp(-out$c0 * distv)) - 0.03), 1e-12)
      }
      assert_scpc_matches_fixture(out, fixtures[[cs$name]], uncond = cs$uncond)
      outputs[[bk]] <- out
    }
    if (length(outputs) == 2L) {
      expect_equal(outputs$rust$scpcstats, outputs$fallback$scpcstats, tolerance = 1e-10)
      expect_equal(outputs$rust$scpccvs, outputs$fallback$scpccvs, tolerance = 1e-10)
      expect_equal(outputs$rust$c0, outputs$fallback$c0, tolerance = 1e-10)
      expect_equal(outputs$rust$cv, outputs$fallback$cv, tolerance = 1e-10)
      expect_identical(outputs$rust$q, outputs$fallback$q)
      expect_identical(outputs$rust$tied, outputs$fallback$tied)
    }
  }

  if (rust_available) {
    # Also check additional uncond/cvs combinations between rust and fallback
    extra_cases <- list(
      list(fit = fit200, data = d200, latlong = FALSE, uncond = FALSE, cvs = FALSE),
      list(fit = fit200, data = d200, latlong = FALSE, uncond = TRUE, cvs = TRUE),
      list(fit = fit200, data = d200, latlong = TRUE, uncond = TRUE, cvs = FALSE),
      list(fit = fit200, data = d200, latlong = TRUE, uncond = TRUE, cvs = TRUE)
    )
    for (cs in extra_cases) {
      if (cs$latlong && !requireNamespace("geodist", quietly = TRUE)) {
        next
      }
      run_pair <- function(bk) {
        if (cs$latlong) {
          scpc(cs$fit, cs$data, lon = "s_1", lat = "s_2", avc = 0.03, uncond = cs$uncond, cvs = cs$cvs, backend = bk)
        } else {
          scpc(cs$fit, cs$data, coords_euclidean = c("s_1", "s_2"), avc = 0.03, uncond = cs$uncond, cvs = cs$cvs, backend = bk)
        }
      }
      out_r <- run_pair("rust")
      out_f <- run_pair("fallback")
      expect_equal(out_r$scpcstats, out_f$scpcstats, tolerance = 1e-10)
      expect_equal(out_r$scpccvs, out_f$scpccvs, tolerance = 1e-10)
      expect_equal(out_r$c0, out_f$c0, tolerance = 1e-10)
      expect_equal(out_r$cv, out_f$cv, tolerance = 1e-10)
      expect_equal(out_r$q, out_f$q)
    }
  }
})

test_that("nthreads = 1 runs Rust on 1 thread and matches multi-threaded and fallback", {
  skip_if_not(is.loaded("wrap__rust_run_exact_scpc", PACKAGE = "scpcR"), "native Rust library not loaded")
  bench_200_path <- find_spur_file("verification/data/bench_200.csv")
  skip_if(is.na(bench_200_path), "bench_200.csv not available")

  d200 <- utils::read.csv(bench_200_path)
  fit200 <- stats::lm(y ~ x1 + x2, data = d200)

  out_1 <- scpc(
    fit200, d200,
    coords_euclidean = c("s_1", "s_2"),
    avc = 0.03, uncond = FALSE, cvs = TRUE,
    backend = "rust", nthreads = 1L
  )
  out_4 <- scpc(
    fit200, d200,
    coords_euclidean = c("s_1", "s_2"),
    avc = 0.03, uncond = FALSE, cvs = TRUE,
    backend = "rust", nthreads = 4L
  )
  out_fb <- scpc(
    fit200, d200,
    coords_euclidean = c("s_1", "s_2"),
    avc = 0.03, uncond = FALSE, cvs = TRUE,
    backend = "fallback"
  )

  expect_identical(out_1$backend, "rust")
  expect_identical(out_4$backend, "rust")
  expect_identical(out_fb$backend, "fallback")
  expect_equal(out_1$scpcstats, out_4$scpcstats, tolerance = 1e-12)
  expect_equal(out_1$scpccvs, out_4$scpccvs, tolerance = 1e-12)
  expect_equal(out_1$scpcstats, out_fb$scpcstats, tolerance = 1e-10)
  expect_equal(out_1$scpccvs, out_fb$scpccvs, tolerance = 1e-10)
})

test_that("SPATIAL_BACKEND and spatial.backend options govern backend selection cleanly", {
  dat <- make_scpc_data(n = 20)
  fit <- stats::lm(y ~ x, data = dat)

  old_opt <- getOption("spatial.backend")
  old_env <- Sys.getenv("SPATIAL_BACKEND", unset = NA_character_)
  on.exit({
    options(spatial.backend = old_opt)
    if (is.na(old_env)) {
      Sys.unsetenv("SPATIAL_BACKEND")
    } else {
      Sys.setenv(SPATIAL_BACKEND = old_env)
    }
  }, add = TRUE)

  options(spatial.backend = NULL)
  Sys.setenv(SPATIAL_BACKEND = "fallback")
  out_env_fb <- scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = 0.1)
  expect_identical(out_env_fb$backend, "fallback")

  if (is.loaded("wrap__rust_run_exact_scpc", PACKAGE = "scpcR")) {
    # Explicit argument overrides env and option
    out_explicit <- scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = 0.1, backend = "rust")
    expect_identical(out_explicit$backend, "rust")

    # Option overrides env
    options(spatial.backend = "rust")
    out_opt <- scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = 0.1)
    expect_identical(out_opt$backend, "rust")

    # Option "fallback" overrides SPATIAL_BACKEND="rust"
    Sys.setenv(SPATIAL_BACKEND = "rust")
    options(spatial.backend = "fallback")
    out_opt_fb <- scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = 0.1)
    expect_identical(out_opt_fb$backend, "fallback")
  }
})

test_that("nthreads resolution honors options, env vars, validation, and _R_CHECK_LIMIT_CORES_", {
  old_opts <- options(scpcR.threads = NULL, mc.cores = NULL)
  old_env_th <- Sys.getenv("SCPCR_NUM_THREADS", unset = NA_character_)
  old_env_lim <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  on.exit({
    options(old_opts)
    if (is.na(old_env_th)) Sys.unsetenv("SCPCR_NUM_THREADS") else Sys.setenv(SCPCR_NUM_THREADS = old_env_th)
    if (is.na(old_env_lim)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_") else Sys.setenv(`_R_CHECK_LIMIT_CORES_` = old_env_lim)
  }, add = TRUE)

  Sys.unsetenv("SCPCR_NUM_THREADS")
  Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
  expect_null(.resolve_scpc_nthreads(NULL))
  expect_identical(.resolve_scpc_nthreads(3), 3L)

  options(mc.cores = 4L)
  expect_identical(.resolve_scpc_nthreads(NULL), 4L)
  Sys.setenv(SCPCR_NUM_THREADS = "3")
  expect_identical(.resolve_scpc_nthreads(NULL), 3L)
  options(scpcR.threads = 5L)
  expect_identical(.resolve_scpc_nthreads(NULL), 5L)

  Sys.setenv(`_R_CHECK_LIMIT_CORES_` = "TRUE")
  expect_identical(.resolve_scpc_nthreads(NULL), 2L)
  expect_identical(.resolve_scpc_nthreads(1L), 1L)
  expect_identical(.resolve_scpc_nthreads(6L), 2L)

  expect_error(.resolve_scpc_nthreads(0), "positive integer")
  expect_error(.resolve_scpc_nthreads(1.5), "positive integer")
  expect_error(.resolve_scpc_nthreads(NA_integer_), "positive integer")

  dat <- make_scpc_data(n = 15)
  fit <- stats::lm(y ~ x, data = dat)
  expect_error(
    scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = 0.1, backend = "bogus"),
    "`backend` must be one of"
  )
  expect_error(
    scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = 0.1, nthreads = 0),
    "positive integer"
  )
})

test_that("unsupported M1 cases fall back in auto and error in explicit rust mode", {
  dat <- make_scpc_data(n = 20)
  dat$lon <- ave(dat$lon, dat$cluster, FUN = mean)
  dat$lat <- ave(dat$lat, dat$cluster, FUN = mean)
  fit <- stats::lm(y ~ x, data = dat)

  out_cl_auto <- scpc(
    fit, dat,
    coords_euclidean = c("lon", "lat"),
    cluster = "cluster",
    avc = 0.1, uncond = TRUE, backend = "auto"
  )
  expect_identical(out_cl_auto$backend, "fallback")

  expect_error(
    scpc(
      fit, dat,
      coords_euclidean = c("lon", "lat"),
      cluster = "cluster",
      avc = 0.1, uncond = TRUE, backend = "rust"
    )
  )

  out_approx_auto <- scpc(
    fit, dat,
    coords_euclidean = c("lon", "lat"),
    avc = 0.1, method = "approx", uncond = TRUE, backend = "auto"
  )
  expect_identical(out_approx_auto$backend, "fallback")
  expect_error(
    scpc(
      fit, dat,
      coords_euclidean = c("lon", "lat"),
      avc = 0.1, method = "approx", uncond = TRUE, backend = "rust"
    )
  )

  if (!is.loaded("wrap__rust_run_exact_scpc", PACKAGE = "scpcR")) {
    expect_error(
      scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = 0.1, backend = "rust"),
      "not loaded"
    )
  }
})

test_that("no-intercept OLS, 1D Euclidean coords, and tied grids behave as expected", {
  skip_if_not(is.loaded("wrap__rust_run_exact_scpc", PACKAGE = "scpcR"), "native Rust library not loaded")
  set.seed(77)
  n <- 40
  d <- data.frame(
    y = stats::rnorm(n),
    x1 = stats::rnorm(n),
    s1 = stats::runif(n),
    s2 = stats::runif(n)
  )
  fit_noint <- stats::lm(y ~ x1 - 1, data = d)
  out_noint_r <- scpc(fit_noint, d, coords_euclidean = c("s1", "s2"), avc = 0.05, uncond = FALSE, cvs = TRUE, backend = "rust")
  out_noint_f <- scpc(fit_noint, d, coords_euclidean = c("s1", "s2"), avc = 0.05, uncond = FALSE, cvs = TRUE, backend = "fallback")
  expect_equal(out_noint_r$scpcstats, out_noint_f$scpcstats, tolerance = 1e-10)
  expect_equal(out_noint_r$scpccvs, out_noint_f$scpccvs, tolerance = 1e-10)

  out_1d_r <- scpc(fit_noint, d, coords_euclidean = "s1", avc = 0.05, uncond = FALSE, backend = "rust")
  out_1d_f <- scpc(fit_noint, d, coords_euclidean = "s1", avc = 0.05, uncond = FALSE, backend = "fallback")
  expect_equal(out_1d_r$scpcstats, out_1d_f$scpcstats, tolerance = 1e-10)

  # Tied regular lattice: auto falls back cleanly, explicit rust stays on rust
  grid <- expand.grid(s1 = 1:4, s2 = 1:4)
  grid$x <- seq_len(nrow(grid)) * 0.25 + stats::rnorm(nrow(grid), sd = 0.05)
  grid$y <- 1 + 0.8 * grid$x + stats::rnorm(nrow(grid), sd = 0.1)
  fit_grid <- stats::lm(y ~ x, data = grid)

  out_grid_auto <- scpc(fit_grid, grid, coords_euclidean = c("s1", "s2"), avc = 0.1, backend = "auto")
  out_grid_rust <- scpc(fit_grid, grid, coords_euclidean = c("s1", "s2"), avc = 0.1, backend = "rust")
  expect_identical(out_grid_auto$backend, "fallback")
  expect_identical(out_grid_rust$backend, "rust")
  expect_true(isTRUE(out_grid_auto$tied))
  expect_true(isTRUE(out_grid_rust$tied))
})

test_that("scale-invariant calibration and co-location floor match defining oracle", {
  bench_200_path <- find_spur_file("verification/data/bench_200.csv")
  skip_if(is.na(bench_200_path), "bench_200.csv not available")

  d30 <- utils::read.csv(bench_200_path)[seq_len(30L), , drop = FALSE]
  fit30 <- stats::lm(y ~ x1 + x2, data = d30)

  backends <- "fallback"
  if (is.loaded("wrap__rust_run_exact_scpc", PACKAGE = "scpcR")) {
    backends <- c("rust", "fallback")
  }

  scales <- c(1, 1e-6, 1e6)
  for (bk in backends) {
    outs_03 <- list()
    for (i in seq_along(scales)) {
      sc <- scales[[i]]
      d_sc <- d30
      d_sc$s_1 <- d30$s_1 * sc
      d_sc$s_2 <- d30$s_2 * sc
      out <- scpc(
        fit30, d_sc,
        coords_euclidean = c("s_1", "s_2"),
        avc = 0.03, uncond = FALSE, cvs = TRUE, backend = bk
      )
      distv <- .lvech(.getdistmat(d_sc[, c("s_1", "s_2"), drop = FALSE], latlong = FALSE))
      expect_lte(abs(mean(exp(-out$c0 * distv)) - 0.03), 1e-12)
      expect_true(isTRUE(out$tied))
      outs_03[[i]] <- out
    }
    expect_equal(outs_03[[2L]]$c0, 6099987.97, tolerance = 1e-6)
    for (i in c(2L, 3L)) {
      sc <- scales[[i]]
      expect_equal(outs_03[[i]]$c0 * sc, outs_03[[1L]]$c0, tolerance = 1e-10)
      expect_equal(outs_03[[i]]$cv, outs_03[[1L]]$cv, tolerance = 1e-10)
      expect_equal(outs_03[[i]]$q, outs_03[[1L]]$q)
    }
  }

  # Untied geometry (d30 at avc = 0.1): full inference and Rust/fallback parity are scale-invariant
  results_by_backend <- list()
  for (bk in backends) {
    outs <- list()
    for (i in seq_along(scales)) {
      sc <- scales[[i]]
      d_sc <- d30
      d_sc$s_1 <- d30$s_1 * sc
      d_sc$s_2 <- d30$s_2 * sc
      out <- scpc(
        fit30, d_sc,
        coords_euclidean = c("s_1", "s_2"),
        avc = 0.1, uncond = FALSE, cvs = TRUE, backend = bk
      )
      distv <- .lvech(.getdistmat(d_sc[, c("s_1", "s_2"), drop = FALSE], latlong = FALSE))
      expect_lte(abs(mean(exp(-out$c0 * distv)) - 0.1), 1e-12)
      expect_false(isTRUE(out$tied))
      outs[[i]] <- out
    }
    base_out <- outs[[1L]]
    for (i in c(2L, 3L)) {
      sc <- scales[[i]]
      expect_equal(outs[[i]]$c0 * sc, base_out$c0, tolerance = 1e-10)
      expect_equal(outs[[i]]$cv, base_out$cv, tolerance = 1e-10)
      expect_equal(outs[[i]]$q, base_out$q)
      expect_equal(outs[[i]]$scpcstats, base_out$scpcstats, tolerance = 1e-10)
      expect_equal(outs[[i]]$scpccvs, base_out$scpccvs, tolerance = 1e-10)
    }
    results_by_backend[[bk]] <- base_out
  }
  if (length(results_by_backend) == 2L) {
    expect_equal(results_by_backend$rust$scpcstats, results_by_backend$fallback$scpcstats, tolerance = 1e-10)
    expect_equal(results_by_backend$rust$scpccvs, results_by_backend$fallback$scpccvs, tolerance = 1e-10)
    expect_equal(results_by_backend$rust$c0, results_by_backend$fallback$c0, tolerance = 1e-10)
  }

  # Co-location floor: 1 zero pair out of 15 pairs -> floor = 1/15
  d_coloc <- data.frame(
    y = c(1.0, 2.0, 1.5, 3.2, 4.1, 2.7),
    x = c(0, 1, 2, 3, 5, 8),
    s1 = c(0, 0, 1, 3, 7, 15),
    s2 = rep(0, 6)
  )
  fit_coloc <- stats::lm(y ~ x, data = d_coloc)
  distv_coloc <- .lvech(.getdistmat(d_coloc[, c("s1", "s2"), drop = FALSE], latlong = FALSE))

  coloc_outs <- list()
  for (bk in backends) {
    out_ok <- scpc(
      fit_coloc, d_coloc,
      coords_euclidean = c("s1", "s2"),
      avc = 0.25, uncond = TRUE, backend = bk
    )
    expect_lte(abs(mean(exp(-out_ok$c0 * distv_coloc)) - 0.25), 1e-12)
    coloc_outs[[bk]] <- out_ok

    expect_error(
      scpc(fit_coloc, d_coloc, coords_euclidean = c("s1", "s2"), avc = 0.05, uncond = TRUE, backend = bk)
    )

    d_ident <- d_coloc
    d_ident$s1 <- rep(1, 6)
    expect_error(
      scpc(fit_coloc, d_ident, coords_euclidean = c("s1", "s2"), avc = 0.5, uncond = TRUE, backend = bk)
    )
  }
  if (length(coloc_outs) == 2L) {
    expect_equal(coloc_outs$rust$scpcstats, coloc_outs$fallback$scpcstats, tolerance = 1e-10)
  }
})

test_that("near-zero t-statistic matches analytic two-point Cauchy tail probability", {
  d_near0 <- data.frame(
    y = c(1 + 1e-6, -1 + 1e-6),
    s = c(0, 1)
  )
  fit_near0 <- stats::lm(y ~ 1, data = d_near0)
  backends <- "fallback"
  if (is.loaded("wrap__rust_run_exact_scpc", PACKAGE = "scpcR")) {
    backends <- c("rust", "fallback")
  }
  for (bk in backends) {
    out <- scpc(
      fit_near0, d_near0,
      coords_euclidean = "s",
      avc = 0.03, uncond = TRUE, backend = bk
    )
    rho <- exp(-out$c0)
    expected_p <- 1 - (2 / pi) * atan(abs(out$scpcstats[1L, "t"]) * sqrt((1 - rho) / (1 + rho)))
    expect_equal(out$q, 1)
    expect_lte(abs(abs(unname(out$scpcstats[1L, "t"])) - 1e-6), 1e-12)
    expect_equal(unname(out$scpcstats[1L, "P>|t|"]), expected_p, tolerance = 1e-10)
    expect_equal(unname(out$scpcstats[1L, "P>|t|"]), 0.99999938, tolerance = 1e-7)
  }
})

test_that("small-unit regressor rank and positive outcome rescaling preserve inference", {
  set.seed(2026)
  p <- sort(stats::runif(50))
  s <- stats::runif(50)
  x_base <- p - 0.5
  x_tiny <- 1e-11 * x_base
  y <- 0.3 + p + stats::rnorm(50, sd = 0.5)
  d_rank <- data.frame(y = y, x = x_tiny, x_base = x_base, p = p, s = s)

  fit_scaled <- stats::lm(y ~ x, data = d_rank)
  fit_base <- stats::lm(y ~ x_base, data = d_rank)

  backends <- "fallback"
  if (is.loaded("wrap__rust_run_exact_scpc", PACKAGE = "scpcR")) {
    backends <- c("rust", "fallback")
  }

  outs_scaled <- list()
  outs_base <- list()
  for (bk in backends) {
    out_sc <- scpc(
      fit_scaled, d_rank,
      coords_euclidean = c("p", "s"),
      avc = 0.03, uncond = FALSE, cvs = TRUE, backend = bk
    )
    out_bs <- scpc(
      fit_base, d_rank,
      coords_euclidean = c("p", "s"),
      avc = 0.03, uncond = FALSE, cvs = TRUE, backend = bk
    )
    expect_equal(out_sc$q, out_bs$q)
    expect_equal(out_sc$cv, out_bs$cv, tolerance = 1e-10)
    expect_equal(unname(out_sc$scpcstats[2L, c("t", "P>|t|")]), unname(out_bs$scpcstats[2L, c("t", "P>|t|")]), tolerance = 1e-10)
    expect_equal(unname(out_sc$scpccvs[2L, ]), unname(out_bs$scpccvs[2L, ]), tolerance = 1e-10)
    outs_scaled[[bk]] <- out_sc
    outs_base[[bk]] <- out_bs
  }
  if (length(outs_scaled) == 2L) {
    expect_equal(outs_scaled$rust$scpcstats, outs_scaled$fallback$scpcstats, tolerance = 1e-10)
    expect_equal(outs_scaled$rust$scpccvs, outs_scaled$fallback$scpccvs, tolerance = 1e-10)
  }

  # Positive outcome rescaling: q, t, P>|t|, cv invariant; Coef, SE, CI proportional
  bench_200_path <- find_spur_file("verification/data/bench_200.csv")
  skip_if(is.na(bench_200_path), "bench_200.csv not available")
  d50 <- utils::read.csv(bench_200_path)[seq_len(50L), , drop = FALSE]
  fit_orig <- stats::lm(y ~ x1 + x2, data = d50)
  fit_tiny <- stats::lm(I(y * 1e-12) ~ x1 + x2, data = d50)
  fit_huge <- stats::lm(I(y * 1e6) ~ x1 + x2, data = d50)

  for (bk in backends) {
    o_orig <- scpc(fit_orig, d50, coords_euclidean = c("s_1", "s_2"), avc = 0.03, uncond = FALSE, cvs = TRUE, backend = bk)
    o_tiny <- scpc(fit_tiny, d50, coords_euclidean = c("s_1", "s_2"), avc = 0.03, uncond = FALSE, cvs = TRUE, backend = bk)
    o_huge <- scpc(fit_huge, d50, coords_euclidean = c("s_1", "s_2"), avc = 0.03, uncond = FALSE, cvs = TRUE, backend = bk)

    for (pair in list(list(out = o_tiny, scale = 1e-12), list(out = o_huge, scale = 1e6))) {
      o_s <- pair$out
      sc <- pair$scale
      expect_equal(o_s$q, o_orig$q)
      expect_equal(o_s$c0, o_orig$c0, tolerance = 1e-10)
      expect_equal(o_s$cv, o_orig$cv, tolerance = 1e-10)
      expect_equal(o_s$scpccvs, o_orig$scpccvs, tolerance = 1e-10)
      expect_equal(o_s$scpcstats[, c("t", "P>|t|")], o_orig$scpcstats[, c("t", "P>|t|")], tolerance = 1e-10)
      expect_equal(o_s$scpcstats[, c("Coef", "Std_Err", "2.5 %", "97.5 %")] / sc,
                   o_orig$scpcstats[, c("Coef", "Std_Err", "2.5 %", "97.5 %")],
                   tolerance = 1e-10)
    }
  }
})

test_that("nthreads = 1 stays fork-safe in parallel::mcparallel after warming threads", {
  skip_if_not(is.loaded("wrap__rust_run_exact_scpc", PACKAGE = "scpcR"), "native Rust library not loaded")
  skip_on_os("windows")
  skip_if_not_installed("parallel")

  bench_200_path <- find_spur_file("verification/data/bench_200.csv")
  skip_if(is.na(bench_200_path), "bench_200.csv not available")
  d60 <- utils::read.csv(bench_200_path)[seq_len(60L), , drop = FALSE]
  fit60 <- stats::lm(y ~ x1 + x2, data = d60)

  warm_1 <- scpc(
    fit60, d60,
    coords_euclidean = c("s_1", "s_2"),
    avc = 0.03, uncond = FALSE, cvs = TRUE,
    backend = "rust", nthreads = 1L
  )
  warm_4 <- scpc(
    fit60, d60,
    coords_euclidean = c("s_1", "s_2"),
    avc = 0.03, uncond = FALSE, cvs = TRUE,
    backend = "rust", nthreads = 4L
  )

  job <- parallel::mcparallel({
    scpc(
      fit60, d60,
      coords_euclidean = c("s_1", "s_2"),
      avc = 0.03, uncond = FALSE, cvs = TRUE,
      backend = "rust", nthreads = 1L
    )
  })
  collected <- NULL
  deadline <- proc.time()[["elapsed"]] + 15
  while (proc.time()[["elapsed"]] < deadline) {
    collected <- parallel::mccollect(job, wait = FALSE, timeout = 0.2)
    if (!is.null(collected)) break
  }
  if (is.null(collected)) {
    tools::pskill(job$pid, tools::SIGKILL)
    parallel::mccollect(job, wait = FALSE)
    fail("mcparallel child timed out after warming Rust thread pools")
  }
  child_out <- collected[[1L]]
  expect_s3_class(child_out, "scpc")
  expect_identical(child_out$backend, "rust")
  expect_equal(child_out$scpcstats, warm_1$scpcstats, tolerance = 1e-12)
  expect_equal(child_out$scpccvs, warm_4$scpccvs, tolerance = 1e-12)
})

test_that("ineligible models route to fallback in auto and error in rust regardless of uncond", {
  dat <- make_scpc_data(n = 25)
  dat$w <- seq(1, 2, length.out = nrow(dat))
  dat$fe <- rep(seq_len(5L), each = 5L)

  fit_weighted <- stats::lm(y ~ x, data = dat, weights = w)
  fit_glm <- stats::glm(y ~ x, data = dat)

  for (uncond_flag in c(TRUE, FALSE)) {
    out_w_auto <- scpc(
      fit_weighted, dat,
      coords_euclidean = c("lon", "lat"),
      avc = 0.1, uncond = uncond_flag, backend = "auto"
    )
    expect_identical(out_w_auto$backend, "fallback")
    expect_error(
      scpc(
        fit_weighted, dat,
        coords_euclidean = c("lon", "lat"),
        avc = 0.1, uncond = uncond_flag, backend = "rust"
      )
    )

    out_glm_auto <- scpc(
      fit_glm, dat,
      coords_euclidean = c("lon", "lat"),
      avc = 0.1, uncond = uncond_flag, backend = "auto"
    )
    expect_identical(out_glm_auto$backend, "fallback")
    expect_error(
      scpc(
        fit_glm, dat,
        coords_euclidean = c("lon", "lat"),
        avc = 0.1, uncond = uncond_flag, backend = "rust"
      )
    )
  }

  if (requireNamespace("fixest", quietly = TRUE)) {
    with_fixest_single_thread({
      fit_fe <- fixest::feols(y ~ x | fe, data = dat)
      for (uncond_flag in c(TRUE, FALSE)) {
        out_fe_auto <- scpc(
          fit_fe, dat,
          coords_euclidean = c("lon", "lat"),
          avc = 0.1, uncond = uncond_flag, backend = "auto"
        )
        expect_identical(out_fe_auto$backend, "fallback")
        expect_error(
          scpc(
            fit_fe, dat,
            coords_euclidean = c("lon", "lat"),
            avc = 0.1, uncond = uncond_flag, backend = "rust"
          )
        )
      }

      if (is.loaded("wrap__rust_run_exact_scpc", PACKAGE = "scpcR")) {
        fit_fe_plain <- fixest::feols(y ~ x, data = dat)
        out_plain_auto <- scpc(
          fit_fe_plain, dat,
          coords_euclidean = c("lon", "lat"),
          avc = 0.1, uncond = FALSE, backend = "auto"
        )
        out_plain_fb <- scpc(
          fit_fe_plain, dat,
          coords_euclidean = c("lon", "lat"),
          avc = 0.1, uncond = FALSE, backend = "fallback"
        )
        expect_identical(out_plain_auto$backend, "rust")
        expect_equal(out_plain_auto$scpcstats, out_plain_fb$scpcstats, tolerance = 1e-10)
      }
    })
  }
})

test_that("forced fallback and ineligible auto calls ignore invalid ambient thread settings", {
  old_opts <- options(scpcR.threads = "bad_option", mc.cores = 0L)
  old_env <- Sys.getenv("SCPCR_NUM_THREADS", unset = NA_character_)
  Sys.setenv(SCPCR_NUM_THREADS = "not_an_int")
  on.exit({
    options(old_opts)
    if (is.na(old_env)) Sys.unsetenv("SCPCR_NUM_THREADS") else Sys.setenv(SCPCR_NUM_THREADS = old_env)
  }, add = TRUE)

  dat <- make_scpc_data(n = 20)
  dat$lon <- ave(dat$lon, dat$cluster, FUN = mean)
  dat$lat <- ave(dat$lat, dat$cluster, FUN = mean)
  fit <- stats::lm(y ~ x, data = dat)

  out_fb <- scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = 0.1, backend = "fallback")
  expect_identical(out_fb$backend, "fallback")

  out_cl <- scpc(fit, dat, coords_euclidean = c("lon", "lat"), cluster = "cluster", avc = 0.1, uncond = TRUE, backend = "auto")
  expect_identical(out_cl$backend, "fallback")

  expect_error(
    scpc(fit, dat, coords_euclidean = c("lon", "lat"), avc = 0.1, backend = "fallback", nthreads = 0),
    "positive integer"
  )
})
