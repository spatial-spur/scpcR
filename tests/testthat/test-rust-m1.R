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

assert_scpc_matches_fixture <- function(out, fixture, tol = 1e-10) {
  expected_stats <- json_rows_to_matrix(fixture$scpcstats)
  rownames(expected_stats) <- unlist(fixture$coef_names)
  colnames(expected_stats) <- c("Coef", "Std_Err", "t", "P>|t|", "2.5 %", "97.5 %")

  expect_equal(unname(out$scpcstats), unname(expected_stats), tolerance = tol)
  expect_equal(rownames(out$scpcstats), rownames(expected_stats))
  expect_equal(out$c0, fixture$c0, tolerance = tol)
  expect_equal(out$cv, fixture$cv, tolerance = tol)
  expect_equal(out$q, fixture$q)

  if (is.null(fixture$scpccvs)) {
    expect_null(out$scpccvs)
  } else {
    expected_cvs <- json_rows_to_matrix(fixture$scpccvs)
    expect_equal(unname(out$scpccvs), unname(expected_cvs), tolerance = tol)
  }
}

test_that("rust and fallback backends match m1_exact_scpc.json within 1e-10", {
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

  for (bk in backends) {
    for (cs in cases) {
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
      assert_scpc_matches_fixture(out, fixtures[[cs$name]], tol = 1e-10)
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
})
