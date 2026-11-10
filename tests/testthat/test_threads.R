test_that("thread control helpers set and retrieve thread configuration correctly", {
  old_th <- scpc_get_num_threads()
  on.exit({
    scpc_set_num_threads("auto")
    options(scpcR.threads = NULL)
    options(mc.cores = NULL)
    Sys.unsetenv("SCPCR_NUM_THREADS")
    Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
  }, add = TRUE)

  # 1. Setting specific thread counts via helper
  scpc_set_num_threads(2)
  expect_equal(scpc_get_num_threads(), 2L)
  expect_equal(getOption("scpcR.threads"), 2L)

  # 2. Restoring to auto
  scpc_set_num_threads("auto")
  expect_null(getOption("scpcR.threads"))

  # 3. Direct options() control
  options(scpcR.threads = 1L)
  expect_equal(scpc_get_num_threads(), 1L)

  # 4. mc.cores fallback
  options(scpcR.threads = NULL)
  options(mc.cores = 3L)
  expect_equal(scpc_get_num_threads(), 3L)

  # 5. CRAN _R_CHECK_LIMIT_CORES_ policy cap at 2
  options(scpcR.threads = NULL)
  options(mc.cores = NULL)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
  expect_equal(scpc_get_num_threads(), 2L)
  Sys.unsetenv("_R_CHECK_LIMIT_CORES_")

  # 6. Invalid inputs throw error
  expect_error(scpc_set_num_threads(0))
  expect_error(scpc_set_num_threads(-2))
  expect_error(scpc_set_num_threads("many"))
})

test_that("scpc(..., nthreads = 1) parameter and single-threaded execution match multi-threaded execution", {
  on.exit(scpc_set_num_threads("auto"), add = TRUE)

  set.seed(42)
  n <- 50
  dat <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    lon = runif(n, -100, -80),
    lat = runif(n, 30, 45)
  )
  fit <- lm(y ~ x, data = dat)

  # Run multi-threaded via auto
  scpc_set_num_threads("auto")
  res_multi <- scpc(fit, data = dat, lon = "lon", lat = "lat", avc = 0.05, cvs = TRUE)

  # Run via nthreads = 1 argument
  res_arg1 <- scpc(fit, data = dat, lon = "lon", lat = "lat", avc = 0.05, cvs = TRUE, nthreads = 1)

  # Run single-threaded via setter
  scpc_set_num_threads(1)
  res_single <- scpc(fit, data = dat, lon = "lon", lat = "lat", avc = 0.05, cvs = TRUE)

  expect_equal(res_arg1$scpcstats, res_multi$scpcstats, tolerance = 1e-14)
  expect_equal(res_arg1$scpccvs, res_multi$scpccvs, tolerance = 1e-14)
  expect_equal(res_single$scpcstats, res_multi$scpcstats, tolerance = 1e-14)
  expect_equal(res_single$scpccvs, res_multi$scpccvs, tolerance = 1e-14)
  expect_equal(res_single$cv, res_multi$cv, tolerance = 1e-14)
})
