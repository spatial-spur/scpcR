test_that("scpc handles fixest IV models without fixed effects", {
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    set.seed(2001)
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
      coord_x = runif(n),
      coord_y = runif(n)
    )

    fit <- fixest::feols(y ~ w | x ~ z, data = dat)
    out <- expect_no_warning(
      scpc(
        fit,
        dat,
        coords_euclidean = c("coord_x", "coord_y"),
        ncoef = 2,
        avc = 0.1,
        uncond = FALSE,
        cvs = TRUE
      )
    )

    expect_s3_class(out, "scpc")
    expect_true(all(is.finite(out$scpcstats)))
    expect_true(all(is.finite(out$scpccvs)))
  })
})

test_that("scpc handles fixest IV models with one-way absorbed fixed effects", {
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    set.seed(2002)
    n_fe <- 20
    t_per_fe <- 5
    n <- n_fe * t_per_fe
    fe <- rep(seq_len(n_fe), each = t_per_fe)
    z <- stats::rnorm(n)
    w <- stats::rnorm(n)
    u <- stats::rnorm(n)
    x <- 0.7 * z + 0.3 * w + u
    y <- 1 + 1.1 * x + 0.4 * w + stats::rnorm(n_fe)[fe] + u
    dat <- data.frame(
      y = y,
      x = x,
      w = w,
      z = z,
      fe = fe,
      coord_x = runif(n),
      coord_y = runif(n)
    )

    fit <- fixest::feols(y ~ w | fe | x ~ z, data = dat)
    out <- expect_no_warning(
      scpc(
        fit,
        dat,
        coords_euclidean = c("coord_x", "coord_y"),
        ncoef = 2,
        avc = 0.1,
        uncond = FALSE,
        cvs = TRUE
      )
    )

    expect_s3_class(out, "scpc")
    expect_true(all(is.finite(out$scpcstats)))
    expect_true(all(is.finite(out$scpccvs)))
  })
})

test_that("scpc handles fixest IV models with two-way absorbed fixed effects", {
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    set.seed(2005)
    n1 <- 5
    n2 <- 4
    reps <- 4
    n <- n1 * n2 * reps
    fe1 <- rep(rep(seq_len(n1), each = n2), each = reps)
    fe2 <- rep(rep(seq_len(n2), times = n1), each = reps)
    z <- stats::rnorm(n)
    w <- stats::rnorm(n)
    u <- stats::rnorm(n)
    x <- 0.7 * z + 0.2 * w + u
    y <- 1 + 1.1 * x + 0.35 * w + stats::rnorm(n1)[fe1] + stats::rnorm(n2)[fe2] + u
    dat <- data.frame(
      y = y,
      x = x,
      w = w,
      z = z,
      fe1 = fe1,
      fe2 = fe2,
      coord_x = runif(n),
      coord_y = runif(n)
    )

    fit <- fixest::feols(y ~ w | fe1 + fe2 | x ~ z, data = dat)
    out <- expect_no_warning(
      scpc(
        fit,
        dat,
        coords_euclidean = c("coord_x", "coord_y"),
        ncoef = 2,
        avc = 0.1,
        uncond = FALSE,
        cvs = TRUE
      )
    )

    expect_s3_class(out, "scpc")
    expect_true(all(is.finite(out$scpcstats)))
    expect_true(all(is.finite(out$scpccvs)))
  })
})

test_that("scpc handles clustered fixest IV models", {
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    set.seed(2003)
    n_cluster <- 40
    cl_size <- 3
    n <- n_cluster * cl_size
    cl <- rep(seq_len(n_cluster), each = cl_size)
    z <- stats::rnorm(n)
    w <- stats::rnorm(n)
    u <- stats::rnorm(n)
    x <- 0.8 * z + 0.2 * w + u
    y <- 1 + 1.15 * x + 0.35 * w + u
    lon_cl <- runif(n_cluster)
    lat_cl <- runif(n_cluster)
    dat <- data.frame(
      y = y,
      x = x,
      w = w,
      z = z,
      cl = cl,
      lon = lon_cl[cl],
      lat = lat_cl[cl]
    )

    fit <- fixest::feols(y ~ w | x ~ z, data = dat)
    out <- expect_no_warning(
      scpc(
        fit,
        dat,
        coords_euclidean = c("lon", "lat"),
        cluster = "cl",
        ncoef = 2,
        avc = 0.1,
        uncond = FALSE,
        cvs = TRUE
      )
    )

    expect_s3_class(out, "scpc")
    expect_true(all(is.finite(out$scpcstats)))
    expect_true(all(is.finite(out$scpccvs)))
  })
})

test_that("scpc handles clustered fixest IV models with absorbed fixed effects", {
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    set.seed(2004)
    n_fe <- 15
    t_per_fe <- 6
    n <- n_fe * t_per_fe
    fe <- rep(seq_len(n_fe), each = t_per_fe)
    cl <- rep(seq_len(n / 2), each = 2)
    z <- stats::rnorm(n)
    w <- stats::rnorm(n)
    u <- stats::rnorm(n)
    x <- 0.75 * z + 0.25 * w + u
    y <- 1 + 1.05 * x + 0.45 * w + stats::rnorm(n_fe)[fe] + u
    lon_cl <- runif(length(unique(cl)))
    lat_cl <- runif(length(unique(cl)))
    dat <- data.frame(
      y = y,
      x = x,
      w = w,
      z = z,
      fe = fe,
      cl = cl,
      lon = lon_cl[cl],
      lat = lat_cl[cl]
    )

    fit <- fixest::feols(y ~ w | fe | x ~ z, data = dat)
    out <- expect_no_warning(
      scpc(
        fit,
        dat,
        coords_euclidean = c("lon", "lat"),
        cluster = "cl",
        ncoef = 2,
        avc = 0.1,
        uncond = FALSE,
        cvs = TRUE
      )
    )

    expect_s3_class(out, "scpc")
    expect_true(all(is.finite(out$scpcstats)))
    expect_true(all(is.finite(out$scpccvs)))
  })
})

test_that("scpc handles clustered fixest IV models in the large-n approximation path", {
  skip_if_not_installed("fixest")

  with_fixest_single_thread({
    set.seed(2006)
    n_cluster <- 40
    cl_size <- 3
    n <- n_cluster * cl_size
    cl <- rep(seq_len(n_cluster), each = cl_size)
    z <- stats::rnorm(n)
    w <- stats::rnorm(n)
    u <- stats::rnorm(n)
    x <- 0.8 * z + 0.2 * w + u
    y <- 1 + 1.15 * x + 0.35 * w + u
    lon_cl <- runif(n_cluster)
    lat_cl <- runif(n_cluster)
    dat <- data.frame(
      y = y,
      x = x,
      w = w,
      z = z,
      cl = cl,
      lon = lon_cl[cl],
      lat = lat_cl[cl]
    )

    fit <- fixest::feols(y ~ w | x ~ z, data = dat)
    out <- expect_no_warning(
      scpc(
        fit,
        dat,
        coords_euclidean = c("lon", "lat"),
        cluster = "cl",
        ncoef = 2,
        avc = 0.1,
        method = "approx",
        large_n_seed = 1,
        uncond = FALSE,
        cvs = TRUE
      )
    )

    expect_s3_class(out, "scpc")
    expect_identical(out$method, "approx")
    expect_equal(out$large_n_seed, 1)
    expect_true(all(is.finite(out$scpcstats)))
    expect_true(all(is.finite(out$scpccvs)))
  })
})
