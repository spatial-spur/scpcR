test_that("native C++ matrix and spatial helpers match pure-R reference implementations", {
  set.seed(42)
  n <- 65
  coords_euc <- cbind(runif(n, -2, 2), runif(n, -2, 2))
  coords_geo <- cbind(runif(n, -120, -70), runif(n, 25, 48))

  # 1. .lvech & .demeanmat
  M_sym <- crossprod(matrix(rnorm(n * n), n, n))
  ref_lvech <- M_sym[lower.tri(M_sym, diag = FALSE)]
  expect_equal(.lvech(M_sym), ref_lvech, tolerance = 1e-15)

  M_rect <- matrix(rnorm(n * 17), n, 17)
  r_dem <- M_rect - rowMeans(M_rect)
  ref_demean <- sweep(r_dem, 2, colMeans(r_dem))
  expect_equal(.demeanmat(M_rect), ref_demean, tolerance = 1e-14)

  # 2. .getdistmat & .getdistvec (Euclidean and Geodesic)
  D_euc <- .getdistmat(coords_euc, latlong = FALSE)
  ref_D_euc <- as.matrix(stats::dist(coords_euc, diag = TRUE, upper = TRUE))
  expect_equal(D_euc, unname(ref_D_euc), tolerance = 1e-14)

  D_geo <- .getdistmat(coords_geo, latlong = TRUE)
  lon_rad <- coords_geo[, 1] * pi / 180
  lat_rad <- coords_geo[, 2] * pi / 180
  ref_D_geo <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      dlon <- 0.5 * (lon_rad[i] - lon_rad[j])
      dlat <- 0.5 * (lat_rad[i] - lat_rad[j])
      ref_D_geo[i, j] <- asin(sqrt(sin(dlat)^2 + cos(lat_rad[i]) * cos(lat_rad[j]) * sin(dlon)^2)) / pi
    }
  }
  expect_equal(D_geo, ref_D_geo, tolerance = 1e-14)

  S1 <- coords_geo[1:30, , drop = FALSE]
  S2 <- coords_geo[31:60, , drop = FALSE]
  ref_dv_geo <- asin(sqrt(
    sin(0.5 * (S1[, 2] - S2[, 2]) * pi / 180)^2 +
      cos(S1[, 2] * pi / 180) * cos(S2[, 2] * pi / 180) * sin(0.5 * (S1[, 1] - S2[, 1]) * pi / 180)^2
  )) / pi
  expect_equal(.getdistvec(S1, S2, latlong = TRUE), ref_dv_geo, tolerance = 1e-14)
  expect_equal(.getdistvec(S1, S2, latlong = FALSE), sqrt(rowSums((S1 - S2)^2)), tolerance = 1e-14)

  # 3. .getavc & .getc0fromavc (including histogram convexity path for m >= 4096)
  coords_med <- cbind(runif(110), runif(110))
  D_med <- .getdistmat(coords_med, latlong = FALSE)
  dv_med <- .lvech(D_med) # length = 5995 >= 4096
  expect_equal(.getavc(3.5, dv_med), mean(exp(-3.5 * dv_med)), tolerance = 1e-14)

  ref_getc0 <- function(dist, avc0) {
    c0 <- c1 <- 10
    while (mean(exp(-c0 * dist)) < avc0) {
      c1 <- c0
      c0 <- 0.5 * c0
    }
    while (mean(exp(-c1 * dist)) > avc0 & c1 < 5000) {
      c0 <- c1
      c1 <- 2 * c1
    }
    repeat {
      c_val <- sqrt(c0 * c1)
      if (mean(exp(-c_val * dist)) > avc0) c0 <- c_val else c1 <- c_val
      if (c1 - c0 < 0.001) break
    }
    c_val
  }
  c0_ref <- ref_getc0(dv_med, 0.03)
  cmax_ref <- ref_getc0(dv_med, 0.00001)
  expect_equal(.getc0fromavc(dv_med, 0.03), c0_ref, tolerance = 1e-12)
  expect_equal(.getc0fromavc(dv_med, 0.00001), cmax_ref, tolerance = 1e-12)
  c_pair <- as.numeric(rcpp_getc0_from_distmat(D_med, c(0.03, 0.00001)))
  expect_equal(c_pair, c(c0_ref, cmax_ref), tolerance = 1e-12)

  # 4. .getOms, rcpp_get_oms_batch, .getrp, .maxrp, .getcv, .setfinalW
  qmax <- 8L
  W <- .getW(D_euc, 5.0, qmax)
  Oms <- .getOms(D_euc, 5.0, 25.0, W, 1.2)
  nc <- .getnc(5.0, 25.0, 1.2)
  ref_Oms <- vector("list", nc)
  c_cur <- 5.0
  for (i in seq_len(nc)) {
    ref_Oms[[i]] <- crossprod(W, exp(-c_cur * D_euc) %*% W)
    c_cur <- c_cur * 1.2
  }
  expect_equal(Oms, ref_Oms, tolerance = 1e-12)

  W2 <- .getW(D_euc, 7.0, qmax)
  batch_Oms <- rcpp_get_oms_batch(D_euc, 5.0, 25.0, list(W, W2), 1.2)
  expect_equal(batch_Oms[[1]], Oms, tolerance = 1e-12)
  expect_equal(batch_Oms[[2]], .getOms(D_euc, 5.0, 25.0, W2, 1.2), tolerance = 1e-12)

  # Vectorized .getcv vs scalar .getcv
  levs <- c(0.32, 0.10, 0.05, 0.01)
  cvs_vec <- .getcv(Oms, 5L, levs)
  cvs_scalar <- vapply(levs, function(lv) .getcv(Oms, 5L, lv), 0.0)
  expect_equal(cvs_vec, cvs_scalar, tolerance = 1e-14)
})

test_that("native C++ large-n helpers match pure-R reference on both small and large n", {
  # Test LCG helpers .raninds and .jumble_s
  ref_next_u <- function(rt) {
    rt <- (64389 * rt + 1) %% 2^32
    list(value = rt / 2^32, state = rt)
  }
  ref_raninds <- function(n, capM, random_t) {
    v <- numeric(capM + 1L)
    nxt <- ref_next_u(random_t)
    random_t <- nxt$state
    j <- floor(n * nxt$value)
    for (i in seq_len(capM + 1L)) {
      v[i] <- j + 1L
      nxt <- ref_next_u(random_t)
      random_t <- nxt$state
      j <- (j + 1 + floor(nxt$value * (n - 1))) %% n
    }
    list(indices = as.integer(v), state = random_t)
  }

  ri_cpp <- .raninds(40L, 5000L, 1)
  ri_ref <- ref_raninds(40L, 5000L, 1)
  expect_identical(ri_cpp$indices, ri_ref$indices)
  expect_equal(ri_cpp$state, ri_ref$state)

  # Test .lnget_Oms on n <= 1200 (pair-table path) and n > 1200 (streaming tile path)
  for (n_test in c(50L, 1300L)) {
    set.seed(100 + n_test)
    S <- cbind(runif(n_test), runif(n_test))
    W <- qr.Q(qr(matrix(rnorm(n_test * 6), n_test, 6)))
    capM <- 20000L
    res_cpp <- .lnget_Oms(S, 4.0, 12.0, W, 1.2, capM = capM, random_t = 7, latlong = FALSE)

    inds_res <- .raninds(n_test, capM, 7)
    inds <- inds_res$indices
    dist <- .getdistvec(
      S[inds[seq_len(capM)], , drop = FALSE],
      S[inds[2:(capM + 1L)], , drop = FALSE],
      latlong = FALSE
    )
    W1 <- W[inds[seq_len(capM)], , drop = FALSE]
    W2 <- W[inds[2:(capM + 1L)], , drop = FALSE]
    nc <- .getnc(4.0, 12.0, 1.2)
    ref_Oms <- vector("list", nc)
    c_cur <- 4.0
    for (i in seq_len(nc)) {
      cd <- exp(-c_cur * dist)
      ref_Oms[[i]] <- diag(ncol(W)) +
        0.5 * (n_test * (n_test - 1) / capM) *
        (crossprod(W1, W2 * cd) + crossprod(W2, W1 * cd))
      c_cur <- c_cur * 1.2
    }
    expect_equal(res_cpp$state, inds_res$state)
    expect_equal(res_cpp$Oms, ref_Oms, tolerance = 1e-11)
  }
})
