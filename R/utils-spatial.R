# ---------------------------------------------------------------------------
# Spatial helpers
# ---------------------------------------------------------------------------

.scpc_spatial_cache <- new.env(parent = emptyenv())

.getavc <- function(c, dist) {
  ## Average pairwise correlation for exponential kernel parameter c
  storage.mode(dist) <- "double"
  .Call("cpp_getavc", as.double(c), dist, PACKAGE = "scpcR")
}

.getdistmat <- function(S, latlong) {
  ## Distance matrix from coordinates.
  ## latlong = TRUE: haversine distances as fractions of Earth's circumference.
  ## latlong = FALSE: Euclidean distances.
  S <- as.matrix(S)
  if (latlong && ncol(S) != 2L) {
    stop("Internal error: geodesic coordinates must have exactly two columns.")
  }
  storage.mode(S) <- "double"
  .Call("cpp_distmat", S, as.logical(latlong), PACKAGE = "scpcR")
}

.getdistvec <- function(S1, S2, latlong) {
  ## Vector of paired distances between two coordinate matrices.
  S1 <- as.matrix(S1)
  S2 <- as.matrix(S2)
  if (!identical(dim(S1), dim(S2))) {
    stop("Internal error: paired distance inputs must have matching dimensions.")
  }
  storage.mode(S1) <- "double"
  storage.mode(S2) <- "double"
  .Call("cpp_distvec", S1, S2, as.logical(latlong), PACKAGE = "scpcR")
}

.getc0fromavc <- function(dist, avc0) {
  ## Solve for exponential kernel parameter c0 such that
  ## mean(exp(-c0 * dist)) == avc0, using bisection.
  storage.mode(dist) <- "double"
  .Call("cpp_getc0_from_avc", dist, as.double(avc0), PACKAGE = "scpcR")
}

#' @importFrom gaussquad legendre.quadrature.rules
.GQ <- gaussquad::legendre.quadrature.rules(40)[[40]]
.GQx <- .GQ[, 1] * 0.5 + 0.5
.GQw <- .GQ[, 2] * 0.5

.getrp <- function(Om, cv) {
  ## Rejection probability for the SCPC test given Omega matrix Om
  ## and normalised critical value cv (= cv/sqrt(q) in paper notation).
  ## Uses Gauss-Legendre quadrature over the characteristic function.
  storage.mode(Om) <- "double"
  .Call("cpp_getrp", Om, as.double(cv), .GQx, .GQw, PACKAGE = "scpcR")
}

.maxrp <- function(Oms, q, cv) {
  ## Largest rejection probability across a list of Omega matrices
  .Call("cpp_maxrp", Oms, as.integer(q), as.double(cv), .GQx, .GQw, PACKAGE = "scpcR")
}

.getcv <- function(Oms, q, level) {
  ## Two-sided critical value that controls size at the given level,
  ## maximised over the Omega grid.
  .Call("cpp_getcv", Oms, as.integer(q), as.double(level), .GQx, .GQw, PACKAGE = "scpcR")
}

.setfinalW <- function(Oms, W, qmax) {
  ## Select the optimal q (number of spatial PCs) that minimises the
  ## expected confidence interval length under i.i.d. errors.
  storage.mode(W) <- "double"
  .Call("cpp_setfinalW", Oms, W, as.integer(qmax), .GQx, .GQw, PACKAGE = "scpcR")
}

.getnc <- function(c0, cmax, cgridfac) {
  ## Number of c-values in the multiplicative grid
  max(2, ceiling(log(cmax / c0) / log(cgridfac)))
}

.getOms <- function(distmat, c0, cmax, W, cgridfac) {
  ## List of Omega(c) = W' Sigma(c) W matrices over the c-grid
  storage.mode(distmat) <- "double"
  storage.mode(W) <- "double"
  .Call(
    "cpp_get_oms",
    distmat,
    as.double(c0),
    as.double(cmax),
    W,
    as.double(cgridfac),
    PACKAGE = "scpcR"
  )
}

.normalize_s <- function(S, latlong) {
  ## This follows the original Princeton large-n branch closely so the
  ## approximation uses the same row order as the reference code.
  S <- as.matrix(S)

  if (!latlong) {
    S <- sweep(S, 2, colMeans(S), `-`)
    rot <- eigen(crossprod(S), symmetric = TRUE)$vectors
    S <- S %*% rot
    if (max(S[, 1]) != max(abs(S[, 1]))) {
      S <- -S
    }
    perm <- do.call(order, as.data.frame(S))
    S <- S[perm, , drop = FALSE]
    S <- sweep(S, 2, apply(S, 2, min), `-`)
    smax <- max(S)
    if (is.finite(smax) && smax > 0) {
      S <- S / smax
    }
  } else {
    S[, 1] <- S[, 1] - mean(S[, 1])
    S[, 1] <- ((S[, 1] + 180) %% 360) - 180
    perm <- order(S[, 2], S[, 1])
    S <- S[perm, , drop = FALSE]
  }

  list(coords = S, perm = perm)
}

.next_u <- function(random_t) {
  random_t <- (64389 * random_t + 1) %% 2^32
  list(value = random_t / 2^32, state = random_t)
}

.jumble_s <- function(S, m, random_t) {
  S <- as.matrix(S)
  storage.mode(S) <- "double"
  .Call("cpp_jumble_s", S, as.integer(m), as.double(random_t), PACKAGE = "scpcR")
}

.ln_subset_evecs <- function(distmat, c0, qmax) {
  distmat <- as.matrix(distmat)
  storage.mode(distmat) <- "double"
  .Call("cpp_ln_subset_evecs", distmat, as.double(c0), as.integer(qmax), PACKAGE = "scpcR")
}

.lnset_wc0 <- function(S, avc0, qmax, minavc, latlong,
                       capN = 20L, m = 1000L, random_t = 1) {
  S <- as.matrix(S)
  storage.mode(S) <- "double"
  .Call(
    "cpp_lnset_wc0",
    S,
    as.double(avc0),
    as.integer(qmax),
    as.double(minavc),
    as.logical(latlong),
    as.integer(capN),
    as.integer(m),
    as.double(random_t),
    PACKAGE = "scpcR"
  )
}

.raninds <- function(n, capM, random_t) {
  .Call("cpp_raninds", as.integer(n), as.integer(capM), as.double(random_t), PACKAGE = "scpcR")
}

.lnget_Oms <- function(S, c0, cmax, W, cgridfac,
                       capM = 1000000L, random_t = 1, latlong = FALSE) {
  S <- as.matrix(S)
  W <- as.matrix(W)
  storage.mode(S) <- "double"
  storage.mode(W) <- "double"
  .Call(
    "cpp_lnget_oms",
    S,
    as.double(c0),
    as.double(cmax),
    W,
    as.double(cgridfac),
    as.integer(capM),
    as.double(random_t),
    as.logical(latlong),
    PACKAGE = "scpcR"
  )
}

.validate_large_n_seed <- function(seed) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed >= 2^32 || seed != floor(seed)) {
    stop("`large_n_seed` must be a single integer-valued number in [0, 2^32).")
  }
  as.numeric(seed)
}

.validate_scpc_method <- function(method) {
  methods <- c("auto", "exact", "approx")
  if (!is.character(method) || length(method) != 1L || !nzchar(method) ||
      !method %in% methods) {
    stop("`method` must be one of \"auto\", \"exact\", or \"approx\".")
  }
  method
}

.resolve_scpc_method <- function(n, method, large_n_threshold = 4500L) {
  ## this is only the old exact-vs-approx routing rule pulled into one
  ## place so the 4500 threshold can be tested without running scpc().
  method <- .validate_scpc_method(method)
  if (identical(method, "auto")) {
    if (n < large_n_threshold) "exact" else "approx"
  } else {
    method
  }
}

.setOmsWfin <- function(coords, avc0, latlong, method = "auto", large_n_seed = 1) {
  ## keep the old large-n algorithm as-is, but let the routing and seed
  ## checks happen up front so the split code stays easier to read.
  method <- .validate_scpc_method(method)
  large_n_seed <- .validate_large_n_seed(large_n_seed)
  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  n <- nrow(coords)

  cgridfac <- 1.2
  minavc <- 0.00001

  large_n_threshold <- 4500L
  large_n_capN <- 20L
  large_n_capM <- 1000000L
  large_n_m <- 1000L
  method_actual <- .resolve_scpc_method(n, method, large_n_threshold = large_n_threshold)

  cached <- .scpc_spatial_cache$last
  if (!is.null(cached) &&
      identical(cached$avc0, as.double(avc0)) &&
      identical(cached$latlong, as.logical(latlong)) &&
      identical(cached$method_actual, method_actual) &&
      (identical(method_actual, "exact") || identical(cached$large_n_seed, large_n_seed)) &&
      identical(dim(cached$coords_in), dim(coords)) &&
      identical(as.vector(cached$coords_in), as.vector(coords))) {
    return(cached$result)
  }

  if (avc0 >= 0.05) {
    qmax <- 10
  } else if (avc0 >= 0.01) {
    qmax <- 20
  } else if (avc0 >= 0.005) {
    qmax <- 60
  } else {
    qmax <- 120
  }

  distmat <- NULL
  coords_use <- coords
  perm <- seq_len(n)
  random_t <- NULL

  if (identical(method_actual, "exact")) {
    distmat <- .getdistmat(coords, latlong)
    correlation_floor <- (sum(distmat == 0) - n) / (n * (n - 1))
    minavc_target <- correlation_floor + (1 - correlation_floor) * minavc
    c_pair <- .Call(
      "cpp_getc0_from_distmat",
      distmat,
      as.double(c(avc0, minavc_target)),
      PACKAGE = "scpcR"
    )
    c0 <- c_pair[[1L]]
    cmax <- c_pair[[2L]]
  } else {
    norm_s <- .normalize_s(coords, latlong)
    coords_use <- norm_s$coords
    perm <- norm_s$perm
  }

  repeat {
    qmax <- min(qmax, n - 1)
    if (identical(method_actual, "exact")) {
      W <- .getW(distmat, c0, qmax)
      Oms <- .getOms(distmat, c0, cmax, W, cgridfac)
      random_t <- NULL
    } else {
      random_t <- large_n_seed
      ln_w <- .lnset_wc0(
        coords_use, avc0, qmax, minavc, latlong,
        capN = large_n_capN, m = large_n_m, random_t = random_t
      )
      W <- ln_w$W
      c0 <- ln_w$c0
      cmax <- ln_w$cmax
      oms_res <- .lnget_Oms(
        coords_use, c0, cmax, W, cgridfac,
        capM = large_n_capM, random_t = ln_w$random_t, latlong = latlong
      )
      Oms <- oms_res$Oms
      random_t <- oms_res$state
    }
    fin <- .setfinalW(Oms, W, qmax)
    if (fin$q < qmax || qmax == n - 1) break
    qmax <- round(qmax + qmax / 2)
  }

  res <- list(
    Wfin = fin$W,
    cvfin = fin$cv,
    Omsfin = Oms,
    c0 = c0,
    cmax = cmax,
    coords = coords_use,
    perm = perm,
    distmat = distmat,
    method = method_actual,
    large_n = identical(method_actual, "approx"),
    random_state = random_t
  )

  .scpc_spatial_cache$last <- list(
    coords_in = coords,
    avc0 = as.double(avc0),
    latlong = as.logical(latlong),
    method_actual = method_actual,
    large_n_seed = large_n_seed,
    result = res
  )

  res
}

.get_conditional_Oms <- function(spc, Wx, latlong, random_state) {
  ## the old code did this inline in scpc(). this helper keeps the same
  ## exact-vs-large-n logic in one place so core.R stays readable.
  if (isTRUE(spc$large_n)) {
    omsx_res <- .lnget_Oms(
      spc$coords,
      spc$c0,
      spc$cmax,
      Wx,
      1.2,
      capM = 1000000L,
      random_t = random_state,
      latlong = latlong
    )
    return(list(Oms = omsx_res$Oms, random_state = omsx_res$state))
  }

  list(
    Oms = .getOms(spc$distmat, spc$c0, spc$cmax, Wx, 1.2),
    random_state = random_state
  )
}
