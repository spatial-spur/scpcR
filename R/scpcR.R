# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.getavc <- function(c, dist) {
  ## Average pairwise correlation for exponential kernel parameter c
  mean(exp(-c * dist))
}

.getdistmat <- function(S, latlong) {
  ## Distance matrix from coordinates.
  ## latlong = TRUE: haversine distances in units of pi * R_earth.
  ## latlong = FALSE: Euclidean distances.
  S <- as.matrix(S)
  if (latlong) {
    if (ncol(S) != 2L) {
      stop("Internal error: geodesic coordinates must have exactly two columns.")
    }
    if (!requireNamespace("geodist", quietly = TRUE))
      stop("Install package 'geodist' for lat/long distances.")
    S <- as.data.frame(S)
    names(S) <- c("lon", "lat")
    D <- geodist::geodist(S, measure = "haversine") / (2 * pi * 6378137)
  } else {
    D <- as.matrix(stats::dist(S, diag = TRUE, upper = TRUE))
  }
  D
}

.getdistvec <- function(S1, S2, latlong) {
  ## Vector of paired distances between two coordinate matrices.
  S1 <- as.matrix(S1)
  S2 <- as.matrix(S2)
  if (!identical(dim(S1), dim(S2))) {
    stop("Internal error: paired distance inputs must have matching dimensions.")
  }
  if (latlong) {
    lon1 <- S1[, 1] * pi / 180
    lat1 <- S1[, 2] * pi / 180
    lon2 <- S2[, 1] * pi / 180
    lat2 <- S2[, 2] * pi / 180
    dlon <- 0.5 * (lon1 - lon2)
    dlat <- 0.5 * (lat1 - lat2)
    asin(sqrt(sin(dlat)^2 + cos(lat1) * cos(lat2) * sin(dlon)^2)) / pi
  } else {
    sqrt(rowSums((S1 - S2)^2))
  }
}

.lvech <- function(mat) {
  ## Lower-triangular elements as a vector
  mat[lower.tri(mat)]
}

.demeanmat <- function(mat) {
  ## Double-demean a matrix (row and column means removed)
  mat <- mat - rowMeans(mat)
  mat <- mat - matrix(colMeans(mat), nrow(mat), ncol(mat), byrow = TRUE)
  mat
}

.getc0fromavc <- function(dist, avc0) {
  ## Solve for exponential kernel parameter c0 such that
  ## mean(exp(-c0 * dist)) == avc0, using bisection.
  c0 <- c1 <- 10
  while (.getavc(c0, dist) < avc0) {
    c1 <- c0
    c0 <- 0.5 * c0
  }
  while (.getavc(c1, dist) > avc0 & c1 < 5000) {
    c0 <- c1
    c1 <- 2 * c1
  }
  repeat {
    c <- sqrt(c0 * c1)
    if (.getavc(c, dist) > avc0) c0 <- c else c1 <- c
    if (c1 - c0 < 0.001) break
  }
  c
}

.getW <- function(distmat, c0, qmax) {
  ## First qmax eigenvectors of demeaned Sigma(c0), prepended with

  ## a normalised constant vector.  All columns have unit length.
  n <- nrow(distmat)
  Sig <- exp(-c0 * distmat)
  Sig_d <- .demeanmat(Sig)
  if (requireNamespace("RSpectra", quietly = TRUE) && qmax < n - 1) {
    eig <- RSpectra::eigs_sym(Sig_d, k = qmax, which = "LM")
    V   <- eig$vectors
  } else {
    V <- eigen(Sig_d, symmetric = TRUE)$vectors[, seq_len(qmax)]
  }
  cbind(rep(1, n) / sqrt(n), V)
}

#' @importFrom gaussquad legendre.quadrature.rules
.GQ  <- gaussquad::legendre.quadrature.rules(40)[[40]]
.GQx <- .GQ[, 1] * 0.5 + 0.5
.GQw <- .GQ[, 2] * 0.5

.getrp <- function(Om, cv) {
  ## Rejection probability for the SCPC test given Omega matrix Om

  ## and normalised critical value cv (= cv/sqrt(q) in paper notation).
  ## Uses Gauss-Legendre quadrature over the characteristic function.
  Omx       <- -cv^2 * Om
  Omx[1, ]  <- Om[1, ]
  evals_raw <- Re(eigen(Omx, only.values = TRUE)$values)
  evals     <- -evals_raw[evals_raw < 0]
  if (!length(evals)) return(0)
  denom <- max(evals_raw)
  if (!is.finite(denom) || denom <= 0) return(0)
  evals <- evals / denom
  tot   <- 0
  for (j in seq_along(.GQx)) {
    u   <- .GQx[j]
    arg <- (1 - u^2) * exp(sum(log1p(evals / (1 - u^2))))
    if (!is.finite(arg) || arg <= 0) next
    tot <- tot + .GQw[j] / sqrt(arg)
  }
  as.numeric(tot * 2 / pi)
}

.maxrp <- function(Oms, q, cv) {
  ## Largest rejection probability across a list of Omega matrices
  rps <- vapply(Oms, function(Om) .getrp(Om[1:(q + 1), 1:(q + 1)], cv), 0.0)
  list(max = max(rps), i = which.max(rps))
}

.getcv <- function(Oms, q, level) {
  ## Two-sided critical value that controls size at the given level,
  ## maximised over the Omega grid.
  rp <- 1
  i  <- 1
  cv0 <- stats::qt(1 - level / 2, df = q) / sqrt(q)
  while (rp > level) {
    cv1 <- cv0
    repeat {
      if (.getrp(Oms[[i]][1:(q + 1), 1:(q + 1)], cv1) > level) {
        cv0 <- cv1
        cv1 <- cv1 + 1 / sqrt(q)
      } else break
    }
    while (cv1 - cv0 > 0.001 / sqrt(q)) {
      cv <- 0.5 * (cv0 + cv1)
      if (.getrp(Oms[[i]][1:(q + 1), 1:(q + 1)], cv) > level) cv0 <- cv else cv1 <- cv
    }
    maxrp <- .maxrp(Oms, q, cv1)
    if (maxrp$i == i) break
    i   <- maxrp$i
    rp  <- maxrp$max
    cv0 <- cv1
  }
  cv1 * sqrt(q)
}

.setfinalW <- function(Oms, W, qmax) {
  ## Select the optimal q (number of spatial PCs) that minimises the
  ## expected confidence interval length under i.i.d. errors.
  cvs     <- lengths <- numeric(qmax)
  for (q in seq_len(qmax)) {
    cvs[q]     <- .getcv(Oms, q, 0.05)
    lengths[q] <- cvs[q] * gamma(0.5 * (q + 1)) / (sqrt(q) * gamma(0.5 * q))
  }
  q_opt <- which.min(lengths)
  list(W = W[, 1:(q_opt + 1), drop = FALSE], cv = cvs[q_opt], q = q_opt)
}

.getnc <- function(c0, cmax, cgridfac) {
  ## Number of c-values in the multiplicative grid
  max(2, ceiling(log(cmax / c0) / log(cgridfac)))
}

.getOms <- function(distmat, c0, cmax, W, cgridfac) {
  ## List of Omega(c) = W' Sigma(c) W matrices over the c-grid
  nc  <- .getnc(c0, cmax, cgridfac)
  Oms <- vector("list", nc)
  c   <- c0
  for (i in seq_len(nc)) {
    Oms[[i]] <- crossprod(W, exp(-c * distmat) %*% W)
    c <- c * cgridfac
  }
  Oms
}

.gettau <- function(y, W) {
  ## SCPC t-statistic
  sqrt(ncol(W) - 1) * crossprod(W[, 1], y) / norm(crossprod(W[, -1], y), type = "2")
}

# ---------------------------------------------------------------------------
# Large-n helpers (Stata approximation branch)
# ---------------------------------------------------------------------------
.normalize_s <- function(S, latlong) {
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
    ## Keep the internal lon/lat order but reproduce Stata's longitude wrap
    ## and lexicographic ordering by latitude, then longitude.
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
  n <- nrow(S)
  for (i in seq_len(m)) {
    nxt <- .next_u(random_t)
    random_t <- nxt$state
    j <- floor(nxt$value * n) + 1L
    tmp <- S[j, ]
    S[j, ] <- S[i, ]
    S[i, ] <- tmp
  }
  list(coords = S, state = random_t)
}

.ln_subset_evecs <- function(distmat, c0, qmax) {
  Sig_d <- .demeanmat(exp(-c0 * distmat))
  n <- nrow(Sig_d)
  if (requireNamespace("RSpectra", quietly = TRUE) && qmax < n - 1) {
    RSpectra::eigs_sym(Sig_d, k = qmax, which = "LM")$vectors
  } else {
    eigen(Sig_d, symmetric = TRUE)$vectors[, seq_len(qmax), drop = FALSE]
  }
}

.lnset_wc0 <- function(S, avc0, qmax, minavc, latlong,
                       capN = 20L, m = 1000L, random_t = 1) {
  n <- nrow(S)
  m <- min(m, n)
  ms <- vector("list", capN)
  block_len <- m * (m - 1) / 2
  distvec <- numeric(capN * block_len)

  r <- S
  for (i in seq_len(capN)) {
    jumbled <- .jumble_s(r, m, random_t)
    r <- jumbled$coords
    random_t <- jumbled$state
    ms[[i]] <- list(
      coords = r[seq_len(m), , drop = FALSE],
      distmat = .getdistmat(r[seq_len(m), , drop = FALSE], latlong)
    )
    idx <- ((i - 1L) * block_len + 1L):(i * block_len)
    distvec[idx] <- .lvech(ms[[i]]$distmat)
  }

  c0 <- .getc0fromavc(distvec, avc0)
  cmax <- .getc0fromavc(distvec, minavc)

  Wall <- matrix(0, n, capN * qmax)
  for (i in seq_len(capN)) {
    W0 <- .ln_subset_evecs(ms[[i]]$distmat, c0, qmax)
    Wx <- matrix(0, n, qmax)
    for (j in seq_len(m)) {
      diff <- sweep(S, 2, ms[[i]]$coords[j, ], `-`)
      v <- exp(-c0 * sqrt(rowSums(diff^2)))
      Wx <- Wx + tcrossprod(v, W0[j, ])
    }
    Wx <- sweep(Wx, 2, colMeans(Wx), `-`)
    norms <- sqrt(colSums(Wx^2))
    norms[!is.finite(norms) | norms == 0] <- 1
    Wx <- sweep(Wx, 2, norms, `/`)
    for (j in seq_len(qmax)) {
      Wall[, (j - 1L) * capN + i] <- Wx[, j]
    }
  }

  W <- matrix(0, n, qmax)
  for (i in seq_len(qmax)) {
    Wx <- Wall[, seq_len(capN * i), drop = FALSE]
    evec <- eigen(crossprod(Wx), symmetric = TRUE)$vectors[, 1, drop = FALSE]
    W[, i] <- as.numeric(Wx %*% evec)
    W[, i] <- W[, i] / sqrt(sum(W[, i]^2))
    Wall <- Wall - W[, i, drop = FALSE] %*% crossprod(W[, i], Wall)
  }

  list(
    W = cbind(rep(1 / sqrt(n), n), W),
    c0 = c0,
    cmax = cmax,
    random_t = random_t
  )
}

.raninds <- function(n, capM, random_t) {
  v <- numeric(capM + 1L)
  nxt <- .next_u(random_t)
  random_t <- nxt$state
  j <- floor(n * nxt$value)

  for (i in seq_len(capM + 1L)) {
    v[i] <- j + 1L
    nxt <- .next_u(random_t)
    random_t <- nxt$state
    j <- (j + 1 + floor(nxt$value * (n - 1))) %% n
  }

  list(indices = as.integer(v), state = random_t)
}

.lnget_Oms <- function(S, c0, cmax, W, cgridfac,
                       capM = 1000000L, random_t = 1, latlong = FALSE) {
  nc <- .getnc(c0, cmax, cgridfac)
  Oms <- vector("list", nc)

  n <- nrow(S)
  inds_res <- .raninds(n, capM, random_t)
  inds <- inds_res$indices
  dist <- .getdistvec(S[inds[seq_len(capM)], , drop = FALSE],
                      S[inds[2:(capM + 1L)], , drop = FALSE],
                      latlong)
  W1 <- W[inds[seq_len(capM)], , drop = FALSE]
  W2 <- W[inds[2:(capM + 1L)], , drop = FALSE]

  c <- c0
  for (i in seq_len(nc)) {
    cd <- exp(-c * dist)
    Oms[[i]] <- diag(ncol(W)) +
      0.5 * (n * (n - 1) / capM) *
      (crossprod(W1, W2 * cd) + crossprod(W2, W1 * cd))
    c <- c * cgridfac
  }

  list(Oms = Oms, state = inds_res$state)
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

# ---------------------------------------------------------------------------
# Orthogonalisation helper (conditional SCPC)
# ---------------------------------------------------------------------------
.orthogonalize_W <- function(W, xj, xjs, model_mat, include_intercept = TRUE, fixef_id = NULL) {
  Wx <- W
  Wx[, 1] <- Wx[, 1] * xj * xjs

  if (ncol(Wx) > 1L) {
    X <- as.matrix(model_mat)
    has_intercept_col <- !is.null(colnames(X)) && "(Intercept)" %in% colnames(X)
    if (isTRUE(include_intercept) && !has_intercept_col) {
      X <- cbind("(Intercept)" = 1, X)
    }
    qrX <- qr(X)
    Rx <- sweep(Wx[, -1, drop = FALSE], 1, xjs, `*`)
    if (!is.null(fixef_id)) {
      Rx <- as.matrix(fixest::demean(Rx, f = fixef_id, nthreads = 1L))
    }
    Rr <- qr.resid(qrX, Rx)
    Wx[, -1] <- sweep(Rr, 1, xj, `*`)
  }
  Wx
}

# ---------------------------------------------------------------------------
# Orthogonalisation helper (conditional SCPC, clustered)
# ---------------------------------------------------------------------------
.orthogonalize_W_cluster <- function(W, cl_vec, xj_indiv, model_mat_indiv,
                                     include_intercept = TRUE) {
  ## Construct the conditional Wx matrix for clustered SCPC following the

  ## Stata set_Wx_cluster algorithm (Mueller & Watson):
  ##   1. Within-cluster normalization of influence directions
  ##   2. Expand cluster-level W to individual observations
  ##   3. Orthogonalize at individual level against regressors
  ##   4. Aggregate back to cluster level
  ##
  ## W:               nclust x (q+1) cluster-level spatial projection matrix
  ## cl_vec:          factor of cluster membership (length n_individual)
  ## xj_indiv:        numeric (length n_individual) influence function directions
  ## model_mat_indiv: n_individual x p model matrix for projection
  nclust <- nrow(W)
  ncol_W <- ncol(W)
  cl_idx <- as.integer(cl_vec)

  ## Within-cluster normalization: xjs_i = xj_i / sqrt(sum_{i in g} xj_i^2)
  xj_sq_sum <- as.numeric(rowsum(xj_indiv^2, cl_vec))
  xjs_indiv <- xj_indiv / sqrt(xj_sq_sum[cl_idx])
  xjs_indiv[!is.finite(xjs_indiv)] <- 0

  ## Expand W from cluster level to individual level
  W_indiv <- W[cl_idx, , drop = FALSE]

  Wx <- matrix(0, nclust, ncol_W)

  ## Column 1: Wx[g,1] = sum_{i in g}(W[g,1] * xjs_i * xj_i)
  Wx[, 1] <- as.numeric(rowsum(W_indiv[, 1] * xjs_indiv * xj_indiv, cl_vec))

  ## Columns > 1: orthogonalize at individual level, aggregate back
  if (ncol_W > 1L) {
    X <- as.matrix(model_mat_indiv)
    has_intercept_col <- !is.null(colnames(X)) && "(Intercept)" %in% colnames(X)
    if (isTRUE(include_intercept) && !has_intercept_col) {
      X <- cbind("(Intercept)" = 1, X)
    }
    qrX <- qr(X)

    for (col in 2:ncol_W) {
      temp <- W_indiv[, col] * xjs_indiv
      resid_col <- qr.resid(qrX, temp)
      Wx[, col] <- as.numeric(rowsum(resid_col * xj_indiv, cl_vec))
    }
  }

  Wx
}

# ---------------------------------------------------------------------------
# Main spatial engine
# ---------------------------------------------------------------------------
.setOmsWfin <- function(coords, avc0, latlong, method = "auto", large_n_seed = 1) {
  n <- nrow(coords)

  cgridfac <- 1.2
  minavc   <- 0.00001

  large_n_threshold <- 4500L
  large_n_capN <- 20L
  large_n_capM <- 1000000L
  large_n_m <- 1000L
  method_actual <- if (identical(method, "auto")) {
    if (n < large_n_threshold) "exact" else "approx"
  } else {
    method
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

  repeat {
    qmax <- min(qmax, n - 1)
    if (identical(method_actual, "exact")) {
      distmat <- .getdistmat(coords, latlong)
      distv <- .lvech(distmat)
      c0 <- .getc0fromavc(distv, avc0)
      cmax <- .getc0fromavc(distv, minavc)
      W <- .getW(distmat, c0, qmax)
      Oms <- .getOms(distmat, c0, cmax, W, cgridfac)
      coords_use <- coords
      perm <- seq_len(n)
      random_t <- NULL
    } else {
      random_t <- large_n_seed
      norm_s <- .normalize_s(coords, latlong)
      coords_use <- norm_s$coords
      perm <- norm_s$perm
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
  list(
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
}

.get_obs_index <- function(model, data) {
  if (inherits(model, "fixest")) {
    if (!requireNamespace("fixest", quietly = TRUE)) {
      stop("Package 'fixest' is required for fixest model objects.")
    }
    obs_index <- fixest::obs(model)
  } else if (inherits(model, "lm")) {
    mf <- stats::model.frame(model)
    rn <- rownames(mf)
    suppressWarnings(obs_index <- as.integer(rn))
    if (anyNA(obs_index)) {
      obs_index <- match(rn, rownames(data))
      if (anyNA(obs_index)) {
        stop("Could not map lm model frame rows back to `data`.")
      }
    }
  } else {
    stop("Model class not supported. Provide an lm or fixest model.")
  }

  if (any(obs_index < 1L | obs_index > nrow(data))) {
    stop("Model observation indices are outside the row range of `data`.")
  }
  obs_index
}

.get_scpc_model_matrix <- function(model) {
  if (inherits(model, "fixest") && isTRUE(model$iv) && isTRUE(model$iv_stage == 2)) {
    mm <- tryCatch(
      utils::getFromNamespace("model.matrix.fixest", "fixest")(model, type = "iv.rhs2"),
      error = function(e) NULL
    )
    if (!is.null(mm)) {
      return(mm)
    }
  }
  stats::model.matrix(model)
}

.has_fixest_fe <- function(model) {
  inherits(model, "fixest") &&
    !is.null(model$fixef_vars) &&
    length(model$fixef_vars) > 0L
}

.get_conditional_projection_setup <- function(model, model_mat, n, uncond) {
  setup <- list(model_mat = as.matrix(model_mat), include_intercept = TRUE, fixef_id = NULL)
  if (isTRUE(uncond) || !inherits(model, "fixest") || !.has_fixest_fe(model)) {
    return(setup)
  }

  mm <- as.matrix(model_mat)
  coef_names <- names(stats::coef(model))
  if (!is.null(colnames(mm)) && all(coef_names %in% colnames(mm))) {
    mm <- mm[, coef_names, drop = FALSE]
  } else if (ncol(mm) == length(coef_names)) {
    colnames(mm) <- coef_names
  } else {
    stop(
      "Could not align demeaned fixest regressors with model coefficients ",
      "(coef count = ", length(coef_names), ", demeaned columns = ", ncol(mm), ")."
    )
  }

  if (nrow(mm) != n) {
    stop("Internal error: fixest conditional regressors have incompatible row count.")
  }

  mm <- tryCatch(
    as.matrix(fixest::demean(mm, f = model$fixef_id, nthreads = 1L)),
    error = function(e) e
  )
  if (inherits(mm, "error")) {
    stop("Could not FE-demean conditional regressors for fixest model: ", conditionMessage(mm))
  }
  list(model_mat = mm, include_intercept = FALSE, fixef_id = model$fixef_id)
}

.resolve_coords_input <- function(data, obs_index, lon, lat, coords_euclidean) {
  use_geodesic <- !is.null(lon) || !is.null(lat)
  use_euclidean <- !is.null(coords_euclidean)

  if (use_geodesic && use_euclidean) {
    stop("Specify either `lon`/`lat` or `coords_euclidean`, not both.")
  }
  if (!use_geodesic && !use_euclidean) {
    stop("Specify coordinates via `lon`/`lat` or `coords_euclidean`.")
  }

  if (use_geodesic) {
    if (is.null(lon) || is.null(lat)) {
      stop("For geodesic coordinates, provide both `lon` and `lat`.")
    }
    if (!is.character(lon) || length(lon) != 1L || !nzchar(lon) ||
        !is.character(lat) || length(lat) != 1L || !nzchar(lat)) {
      stop("`lon` and `lat` must each be a single column name.")
    }
    miss <- setdiff(c(lon, lat), names(data))
    if (length(miss) > 0) {
      stop("Coordinate variables not found in data: ", paste(miss, collapse = ", "))
    }
    coords <- data[obs_index, c(lon, lat), drop = FALSE]
    if (!all(vapply(coords, is.numeric, logical(1)))) {
      stop("`lon` and `lat` must reference numeric columns.")
    }
    if (any(!is.finite(as.matrix(coords)))) {
      stop("Geodesic coordinates must be finite.")
    }
    if (any(coords[[lon]] < -180 | coords[[lon]] > 180)) {
      stop("Longitude values must be in [-180, 180].")
    }
    if (any(coords[[lat]] < -90 | coords[[lat]] > 90)) {
      stop("Latitude values must be in [-90, 90].")
    }
    return(list(coords = as.matrix(coords), latlong = TRUE))
  }

  if (!is.character(coords_euclidean) || length(coords_euclidean) < 1L) {
    stop("`coords_euclidean` must be a character vector with at least one column name.")
  }
  miss <- setdiff(coords_euclidean, names(data))
  if (length(miss) > 0) {
    stop("Coordinate variables not found in data: ", paste(miss, collapse = ", "))
  }
  coords <- data[obs_index, coords_euclidean, drop = FALSE]
  if (!all(vapply(coords, is.numeric, logical(1)))) {
    stop("`coords_euclidean` columns must be numeric.")
  }
  if (any(!is.finite(as.matrix(coords)))) {
    stop("Euclidean coordinates must be finite.")
  }
  list(coords = as.matrix(coords), latlong = FALSE)
}

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

#' Spatial Correlation-Robust Inference (SCPC)
#'
#' Compute spatial correlation-robust inference for regression coefficients
#' using the SCPC method of Mueller and Watson (2022, 2023).
#'
#' @param model Fitted model object.
#'   Currently supported: \code{\link{lm}} and \pkg{fixest}
#'   (\code{\link[fixest]{feols}}) objects, including IV models.
#' @param data Data frame used to fit \code{model}.  Must contain the
#'   coordinate columns referenced by \code{lon}/\code{lat} or
#'   \code{coords_euclidean}.
#' @param lon Character string naming the longitude column in \code{data}.
#'   Must be supplied together with \code{lat}.
#' @param lat Character string naming the latitude column in \code{data}.
#'   Must be supplied together with \code{lon}.
#' @param coords_euclidean Character vector of one or more column names in
#'   \code{data} for Euclidean coordinates.  Supply either
#'   \code{lon}/\code{lat} or \code{coords_euclidean}, not both.
#' @param cluster Optional character string naming a clustering variable
#'   in \code{data}.  When clustering is used, coordinates are taken from
#'   the first observation in each cluster; they should be constant
#'   within clusters.
#' @param ncoef Integer; number of coefficients to report.  Default
#'   \code{NULL} reports all coefficients.
#' @param avc Numeric; upper bound on the average pairwise correlation
#'   for which size is controlled.  Must be in \code{(0.001, 0.99)}.
#'   Default is 0.03.
#' @param method Character string selecting the computational method.
#'   \code{"auto"} (default) chooses \code{"exact"} for smaller problems
#'   and \code{"approx"} for larger ones (currently at \eqn{n < 4500}
#'   versus \eqn{n >= 4500}). \code{"exact"} computes the full pairwise
#'   distance matrix and full SCPC objects. \code{"approx"} uses a
#'   randomized approximation that avoids forming the full dense
#'   distance and kernel matrices.
#' @param large_n_seed Numeric; integer-valued seed used by the
#'   randomized approximation when \code{method = "approx"} or when
#'   \code{method = "auto"} selects it. Ignored when
#'   \code{method = "exact"}. Default is \code{1}.
#' @param uncond Logical; if \code{TRUE}, report unconditional critical
#'   values only (skip the conditional adjustment of Mueller and Watson
#'   2023).  Default is \code{FALSE}.
#' @param cvs Logical; if \code{TRUE}, include per-coefficient critical
#'   values at the 32\%, 10\%, 5\%, and 1\% levels.  Default is
#'   \code{FALSE}.
#'
#' @return An object of class \code{"scpc"} with components:
#' \describe{
#'   \item{\code{scpcstats}}{Matrix of coefficient estimates, SCPC standard
#'     errors, t-statistics, p-values, and 95\% confidence interval bounds.}
#'   \item{\code{scpccvs}}{Matrix of two-sided critical values at the 32\%,
#'     10\%, 5\%, and 1\% levels, or \code{NULL} when \code{cvs = FALSE}.}
#'   \item{\code{W}}{The spatial projection matrix (n x (q+1)).}
#'   \item{\code{avc}}{The average correlation bound used.}
#'   \item{\code{c0}}{Exponential kernel parameter corresponding to
#'     \code{avc}.}
#'   \item{\code{cv}}{The unconditional 5\% critical value.}
#'   \item{\code{q}}{Number of spatial principal components selected.}
#'   \item{\code{method}}{Spatial algorithm actually used:
#'     \code{"exact"} or \code{"approx"}.}
#'   \item{\code{large_n_seed}}{Seed used for the large-\eqn{n}
#'     approximation branch.}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @references
#' Mueller, U. K. and Watson, M. W. (2022).
#' \dQuote{Spatial Correlation Robust Inference.}
#' \emph{Econometrica}, 90(6), 2901--2935.
#' \doi{10.3982/ECTA19465}
#'
#' Mueller, U. K. and Watson, M. W. (2023).
#' \dQuote{Spatial Correlation Robust Inference in Linear Regression and
#' Panel Models.}
#' \emph{Journal of Business & Economic Statistics}, 41(4), 1050--1064.
#' \doi{10.1080/07350015.2022.2127737}
#'
#' @examples
#' set.seed(42)
#' n <- 60
#' d <- data.frame(
#'   y   = rnorm(n),
#'   x   = rnorm(n),
#'   lon = runif(n, -100, -80),
#'   lat = runif(n, 30, 45)
#' )
#' fit <- lm(y ~ x, data = d)
#'
#' # Euclidean coordinates, unconditional
#' out <- scpc(fit, data = d, coords_euclidean = c("lon", "lat"),
#'             avc = 0.1, uncond = TRUE)
#' out
#'
#' # Extract coefficients and confidence intervals
#' coef(out)
#' confint(out)
#'
#' @export
scpc <- function(model,
                 data,
                 lon      = NULL,
                 lat      = NULL,
                 coords_euclidean = NULL,
                 cluster  = NULL,
                 ncoef    = NULL,
                 avc      = 0.03,
                 method   = "auto",
                 large_n_seed = 1,
                 uncond   = FALSE,
                 cvs      = FALSE) {

  if (avc <= 0.001 || avc >= 0.99)
    stop("Option avc() must be in (0.001, 0.99).")
  method <- .validate_scpc_method(method)
  large_n_seed <- .validate_large_n_seed(large_n_seed)

  ## 1. Influence functions --------------------------------------------------
  S    <- sandwich::estfun(model)
  model_mat <- .get_scpc_model_matrix(model)
  n    <- nrow(S); p <- ncol(S); neff <- n
  if (nrow(model_mat) != n) {
    stop("Model matrix row count does not match score matrix row count.")
  }

  cond_setup <- .get_conditional_projection_setup(model, model_mat, n = n, uncond = uncond)
  model_mat_cond <- cond_setup$model_mat
  cond_include_intercept <- cond_setup$include_intercept
  cond_fixef_id <- cond_setup$fixef_id
  if (nrow(model_mat_cond) != n || ncol(model_mat_cond) != p) {
    stop(
      "Conditional projection matrix dimensions are incompatible with model coefficients ",
      "(rows = ", nrow(model_mat_cond), ", cols = ", ncol(model_mat_cond),
      ", expected ", n, " x ", p, ")."
    )
  }

  obs_index <- .get_obs_index(model, data)
  coord_info <- .resolve_coords_input(data, obs_index, lon, lat, coords_euclidean)
  coords  <- coord_info$coords
  latlong <- coord_info$latlong

  if (!is.null(cluster)) {
    if (!is.character(cluster) || length(cluster) != 1L || !nzchar(cluster)) {
      stop("`cluster` must be a single column name.")
    }
    if (!cluster %in% names(data)) {
      stop("Cluster variable not found in data: ", cluster)
    }
    if (!is.null(cond_fixef_id) && !isTRUE(uncond)) {
      stop("Conditional SCPC with absorbed fixed effects and external clustering is not yet implemented; use `uncond = TRUE`.")
    }
    cl_vec <- factor(data[[cluster]][obs_index])

    ## Warn if coordinates vary within clusters
    coord_by_cl <- rowsum(coords, cl_vec)
    n_per_cl    <- as.numeric(table(cl_vec))
    coord_means <- coord_by_cl / n_per_cl
    coord_first <- coords[match(levels(cl_vec), cl_vec), , drop = FALSE]
    if (max(abs(coord_means - coord_first)) > 1e-8) {
      warning("Coordinates vary within clusters. The first observation's ",
              "coordinates are used for each cluster; consider averaging ",
              "coordinates within clusters before calling scpc().")
    }

    S              <- rowsum(S, cl_vec)
    coords         <- coord_first
    neff           <- nrow(S)
  }

  if (ncol(coords) == 1) coords <- cbind(coords, 0)

  ## 2. Spatial kernel (unconditional) ---------------------------------------
  spc  <- .setOmsWfin(
    coords, avc, latlong,
    method = method,
    large_n_seed = large_n_seed
  )
  Wfin <- spc$Wfin; cvfin <- spc$cvfin; Omsfin <- spc$Omsfin
  perm <- spc$perm
  q    <- ncol(Wfin) - 1

  ## 3. Bread ----------------------------------------------------------------
  bread_inv <- sandwich::bread(model) / n

  ## 4. Loop over coefficients ------------------------------------------------
  k_use <- if (is.null(ncoef)) p else min(ncoef, p)
  out   <- matrix(NA_real_, k_use, 6)
  levs  <- c(0.32, 0.10, 0.05, 0.01)
  cvs_mat    <- if (cvs) matrix(NA_real_, k_use, 4) else NULL
  cvs_uncond <- if (cvs) vapply(levs, function(lv) .getcv(Omsfin, q, lv), 0.0) else NULL
  coef_names <- names(stats::coef(model))[seq_len(k_use)]
  rownames(out) <- coef_names
  colnames(out) <- c("Coef", "Std_Err", "t", "P>|t|", "2.5 %", "97.5 %")
  large_n_random_state <- spc$random_state

  for (j in seq_len(k_use)) {
    wj  <- as.numeric(neff * bread_inv[j, ] %*% t(S)) + stats::coef(model)[j]
    wj_perm <- wj[perm]

    ## unconditional statistic -----------------------------------------------
    tau_u <- as.numeric(sqrt(q) * crossprod(Wfin[, 1], wj_perm) /
      sqrt(sum((t(Wfin[, -1]) %*% wj_perm)^2)))
    SE    <- as.numeric(sqrt(sum((t(Wfin[, -1]) %*% wj_perm)^2)) / (sqrt(q) * sqrt(neff)))
    p_u   <- .maxrp(Omsfin, q, abs(tau_u) / sqrt(q))$max

    ## conditional branch (skip when uncond = TRUE) --------------------------
    if (!uncond) {
      if (!is.null(cluster)) {
        ## Clustered conditional: individual-level orthogonalization
        cl_vec_scpc <- if (identical(perm, seq_len(length(perm)))) {
          cl_vec
        } else {
          factor(as.character(cl_vec), levels = levels(cl_vec)[perm])
        }
        xj_indiv <- as.numeric(neff * bread_inv[j, ] %*% t(model_mat_cond))
        Wx <- .orthogonalize_W_cluster(
          Wfin, cl_vec_scpc, xj_indiv, model_mat_cond,
          include_intercept = cond_include_intercept
        )
      } else {
        ## Non-clustered conditional
        xj  <- as.numeric(neff * bread_inv[j, ] %*% t(model_mat_cond))[perm]
        xjs <- sign(xj)
        Wx  <- .orthogonalize_W(
          Wfin, xj, xjs, model_mat_cond[perm, , drop = FALSE],
          include_intercept = cond_include_intercept,
          fixef_id = cond_fixef_id
        )
      }
      if (spc$large_n) {
        omsx_res <- .lnget_Oms(
          spc$coords, spc$c0, spc$cmax, Wx, 1.2,
          capM = 1000000L, random_t = large_n_random_state, latlong = latlong
        )
        Omsx <- omsx_res$Oms
        large_n_random_state <- omsx_res$state
      } else {
        Omsx <- .getOms(spc$distmat, spc$c0, spc$cmax, Wx, 1.2)
      }
      p_c     <- .maxrp(Omsx, q, abs(tau_u) / sqrt(q))$max
      cvx     <- .getcv(Omsx, q, 0.05)
      p_final <- max(p_u, p_c)
      cv      <- max(cvfin, cvx)
    } else {
      p_final <- p_u
      cv      <- cvfin
    }

    out[j, ] <- c(stats::coef(model)[j], SE, tau_u, p_final,
                  stats::coef(model)[j] - cv * SE,
                  stats::coef(model)[j] + cv * SE)

    if (cvs) {
      cvs_vec <- cvs_uncond
      if (!uncond) {
        cvs_cond <- vapply(levs, function(lv) .getcv(Omsx, q, lv), 0.0)
        cvs_vec  <- pmax(cvs_vec, cvs_cond)
      }
      cvs_mat[j, ] <- cvs_vec
    }
  }

  ## Label the cvs matrix ----------------------------------------------------
  if (!is.null(cvs_mat)) {
    rownames(cvs_mat) <- coef_names
    colnames(cvs_mat) <- c("32%", "10%", "5%", "1%")
  }

  structure(list(
    scpcstats = out,
    scpccvs   = cvs_mat,
    W         = Wfin,
    avc       = avc,
    c0        = spc$c0,
    cv        = cvfin,
    q         = q,
    method    = spc$method,
    large_n_seed = large_n_seed,
    call      = match.call()
  ), class = "scpc")
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' @rdname scpc
#' @param x An object of class \code{"scpc"}.
#' @param \dots Further arguments (currently unused).
#' @export
print.scpc <- function(x, ...) {
  cat("\nSCPC Inference (ncoef =", nrow(x$scpcstats), ", q =", x$q, ")\n\n")
  stats::printCoefmat(x$scpcstats[, 1:4], P.values = TRUE, has.Pvalue = TRUE)
  if (!is.null(x$scpccvs)) {
    cat("\nTwo-sided critical values:\n")
    print(x$scpccvs)
  }
  invisible(x)
}

#' @rdname scpc
#' @param object An object of class \code{"scpc"}.
#' @export
summary.scpc <- function(object, ...) {
  cat("\nSCPC Inference (ncoef =", nrow(object$scpcstats), ", q =", object$q,
      ", avc =", object$avc, ")\n\n")
  stats::printCoefmat(object$scpcstats[, 1:4, drop = FALSE],
                      P.values = TRUE, has.Pvalue = TRUE,
                      signif.stars = TRUE)
  cat("\n95% Confidence Intervals:\n")
  print(object$scpcstats[, 5:6, drop = FALSE])
  if (!is.null(object$scpccvs)) {
    cat("\nTwo-sided critical values:\n")
    print(object$scpccvs)
  }
  invisible(object)
}

#' @rdname scpc
#' @export
coef.scpc <- function(object, ...) {
  out <- object$scpcstats[, "Coef"]
  names(out) <- rownames(object$scpcstats)
  out
}

#' @rdname scpc
#' @param parm Character vector of coefficient names, or numeric indices.
#'   Default is all coefficients.
#' @param level Confidence level.  Only 0.95 is available unless
#'   \code{cvs = TRUE} was used, in which case 0.68, 0.90, 0.95, and 0.99
#'   are supported.
#' @export
confint.scpc <- function(object, parm = NULL, level = 0.95, ...) {
  st  <- object$scpcstats
  nms <- rownames(st)

  if (is.null(parm)) {
    idx <- seq_len(nrow(st))
  } else if (is.character(parm)) {
    idx <- match(parm, nms)
    if (anyNA(idx)) stop("Unknown coefficient(s): ", paste(parm[is.na(idx)], collapse = ", "))
  } else {
    idx <- parm
  }

  avail <- c(0.68, 0.90, 0.95, 0.99)
  avail_idx <- match(TRUE, abs(avail - level) < 1e-8)

  if (abs(level - 0.95) < 1e-8) {
    ci <- st[idx, c("2.5 %", "97.5 %"), drop = FALSE]
  } else if (!is.null(object$scpccvs) && !is.na(avail_idx)) {
    cv_vals <- object$scpccvs[idx, avail_idx]
    se_vals <- st[idx, "Std_Err"]
    co_vals <- st[idx, "Coef"]
    ci <- cbind(co_vals - cv_vals * se_vals, co_vals + cv_vals * se_vals)
  } else {
    avail_msg <- "0.95"
    if (!is.null(object$scpccvs)) avail_msg <- paste(c(avail_msg, "0.68, 0.90, 0.99"), collapse = ", ")
    stop("Confidence level ", level, " is not available. Available levels: ", avail_msg)
  }

  pct <- format(100 * c((1 - level) / 2, 1 - (1 - level) / 2), digits = 3)
  colnames(ci) <- paste0(pct, " %")
  rownames(ci) <- nms[idx]
  ci
}
