# ---------------------------------------------------------------------------
# Spatial helpers
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
    if (!requireNamespace("geodist", quietly = TRUE)) {
      stop("Install package 'geodist' for lat/long distances.")
    }
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

#' @importFrom gaussquad legendre.quadrature.rules
.GQ <- gaussquad::legendre.quadrature.rules(40)[[40]]
.GQx <- .GQ[, 1] * 0.5 + 0.5
.GQw <- .GQ[, 2] * 0.5

.getrp <- function(Om, cv) {
  ## Rejection probability for the SCPC test given Omega matrix Om
  ## and normalised critical value cv (= cv/sqrt(q) in paper notation).
  ## Uses Gauss-Legendre quadrature over the characteristic function.
  Omx <- -cv^2 * Om
  Omx[1, ] <- Om[1, ]
  evals_raw <- Re(eigen(Omx, only.values = TRUE)$values)
  evals <- -evals_raw[evals_raw < 0]
  if (!length(evals)) return(0)
  denom <- max(evals_raw)
  if (!is.finite(denom) || denom <= 0) return(0)
  evals <- evals / denom
  tot <- 0
  for (j in seq_along(.GQx)) {
    u <- .GQx[j]
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
  i <- 1
  cv0 <- stats::qt(1 - level / 2, df = q) / sqrt(q)
  while (rp > level) {
    cv1 <- cv0
    repeat {
      if (.getrp(Oms[[i]][1:(q + 1), 1:(q + 1)], cv1) > level) {
        cv0 <- cv1
        cv1 <- cv1 + 1 / sqrt(q)
      } else {
        break
      }
    }
    while (cv1 - cv0 > 0.001 / sqrt(q)) {
      cv <- 0.5 * (cv0 + cv1)
      if (.getrp(Oms[[i]][1:(q + 1), 1:(q + 1)], cv) > level) cv0 <- cv else cv1 <- cv
    }
    maxrp <- .maxrp(Oms, q, cv1)
    if (maxrp$i == i) break
    i <- maxrp$i
    rp <- maxrp$max
    cv0 <- cv1
  }
  cv1 * sqrt(q)
}

.setfinalW <- function(Oms, W, qmax) {
  ## Select the optimal q (number of spatial PCs) that minimises the
  ## expected confidence interval length under i.i.d. errors.
  cvs <- lengths <- numeric(qmax)
  for (q in seq_len(qmax)) {
    cvs[q] <- .getcv(Oms, q, 0.05)
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
  nc <- .getnc(c0, cmax, cgridfac)
  Oms <- vector("list", nc)
  c <- c0
  for (i in seq_len(nc)) {
    Oms[[i]] <- crossprod(W, exp(-c * distmat) %*% W)
    c <- c * cgridfac
  }
  Oms
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
  dist <- .getdistvec(
    S[inds[seq_len(capM)], , drop = FALSE],
    S[inds[2:(capM + 1L)], , drop = FALSE],
    latlong
  )
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
  n <- nrow(coords)

  cgridfac <- 1.2
  minavc <- 0.00001

  large_n_threshold <- 4500L
  large_n_capN <- 20L
  large_n_capM <- 1000000L
  large_n_m <- 1000L
  method_actual <- .resolve_scpc_method(n, method, large_n_threshold = large_n_threshold)

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
