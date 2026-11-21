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
    distances <- stats::dist(S)
    D <- as.matrix(distances)
    if (!all(is.finite(distances) & distances >= 1e-150)) {
      for (i in seq_len(nrow(S) - 1L)) {
        j <- (i + 1L):nrow(S)
        j <- j[!is.finite(D[i, j]) | D[i, j] < 1e-150]
        if (!length(j)) next
        base <- S[rep(i, length(j)), , drop = FALSE]
        D[i, j] <- D[j, i] <- .getdistvec(S[j, , drop = FALSE], base, latlong = FALSE)
      }
    }
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
    a <- sin(dlat)^2 + cos(lat1) * cos(lat2) * sin(dlon)^2
    asin(sqrt(pmin(pmax(a, 0), 1))) / pi
  } else {
    delta <- S1 - S2
    d <- sqrt(rowSums(delta^2))
    fix <- !is.finite(d) | d < 1e-150
    if (any(fix)) {
      small <- delta[fix, , drop = FALSE]
      scale <- apply(abs(small), 1L, max)
      if (any(!is.finite(scale))) stop("Pairwise distance exceeds the floating-point range.")
      d[fix] <- scale * sqrt(rowSums((small / ifelse(scale > 0, scale, 1))^2))
      if (any(!is.finite(d))) stop("Pairwise distance exceeds the floating-point range.")
    }
    d
  }
}

.getc0fromavc <- function(dist, avc0) {
  ## Solve for exponential kernel parameter c0 such that
  ## mean(exp(-c0 * dist)) == avc0, using a scale-invariant bracket
  ## and Illinois / bisection root solver.
  dist <- as.numeric(dist)
  n_dist <- length(dist)
  if (n_dist == 0L) {
    stop("At least two locations are required.")
  }
  if (any(!is.finite(dist) | dist < 0)) {
    stop("Distances must be finite and nonnegative.")
  }
  floor_corr <- mean(dist == 0)
  pos <- dist[dist > 0]
  if (length(pos) == 0L) {
    stop("At least two distinct locations are required.")
  }
  if (!is.numeric(avc0) || length(avc0) != 1L || !is.finite(avc0) || avc0 <= floor_corr || avc0 >= 1.0) {
    stop("Target correlation must lie above the co-located-pair fraction and below one.")
  }

  mean_pos <- mean(pos)
  target_pos <- (avc0 - floor_corr) / (1 - floor_corr)
  near_one <- avc0 > 0.9999
  eval_res <- function(c) {
    if (near_one) {
      (1 - avc0) + sum(expm1(-c * pos)) / n_dist
    } else {
      mean(exp(-c * dist)) - avc0
    }
  }

  raw_init <- if (is.finite(mean_pos) && mean_pos > 0) {
    cand <- -log(target_pos) / mean_pos
    if (is.finite(cand) && cand > 0) cand else NA_real_
  } else {
    NA_real_
  }
  if (is.na(raw_init)) {
    max_d <- max(pos)
    raw_init <- min(max(.Machine$double.xmin, 1 / max_d), .Machine$double.xmax * 0.5)
  }

  c0 <- raw_init
  c1 <- raw_init
  f0 <- eval_res(c0)
  if (f0 == 0) return(c0)
  f1 <- f0

  if (f0 < 0) {
    c1 <- c0
    f1 <- f0
    for (iter in seq_len(1024L)) {
      c0 <- 0.5 * c0
      if (c0 == 0 || !is.finite(c0)) {
        stop("Decay calibration lower bracket underflowed; distances may be degenerate.")
      }
      f0 <- eval_res(c0)
      if (f0 == 0) return(c0)
      if (f0 > 0) break
      c1 <- c0
      f1 <- f0
    }
    if (f0 < 0) {
      stop("Decay calibration lower bracket did not enclose root.")
    }
  } else {
    for (iter in seq_len(1024L)) {
      c0 <- c1
      f0 <- f1
      if (c1 >= .Machine$double.xmax * 0.5) {
        c1 <- .Machine$double.xmax
        f1 <- eval_res(c1)
        break
      }
      c1 <- 2 * c1
      f1 <- eval_res(c1)
      if (f1 == 0) return(c1)
      if (f1 < 0) break
    }
    if (f1 > 0) {
      stop("Spatial decay calibration failed to bracket upper decay scale.")
    }
  }

  scale_gap <- max(min(avc0 - floor_corr, 1 - avc0), 1e-15)
  residual_tol <- 1e-12 * scale_gap
  side <- 0L
  cmid <- c0 + 0.5 * (c1 - c0)
  for (iter in seq_len(200L)) {
    fallback_mid <- if (c0 > 0 && c1 / c0 > 4) {
      geom <- sqrt(c0 * c1)
      if (is.finite(geom) && geom > c0 && geom < c1) geom else c0 + 0.5 * (c1 - c0)
    } else {
      c0 + 0.5 * (c1 - c0)
    }
    if (iter %% 4L == 0L && c0 > 0 && c1 / c0 > 4) {
      cmid <- fallback_mid
    } else {
      cmid <- c0 - f0 * ((c1 - c0) / (f1 - f0))
      if (!is.finite(cmid) || cmid <= c0 || cmid >= c1) {
        cmid <- fallback_mid
      }
    }
    fmid <- eval_res(cmid)
    if (abs(fmid) <= residual_tol) {
      return(cmid)
    }
    if (fmid > 0) {
      c0 <- cmid
      f0 <- fmid
      if (side == -1L) f1 <- 0.5 * f1
      side <- -1L
    } else {
      c1 <- cmid
      f1 <- fmid
      if (side == 1L) f0 <- 0.5 * f0
      side <- 1L
    }
    c_est <- c0 + 0.5 * (c1 - c0)
    if ((c1 - c0) <= 1e-12 * max(c_est, .Machine$double.xmin)) {
      return(c_est)
    }
  }
  c_est <- c0 + 0.5 * (c1 - c0)
  if (abs(eval_res(c_est)) <= 1e-8 * scale_gap || (c1 - c0) <= 1e-12 * max(c_est, .Machine$double.xmin)) {
    return(c_est)
  }
  stop("Spatial decay calibration did not converge.")
}

.angular_deficit <- function(sine2, denominator, normalized_eigenvalues) {
  product <- rep.int(1.0, length(sine2))
  deficit <- rep.int(0.0, length(sine2))
  for (lam in normalized_eigenvalues) {
    extra <- lam * denominator
    inv <- 1.0 / (sine2 + extra)
    deficit <- deficit + product * extra * inv
    product <- product * sine2 * inv
  }
  deficit / ((1.0 + sqrt(product)) * denominator)
}

.angular_integrand <- function(cosine2, eigenvalues) {
  if (cosine2 <= 0) return(0.0)
  large <- eigenvalues > cosine2
  log_terms <- log1p(eigenvalues[!large] / cosine2)
  large_evals <- eigenvalues[large]
  exp(-0.5 * (sum(log_terms) + sum(log(large_evals) - log(cosine2) + log1p(cosine2 / large_evals))))
}

#' @importFrom gaussquad legendre.quadrature.rules
.GQ <- gaussquad::legendre.quadrature.rules(40)[[40]]
.GQx <- .GQ[, 1] * 0.5 + 0.5
.GQw <- .GQ[, 2] * 0.5

.bivariate_rejection <- function(a, b, d, cv) {
  # The signed product (X - cv Y)(X + cv Y) has an exact Gaussian quadrant probability.
  tol <- 64 * .Machine$double.eps
  magnitude <- max(abs(c(a, b, d)))
  if (magnitude == 0) return(0)
  v <- c(a, d) / magnitude + tol
  if (min(v) < 0 || abs(b) / magnitude > prod(sqrt(v))) stop("SCPC covariance is not positive semidefinite.")
  a <- max(0, a)
  d <- max(0, d)
  if (min(a, d) == 0) return(as.numeric(a > 0))
  # Exact products in scaled units distinguish a tiny determinant from a tie.
  e <- pmax(-1022, pmin(1023, floor(log2(c(a, d)))))
  eb <- floor(sum(e) / 2)
  an <- a / 2^e[1L]
  dn <- (d / 2^e[2L]) * 2^(sum(e) - 2*eb)
  bn <- b / 2^eb
  ad <- an * dn
  bb <- bn * bn
  det <- if (bb > 2*ad) 0 else max(0, (ad - bb) +
    (.scpc_product_tail(an,dn,ad) - .scpc_product_tail(bn,bn,bb)))
  t <- cv * 2^(eb - e[1L])
  if (t > 2^500) {
    spread <- 2 * sqrt(det) / t
    difference <- dn
  } else {
    td <- t * dn
    ttd <- t * td
    spread <- 2 * t * sqrt(det)
    difference <- (ttd - an) + (t * .scpc_product_tail(t,dn,td) + .scpc_product_tail(td,t,ttd))
  }
  if (spread == 0) return(as.numeric(difference < 0))
  atan2(spread, difference) / pi
}

# Dekker's error-free product. Inputs have already been scaled by a power of two.
.scpc_product_tail <- function(x, y, p) {
  sx <- 134217729 * x
  sy <- 134217729 * y
  xh <- sx - (sx - x)
  yh <- sy - (sy - y)
  xl <- x - xh
  yl <- y - yh
  xl * yl - (((p - xh * yh) - xl * yh) - xh * yl)
}

.one_factor_rejection <- function(Om, cv) {
  if (Om[1L, 1L] <= 0 || Om[1L, 1L] * Om[2L, 2L] != Om[2L, 1L]^2) return(NULL)
  scale <- 2^max(-1022, min(1023, floor(log2(max(diag(Om))))))
  A <- Om / scale
  if (any(A * scale != Om) || any(Om != t(Om))) return(NULL)
  # Rank one is certified by vanishing minors. Keep every representable positive direction.
  a <- A[1L, 1L]
  rest <- A[-1L, -1L, drop = FALSE]
  left <- A[-1L, 1L]
  right <- rep(left, each = length(left))
  p <- a * rest
  if (any(p != left * right) || any(p == 0 & (rest != 0 | (left != 0 & right != 0))) ||
      any(p != 0 & abs(p) < 2^-918)) return(NULL)
  if (any(.scpc_product_tail(a, rest, p) != .scpc_product_tail(left, right, p))) return(NULL)
  denominator <- sum(diag(A)[-1L])
  as.numeric(if (cv > 1) (a / cv) / cv > denominator else a > (denominator * cv) * cv)
}

.rejection_eigenvalues <- function(Om, cv, tolerance) {
  # Factor Om = L L' before forming the symmetric weighted quadratic form.
  v <- diag(Om)
  if (min(v) < -tolerance) stop("SCPC covariance is not positive semidefinite.")
  sd <- sqrt(pmax(0, c(v[1L], rep(max(v[-1L]), length(v) - 1L))))
  zero <- sd == 0
  divisor <- sd
  divisor[zero] <- 1
  corr <- sweep(Om / divisor, 2L, divisor, "/")
  corr[zero, ] <- 0
  corr[, zero] <- 0
  if (any(abs(Om[zero, , drop = FALSE]) > tolerance) || any(!is.finite(corr))) {
    scale <- max(0, v)
    if (scale == 0) stop("SCPC covariance is not positive semidefinite.")
    corr <- Om / scale
    sd[] <- sqrt(scale)
  }
  L <- tryCatch(t(chol(corr)), error = function(e) NULL)
  if (is.null(L)) {
    ev <- eigen(corr, symmetric = TRUE)
    if (min(ev$values) < -32 * length(v) * .Machine$double.eps) {
      scale <- max(0, v)
      if (scale == 0) stop("SCPC covariance is not positive semidefinite.")
      ev <- eigen(Om / scale, symmetric = TRUE)
      if (min(ev$values) < -32 * length(v) * .Machine$double.eps) stop("SCPC covariance is not positive semidefinite.")
      sd[] <- sqrt(scale)
    }
    keep <- ev$values > 0
    if (!any(keep)) return(0)
    L <- sweep(ev$vectors[, keep, drop = FALSE], 2L, sqrt(ev$values[keep]), "*")
  }
  if (cv > 1) sd[1L] <- sd[1L] / cv else sd[-1L] <- cv * sd[-1L]
  if (max(sd) == 0) return(0)
  W <- (sd / max(sd)) * L
  eigen(crossprod(W[1L, , drop = FALSE]) - crossprod(W[-1L, , drop = FALSE]),
        symmetric = TRUE, only.values = TRUE)$values
}

.getrp <- function(Om, cv) {
  ## Rejection probability for the SCPC test given Omega matrix Om
  ## and normalised critical value cv (= cv/sqrt(q) in paper notation).
  if (!is.finite(cv)) {
    stop("SCPC critical value must be finite.")
  }
  if (cv < 0) return(1.0)
  if (!is.matrix(Om) || nrow(Om) != ncol(Om) || nrow(Om) < 2L || any(!is.finite(Om))) {
    stop("SCPC covariance must be finite and square with at least two rows.")
  }
  tolerance <- 32 * nrow(Om) * .Machine$double.eps * max(0, diag(Om))
  if (any(abs(Om - t(Om)) > tolerance)) stop("SCPC covariance is not symmetric.")
  if (nrow(Om) > 2L) {
    active <- which(diag(Om)[-1L] != 0) + 1L
    if (length(active) <= 1L) {
      j <- if (length(active)) active else 2L
      keep <- c(1L, j)
      if (all(Om[-keep, ] == 0) && all(Om[, -keep] == 0)) Om <- Om[keep, keep]
    }
  }
  if (nrow(Om) == 2L) {
    b <- Om[1L, 2L] + 0.5 * (Om[2L, 1L] - Om[1L, 2L])
    return(.bivariate_rejection(Om[1L, 1L], b, Om[2L, 2L], cv))
  }
  atom <- .one_factor_rejection(Om, cv)
  if (!is.null(atom)) return(atom)
  evals_raw <- .rejection_eigenvalues(Om, cv, tolerance)
  denom <- max(evals_raw)
  if (!is.finite(denom) || denom <= 0) return(0.0)
  evals <- -evals_raw[evals_raw < 0]
  if (length(evals) == 0L) return(1.0)
  evals <- evals / denom
  if (length(evals) == 1L) {
    return(as.numeric((2 / pi) * atan2(1, sqrt(evals[1L]))))
  }
  if (max(evals) - min(evals) <= 1e-12 * min(evals)) {
    q_eff <- length(evals)
    return(as.numeric(2 * stats::pt(-sqrt(q_eff * evals[1L]), df = q_eff)))
  }
  scale2 <- sum(evals)
  if (scale2 == 0) return(1.0)
  complement <- (scale2 < 0.0625)
  if (complement) evals <- evals / scale2

  eval_panel <- function(left, right) {
    angle <- left + (right - left) * .GQx
    cosine2 <- cos(angle)^2
    if (complement) {
      sine2 <- sin(angle)^2
      vals <- .angular_deficit(sine2, cosine2 + scale2 * sine2, evals)
    } else {
      product <- rep.int(1.0, length(angle))
      for (lam in evals) {
        product <- product * (cosine2 / (cosine2 + lam))
      }
      vals <- sqrt(product)
      tiny <- product < .Machine$double.xmin
      if (any(tiny)) {
        vals[tiny] <- vapply(cosine2[tiny], .angular_integrand, 0.0, eigenvalues = evals)
      }
    }
    sum(.GQw * vals) * (right - left)
  }

  end <- pi / 2
  whole <- eval_panel(0, end)
  base_tol <- if (abs(whole) < 1e-8) 1e-10 * abs(whole) else 1e-11
  pending <- list(list(left = 0, right = end, whole = whole, tol = base_tol, depth = 0L))
  total <- 0.0
  while (length(pending) > 0L) {
    cur <- pending[[length(pending)]]
    pending[[length(pending)]] <- NULL
    mid <- 0.5 * (cur$left + cur$right)
    left_val <- eval_panel(cur$left, mid)
    right_val <- eval_panel(mid, cur$right)
    if (cur$depth >= 20L || abs(left_val + right_val - cur$whole) <= max(cur$tol, 1e-11 * abs(cur$whole))) {
      total <- total + left_val + right_val
    } else {
      pending[[length(pending) + 1L]] <- list(
        left = mid, right = cur$right, whole = right_val, tol = cur$tol * 0.5, depth = cur$depth + 1L
      )
      pending[[length(pending) + 1L]] <- list(
        left = cur$left, right = mid, whole = left_val, tol = cur$tol * 0.5, depth = cur$depth + 1L
      )
    }
  }
  prob <- if (complement) 1.0 - sqrt(scale2) * total * (2 / pi) else total * (2 / pi)
  as.numeric(max(0.0, min(1.0, prob)))
}

.maxrp <- function(Oms, q, cv) {
  ## Largest rejection probability across a list of Omega matrices
  rps <- vapply(Oms, function(Om) .getrp(Om[1:(q + 1), 1:(q + 1), drop = FALSE], cv), 0.0)
  list(max = max(rps), i = which.max(rps))
}

.getcv <- function(Oms, q, level) {
  ## Two-sided critical value that controls size at the given level,
  ## maximised over the Omega grid.
  step <- 1 / sqrt(q)
  tol <- 0.001 * step
  rp <- 1
  i <- 1L
  cv0 <- stats::qt(1 - level / 2, df = q) * step
  for (outer_iter in seq_len(length(Oms) + 4L)) {
    if (rp <= level) break
    Om_i <- Oms[[i]][1:(q + 1), 1:(q + 1), drop = FALSE]
    cv1 <- cv0
    for (iter in seq_len(1024L)) {
      if (.getrp(Om_i, cv1) > level) {
        cv0 <- cv1
        cv1 <- cv1 + if (iter <= 32L) step else cv1
      } else {
        break
      }
    }
    if (.getrp(Om_i, cv1) > level) {
      stop("SCPC critical value could not be bracketed.")
    }
    while (cv1 - cv0 > tol) {
      cv <- 0.5 * (cv0 + cv1)
      if (cv == cv0 || cv == cv1) break
      if (.getrp(Om_i, cv) > level) cv0 <- cv else cv1 <- cv
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
    j <- i + floor(nxt$value * (n - i + 1L))
    if (j != i) {
      tmp <- S[j, ]
      S[j, ] <- S[i, ]
      S[i, ] <- tmp
    }
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

  floor_corr <- mean(distvec == 0)
  minavc_eff <- floor_corr + (1 - floor_corr) * minavc
  c0 <- .getc0fromavc(distvec, avc0)
  cmax <- .getc0fromavc(distvec, minavc_eff)

  Wall <- matrix(0, n, capN * qmax)
  for (i in seq_len(capN)) {
    W0 <- .ln_subset_evecs(ms[[i]]$distmat, c0, qmax)
    Wx <- matrix(0, n, qmax)
    for (j in seq_len(m)) {
      d_cross <- .getdistvec(
        S,
        matrix(ms[[i]]$coords[j, ], nrow = n, ncol = ncol(S), byrow = TRUE),
        latlong
      )
      v <- exp(-c0 * d_cross)
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
  min_norm <- sqrt(.Machine$double.eps)
  for (i in seq_len(qmax)) {
    Wx <- Wall[, seq_len(capN * i), drop = FALSE]
    evec <- eigen(crossprod(Wx), symmetric = TRUE)$vectors[, 1, drop = FALSE]
    W[, i] <- as.numeric(Wx %*% evec)
    norm_i <- sqrt(sum(W[, i]^2))
    if (!is.finite(norm_i) || norm_i <= min_norm) {
      stop("degenerate principal direction in Nystrom landmark consolidation", call. = FALSE)
    }
    W[, i] <- W[, i] / norm_i
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

.check_sampled_covariance <- function(omega) {
  rel_tol <- 32 * nrow(omega) * .Machine$double.eps
  scale <- max(abs(omega))
  if (!is.finite(scale) ||
      (scale > 0 && min(eigen(omega / scale, symmetric = TRUE, only.values = TRUE)$values) < -rel_tol)) {
    stop(errorCondition(
      'Invalid sampled covariance in approximate SCPC. Use method = "exact" if feasible.',
      class = "spatial_invalid_sampled_covariance"
    ))
  }
}

.scpc_pair_sample <- function(S, capM, random_t, latlong) {
  draw <- .raninds(nrow(S), capM, random_t)
  idx <- draw$indices
  draw$dist <- .getdistvec(S[idx[seq_len(capM)], , drop = FALSE],
                           S[idx[seq_len(capM) + 1L], , drop = FALSE], latlong)
  draw
}

.lnget_Oms <- function(S, c0, cmax, W, cgridfac,
                       capM = 1000000L, random_t = 1, latlong = FALSE,
                       pairs = .scpc_pair_sample(S, capM, random_t, latlong)) {
  nc <- .getnc(c0, cmax, cgridfac)
  Oms <- vector("list", nc)
  n <- nrow(S)
  inds <- pairs$indices
  dist <- pairs$dist
  capM <- length(dist)
  W1 <- W[inds[seq_len(capM)], , drop = FALSE]
  W2 <- W[inds[2:(capM + 1L)], , drop = FALSE]

  gram <- crossprod(W)
  c <- c0
  for (i in seq_len(nc)) {
    cd <- exp(-c * dist)
    Oms[[i]] <- gram +
      0.5 * (n * (n - 1) / capM) *
      (crossprod(W1, W2 * cd) + crossprod(W2, W1 * cd))
    .check_sampled_covariance(Oms[[i]])
    c <- c * cgridfac
  }

  list(Oms = Oms, state = pairs$state)
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
  if (n < 2L) {
    stop("At least two locations are required.")
  }

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
  tied <- FALSE

  if (identical(method_actual, "exact")) {
    distmat <- .getdistmat(coords, latlong)
    distv <- .lvech(distmat)
    floor_corr <- mean(distv == 0)
    minavc_eff <- floor_corr + (1 - floor_corr) * minavc
    c0 <- .getc0fromavc(distv, avc0)
    cmax <- .getc0fromavc(distv, minavc_eff)
    coords_use <- coords
    perm <- seq_len(n)
    random_t <- NULL
  } else {
    norm_s <- .normalize_s(coords, latlong)
    coords_use <- norm_s$coords
    perm <- norm_s$perm
  }

  q_limit <- if (identical(method_actual, "approx")) min(large_n_m, n) - 1L else n - 1L
  repeat {
    qmax <- min(qmax, q_limit)
    if (identical(method_actual, "exact")) {
      W <- .getW(distmat, c0, qmax)
      tied <- isTRUE(attr(W, "tied"))
      attr(W, "tied") <- NULL
      Oms <- .getOms(distmat, c0, cmax, W, cgridfac)
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
      tied <- FALSE
    }
    fin <- .setfinalW(Oms, W, qmax)
    if (fin$q < qmax || qmax >= q_limit) break
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
    random_state = random_t,
    tied = isTRUE(tied)
  )
}

.get_conditional_Oms <- function(spc, Wx, pairs) {
  if (isTRUE(spc$large_n)) {
    return(.lnget_Oms(spc$coords, spc$c0, spc$cmax, Wx, 1.2, pairs = pairs)$Oms)
  }
  .getOms(spc$distmat, spc$c0, spc$cmax, Wx, 1.2)
}
