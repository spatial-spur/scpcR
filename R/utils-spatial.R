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

.setOmsWfin <- function(distmat, avc0) {
  n <- nrow(distmat)
  distv <- .lvech(distmat)

  cgridfac <- 1.2
  minavc <- 0.00001

  if (avc0 >= 0.05) {
    qmax <- 10
  } else if (avc0 >= 0.01) {
    qmax <- 20
  } else if (avc0 >= 0.005) {
    qmax <- 60
  } else {
    qmax <- 120
  }

  c0 <- .getc0fromavc(distv, avc0)
  cmax <- .getc0fromavc(distv, minavc)

  repeat {
    qmax <- min(qmax, n - 1)
    W <- .getW(distmat, c0, qmax)
    Oms <- .getOms(distmat, c0, cmax, W, cgridfac)
    fin <- .setfinalW(Oms, W, qmax)
    if (fin$q < qmax || qmax == n - 1) break
    qmax <- round(qmax + qmax / 2)
  }
  list(
    Wfin = fin$W,
    cvfin = fin$cv,
    Omsfin = Oms,
    c0 = c0,
    cmax = cmax
  )
}
