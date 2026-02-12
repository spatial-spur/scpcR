
.getavc <- function(c, dist) {
  ## computes average correlation given vector of pairwise distances
  mean(exp(-c * dist))
}

.getdistmat <- function(S, latlong) {
  ## computes matrix of distances from locations
  ## if external latlong=0, distances computed as norm, otherwise latitude / longitude 
  ## and haversine formula for sphere with radius 1/Pi
  S <- as.matrix(S)
  if(latlong) {
    if(ncol(S) != 2L) {
      stop("Internal error: geodesic coordinates must have exactly two columns.")
    }
    if(!requireNamespace("geodist", quietly = TRUE))
      stop("Install package 'geodist' for lat/long distances.")
    S <- as.data.frame(S)
    names(S) <- c("lon", "lat")
    D <- geodist::geodist(S, measure = "haversine") / (2 * pi * 6378137)  # units pi*R_earth
  } else {
    D <- as.matrix(stats::dist(S, diag = TRUE, upper = TRUE))  # Euclidean
  }
  D
}

.lvech <- function(mat) {
  ## takes lower triangular part of matrix and puts it into vector
  mat[lower.tri(mat)]
}

.demeanmat <- function(mat) {
  ## demeans matrix
  mat <- mat - rowMeans(mat)
  mat <- mat - matrix(colMeans(mat), nrow(mat), ncol(mat), byrow = TRUE)
  return(mat)
}

.getc0fromavc <- function(dist, avc0) {
  ## solves for c0 from vector of pairwise distances given avc0
  c0 <- c1 <- 10
  while(.getavc(c0, dist) < avc0) {
    c1 <- c0
    c0 <- 0.5 * c0
  }
  while(.getavc(c1,dist) > avc0 & c1 < 5000){
    c0 <- c1
    c1 <- 2 * c1
  }
  repeat{
    c <- sqrt(c0 * c1)
    if(.getavc(c, dist) > avc0) c0 <- c else c1 <- c
    if(c1 - c0 < 0.001) break
  }
  return(c)
}

.getW <- function(distmat, c0, qmax) {
  ## computes first qmax eigenvectors of demeaned Sigma(c0) and stores them in W, 
  ## along with column vector of ones
  ## all columns of W are normalized to length 1
  n <- nrow(distmat)
  Sig <- exp(-c0 * distmat)
  Sig_d <- .demeanmat(Sig)
  if(requireNamespace("RSpectra", quietly = TRUE) && qmax < n - 1) {
    eig <- RSpectra::eigs_sym(Sig_d, k = qmax, which = "LM")
    V   <- eig$vectors
  } else {
    V <- eigen(Sig_d, symmetric = TRUE)$vectors[, seq_len(qmax)]
  }
  return(cbind(rep(1, n) / sqrt(n), V))
}

.GQ <- gaussquad::legendre.quadrature.rules(40)[[40]]  # 40-point table
.GQx <- .GQ[, 1] * 0.5 + 0.5  # nodes, transformed to [0,1]
.GQw <- .GQ[, 2] * 0.5        # weights, transformed to [0,1]

.getrp <- function(Om, cv) {
  ## computes rejection probability given Om and cv; cv here is cv/sqrt(q) in notation of paper
  Omx        <- -cv^2 * Om
  Omx[1,  ]  <- Om[1, ]               # restore first row
  evals_raw  <- Re(eigen(Omx, only.values = TRUE)$values)
  evals      <- -evals_raw[evals_raw < 0]
  if(!length(evals)) return(0)
  denom <- max(evals_raw)
  if(!is.finite(denom) || denom <= 0) return(0)
  evals <- evals / denom
  tot        <- 0
  for(j in seq_along(.GQx)) {
    u   <- .GQx[j]
    den <- sqrt((1 - u^2) * exp(sum(log1p(evals / (1 - u^2)))))
    tot <- tot + .GQw[j] / den
  }
  return(as.numeric(tot * 2 / pi))              # final probability
}

.maxrp <- function(Oms, q, cv) {
  ## computes largest rejection probability given vector of Om matrices stored in Oms given q and cv
  rps <- vapply(Oms, function(Om) .getrp(Om[1:(q + 1), 1:(q + 1)], cv), 0.0)
  return(list(max = max(rps), i = which.max(rps)))  # max rp and index
}

.getcv <- function(Oms, q, level) {
  ## computes (two-sided) critical value from q and Oms of given level
  rp <- 1
  i <- 1
  cv0 <- stats::qt(1 - level/2, df = q) / sqrt(q)  
  while(rp>level) {
    cv1 <- cv0
    repeat {
      if(.getrp(Oms[[i]][1:(q + 1), 1:(q + 1)], cv1) > level) {
        cv0 <- cv1
        cv1 <- cv1 + 1 / sqrt(q)
      } else break
    }
    while(cv1 - cv0 > 0.001 / sqrt(q)) {
      cv <- 0.5 * (cv0 + cv1)
      if(.getrp(Oms[[i]][1:(q + 1), 1:(q + 1)], cv) > level) cv0 <- cv else cv1 <- cv
    }
    maxrp <- .maxrp(Oms, q, cv1)
    if (maxrp$i == i) break
    i <- maxrp$i
    rp <- maxrp$max
    cv0 <- cv1
  }
  return(cv1 * sqrt(q))
}

.setfinalW <- function(Oms, W, qmax) {
  ## solves for optimal q and cv from Oms, stores results in W and cv
  cvs     <- lengths <- numeric(qmax)
  for(q in seq_len(qmax)) {
    cvs[q]     <- .getcv(Oms, q, 0.05)
    lengths[q] <- cvs[q] * gamma(0.5 * (q + 1)) / (sqrt(q) * gamma(0.5 * q))
  }
  q_opt <- which.min(lengths)
  list(W = W[, 1:(q_opt + 1), drop = FALSE], cv = cvs[q_opt], q = q_opt)
}

.getnc <- function(c0, cmax, cgridfac) {
  ## computes number of c-values in grid so that largest c is at least cmax
  max(2, ceiling(log(cmax/c0) / log(cgridfac)))
}

.getOms <- function(distmat, c0, cmax, W, cgridfac) {
  ## computes vector of qmax x qmax Om(c) matrices from distmat, c0 and W
  nc  <- .getnc(c0, cmax, cgridfac)
  Oms <- vector("list", nc)
  c   <- c0
  for(i in seq_len(nc)) {
    Oms[[i]] <- crossprod(W, exp(-c * distmat) %*% W)  # t(W) Sigma(c) W
    c <- c * cgridfac
  }
  return(Oms)
}

.gettau <- function(y, W) {
  ## computes SCPC t-stat
  return(sqrt(ncol(W)-1) * crossprod(W[, 1], y) / norm(crossprod(W[,-1], y), type="2"))
}


# Orthogonalisation helper (needed for conditional version)
.orthogonalize_W <- function(W, xj, xjs, model_mat) {
  Wx <- W
  Wx[, 1] <- Wx[, 1] * xj * xjs

  if (ncol(Wx) > 1L) {
    X <- cbind("(Intercept)" = 1, as.matrix(model_mat))
    qrX <- qr(X)
    Rx <- sweep(Wx[, -1, drop = FALSE], 1, xjs, `*`)
    Rr <- qr.resid(qrX, Rx)
    Wx[, -1] <- sweep(Rr, 1, xj, `*`)
  }
  Wx
}

# ---------------------------------------------------------------------------
# 3. Main spatial engine (Mata setOmsWfin)
# ---------------------------------------------------------------------------
.setOmsWfin <- function(distmat, avc0) {
  
  n <- nrow(distmat)
  distv <- .lvech(distmat)
  
  cgridfac <- 1.2 # factor in c-grid for size control
  minavc   <- 0.00001 # minimal avc value for which size control is checked
  
  if(avc0>=0.05){
    qmax <- 10
  }else if(avc0>=0.01){
    qmax <- 20
  }else if(avc0>=0.005){
    qmax <- 60
  }else{
    qmax <- 120
  }

  c0    <- .getc0fromavc(distv, avc0)
  cmax  <- .getc0fromavc(distv, minavc)
  
  repeat {
    qmax <- min(qmax, n - 1)
    W   <- .getW(distmat, c0, qmax)
    Oms <- .getOms(distmat, c0, cmax, W, cgridfac)
    fin <- .setfinalW(Oms, W, qmax)
    if(fin$q < qmax || qmax == n - 1) break
    qmax <- round(qmax + qmax / 2)      # mimic ado's escalation
  }
  list(Wfin = fin$W, cvfin = fin$cv, Omsfin = Oms,
       c0 = c0, cmax = cmax)
}

# .setOmsWfin <- function(S, avc0, latlong) {
#   ## S   = coordinates matrix (n x d)
#   n     <- nrow(S)
#   distmat <- .getdistmat(S, latlong)
#   distv <- .lvech(distmat)
#   
#   cgridfac <- 1.2 # factor in c-grid for size control
#   minavc   <- 0.00001 # minimal avc value for which size control is checked
#   
#   if(avc0>=0.05){
#     qmax <- 10
#   }else if(avc0>=0.01){
#     qmax <- 20
#   }else if(avc0>=0.005){
#     qmax <- 60
#   }else{
#     qmax <- 120
#   }
#   
#   c0    <- .getc0fromavc(distv, avc0)
#   cmax  <- .getc0fromavc(distv, minavc)
#   
#   counter <- 0
#   repeat {
#     print(counter)
#     counter <- counter + 1
#     qmax <- min(qmax, n - 1)
#     W   <- .getW(distmat, c0, qmax)
#     Oms <- .getOms(distmat, c0, cmax, W, cgridfac)
#     fin <- .setfinalW(Oms, W, qmax)
#     if(fin$q < qmax || qmax == n - 1) break
#     qmax <- round(qmax + qmax / 2)      # mimic ado's escalation
#   }
#   list(Wfin = fin$W, cvfin = fin$cv, Omsfin = Oms, c0 = c0, cmax = cmax)
# }

.get_obs_index <- function(model, data) {
  if ("fixest" %in% class(model)) {
    if(!requireNamespace("fixest", quietly = TRUE)) {
      stop("Package 'fixest' is required for fixest model objects.")
    }
    obs_index <- fixest::obs(model)
  } else if ("lm" %in% class(model)) {
    mf <- stats::model.frame(model)
    rn <- rownames(mf)
    suppressWarnings(obs_index <- as.integer(rn))
    if(anyNA(obs_index)) {
      obs_index <- match(rn, rownames(data))
      if(anyNA(obs_index)) {
        stop("Could not map lm model frame rows back to `data`.")
      }
    }
  } else {
    stop("Model class not supported. Provide an lm or fixest model.")
  }

  if(any(obs_index < 1L | obs_index > nrow(data))) {
    stop("Model observation indices are outside the row range of `data`.")
  }
  obs_index
}

.resolve_coords_input <- function(data, obs_index, lon, lat, coord_euclidean) {
  use_geodesic <- !is.null(lon) || !is.null(lat)
  use_euclidean <- !is.null(coord_euclidean)

  if(use_geodesic && use_euclidean) {
    stop("Specify either `lon`/`lat` or `coord_euclidean`, not both.")
  }
  if(!use_geodesic && !use_euclidean) {
    stop("Specify coordinates via `lon`/`lat` or `coord_euclidean`.")
  }

  if(use_geodesic) {
    if(is.null(lon) || is.null(lat)) {
      stop("For geodesic coordinates, provide both `lon` and `lat`.")
    }
    if(!is.character(lon) || length(lon) != 1L || !nzchar(lon) ||
       !is.character(lat) || length(lat) != 1L || !nzchar(lat)) {
      stop("`lon` and `lat` must each be a single column name.")
    }
    miss <- setdiff(c(lon, lat), names(data))
    if(length(miss) > 0) {
      stop("Coordinate variables not found in data: ", paste(miss, collapse = ", "))
    }
    coords <- data[obs_index, c(lon, lat), drop = FALSE]
    if(!all(vapply(coords, is.numeric, logical(1)))) {
      stop("`lon` and `lat` must reference numeric columns.")
    }
    if(any(!is.finite(as.matrix(coords)))) {
      stop("Geodesic coordinates must be finite.")
    }
    if(any(coords[[lon]] < -180 | coords[[lon]] > 180)) {
      stop("Longitude values must be in [-180, 180].")
    }
    if(any(coords[[lat]] < -90 | coords[[lat]] > 90)) {
      stop("Latitude values must be in [-90, 90].")
    }
    return(list(coords = as.matrix(coords), latlong = TRUE))
  }

  if(!is.character(coord_euclidean) || length(coord_euclidean) < 1L) {
    stop("`coord_euclidean` must be a character vector with at least one column name.")
  }
  miss <- setdiff(coord_euclidean, names(data))
  if(length(miss) > 0) {
    stop("Coordinate variables not found in data: ", paste(miss, collapse = ", "))
  }
  coords <- data[obs_index, coord_euclidean, drop = FALSE]
  if(!all(vapply(coords, is.numeric, logical(1)))) {
    stop("`coord_euclidean` columns must be numeric.")
  }
  if(any(!is.finite(as.matrix(coords)))) {
    stop("Euclidean coordinates must be finite.")
  }
  list(coords = as.matrix(coords), latlong = FALSE)
}

# ---------------------------------------------------------------------------
# 4. Public API - scpc()
# ---------------------------------------------------------------------------
scpc <- function(model,
                 data,
                 lon      = NULL,
                 lat      = NULL,
                 coord_euclidean = NULL,
                 cluster  = NULL,
                 k        = 10,
                 avc      = 0.03,
                 uncond   = FALSE,
                 cvs      = FALSE) {
  
  if(avc <= 0.001 || avc >= 0.99)
    stop("Option avc() must be in (0.001, 0.99).")
  
  ## 1. Influence functions ------------------------------------------
  S    <- sandwich::estfun(model)
  model_mat <- stats::model.matrix(model)
  n    <- nrow(S); p <- ncol(S); neff <- n
  obs_index <- .get_obs_index(model, data)
  coord_info <- .resolve_coords_input(data, obs_index, lon, lat, coord_euclidean)
  coords <- coord_info$coords
  latlong <- coord_info$latlong

  if(!is.null(cluster)) {
    if(length(cluster) == nrow(data)) {
      cluster <- cluster[obs_index]
    }
    cluster <- factor(cluster)
    if(length(cluster) != n) stop("cluster length mismatch")
    S       <- rowsum(S, cluster)
    model_mat <- rowsum(model_mat, cluster)
    coords  <- coords[match(unique(cluster), cluster), ]
    neff    <- nrow(S)
  }
  
  if(ncol(coords) == 1) coords <- cbind(coords, 0)
  
  ## 2. Spatial kernel (unconditional) -------------------------------
  D        <- .getdistmat(coords, latlong)
  spc      <- .setOmsWfin(D, avc)
  Wfin     <- spc$Wfin; cvfin <- spc$cvfin; Omsfin <- spc$Omsfin
  q        <- ncol(Wfin) - 1

  ## 3. Bread ---------------------------------------------------------
  # Match Stata's e(V_modelbased) scaling, which is tied to the model sample
  # size (not the clustered effective sample size).
  bread_inv <- sandwich::bread(model) / n

  ## 4. Loop over k coefficients -------------------------------------
  k_use <- min(k, p)
  out   <- matrix(NA_real_, k_use, 6)
  cvs_mat <- if(cvs) matrix(NA_real_, k_use, 4) else NULL
  levs <- c(0.32, 0.10, 0.05, 0.01)
  cvs_uncond <- if (cvs) vapply(levs, function(lv) .getcv(Omsfin, q, lv), 0.0) else NULL
  rownames(out) <- names(stats::coef(model))[seq_len(k_use)]
  colnames(out) <- c("Coef", "Std_Err", "t", "P>|t|", "CI_low", "CI_high")
  
  for(j in seq_len(k_use)) {
    ## coefficient-specific score & pseudo-outcome
    xj      <- as.numeric(neff * bread_inv[j, ] %*% t(model_mat))      # scpc_x
    #xj_norm <- xj / sqrt(sum(xj^2))                           # scpc_xs
    xjs <- sign(xj)
    wj      <- as.numeric(neff * bread_inv[j, ] %*% t(S)) + stats::coef(model)[j]  # scpc w-vector

    ## unconditional statistic --------------------------------------
    tau_u  <- as.numeric(sqrt(q) * crossprod(Wfin[, 1], wj) /
      sqrt(sum((t(Wfin[, -1]) %*% wj)^2)))
    SE     <- as.numeric(sqrt(sum((t(Wfin[, -1]) %*% wj)^2)) / (sqrt(q) * sqrt(neff)))
    p_u    <- .maxrp(Omsfin, q, abs(tau_u) / sqrt(q))$max
    
    ## conditional branch (skip when uncond=TRUE) --------------------
    if(!uncond) {
      Wx      <- .orthogonalize_W(Wfin, xj, xjs, model_mat)
      Omsx    <- .getOms(D, spc$c0, spc$cmax, Wx, 1.2)
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
    
    if(cvs) {
      cvs_vec <- cvs_uncond
      if(!uncond) {
        cvs_cond <- vapply(levs, function(lv) .getcv(Omsx, q, lv), 0.0)
        cvs_vec <- pmax(cvs_vec, cvs_cond)
      }
      cvs_mat[j, ] <- cvs_vec
    }
  }
  
  ## 5. Critical-value table (32/10/5/1 %) when requested ------------
  # cvs_mat <- if(cvs) {
  #   levs <- c(0.32, 0.10, 0.05, 0.01)
  #   res  <- sapply(levs, function(lv) .getcv(Omsfin, q, lv))
  #   #dimnames(res) <- list("Two-Sided", paste0(levs * 100, "%"))
  #   res
  # } else NULL
  
  structure(list(
    scpcstats = out,
    scpccvs   = cvs_mat,
    W         = Wfin,
    avc       = avc,
    c0        = spc$c0,
    cv        = cvfin,
    q         = q,
    call      = match.call())
    , class = "scpc")
}

print.scpc <- function(x, ...) {
  cat("\nSCPC Inference (k =", nrow(x$scpcstats), ", q =", x$q, ")\n\n")
  stats::printCoefmat(x$scpcstats[, 1:4], P.values = TRUE, has.Pvalue = TRUE)
  if(!is.null(x$scpccvs)) {
    cat("\nTwo-sided critical values:\n")
    print(x$scpccvs)
  }
  invisible(x)
}
