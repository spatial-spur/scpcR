
.getavc <- function(c, dist) {
  ## computes average correlation given vector of pairwise distances
  mean(exp(-c * dist))
}

.getdistmat <- function(S, latlong) {
  ## computes matrix of distances from locations
  ## if external latlong=0, distances computed as norm, otherwise latitude / longitude 
  ## and haversine formula for sphere with radius 1/Pi
  S <- as.matrix(S)
  n <- nrow(S)
  if(latlong) {
    if(!requireNamespace("geodist", quietly = TRUE))
      stop("Install package 'geodist' for lat/long distances.")
    D <- geodist::geodist(S, measure="haversine") / (2 * pi * 6378137)      #  → units πR⊕
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

.GQ <- gaussquad::legendre.quadrature.rules(40)[[40]]  # 40‑point table
.GQx <- .GQ[, 1] * 0.5 + 0.5  # nodes, transformed to [0,1]
.GQw <- .GQ[, 2] * 0.5        # weights, transformed to [0,1]

.getrp <- function(Om, cv) {
  ## computes rejection probability given Om and cv; cv here is cv/sqrt(q) in notation of paper
  Omx        <- -cv^2 * Om
  Omx[1,  ]  <- Om[1, ]               # restore first row
  evals      <- Re(eigen(Omx, only.values = TRUE)$values)
  evals      <- -evals[evals < 0]/max(evals) 
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
    Oms[[i]] <- crossprod(W, exp(-c * distmat) %*% W)  # t(W) Σ(c) W
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
  for(j in 2:ncol(Wx)) {
    rx <- Wx[, j] * xjs
    r <- stats::residuals(stats::lm(rx ~ model_mat))
    Wx[, j] <- r * xj
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
    qmax <- round(qmax + qmax/2)      # mimic ado’s escalation
  }
  list(Wfin = fin$W, cvfin = fin$cv, Omsfin = Oms,
       c0 = c0, cmax = cmax)
}

# .setOmsWfin <- function(S, avc0, latlong) {
#   ## S   = coordinates matrix (n × d)
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
#     qmax <- round(qmax + qmax/2)      # mimic ado’s escalation
#   }
#   list(Wfin = fin$W, cvfin = fin$cv, Omsfin = Oms, c0 = c0, cmax = cmax)
# }

.get_coords <- function(model, data, coords_var) {
  # Check if coords_var exists in data
  if (!all(coords_var %in% names(data))) {
    stop("Coordinate variables not found in data.")
  }
  
  # Determine the estimation sample used in the model
  if ("fixest" %in% class(model)) {
    # fixest: use obs()
    obs_index <- fixest::obs(model)
    data_used <- data[obs_index, , drop = FALSE]
  } else if ("lm" %in% class(model)) {
    # lm: use rownames from model.frame
    mf <- model.frame(model)
    obs_index <- as.numeric(rownames(mf))
    data_used <- data[obs_index, , drop = FALSE]
  } else {
    stop("Model class not supported. Provide an lm or fixest model.")
  }
  
  # Get coordinates for estimation sample
  coords <- data_used[, coords_var, drop = FALSE]
}

# ---------------------------------------------------------------------------
# 4. Public API  –  scpc()
# ---------------------------------------------------------------------------
scpc <- function(model,
                 data,
                 coords_var,
                 cluster  = NULL,
                 k        = 10,
                 avc      = 0.03,
                 latlong  = FALSE,
                 uncond   = FALSE,
                 cvs      = FALSE) {
  
  if(avc <= 0.001 || avc >= 0.99)
    stop("Option avc() must be in (0.001, 0.99).")
  
  ## 1. Influence functions ------------------------------------------
  S    <- sandwich::estfun(model)
  model_mat <- model.matrix(model)
  n    <- nrow(S); p <- ncol(S); neff <- n

  if(!is.null(cluster)) {
    cluster <- factor(cluster)
    if(length(cluster) != n) stop("cluster length mismatch")
    S       <- rowsum(S, cluster)
    model_mat <- rowsum(model_mat, cluster)
    coords  <- coords[match(unique(cluster), cluster), ]
    neff    <- nrow(S)
  }
  
  coords <- as.matrix(.get_coords(model, data, coords_var))
  if(ncol(coords) == 1) coords <- cbind(coords, 0)
  
  ## 2. Spatial kernel (unconditional) -------------------------------
  D        <- .getdistmat(coords, latlong)
  spc      <- .setOmsWfin(D, avc)
  Wfin     <- spc$Wfin; cvfin <- spc$cvfin; Omsfin <- spc$Omsfin
  q        <- ncol(Wfin) - 1

  ## 3. Bread ---------------------------------------------------------
  bread_inv <- sandwich::bread(model) / neff

  ## 4. Loop over k coefficients -------------------------------------
  k_use <- min(k, p)
  out   <- matrix(NA_real_, k_use, 6)
  cvs_mat <- if(cvs) matrix(NA_real_, k_use, 4) else NULL
  rownames(out) <- names(stats::coef(model))[seq_len(k_use)]
  colnames(out) <- c("Coef", "Std_Err", "t", "P>|t|", "CI_low", "CI_high")
  
  for(j in seq_len(k_use)) {
    ## coefficient‑specific score & pseudo‑outcome
    xj      <- as.numeric(neff * bread_inv[j, ] %*% t(model_mat))      # scpc_x
    #xj_norm <- xj / sqrt(sum(xj^2))                           # scpc_xs
    xjs <- sign(xj)
    wj      <- as.numeric(neff * bread_inv[j, ] %*% t(S)) + stats::coef(model)[j]       # scpc w‑vector

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
    
    cvs_mat[j, ] <- if(cvs) {
      levs <- c(0.32, 0.10, 0.05, 0.01)
      cvs_vec <- sapply(levs, function(lv) .getcv(Omsfin, q, lv))
      if(!uncond) {
        for(i in seq_along(cvs_vec)) {
          cvs_vec[i] <- max(cvs_vec[i], .getcv(Omsx, q, levs[i]))
        }
      }
      cvs_vec
    } else NULL
  }
  
  ## 5. Critical‑value table (32/10/5/1 %) when requested ------------
  # cvs_mat <- if(cvs) {
  #   levs <- c(0.32, 0.10, 0.05, 0.01)
  #   res  <- sapply(levs, function(lv) .getcv(Omsfin, q, lv))
  #   #dimnames(res) <- list("Two‑Sided", paste0(levs * 100, "%"))
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
  printCoefmat(x$scpcstats[,1:4], P.values = TRUE, has.Pvalue = TRUE)
  if(!is.null(x$scpccvs)) {
    cat("\nTwo‑sided critical values:\n")
    print(x$scpccvs)
  }
  invisible(x)
}
