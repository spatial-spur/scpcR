# ---------------------------------------------------------------------------
# Matrix helpers
# ---------------------------------------------------------------------------

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

.getW <- function(distmat, c0, qmax) {
  ## First qmax eigenvectors of demeaned Sigma(c0), prepended with
  ## a normalised constant vector.  All columns have unit length.
  n <- nrow(distmat)
  Sig <- exp(-c0 * distmat)
  Sig_d <- .demeanmat(Sig)
  Sig_centered <- Sig_d - (n + 1) / n
  k_check <- min(qmax + 1L, n - 1L)
  if (requireNamespace("RSpectra", quietly = TRUE) && k_check < n - 1L) {
    eig <- RSpectra::eigs_sym(Sig_centered, k = k_check, which = "LA")
    vals <- Re(eig$values)
    ord <- order(vals, decreasing = TRUE)
    vals_check <- vals[ord[seq_len(k_check)]]
    V <- eig$vectors[, ord[seq_len(qmax)], drop = FALSE]
  } else {
    eig <- eigen(Sig_centered, symmetric = TRUE)
    vals <- Re(eig$values)
    ord <- order(vals, decreasing = TRUE)
    vals_check <- vals[ord[seq_len(k_check)]]
    V <- eig$vectors[, ord[seq_len(qmax)], drop = FALSE]
  }
  tied <- FALSE
  if (length(vals_check) >= 2L) {
    tol <- 1.6e-11 * max(abs(vals_check[1L]), 1.0)
    tied <- any(vals_check[-length(vals_check)] - vals_check[-1L] <= tol)
  }
  W <- cbind(rep(1, n) / sqrt(n), V)
  attr(W, "tied") <- tied
  W
}

.gettau <- function(y, W) {
  ## SCPC t-statistic
  sqrt(ncol(W) - 1) * crossprod(W[, 1], y) / norm(crossprod(W[, -1], y), type = "2")
}

# ---------------------------------------------------------------------------
# Residualisation helpers
# ---------------------------------------------------------------------------
.scale_design_columns <- function(X) {
  if (ncol(X) == 0L) return(X)
  sweep(X, 2L, .design_column_scales(X), "/")
}

.design_column_scales <- function(X) {
  scales <- apply(abs(X), 2L, max)
  scales[scales == 0] <- 1
  scales
}

.scpc_score_sign <- function(x) {
  out <- sign(x)
  out[out == 0] <- 1
  out
}

.coerce_numeric_matrix <- function(A, arg) {
  vec_input <- is.null(dim(A))
  A <- as.matrix(A)

  if (nrow(A) == 0L) {
    stop(arg, " must have at least one row.")
  }
  if (any(!is.finite(A))) {
    stop(arg, " contains non-finite values.")
  }

  storage.mode(A) <- "double"
  list(mat = A, vec_input = vec_input)
}

.demean_for_scpc <- function(A, fixef_id, arg, restore_scale = TRUE) {
  A <- as.matrix(A)
  if (is.null(fixef_id) || ncol(A) == 0L) {
    return(A)
  }
  if (!requireNamespace("fixest", quietly = TRUE)) {
    stop("Package 'fixest' is required for fixed-effect demeaning.")
  }

  # A fixed-effect shift leaves the within residual unchanged and keeps small
  # within-group differences from being rounded away during column scaling.
  group <- if (is.list(fixef_id)) fixef_id[[1L]] else as.matrix(fixef_id)[, 1L]
  shifted <- A - A[match(group, group), , drop = FALSE]
  if (all(is.finite(shifted))) A <- shifted
  scales <- .design_column_scales(A)
  out <- tryCatch(
    fixest::demean(sweep(A, 2L, scales, "/"), f = fixef_id, nthreads = 1L, tol = 1e-12),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    stop("Could not FE-demean ", arg, ": ", conditionMessage(out))
  }

  out <- as.matrix(out)
  if (nrow(out) != nrow(A) || ncol(out) != ncol(A)) {
    stop(arg, " changed dimensions during fixed-effect demeaning.")
  }

  if (restore_scale) sweep(out, 2L, scales, "*") else out
}

.projection_design <- function(X, fixef_id, arg) {
  if (is.null(fixef_id) || ncol(X) == 0L) return(.scale_design_columns(X))
  X <- .demean_for_scpc(X, fixef_id, arg, restore_scale = FALSE)
  # An absorbed column can retain demeaning roundoff. Discard it before the
  # second scaling so the remnant cannot become a usable design direction.
  absorbed <- apply(abs(X), 2L, max) <= 1e-12
  X[, absorbed] <- 0
  .scale_design_columns(X)
}

.restore_residualizer_shape <- function(A, vec_input) {
  if (isTRUE(vec_input)) {
    return(as.numeric(A[, 1L]))
  }
  A
}

.make_ols_residualizer <- function(X, include_intercept = TRUE, fixef_id = NULL, tol = 1e-10) {
  X <- .coerce_numeric_matrix(X, "OLS residualizer X")$mat
  has_intercept <- !is.null(colnames(X)) && "(Intercept)" %in% colnames(X)
  if (isTRUE(include_intercept) && is.null(fixef_id) && !has_intercept) {
    X <- cbind("(Intercept)" = 1, X)
  }

  X <- .projection_design(X, fixef_id, "OLS residualizer X")
  qrX <- qr(X, tol = tol)
  if (qrX$rank < ncol(X)) {
    stop("OLS residualizer is rank deficient.")
  }

  function(Y) {
    y_in <- .coerce_numeric_matrix(Y, "OLS residualizer Y")
    Y <- y_in$mat
    if (nrow(Y) != nrow(X)) {
      stop("OLS residualizer received Y with an incompatible row count.")
    }
    if (ncol(Y) == 0L) {
      return(.restore_residualizer_shape(Y, y_in$vec_input))
    }

    Y <- .demean_for_scpc(Y, fixef_id, "OLS residualizer Y")
    .restore_residualizer_shape(qr.resid(qrX, Y), y_in$vec_input)
  }
}

.build_iv_projector <- function(X, Z, fixef_id = NULL, tol = 1e-10) {
  X <- .coerce_numeric_matrix(X, "IV residualizer X")$mat
  Z <- .coerce_numeric_matrix(Z, "IV residualizer Z")$mat
  if (nrow(X) != nrow(Z)) {
    stop("IV residualizer requires X and Z to have the same number of rows.")
  }

  X <- .projection_design(X, fixef_id, "IV residualizer X")
  Z <- .projection_design(Z, fixef_id, "IV residualizer Z")
  qrZ <- qr(Z, tol = tol)
  p <- ncol(X)
  if (qrZ$rank < p) stop("IV residualizer is rank deficient: too few independent instruments.")
  Qz <- qr.Q(qrZ)[, seq_len(qrZ$rank), drop = FALSE]

  # H = X (Q_Z' X)^+ Q_Z'. Use the same scaled decomposition to
  # check identification and construct the two factors, without normal equations.
  design <- svd(crossprod(Qz, X), nu = p, nv = p)
  if (design$d[p] <= tol * design$d[1L] || design$d[1L] == 0) {
    stop("IV residualizer is rank deficient: projected design is not full rank.")
  }
  V <- Qz %*% (sweep(design$u, 2L, design$d, "/") %*% t(design$v))
  list(U = X, V = V)
}

.make_iv_residualizer <- function(X, Z, fixef_id = NULL, tol = 1e-10, transpose = FALSE) {
  proj <- .build_iv_projector(X, Z, fixef_id, tol)
  left <- if (isTRUE(transpose)) proj$V else proj$U
  right <- if (isTRUE(transpose)) proj$U else proj$V

  function(Y) {
    y_in <- .coerce_numeric_matrix(Y, "IV residualizer Y")
    Y <- y_in$mat
    if (nrow(Y) != nrow(left)) {
      stop("IV residualizer received Y with an incompatible row count.")
    }
    if (ncol(Y) == 0L) {
      return(.restore_residualizer_shape(Y, y_in$vec_input))
    }

    Y <- .demean_for_scpc(Y, fixef_id, "IV residualizer Y")
    .restore_residualizer_shape(Y - left %*% crossprod(right, Y), y_in$vec_input)
  }
}

# ---------------------------------------------------------------------------
# Orthogonalisation helper (conditional SCPC)
# ---------------------------------------------------------------------------
.orthogonalize_W <- function(W, xj, xjs, model_mat, include_intercept = TRUE, fixef_id = NULL) {
  Wx <- W
  Wx[, 1] <- Wx[, 1] * xj * xjs

  if (ncol(Wx) > 1L) {
    residualize <- .make_ols_residualizer(
      model_mat,
      include_intercept = include_intercept,
      fixef_id = fixef_id
    )
    ## Conditional OLS weights are S_x M_V D_x B_q.
    Rx <- sweep(Wx[, -1, drop = FALSE], 1, xj, `*`)
    Rr <- residualize(Rx)
    Wx[, -1] <- sweep(Rr, 1, xjs, `*`)
  }
  Wx
}

.orthogonalize_W_iv <- function(Wfin, xj, xjs, residualize) {
  if (length(xj) != nrow(Wfin)) {
    stop("Conditional IV projection received xj with an incompatible length.")
  }
  if (length(xjs) != nrow(Wfin)) {
    stop("Conditional IV projection received xjs with an incompatible length.")
  }
  if (any(!is.finite(xj)) || any(!is.finite(xjs))) {
    stop("Conditional IV projection received non-finite xj or xjs values.")
  }

  Wx <- Wfin
  Wx[, 1] <- Wfin[, 1] * xj * xjs

  if (ncol(Wfin) > 1L) {
    Y <- sweep(Wfin[, -1, drop = FALSE], 1, xj, `*`)
    R <- residualize(Y)
    Wx[, -1] <- sweep(R, 1, xjs, `*`)
  }

  if (!all(is.finite(Wx))) {
    stop("Conditional IV projection produced non-finite W values.")
  }

  Wx
}

# ---------------------------------------------------------------------------
# Orthogonalisation helper (conditional SCPC, clustered)
# ---------------------------------------------------------------------------
.orthogonalize_W_cluster <- function(W, cl_vec, xj_indiv, model_mat_indiv,
                                     include_intercept = TRUE, fixef_id = NULL) {
  ## Clustered OLS denominator weights are Xi_s' M_V Xi B_q.
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
    residualize <- .make_ols_residualizer(
      model_mat_indiv,
      include_intercept = include_intercept,
      fixef_id = fixef_id
    )
    Y <- sweep(W_indiv[, -1, drop = FALSE], 1, xj_indiv, `*`)
    R <- residualize(Y)
    Wx[, -1] <- rowsum(sweep(R, 1, xjs_indiv, `*`), cl_vec)
  }

  Wx
}

.orthogonalize_W_cluster_iv <- function(Wfin, cl_vec, xj_indiv, residualize) {
  nclust <- nrow(Wfin)
  ncol_W <- ncol(Wfin)
  cl_idx <- as.integer(cl_vec)

  if (length(xj_indiv) != length(cl_vec)) {
    stop("Clustered IV projection received xj_indiv with an incompatible length.")
  }
  if (any(!is.finite(xj_indiv))) {
    stop("Clustered IV projection received non-finite xj_indiv values.")
  }

  xj_sq_sum <- as.numeric(rowsum(xj_indiv^2, cl_vec))
  xjs_indiv <- xj_indiv / sqrt(xj_sq_sum[cl_idx])
  xjs_indiv[!is.finite(xjs_indiv)] <- 0

  W_indiv <- Wfin[cl_idx, , drop = FALSE]
  Wx <- matrix(0, nclust, ncol_W)

  Wx[, 1] <- as.numeric(rowsum(W_indiv[, 1] * xjs_indiv * xj_indiv, cl_vec))

  if (ncol_W > 1L) {
    Y <- sweep(W_indiv[, -1, drop = FALSE], 1, xj_indiv, `*`)
    R <- residualize(Y)
    Wx[, -1] <- rowsum(sweep(R, 1, xjs_indiv, `*`), cl_vec)
  }

  if (!all(is.finite(Wx))) {
    stop("Clustered conditional IV projection produced non-finite W values.")
  }

  Wx
}
