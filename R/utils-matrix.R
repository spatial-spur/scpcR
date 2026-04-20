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
  if (requireNamespace("RSpectra", quietly = TRUE) && qmax < n - 1) {
    eig <- RSpectra::eigs_sym(Sig_d, k = qmax, which = "LM")
    V <- eig$vectors
  } else {
    V <- eigen(Sig_d, symmetric = TRUE)$vectors[, seq_len(qmax)]
  }
  cbind(rep(1, n) / sqrt(n), V)
}

.gettau <- function(y, W) {
  ## SCPC t-statistic
  sqrt(ncol(W) - 1) * crossprod(W[, 1], y) / norm(crossprod(W[, -1], y), type = "2")
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

.orthogonalize_W_fixest_iv <- function(W, xj, xjs, template_model, template_data) {
  Wx <- W
  Wx[, 1] <- Wx[, 1] * xj * xjs

  if (ncol(Wx) > 1L) {
    aux_data <- template_data
    aux_var <- "scpc_rx"
    while (aux_var %in% names(aux_data)) {
      aux_var <- paste0(aux_var, "_")
    }
    aux_formula <- stats::as.formula(
      paste(aux_var, "~ ."),
      env = environment(stats::formula(template_model))
    )

    for (col in 2:ncol(Wx)) {
      aux_data[[aux_var]] <- W[, col] * xjs
      aux_fit <- update(
        template_model,
        aux_formula,
        data = aux_data,
        nthreads = 1L,
        notes = FALSE
      )
      Wx[, col] <- stats::residuals(aux_fit) * xj
    }
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
