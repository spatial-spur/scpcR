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

.col_norms_safe <- function(M) {
  M <- as.matrix(M)
  if (length(M) == 0L) return(numeric(ncol(M)))
  scales <- apply(abs(M), 2L, max)
  out <- numeric(ncol(M))
  nz <- is.finite(scales) & scales > 0
  if (any(nz)) {
    scaled <- sweep(M[, nz, drop = FALSE], 2L, scales[nz], "/")
    out[nz] <- scales[nz] * sqrt(colSums(scaled^2))
  }
  out[!is.finite(scales)] <- scales[!is.finite(scales)]
  out
}

.group_max <- function(x, group) {
  as.numeric(tapply(x, group, max, default = 0))
}

.scaled_cluster_norm <- function(x, group) {
  scale <- .group_max(abs(x), group)
  row_scale <- scale[as.integer(group)]
  row_scale[row_scale == 0] <- 1
  scaled <- x / row_scale
  norm <- numeric(nlevels(group))
  norm[tabulate(as.integer(group), nlevels(group)) > 0L] <- sqrt(as.numeric(rowsum(scaled^2, group)))
  list(scale = scale, scaled = scaled, norm = norm)
}

.normalize_cluster_direction <- function(xj_indiv, cl_vec) {
  parts <- .scaled_cluster_norm(xj_indiv, cl_vec)
  norm <- parts$norm[as.integer(cl_vec)]
  norm[norm == 0] <- 1
  parts$scaled / norm
}

.cluster_l2_norms <- function(xj_indiv, cl_vec) {
  parts <- .scaled_cluster_norm(xj_indiv, cl_vec)
  parts$scale * parts$norm
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

.scpc_spacing <- function(mag) {
  exponent <- floor(log2(pmax(mag, .Machine$double.xmin)))
  power <- 2^pmin(exponent, .Machine$double.max.exp - 1)
  # log2 can round a value just below a power of two up to that power.
  high <- mag >= .Machine$double.xmin & power > mag
  power[high] <- power[high] / 2
  power * .Machine$double.eps
}

.demean_one_way <- function(A, ids, drop_absorbed) {
  group <- factor(ids)
  index <- as.integer(group)
  first <- match(seq_len(nlevels(group)), index)
  for (j in seq_len(ncol(A))) {
    x <- A[, j]
    largest <- .group_max(abs(x), group)
    anchor <- x[first]
    delta <- .group_max(abs(x - anchor[index]), group)
    shift <- delta < largest
    anchor[!shift] <- 0
    scale <- ifelse(shift, delta, largest)
    scale[scale == 0] <- 1
    scaled <- (x - anchor[index]) / scale[index]
    center <- as.numeric(tapply(scaled, group, mean))
    residual <- scaled - center[index]
    if (drop_absorbed) {
      drop <- 2 * .group_max(abs(residual), group) < .scpc_spacing(largest) / scale
      residual[drop[index]] <- 0
    }
    A[, j] <- residual * scale[index]
  }
  if (any(!is.finite(A))) stop("Fixed-effect demeaning produced a non-finite entry.")
  A
}

.demean_for_scpc <- function(A, fixef_id, arg, restore_scale = TRUE, drop_absorbed = FALSE) {
  A <- as.matrix(A)
  if (is.null(fixef_id) || ncol(A) == 0L) {
    return(A)
  }
  if (!requireNamespace("fixest", quietly = TRUE)) {
    stop("Package 'fixest' is required for fixed-effect demeaning.")
  }

  groups <- if (is.list(fixef_id)) fixef_id else as.list(as.data.frame(fixef_id))
  if (length(groups) == 1L) {
    out <- .demean_one_way(A, groups[[1L]], drop_absorbed)
    return(if (restore_scale) out else .scale_design_columns(out))
  }
  original_scales <- apply(abs(A), 2L, max)
  # Every anchor shift is in the FE span. Apply shifts cumulatively so
  # large offsets in multiple fixed-effect dimensions are removed up front.
  shifted <- A
  best <- original_scales
  for (group in groups) {
    candidate <- shifted - shifted[match(group, group), , drop = FALSE]
    sizes <- apply(abs(candidate), 2L, max)
    take <- is.finite(sizes) & sizes < best
    shifted[, take] <- candidate[, take, drop = FALSE]
    best[take] <- sizes[take]
  }
  scales <- .design_column_scales(shifted)
  demean <- function(values) {
    tryCatch(
      as.matrix(fixest::demean(values, f = groups, nthreads = 1L, tol = 1e-14)),
      error = function(e) stop("Could not FE-demean ", arg, ": ", conditionMessage(e))
    )
  }
  out <- demean(sweep(shifted, 2L, scales, "/"))

  if (nrow(out) != nrow(A) || ncol(out) != ncol(A)) {
    stop(arg, " changed dimensions during fixed-effect demeaning.")
  }

  # A difference smaller than one input spacing can be introduced by
  # rounding the absorbed levels. Retain contrasts whose demeaned
  # amplitude reaches half an input spacing (1 ULP peak-to-peak).
  noise <- 0.5 * (.scpc_spacing(original_scales) / scales)
  size <- apply(abs(out), 2L, max)
  group_codes <- lapply(groups, function(group) match(group, unique(group)))
  counts <- lapply(group_codes, tabulate)
  needs_refinement <- function() {
    bad <- rep(FALSE, ncol(out))
    for (i in seq_along(group_codes)) {
      means <- rowsum(out, group_codes[[i]], reorder = FALSE) / counts[[i]]
      bad <- bad | (apply(abs(means), 2L, max) > 1e-11 * size)
    }
    bad & size >= noise
  }
  # fixest stops on an iteration tolerance. Verify the group moments and
  # reproject any residual that needs it at its own scale.
  refine <- needs_refinement()
  for (step in seq_len(6L)) {
    if (!any(refine)) break
    local_scale <- size[refine]
    out[, refine] <- sweep(
      demean(sweep(out[, refine, drop = FALSE], 2L, local_scale, "/")), 2L, local_scale, "*"
    )
    size[refine] <- apply(abs(out[, refine, drop = FALSE]), 2L, max)
    refine <- needs_refinement()
  }
  if (any(refine)) {
    stop("Could not FE-demean ", arg, " to the required group-moment tolerance.")
  }
  if (drop_absorbed) out[, size < noise] <- 0
  if (restore_scale) sweep(out, 2L, scales, "*") else out
}

.projection_design <- function(X, fixef_id, arg) {
  if (is.null(fixef_id) || ncol(X) == 0L) return(.scale_design_columns(X))
  .scale_design_columns(.demean_for_scpc(
    X, fixef_id, arg, restore_scale = FALSE, drop_absorbed = TRUE
  ))
}

.restore_residualizer_shape <- function(A, vec_input) {
  if (isTRUE(vec_input)) {
    return(as.numeric(A[, 1L]))
  }
  A
}

.scpc_svd_basis <- function(design, tol, require_full_rank) {
  if (ncol(design) == 0L) return(matrix(0, nrow(design), 0L))
  decomp <- svd(design, nv = 0L)
  keep <- decomp$d > tol * decomp$d[1L]
  if (require_full_rank && sum(keep) < ncol(design)) {
    stop("OLS residualizer is rank deficient.", call. = FALSE)
  }
  decomp$u[, keep, drop = FALSE]
}

.make_ols_residualizer <- function(X, include_intercept = TRUE, fixef_id = NULL, tol = 1e-10,
                                  design_within = FALSE) {
  X <- .coerce_numeric_matrix(X, "OLS residualizer X")$mat
  has_intercept <- !is.null(colnames(X)) && "(Intercept)" %in% colnames(X)
  if (isTRUE(include_intercept) && is.null(fixef_id) && !has_intercept) {
    X <- cbind("(Intercept)" = 1, X)
  }

  X <- .projection_design(X, if (isTRUE(design_within)) NULL else fixef_id, "OLS residualizer X")
  Q <- .scpc_svd_basis(X, tol, require_full_rank = TRUE)

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
    .restore_residualizer_shape(Y - Q %*% crossprod(Q, Y), y_in$vec_input)
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
  Qz <- .scpc_svd_basis(Z, tol, require_full_rank = FALSE)
  p <- ncol(X)
  if (ncol(Qz) < p) stop("IV residualizer is rank deficient: too few independent instruments.")

  # H = X (Q_Z' X)^+ Q_Z'. Use the same scaled decomposition to
  # check identification and construct the two factors, without normal equations.
  design <- svd(crossprod(Qz, X), nu = p, nv = p)
  if (design$d[p] <= tol * max(1, design$d[1L])) {
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
.orthogonalize_W <- function(W, xj, xjs, model_mat, include_intercept = TRUE, fixef_id = NULL,
                             design_within = FALSE) {
  residualize <- .make_ols_residualizer(model_mat, include_intercept, fixef_id,
                                       design_within = design_within)
  .conditional_weights(W, xj, xjs, residualize)
}

.conditional_weights <- function(Wfin, xj, xjs, residualize) {
  if (length(xj) != nrow(Wfin)) {
    stop("Conditional projection received xj with an incompatible length.")
  }
  if (length(xjs) != nrow(Wfin)) {
    stop("Conditional projection received xjs with an incompatible length.")
  }
  if (any(!is.finite(xj)) || any(!is.finite(xjs))) {
    stop("Conditional projection received non-finite xj or xjs values.")
  }

  Wx <- Wfin
  Wx[, 1] <- Wfin[, 1] * xj * xjs

  if (ncol(Wfin) > 1L) {
    Y <- sweep(Wfin[, -1, drop = FALSE], 1, xj, `*`)
    R <- residualize(Y)
    Wx[, -1] <- sweep(R, 1, xjs, `*`)
  }

  if (!all(is.finite(Wx))) {
    stop("Conditional projection produced non-finite W values.")
  }

  Wx
}

# ---------------------------------------------------------------------------
# Orthogonalisation helper (conditional SCPC, clustered)
# ---------------------------------------------------------------------------
.orthogonalize_W_cluster <- function(W, cl_vec, xj_indiv, model_mat_indiv,
                                     include_intercept = TRUE, fixef_id = NULL,
                                     design_within = FALSE) {
  residualize <- .make_ols_residualizer(model_mat_indiv, include_intercept, fixef_id,
                                       design_within = design_within)
  .conditional_cluster_weights(W, cl_vec, xj_indiv, residualize)
}

.conditional_cluster_weights <- function(Wfin, cl_vec, xj_indiv, residualize,
                                         xjs_indiv = .normalize_cluster_direction(xj_indiv, cl_vec)) {
  nclust <- nrow(Wfin)
  ncol_W <- ncol(Wfin)
  cl_idx <- as.integer(cl_vec)

  if (length(xj_indiv) != length(cl_vec)) {
    stop("Clustered projection received xj_indiv with an incompatible length.")
  }
  if (any(!is.finite(xj_indiv))) {
    stop("Clustered projection received non-finite xj_indiv values.")
  }

  W_indiv <- Wfin[cl_idx, , drop = FALSE]
  Wx <- matrix(0, nclust, ncol_W)

  Wx[, 1] <- as.numeric(rowsum(W_indiv[, 1] * xjs_indiv * xj_indiv, cl_vec))

  if (ncol_W > 1L) {
    Y <- sweep(W_indiv[, -1, drop = FALSE], 1, xj_indiv, `*`)
    R <- residualize(Y)
    Wx[, -1] <- rowsum(sweep(R, 1, xjs_indiv, `*`), cl_vec)
  }

  if (!all(is.finite(Wx))) {
    stop("Clustered conditional projection produced non-finite W values.")
  }

  Wx
}
