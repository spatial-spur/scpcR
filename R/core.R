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
#'   in \code{data}. Coordinates must be constant within clusters.
#' @param ncoef Integer; number of coefficients to report.  Default
#'   \code{NULL} reports all coefficients.
#' @param avc Numeric; upper bound on the average pairwise correlation
#'   for which size is controlled.  Must be in \code{(0.001, 0.99)}.
#'   Default is 0.03.
#' @param method Character string selecting the spatial algorithm.
#'   \code{"auto"} chooses \code{"exact"} for smaller problems and
#'   \code{"approx"} once the large-\eqn{n} threshold is reached.
#'   \code{"exact"} always uses the full distance matrix.
#'   \code{"approx"} uses the original large-\eqn{n} approximation branch.
#' @param large_n_seed Numeric; integer-valued seed used by the
#'   large-\eqn{n} approximation branch.  Ignored when
#'   \code{method = "exact"}. Default is 1.
#' @param uncond Logical; if \code{TRUE}, report unconditional critical
#'   values only (skip the conditional adjustment of Mueller and Watson
#'   2023).  Default is \code{FALSE}.
#' @param cvs Logical; if \code{TRUE}, include per-coefficient critical
#'   values at the 32\%, 10\%, 5\%, and 1\% levels.  Default is
#'   \code{FALSE}.
#' @param backend Character string selecting the computation backend:
#'   \code{"auto"} (default, uses native Rust when available and eligible),
#'   \code{"rust"} (requires the native Rust backend), or
#'   \code{"fallback"} (pure-R implementation). When omitted, checks
#'   \code{getOption("spatial.backend")} and \code{Sys.getenv("SPATIAL_BACKEND")}.
#' @param nthreads Optional positive integer thread count for parallel
#'   Rust execution, or \code{NULL} (default) to check
#'   \code{getOption("scpcR.threads")}, \code{Sys.getenv("SCPCR_NUM_THREADS")},
#'   and \code{getOption("mc.cores")}. Note that \code{nthreads = 1} executes
#'   Rust on a single thread without falling back to R.
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
#'   \item{\code{method}}{Spatial algorithm used:
#'     \code{"exact"} or \code{"approx"}.}
#'   \item{\code{large_n_seed}}{Seed used by the large-\eqn{n}
#'     approximation branch.}
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{backend}}{Backend used: \code{"rust"} or \code{"fallback"}.}
#'   \item{\code{fallback_reason}}{Reason for selecting the fallback, if used.}
#'   \item{\code{tied}}{Logical; \code{TRUE} when leading spatial eigenvalues
#'     are numerically tied.}
#' }
#' @details Weighted regression fits are unsupported. An indefinite sampled
#'   covariance causes approximate SCPC to report an error. Use exact SCPC
#'   when the sample size permits.
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
#' @useDynLib scpcR, .registration = TRUE
#' @export
scpc <- function(model,
                 data,
                 lon = NULL,
                 lat = NULL,
                 coords_euclidean = NULL,
                 cluster = NULL,
                 ncoef = NULL,
                 avc = 0.03,
                 method = "auto",
                 large_n_seed = 1,
                 uncond = FALSE,
                 cvs = FALSE,
                 backend = c("auto", "rust", "fallback"),
                 nthreads = NULL) {

  backend_resolved <- .resolve_scpc_backend(backend, backend_missing = missing(backend))
  if (!is.null(nthreads)) {
    .parse_positive_int_threads(nthreads, "`nthreads`")
  }

  if (!is.numeric(avc) || length(avc) != 1L || is.na(avc) || !is.finite(avc)) {
    stop("`avc` must be a single finite numeric value.")
  }
  if (avc <= 0.001 || avc >= 0.99) {
    stop("`avc` must lie in (0.001, 0.99).")
  }
  method <- .validate_scpc_method(method)
  large_n_seed <- .validate_large_n_seed(large_n_seed)
  if (
    !is.null(ncoef) &&
    (
      !is.numeric(ncoef) ||
      length(ncoef) != 1L ||
      is.na(ncoef) ||
      !is.finite(ncoef) ||
      ncoef <= 0 ||
      abs(ncoef - round(ncoef)) > 0
    )
  ) {
    stop("`ncoef` must be NULL or a single positive integer.")
  }
  if (!is.logical(uncond) || length(uncond) != 1L || is.na(uncond)) {
    stop("`uncond` must be a single TRUE/FALSE value.")
  }
  if (!is.logical(cvs) || length(cvs) != 1L || is.na(cvs)) {
    stop("`cvs` must be a single TRUE/FALSE value.")
  }
  if (!is.null(cluster) && (!is.character(cluster) || length(cluster) != 1L || !nzchar(cluster))) {
    stop("`cluster` must be a single non-empty column name.")
  }

  # The conditional benchmark is defined for unweighted fits.
  if (.has_regression_weights(model, data)) {
    stop("SCPC does not support models fitted with regression weights.", call. = FALSE)
  }

  if (inherits(model, "glm")) {
    fam <- model$family
    if (is.null(fam) || !identical(fam$family, "gaussian") || !identical(fam$link, "identity")) {
      stop("SCPC only supports Gaussian identity-link GLMs.", call. = FALSE)
    }
  }

  ## 1. Influence functions --------------------------------------------------
  coef_vec <- stats::coef(model)
  if (length(coef_vec) == 0L) {
    stop("SCPC requires at least one estimated model coefficient.", call. = FALSE)
  }
  S <- sandwich::estfun(model)
  raw_coef_names <- names(coef_vec)
  is_fixest_iv <- .is_fixest_iv_second_stage(model)
  iv_cond_design <- NULL
  if (is_fixest_iv && !uncond) {
    iv_cond_design <- .get_fixest_iv_design(model)
    model_mat <- iv_cond_design$model_mat
  } else {
    model_mat <- .get_scpc_model_matrix(model)
  }
  if (inherits(model[["na.action"]], "exclude") && nrow(S) > nrow(model_mat)) {
    na_idx <- as.integer(model[["na.action"]])
    if (length(na_idx) > 0L && nrow(S) - length(na_idx) == nrow(model_mat)) {
      S <- S[-na_idx, , drop = FALSE]
    }
  }
  n <- nrow(S)
  p <- ncol(S)
  neff <- n
  if (nrow(model_mat) != n) {
    stop("Model matrix row count does not match score matrix row count.")
  }
  if (ncol(model_mat) != p) {
    stop("Model matrix column count does not match the coefficient vector.")
  }

  obs_index <- .get_obs_index(model, data)

  # Keep the model matrices in estimation-sample order. The selected backend
  # prepares its own numerical design after dispatch.
  cond_fixef_id <- if (!uncond && .has_fixest_fe(model)) .get_fixest_fe_spec(model) else NULL

  coord_info <- .resolve_coords_input(data, obs_index, lon, lat, coords_euclidean)
  coords <- coord_info$coords
  latlong <- coord_info$latlong

  cl_vec <- NULL
  if (!is.null(cluster)) {
    if (!cluster %in% names(data)) {
      stop("Cluster variable not found in data: ", cluster)
    }
    cl_vec <- factor(data[[cluster]][obs_index])
    if (anyNA(cl_vec)) stop("Cluster IDs must be present for every observation in the estimation sample.")

    ## Coordinates must be constant within clusters.
    coord_first <- coords[match(levels(cl_vec), cl_vec), , drop = FALSE]
    if (any(coords != coord_first[as.integer(cl_vec), , drop = FALSE])) {
      stop(
        "Coordinates vary within clusters. Provide coordinates that are constant ",
        "within each cluster before calling `scpc()`."
      )
    }

    coords <- coord_first
    neff <- nrow(coords)
  }

  if (ncol(coords) == 1) coords <- cbind(coords, 0)

  ## 2. Bread & coefficient subset ------------------------------------------
  bread_inv <- sandwich::bread(model) / n
  k_use <- if (is.null(ncoef)) p else min(ncoef, p)
  coef_names <- raw_coef_names[seq_len(k_use)]

  ## 3. Native Rust dispatch (exact unweighted OLS, IV, and 1-way FE)
  fe_requires_fallback <- !is.null(cond_fixef_id) &&
    (length(model$fixef_vars) != 1L || inherits(cond_fixef_id, "scpc_varying_slopes"))
  fe_id_0 <- NULL
  if (!is.null(cond_fixef_id) && !fe_requires_fallback) {
    fe_raw <- if (is.list(cond_fixef_id)) {
      cond_fixef_id[[1L]]
    } else if (!is.null(dim(cond_fixef_id))) {
      cond_fixef_id[, 1L]
    } else {
      cond_fixef_id
    }
    fe_id_0 <- as.integer(factor(fe_raw)) - 1L
    storage.mode(fe_id_0) <- "integer"
  }

  method_actual <- if (identical(method, "approx") || (identical(method, "auto") && neff >= 4500L)) {
    "approx"
  } else {
    "exact"
  }
  rust_eligible <- (method_actual %in% c("exact", "approx")) &&
    !fe_requires_fallback &&
    !inherits(model, "glm") &&
    (inherits(model, "lm") || inherits(model, "fixest")) &&
    k_use >= 1L
  rust_loaded <- !identical(backend_resolved, "fallback") &&
    is.loaded("wrap__rust_run_model_scpc", PACKAGE = "scpcR")
  rust_forked <- rust_loaded && .scpc_accelerate_state()

  if (identical(backend_resolved, "rust")) {
    if (!rust_loaded) {
      stop("`backend = \"rust\"` was requested, but the native Rust backend is not loaded in `scpcR`.")
    }
    if (rust_forked) {
      stop("Accelerate was used before this process forked. Use a fresh process or backend = 'fallback'.",
           call. = FALSE)
    }
    if (!rust_eligible) {
      stop(
        "`backend = \"rust\"` currently supports unweighted OLS, IV, and 1-way FE models only ",
        "(no weights, no multi-way absorbed fixed effects)."
      )
    }
  }

  fallback_reason <- if (identical(backend_resolved, "fallback")) {
    "requested"
  } else if (rust_forked) {
    "forked_process"
  } else if (!rust_loaded) {
    "native_unavailable"
  } else {
    "unsupported_model"
  }
  native_tied <- FALSE
  if (rust_eligible && rust_loaded && !rust_forked) {
    nthreads_resolved <- .resolve_scpc_nthreads(nthreads)
    rust_tied_fallback <- FALSE
    rust_res <- tryCatch(
      .run_scpc_rust_model(
        coords, latlong, avc, cvs, coef_vec[seq_len(k_use)],
        S, bread_inv[seq_len(k_use), , drop = FALSE],
        design = if (uncond) NULL else if (is_fixest_iv) iv_cond_design$X else model_mat,
        instruments = if (!uncond && is_fixest_iv) iv_cond_design$Z else NULL,
        direction_design = if (!uncond && is_fixest_iv) model_mat else NULL,
        cluster_id = if (is.null(cl_vec)) NULL else as.integer(cl_vec) - 1L,
        fe_id = fe_id_0, nthreads = nthreads_resolved,
        seed = if (identical(method_actual, "approx")) as.double(large_n_seed) else NULL
      ),
      error = function(e) {
        if (identical(method_actual, "approx") &&
            identical(backend_resolved, "auto") &&
            inherits(e, "spatial_tied_eigenspace")) {
          rust_tied_fallback <<- TRUE
          return(NULL)
        }
        stop(e)
      }
    )
    native_tied <- rust_tied_fallback || isTRUE(rust_res$tied)
    if (native_tied) fallback_reason <- "tied_eigenspace"
    if (!isTRUE(rust_tied_fallback) && !is.null(rust_res) &&
        (!isTRUE(rust_res$tied) || identical(backend_resolved, "rust"))) {
      out <- rust_res$scpcstats
      rownames(out) <- coef_names
      colnames(out) <- c("Coef", "Std_Err", "t", "P>|t|", "2.5 %", "97.5 %")

      cvs_mat <- rust_res$scpccvs
      if (!is.null(cvs_mat)) {
        rownames(cvs_mat) <- coef_names
        colnames(cvs_mat) <- c("32%", "10%", "5%", "1%")
      }

      return(structure(list(
        scpcstats = out,
        scpccvs = cvs_mat,
        W = rust_res$wfin,
        avc = avc,
        c0 = rust_res$c0,
        cv = rust_res$cvfin,
        q = as.numeric(rust_res$q),
        method = method_actual,
        large_n_seed = large_n_seed,
        call = match.call(),
        backend = "rust",
        fallback_reason = NULL,
        tied = isTRUE(rust_res$tied)
      ), class = "scpc"))
    }
  }

  if (!is.null(cluster)) S <- rowsum(S, cl_vec)

  ## 4. Spatial kernel (pure-R fallback, unconditional) ----------------------
  spc <- .setOmsWfin(
    coords,
    avc0 = avc,
    latlong = latlong,
    method = method,
    large_n_seed = large_n_seed
  )
  Wfin <- spc$Wfin
  cvfin <- spc$cvfin
  Omsfin <- spc$Omsfin
  perm <- spc$perm
  q <- ncol(Wfin) - 1
  # Common random pairs make cutoffs invariant to coefficient display order.
  conditional_pairs <- if (!uncond && spc$large_n) {
    .scpc_pair_sample(spc$coords, 1000000L, spc$random_state, latlong)
  } else NULL

  # Clustered projectors act in observation order. Unclustered Nyström
  # weights use the sampled site ordering, including for fixed effects.
  if (!uncond) {
    model_mat_cond <- model_mat
    if (!is.null(cond_fixef_id)) {
      model_mat_cond <- .demean_for_scpc(
        model_mat, cond_fixef_id, "fixest conditional regressors",
        drop_absorbed = !is_fixest_iv
      )
    }
    obs_order <- if (is.null(cluster)) perm else seq_len(n)
    fe_ordered <- .permute_fixef_id(cond_fixef_id, obs_order)
    residualize <- if (is_fixest_iv) {
      .make_iv_residualizer(iv_cond_design$X[obs_order, , drop = FALSE],
                            iv_cond_design$Z[obs_order, , drop = FALSE],
                            fixef_id = fe_ordered, transpose = TRUE)
    } else {
      # The fitted model matrix already contains any intercept.
      .make_ols_residualizer(model_mat_cond[obs_order, , drop = FALSE],
                             include_intercept = FALSE,
                             fixef_id = fe_ordered, design_within = TRUE)
    }
    if (!is.null(cluster)) {
      cl_vec_scpc <- factor(as.character(cl_vec), levels = levels(cl_vec)[perm])
    }
  }

  ## 5. Loop over coefficients -----------------------------------------------
  out <- matrix(NA_real_, k_use, 6)
  cvs_mat <- if (cvs) matrix(NA_real_, k_use, 4) else NULL
  cvs_uncond <- if (cvs) {
    c(
      .getcv(Omsfin, q, 0.32),
      .getcv(Omsfin, q, 0.10),
      cvfin,
      .getcv(Omsfin, q, 0.01)
    )
  } else {
    NULL
  }
  rownames(out) <- coef_names
  colnames(out) <- c("Coef", "Std_Err", "t", "P>|t|", "2.5 %", "97.5 %")

  for (j in seq_len(k_use)) {
    coef_j <- coef_vec[[j]]
    wj <- as.numeric(neff * bread_inv[j, ] %*% t(S)) + coef_j
    wj_perm <- wj[perm]

    ## unconditional statistic -----------------------------------------------
    wj_scale <- max(abs(wj_perm))
    if (!is.finite(wj_scale) || wj_scale == 0) {
      stop("SCPC variance projection is degenerate for coefficient: ", coef_names[[j]])
    }
    wj_perm <- wj_perm / wj_scale
    wj_norm <- .col_norms_safe(wj_perm)
    denom <- .col_norms_safe(crossprod(Wfin[, -1, drop = FALSE], wj_perm))
    if (!is.finite(denom) || !is.finite(wj_norm) || wj_norm == 0 || denom <= 1e-12 * wj_norm) {
      stop("SCPC variance projection is degenerate for coefficient: ", coef_names[[j]])
    }
    tau_u <- as.numeric(sqrt(q) * crossprod(Wfin[, 1], wj_perm) / denom)
    SE <- as.numeric((denom / (sqrt(q) * sqrt(neff))) * wj_scale)
    p_u <- .maxrp(Omsfin, q, abs(tau_u) / sqrt(q))$max

    ## conditional branch (skip when uncond = TRUE) --------------------------
    if (!uncond) {
      xj <- as.numeric(neff * (model_mat_cond %*% t(bread_inv[j, , drop = FALSE])))
      if (any(!is.finite(xj))) stop("Conditional score direction contains non-finite values.")
      xjs <- if (is.null(cluster)) .scpc_score_sign(xj) else .normalize_cluster_direction(xj, cl_vec_scpc)
      direction_scale <- max(abs(xj))
      if (direction_scale > 0) xj <- xj / direction_scale
      if (is.null(cluster)) {
        xj <- xj[perm]
        pre_max <- max(.col_norms_safe(sweep(Wfin[, -1, drop = FALSE], 1L, xj, `*`)))
        Wx <- .conditional_weights(Wfin, xj, xjs[perm], residualize)
      } else {
        cl_norm <- .cluster_l2_norms(xj, cl_vec_scpc)
        pre_max <- max(.col_norms_safe(sweep(Wfin[, -1, drop = FALSE], 1L, cl_norm, `*`)))
        Wx <- .conditional_cluster_weights(Wfin, cl_vec_scpc, xj, residualize, xjs)
      }

      if (!all(is.finite(Wx))) {
        stop("Conditional projection produced non-finite W values for coefficient: ", coef_names[[j]])
      }
      if (nrow(Wx) != nrow(Wfin) || ncol(Wx) != ncol(Wfin)) {
        stop("Conditional projection produced W with incompatible dimensions for coefficient: ", coef_names[[j]])
      }
      if (ncol(Wx) > 1L) {
        col_norms <- .col_norms_safe(Wx[, -1, drop = FALSE])
        if (!is.finite(pre_max) || pre_max == 0 || any(!is.finite(col_norms)) || all(col_norms <= 1e-12 * pre_max)) {
          stop("Conditional projection collapsed all non-constant W columns for coefficient: ", coef_names[[j]])
        }
      }

      # The Gaussian rejection event is invariant to one common scale of W.
      Wx <- Wx / max(abs(Wx))
      Omsx <- .get_conditional_Oms(spc, Wx, conditional_pairs)
      p_c <- .maxrp(Omsx, q, abs(tau_u) / sqrt(q))$max
      cvx <- .getcv(Omsx, q, 0.05)
      p_final <- max(p_u, p_c)
      cv <- max(cvfin, cvx)
    } else {
      p_final <- p_u
      cv <- cvfin
    }

    out[j, ] <- c(
      coef_j, SE, tau_u, p_final,
      coef_j - cv * SE,
      coef_j + cv * SE
    )

    if (cvs) {
      cvs_vec <- cvs_uncond
      if (!uncond) {
        cvs_cond <- c(
          .getcv(Omsx, q, 0.32),
          .getcv(Omsx, q, 0.10),
          cvx,
          .getcv(Omsx, q, 0.01)
        )
        cvs_vec <- pmax(cvs_vec, cvs_cond)
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
    scpccvs = cvs_mat,
    W = Wfin,
    avc = avc,
    c0 = spc$c0,
    cv = cvfin,
    q = q,
    method = spc$method,
    large_n_seed = large_n_seed,
    call = match.call(),
    backend = "fallback",
    tied = native_tied || isTRUE(spc$tied),
    fallback_reason = fallback_reason
  ), class = "scpc")
}

# ---------------------------------------------------------------------------
# Backend and thread configuration helpers
# ---------------------------------------------------------------------------

.validate_scpc_backend_value <- function(val, source = "`backend`") {
  valid <- c("auto", "rust", "fallback")
  if (!is.character(val) || length(val) != 1L || is.na(val) || !nzchar(val) || !val %in% valid) {
    stop(source, " must be one of \"auto\", \"rust\", or \"fallback\".", call. = FALSE)
  }
  val
}

.resolve_scpc_backend <- function(backend, backend_missing) {
  if (!backend_missing && !is.null(backend) && !identical(backend, c("auto", "rust", "fallback"))) {
    return(.validate_scpc_backend_value(backend, "`backend`"))
  }

  opt <- getOption("spatial.backend")
  if (!is.null(opt)) {
    return(.validate_scpc_backend_value(opt, "Option `spatial.backend`"))
  }

  env <- Sys.getenv("SPATIAL_BACKEND", unset = "")
  if (nzchar(env)) {
    return(.validate_scpc_backend_value(env, "Environment variable `SPATIAL_BACKEND`"))
  }

  "auto"
}

.parse_positive_int_threads <- function(val, source = "`nthreads`") {
  if (is.character(val)) {
    if (length(val) != 1L || is.na(val) || !nzchar(trimws(val)) || !grepl("^[0-9]+$", trimws(val))) {
      stop(source, " must be NULL or a single positive integer.", call. = FALSE)
    }
    val <- suppressWarnings(as.numeric(trimws(val)))
  }
  if (
    !is.numeric(val) ||
    length(val) != 1L ||
    is.na(val) ||
    !is.finite(val) ||
    val < 1 ||
    abs(val - round(val)) > 0 ||
    val > .Machine$integer.max
  ) {
    stop(source, " must be NULL or a single positive integer.", call. = FALSE)
  }
  as.integer(round(val))
}

.scpc_runtime_env <- new.env(parent = emptyenv())
.scpc_runtime_env$init_pid <- NULL
.scpc_runtime_env$accelerate_pid <- NULL

#' Internal Accelerate fork-state helper
#'
#' @param mark Logical; if `TRUE`, record the current process ID as having
#'   initialized Apple Accelerate.
#' @return Logical scalar; `TRUE` if Apple Accelerate was initialized in a
#'   parent process before `fork()`.
#' @keywords internal
#' @export
.scpc_accelerate_state <- function(mark = FALSE) {
  if (isTRUE(mark) && is.null(.scpc_runtime_env$accelerate_pid)) {
    .scpc_runtime_env$accelerate_pid <- Sys.getpid()
  }
  own <- is.loaded("wrap__rust_build_info", PACKAGE = "scpcR") &&
    isTRUE(.Call("wrap__rust_build_info", PACKAGE = "scpcR")$inherited_accelerate_state)
  owner <- .scpc_runtime_env$accelerate_pid
  own || (!is.null(owner) && !identical(owner, Sys.getpid()))
}

.mark_scpc_accelerate_used <- function() {
  if (is.loaded("wrap__rust_build_info", PACKAGE = "scpcR") &&
      identical(.Call("wrap__rust_build_info", PACKAGE = "scpcR")$linalg_provider, "accelerate")) {
    .scpc_accelerate_state(mark = TRUE)
  }
}

.onLoad <- function(libname, pkgname) {
  .scpc_runtime_env$init_pid <- Sys.getpid()
}

.is_scpc_forked_child <- function() {
  init <- .scpc_runtime_env$init_pid
  if (is.null(init)) {
    .scpc_runtime_env$init_pid <- Sys.getpid()
    return(FALSE)
  }
  !identical(Sys.getpid(), init)
}

.is_r_check_limit_cores <- function() {
  val <- trimws(Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = ""))
  nzchar(val) && !tolower(val) %in% c("false", "0", "no")
}

.resolve_scpc_nthreads <- function(nthreads = NULL) {
  resolved <- NULL
  if (!is.null(nthreads)) {
    resolved <- .parse_positive_int_threads(nthreads, "`nthreads`")
  } else if (.is_scpc_forked_child()) {
    return(1L)
  } else {
    opt_threads <- getOption("scpcR.threads")
    if (!is.null(opt_threads)) {
      resolved <- .parse_positive_int_threads(opt_threads, "Option `scpcR.threads`")
    } else {
      env_threads <- Sys.getenv("SCPCR_NUM_THREADS", unset = "")
      if (nzchar(trimws(env_threads))) {
        resolved <- .parse_positive_int_threads(env_threads, "Environment variable `SCPCR_NUM_THREADS`")
      } else {
        opt_mc <- getOption("mc.cores")
        if (!is.null(opt_mc)) {
          resolved <- .parse_positive_int_threads(opt_mc, "Option `mc.cores`")
        }
      }
    }
  }

  if (.is_r_check_limit_cores()) {
    if (is.null(resolved) || resolved > 2L) {
      resolved <- 2L
    }
  }

  resolved
}

.run_scpc_rust_model <- function(coords, latlong, avc, cvs, coef, scores, bread,
                                 design, instruments, direction_design,
                                 cluster_id, fe_id, nthreads, seed) {
  as_double_matrix <- function(x) {
    if (is.null(x)) return(NULL)
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    x
  }
  .mark_scpc_accelerate_used()
  res <- .Call(
    "wrap__rust_run_model_scpc", as_double_matrix(coords), isTRUE(latlong),
    as.double(avc), isTRUE(cvs), as.double(coef), as_double_matrix(scores),
    as_double_matrix(bread), as_double_matrix(design), as_double_matrix(instruments),
    as_double_matrix(direction_design), cluster_id, fe_id, nthreads, seed,
    PACKAGE = "scpcR"
  )
  if (inherits(res, "extendr_error")) {
    if (inherits(res$value, "condition")) stop(res$value)
    stop(res$value, call. = FALSE)
  }
  res
}
