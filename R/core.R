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
  if (.has_regression_weights(model)) {
    stop("SCPC does not support models fitted with regression weights.", call. = FALSE)
  }

  ## 1. Influence functions --------------------------------------------------
  S <- sandwich::estfun(model)
  coef_vec <- stats::coef(model)
  raw_coef_names <- names(coef_vec)
  is_fixest_iv <- .is_fixest_iv_second_stage(model)
  iv_cond_design <- NULL
  if (is_fixest_iv) {
    iv_cond_design <- .get_fixest_iv_design(model)
    model_mat <- iv_cond_design$model_mat
  } else {
    model_mat <- .get_scpc_model_matrix(model)
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

  cond_include_intercept <- TRUE
  cond_fixef_id <- NULL
  model_mat_cond <- model_mat
  iv_residualizer <- NULL
  if (!isTRUE(uncond) && is_fixest_iv) {
    cond_fixef_id <- iv_cond_design$fixef_id
    if (!is.null(cond_fixef_id)) {
      model_mat_cond <- .demean_for_scpc(
        iv_cond_design$model_mat,
        cond_fixef_id,
        "fixest IV second-stage model matrix"
      )
    }
    iv_residualizer <- .make_iv_residualizer(
      iv_cond_design$X,
      iv_cond_design$Z,
      fixef_id = cond_fixef_id,
      transpose = TRUE
    )
    if (nrow(model_mat_cond) != n || ncol(model_mat_cond) != p) {
      stop(
        "Conditional projection matrix dimensions are incompatible with model coefficients ",
        "(rows = ", nrow(model_mat_cond), ", cols = ", ncol(model_mat_cond),
        ", expected ", n, " x ", p, ")."
      )
    }
  } else if (!isTRUE(uncond)) {
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
  }

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
    coord_by_cl <- rowsum(coords, cl_vec)
    n_per_cl <- as.numeric(table(cl_vec))
    coord_means <- coord_by_cl / n_per_cl
    coord_first <- coords[match(levels(cl_vec), cl_vec), , drop = FALSE]
    if (max(abs(coord_means - coord_first)) > 1e-8) {
      stop(
        "Coordinates vary within clusters. Provide coordinates that are constant ",
        "within each cluster before calling `scpc()`."
      )
    }

    S <- rowsum(S, cl_vec)
    coords <- coord_first
    neff <- nrow(S)
  }

  if (ncol(coords) == 1) coords <- cbind(coords, 0)

  ## 2. Bread & coefficient subset ------------------------------------------
  bread_inv <- sandwich::bread(model) / n
  k_use <- if (is.null(ncoef)) p else min(ncoef, p)
  coef_names <- raw_coef_names[seq_len(k_use)]

  ## 3. Native Rust dispatch (exact unweighted OLS, IV, and 1-way FE)
  fe_is_multiway <- !isTRUE(uncond) &&
    !is.null(cond_fixef_id) &&
    (length(model$fixef_vars) != 1L || !is.null(model$fixef_terms))
  fe_id_0 <- NULL
  if (!isTRUE(uncond) && !is.null(cond_fixef_id) && !fe_is_multiway) {
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
    !fe_is_multiway &&
    !inherits(model, "glm") &&
    (inherits(model, "lm") || inherits(model, "fixest")) &&
    k_use >= 1L
  rust_sym <- if (identical(method_actual, "approx")) {
    "wrap__rust_run_approx_scpc"
  } else if (!isTRUE(uncond) && (!is.null(fe_id_0) || isTRUE(is_fixest_iv))) {
    "wrap__rust_run_exact_scpc_full"
  } else if (is.null(cluster)) {
    "wrap__rust_run_exact_scpc"
  } else {
    "wrap__rust_run_exact_scpc_clustered"
  }
  rust_loaded <- is.loaded(rust_sym, PACKAGE = "scpcR")

  if (identical(backend_resolved, "rust")) {
    if (!rust_loaded) {
      stop("`backend = \"rust\"` was requested, but the native Rust backend is not loaded in `scpcR`.")
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
  } else if (!rust_loaded) {
    "native_unavailable"
  } else {
    "unsupported_model"
  }
  native_tied <- FALSE
  if (!identical(backend_resolved, "fallback") && rust_eligible && rust_loaded) {
    nthreads_resolved <- .resolve_scpc_nthreads(nthreads)
    rust_tied_fallback <- FALSE
    rust_res <- tryCatch(
      .run_scpc_rust_m1(
        coords = coords,
        latlong = latlong,
        avc = avc,
        uncond = uncond,
        cvs = cvs,
        coef_vec = coef_vec,
        S = S,
        bread_inv = bread_inv,
        model_mat_cond = model_mat_cond,
        cond_include_intercept = cond_include_intercept,
        k_use = k_use,
        nthreads = nthreads_resolved,
        cl_vec = cl_vec,
        fe_id_0 = fe_id_0,
        is_fixest_iv = is_fixest_iv,
        iv_cond_design = iv_cond_design,
        cond_fixef_id = cond_fixef_id,
        method_actual = method_actual,
        large_n_seed = large_n_seed
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
  large_n_random_state <- spc$random_state

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
    wj_norm <- sqrt(sum(wj_perm^2))
    denom <- sqrt(sum((t(Wfin[, -1, drop = FALSE]) %*% wj_perm)^2))
    if (!is.finite(denom) || !is.finite(wj_norm) || wj_norm == 0 || denom <= 1e-12 * wj_norm) {
      stop("SCPC variance projection is degenerate for coefficient: ", coef_names[[j]])
    }
    tau_u <- as.numeric(sqrt(q) * crossprod(Wfin[, 1], wj_perm) / denom)
    SE <- as.numeric(denom / (sqrt(q) * sqrt(neff)))
    p_u <- .maxrp(Omsfin, q, abs(tau_u) / sqrt(q))$max

    ## conditional branch (skip when uncond = TRUE) --------------------------
    if (!uncond) {
      cl_vec_scpc <- NULL
      if (!is.null(cluster)) {
        cl_vec_scpc <- if (!identical(perm, seq_len(length(perm)))) {
          factor(as.character(cl_vec), levels = levels(cl_vec)[perm])
        } else {
          cl_vec
        }
      }

      if (!is.null(cluster)) {
        xj_indiv <- as.numeric(neff * bread_inv[j, ] %*% t(model_mat_cond))
        cl_sq <- as.numeric(rowsum(xj_indiv^2, cl_vec_scpc))
        pre_max <- if (ncol(Wfin) > 1L) {
          max(sqrt(colSums(sweep(Wfin[, -1, drop = FALSE], 1, sqrt(cl_sq), `*`)^2)))
        } else {
          0
        }
        if (is_fixest_iv) {
          Wx <- .orthogonalize_W_cluster_iv(
            Wfin = Wfin,
            cl_vec = cl_vec_scpc,
            xj_indiv = xj_indiv,
            residualize = iv_residualizer
          )
        } else {
          ## Clustered conditional: individual-level orthogonalization
          Wx <- .orthogonalize_W_cluster(
            Wfin, cl_vec_scpc, xj_indiv, model_mat_cond,
            include_intercept = cond_include_intercept,
            fixef_id = cond_fixef_id
          )
        }
      } else if (is_fixest_iv) {
        xj_raw <- as.numeric(neff * bread_inv[j, ] %*% t(model_mat_cond))
        if (isTRUE(spc$large_n)) {
          residualize <- .make_iv_residualizer(
            iv_cond_design$X[perm, , drop = FALSE],
            iv_cond_design$Z[perm, , drop = FALSE],
            fixef_id = .permute_fixef_id(iv_cond_design$fixef_id, perm),
            transpose = TRUE
          )
          xj <- xj_raw[perm]
        } else {
          residualize <- iv_residualizer
          xj <- xj_raw
        }
        xjs <- .scpc_score_sign(xj)
        Rx <- sweep(Wfin[, -1, drop = FALSE], 1, xj, `*`)
        pre_max <- if (ncol(Wfin) > 1L) max(sqrt(colSums(Rx^2))) else 0
        Wx <- .orthogonalize_W_iv(
          Wfin = Wfin,
          xj = xj,
          xjs = xjs,
          residualize = residualize
        )
      } else {
        ## Non-clustered conditional
        xj <- as.numeric(neff * bread_inv[j, ] %*% t(model_mat_cond))
        model_mat_cond_j <- model_mat_cond
        if (isTRUE(spc$large_n)) {
          xj <- xj[perm]
          model_mat_cond_j <- model_mat_cond[perm, , drop = FALSE]
          cond_fixef_id_j <- .permute_fixef_id(cond_fixef_id, perm)
        } else {
          cond_fixef_id_j <- cond_fixef_id
        }
        xjs <- .scpc_score_sign(xj)
        Rx <- sweep(Wfin[, -1, drop = FALSE], 1, xj, `*`)
        pre_max <- if (ncol(Wfin) > 1L) max(sqrt(colSums(Rx^2))) else 0
        Wx <- .orthogonalize_W(
          Wfin, xj, xjs, model_mat_cond_j,
          include_intercept = cond_include_intercept,
          fixef_id = cond_fixef_id_j
        )
      }

      if (!all(is.finite(Wx))) {
        stop("Conditional projection produced non-finite W values for coefficient: ", coef_names[[j]])
      }
      if (nrow(Wx) != nrow(Wfin) || ncol(Wx) != ncol(Wfin)) {
        stop("Conditional projection produced W with incompatible dimensions for coefficient: ", coef_names[[j]])
      }
      if (ncol(Wx) > 1L) {
        col_norms <- sqrt(colSums(Wx[, -1, drop = FALSE]^2))
        if (!is.finite(pre_max) || pre_max == 0 || any(!is.finite(col_norms)) || all(col_norms <= 1e-12 * pre_max)) {
          stop("Conditional projection collapsed all non-constant W columns for coefficient: ", coef_names[[j]])
        }
      }

      omsx_res <- .get_conditional_Oms(
        spc,
        Wx,
        latlong = latlong,
        random_state = large_n_random_state
      )
      Omsx <- omsx_res$Oms
      large_n_random_state <- omsx_res$random_state
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

.is_r_check_limit_cores <- function() {
  val <- trimws(Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = ""))
  nzchar(val) && !tolower(val) %in% c("false", "0", "no")
}

.resolve_scpc_nthreads <- function(nthreads = NULL) {
  resolved <- NULL
  if (!is.null(nthreads)) {
    resolved <- .parse_positive_int_threads(nthreads, "`nthreads`")
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

.build_ols_orthonormal_basis <- function(model_mat_cond, include_intercept = TRUE, tol = 1e-10) {
  X <- .coerce_numeric_matrix(model_mat_cond, "OLS residualizer X")$mat
  has_intercept <- !is.null(colnames(X)) && "(Intercept)" %in% colnames(X)
  if (isTRUE(include_intercept) && !has_intercept) {
    X <- cbind("(Intercept)" = 1, X)
  }
  if (ncol(X) == 0L) {
    return(matrix(0, nrow = nrow(X), ncol = 0L))
  }

  qrX <- qr(.scale_design_columns(X), tol = tol)
  if (qrX$rank < ncol(X)) {
    stop("OLS residualizer is rank deficient.", call. = FALSE)
  }
  if (qrX$rank == 0L) {
    return(matrix(0, nrow = nrow(X), ncol = 0L))
  }

  Q_x <- qr.Q(qrX)[, seq_len(qrX$rank), drop = FALSE]
  if (is.null(dim(Q_x))) {
    Q_x <- matrix(Q_x, nrow = nrow(X), ncol = qrX$rank)
  }
  storage.mode(Q_x) <- "double"
  Q_x
}

.run_scpc_rust_m1 <- function(coords,
                              latlong,
                              avc,
                              uncond,
                              cvs,
                              coef_vec,
                              S,
                              bread_inv,
                              model_mat_cond,
                              cond_include_intercept,
                              k_use,
                              nthreads,
                              cl_vec = NULL,
                              fe_id_0 = NULL,
                              is_fixest_iv = FALSE,
                              iv_cond_design = NULL,
                              cond_fixef_id = NULL,
                              method_actual = "exact",
                              large_n_seed = 1) {
  n <- nrow(S)
  coords_mat <- as.matrix(coords)
  storage.mode(coords_mat) <- "double"

  coef_use <- as.double(coef_vec[seq_len(k_use)])
  bread_sub <- bread_inv[seq_len(k_use), , drop = FALSE]
  wj_mat <- sweep(n * (S %*% t(bread_sub)), 2L, coef_use, "+")
  storage.mode(wj_mat) <- "double"

  cl_drop <- if (is.null(cl_vec)) NULL else droplevels(cl_vec)
  cl_id_0 <- if (is.null(cl_drop)) {
    NULL
  } else {
    id <- as.integer(cl_drop) - 1L
    storage.mode(id) <- "integer"
    id
  }

  if (isTRUE(uncond)) {
    rust_xj <- NULL
    rust_xjs <- NULL
    proj_u <- NULL
    proj_v <- NULL
    fe_id_0 <- NULL
  } else {
    xj_mat <- n * (model_mat_cond %*% t(bread_sub))
    storage.mode(xj_mat) <- "double"
    if (is.null(cl_drop)) {
      xjs_mat <- .scpc_score_sign(xj_mat)
    } else {
      cl_sq <- as.matrix(rowsum(xj_mat^2, cl_drop))
      cl_norm <- sqrt(cl_sq)[as.integer(cl_drop), , drop = FALSE]
      xjs_mat <- ifelse(cl_norm > 0, xj_mat / cl_norm, 0)
    }
    storage.mode(xjs_mat) <- "double"

    if (isTRUE(is_fixest_iv)) {
      iv_proj <- .build_iv_projector(
        iv_cond_design$X,
        iv_cond_design$Z,
        fixef_id = cond_fixef_id
      )
      # Denominator columns are Xi_s' M_IV' Xi R, with M_IV = I - U V'.
      proj_u <- iv_proj$V
      proj_v <- iv_proj$U
      rust_xj <- xj_mat
      rust_xjs <- xjs_mat
    } else {
      Q_x <- .build_ols_orthonormal_basis(
        model_mat_cond,
        include_intercept = cond_include_intercept
      )
      proj_u <- Q_x
      proj_v <- Q_x
      rust_xj <- xj_mat
      rust_xjs <- xjs_mat
    }
  }

  if (identical(method_actual, "approx")) {
    norm_s <- .normalize_s(coords_mat, isTRUE(latlong))
    coords_mat <- as.matrix(norm_s$coords)
    storage.mode(coords_mat) <- "double"
    perm <- norm_s$perm
    wj_mat <- wj_mat[perm, , drop = FALSE]
    if (!isTRUE(uncond)) {
      if (!is.null(cl_drop)) {
        cl_perm <- factor(as.character(cl_drop), levels = levels(cl_drop)[perm])
        cl_id_0 <- as.integer(cl_perm) - 1L
        storage.mode(cl_id_0) <- "integer"
      } else {
        rust_xj <- rust_xj[perm, , drop = FALSE]
        rust_xjs <- rust_xjs[perm, , drop = FALSE]
        proj_u <- proj_u[perm, , drop = FALSE]
        proj_v <- proj_v[perm, , drop = FALSE]
        if (!is.null(fe_id_0)) {
          fe_id_0 <- fe_id_0[perm]
        }
      }
    }
  }

  res <- if (identical(method_actual, "approx")) {
    .Call(
      "wrap__rust_run_approx_scpc",
      coords_mat,
      isTRUE(latlong),
      as.double(avc),
      isTRUE(uncond),
      isTRUE(cvs),
      coef_use,
      wj_mat,
      rust_xj,
      rust_xjs,
      proj_u,
      proj_v,
      0.001,
      nthreads,
      cl_id_0,
      fe_id_0,
      as.double(large_n_seed),
      1000L,
      20L,
      1000000L,
      PACKAGE = "scpcR"
    )
  } else if (!isTRUE(uncond) && (!is.null(fe_id_0) || isTRUE(is_fixest_iv))) {
    .Call(
      "wrap__rust_run_exact_scpc_full",
      coords_mat,
      isTRUE(latlong),
      as.double(avc),
      isTRUE(uncond),
      isTRUE(cvs),
      coef_use,
      wj_mat,
      rust_xj,
      rust_xjs,
      proj_u,
      proj_v,
      0.001,
      nthreads,
      cl_id_0,
      fe_id_0,
      PACKAGE = "scpcR"
    )
  } else if (is.null(cl_id_0)) {
    .Call(
      "wrap__rust_run_exact_scpc",
      coords_mat,
      isTRUE(latlong),
      as.double(avc),
      isTRUE(uncond),
      isTRUE(cvs),
      coef_use,
      wj_mat,
      rust_xj,
      rust_xjs,
      proj_u,
      proj_v,
      0.001,
      nthreads,
      PACKAGE = "scpcR"
    )
  } else {
    .Call(
      "wrap__rust_run_exact_scpc_clustered",
      coords_mat,
      isTRUE(latlong),
      as.double(avc),
      isTRUE(uncond),
      isTRUE(cvs),
      coef_use,
      wj_mat,
      rust_xj,
      rust_xjs,
      proj_u,
      proj_v,
      0.001,
      nthreads,
      cl_id_0,
      PACKAGE = "scpcR"
    )
  }

  if (inherits(res, "extendr_error")) {
    if (inherits(res$value, "condition")) stop(res$value)
    stop(res$value, call. = FALSE)
  }
  res
}
