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
#' @param nthreads Optional integer; number of threads to use for parallel
#'   computations. Defaults to \code{NULL}, which respects package options,
#'   environment variables, or automatic detection.
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
#' }
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
                 nthreads = NULL) {

  if (!is.null(nthreads)) {
    nth <- .get_scpc_threads(nthreads)
    old_opt <- options(scpcR.threads = nth)
    on.exit(options(old_opt), add = TRUE)
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
      fixef_id = cond_fixef_id
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

  if (!is.null(cluster)) {
    if (!cluster %in% names(data)) {
      stop("Cluster variable not found in data: ", cluster)
    }
    if (!is_fixest_iv && !is.null(cond_fixef_id) && !isTRUE(uncond)) {
      stop("Conditional SCPC with absorbed fixed effects and external clustering is not yet implemented; use `uncond = TRUE`.")
    }
    cl_vec <- factor(data[[cluster]][obs_index])

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

  ## 2. Spatial kernel (unconditional) ---------------------------------------
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

  ## 3. Bread ----------------------------------------------------------------
  bread_inv <- sandwich::bread(model) / n

  ## 4. Loop over coefficients ------------------------------------------------
  k_use <- if (is.null(ncoef)) p else min(ncoef, p)
  out <- matrix(NA_real_, k_use, 6)
  levs <- c(0.32, 0.10, 0.05, 0.01)
  cvs_mat <- if (cvs) matrix(NA_real_, k_use, 4) else NULL
  cvs_uncond <- if (cvs) .getcv(Omsfin, q, levs) else NULL
  coef_names <- raw_coef_names[seq_len(k_use)]
  rownames(out) <- coef_names
  colnames(out) <- c("Coef", "Std_Err", "t", "P>|t|", "2.5 %", "97.5 %")

  bread_sub <- neff * bread_inv[seq_len(k_use), , drop = FALSE]
  W_mat <- tcrossprod(S, bread_sub)
  X_mat <- if (!uncond) tcrossprod(model_mat_cond, bread_sub) else NULL

  cl_vec_scpc <- NULL
  ols_residualizer <- NULL
  iv_residualizer_use <- iv_residualizer
  model_mat_cond_use <- model_mat_cond
  cond_fixef_id_use <- cond_fixef_id

  if (!uncond) {
    if (!is.null(cluster)) {
      cl_vec_scpc <- if (!identical(perm, seq_len(length(perm)))) {
        factor(as.character(cl_vec), levels = levels(cl_vec)[perm])
      } else {
        cl_vec
      }
      if (!is_fixest_iv && ncol(Wfin) > 1L) {
        ols_residualizer <- .make_ols_residualizer(
          model_mat_cond,
          include_intercept = cond_include_intercept
        )
      }
    } else if (is_fixest_iv) {
      if (isTRUE(spc$large_n)) {
        iv_residualizer_use <- .make_iv_residualizer(
          iv_cond_design$X[perm, , drop = FALSE],
          iv_cond_design$Z[perm, , drop = FALSE],
          fixef_id = .permute_fixef_id(iv_cond_design$fixef_id, perm)
        )
      }
    } else {
      if (isTRUE(spc$large_n)) {
        model_mat_cond_use <- model_mat_cond[perm, , drop = FALSE]
        cond_fixef_id_use <- .permute_fixef_id(cond_fixef_id, perm)
      }
      if (ncol(Wfin) > 1L) {
        ols_residualizer <- .make_ols_residualizer(
          model_mat_cond_use,
          include_intercept = cond_include_intercept,
          fixef_id = cond_fixef_id_use
        )
      }
    }
  }

  Wfin_1 <- Wfin[, 1]
  Wfin_rest <- Wfin[, -1, drop = FALSE]
  sqrt_q <- sqrt(q)
  sqrt_neff <- sqrt(neff)

  tau_u_vec <- numeric(k_use)
  se_vec <- numeric(k_use)
  p_u_vec <- numeric(k_use)
  Wx_list <- if (!uncond) vector("list", k_use) else NULL

  for (j in seq_len(k_use)) {
    coef_j <- coef_vec[[j]]
    wj <- W_mat[, j] + coef_j
    wj_perm <- wj[perm]

    ## unconditional statistic -----------------------------------------------
    denom <- sqrt(sum(crossprod(Wfin_rest, wj_perm)^2))
    if (!is.finite(denom) || denom <= sqrt(.Machine$double.eps)) {
      stop("SCPC variance projection is degenerate for coefficient: ", coef_names[[j]])
    }
    tau_u <- as.numeric(sqrt_q * sum(Wfin_1 * wj_perm) / denom)
    SE <- as.numeric(denom / (sqrt_q * sqrt_neff))
    p_u <- .maxrp(Omsfin, q, abs(tau_u) / sqrt_q)$max

    tau_u_vec[j] <- tau_u
    se_vec[j] <- SE
    p_u_vec[j] <- p_u

    ## conditional branch (skip when uncond = TRUE) --------------------------
    if (!uncond) {
      if (!is.null(cluster)) {
        xj_indiv <- X_mat[, j]
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
            residualize = ols_residualizer
          )
        }
      } else if (is_fixest_iv) {
        xj_raw <- X_mat[, j]
        xj <- if (isTRUE(spc$large_n)) xj_raw[perm] else xj_raw
        xjs <- sign(xj)
        Wx <- .orthogonalize_W_iv(
          Wfin = Wfin,
          xj = xj,
          xjs = xjs,
          residualize = iv_residualizer_use
        )
      } else {
        ## Non-clustered conditional
        xj <- X_mat[, j]
        if (isTRUE(spc$large_n)) {
          xj <- xj[perm]
        }
        xjs <- sign(xj)
        Wx <- .orthogonalize_W(
          Wfin, xj, xjs, model_mat_cond_use,
          include_intercept = cond_include_intercept,
          fixef_id = cond_fixef_id_use,
          residualize = ols_residualizer
        )
      }

      if (!all(is.finite(Wx))) {
        stop("Conditional projection produced non-finite W values for coefficient: ", coef_names[[j]])
      }
      if (nrow(Wx) != nrow(Wfin) || ncol(Wx) != ncol(Wfin)) {
        stop("Conditional projection produced W with incompatible dimensions for coefficient: ", coef_names[[j]])
      }
      if (
        ncol(Wx) > 1L &&
        all(colSums(abs(Wx[, -1, drop = FALSE])) <= sqrt(.Machine$double.eps))
      ) {
        stop("Conditional projection collapsed all non-constant W columns for coefficient: ", coef_names[[j]])
      }

      storage.mode(Wx) <- "double"
      Wx_list[[j]] <- Wx
    }
  }

  Omsx_batch <- NULL
  if (!uncond && !isTRUE(spc$large_n)) {
    storage.mode(spc$distmat) <- "double"
    Omsx_batch <- rcpp_get_oms_batch(
      spc$distmat,
      as.double(spc$c0),
      as.double(spc$cmax),
      Wx_list,
      1.2,
      as.integer(.get_scpc_threads())
    )
  }

  for (j in seq_len(k_use)) {
    coef_j <- coef_vec[[j]]
    SE <- se_vec[j]
    tau_u <- tau_u_vec[j]
    p_u <- p_u_vec[j]

    if (!uncond) {
      if (isTRUE(spc$large_n)) {
        omsx_res <- .get_conditional_Oms(
          spc,
          Wx_list[[j]],
          latlong = latlong,
          random_state = large_n_random_state
        )
        Omsx <- omsx_res$Oms
        large_n_random_state <- omsx_res$random_state
      } else {
        Omsx <- Omsx_batch[[j]]
      }
      p_c <- .maxrp(Omsx, q, abs(tau_u) / sqrt_q)$max
      if (cvs) {
        cvs_cond <- .getcv(Omsx, q, levs)
        cvx <- cvs_cond[[3L]]
      } else {
        cvx <- .getcv(Omsx, q, 0.05)
      }
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
    call = match.call()
  ), class = "scpc")
}
