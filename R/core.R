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
#' @param uncond Logical; if \code{TRUE}, report unconditional critical
#'   values only (skip the conditional adjustment of Mueller and Watson
#'   2023).  Default is \code{FALSE}.
#' @param cvs Logical; if \code{TRUE}, include per-coefficient critical
#'   values at the 32\%, 10\%, 5\%, and 1\% levels.  Default is
#'   \code{FALSE}.
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
                 uncond = FALSE,
                 cvs = FALSE) {

  if (!is.numeric(avc) || length(avc) != 1L || is.na(avc) || !is.finite(avc)) {
    stop("`avc` must be a single finite numeric value.")
  }
  if (avc <= 0.001 || avc >= 0.99) {
    stop("`avc` must lie in (0.001, 0.99).")
  }
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
  D <- .getdistmat(coords, latlong)
  spc <- .setOmsWfin(D, avc)
  Wfin <- spc$Wfin
  cvfin <- spc$cvfin
  Omsfin <- spc$Omsfin
  q <- ncol(Wfin) - 1

  ## 3. Bread ----------------------------------------------------------------
  bread_inv <- sandwich::bread(model) / n

  ## 4. Loop over coefficients ------------------------------------------------
  k_use <- if (is.null(ncoef)) p else min(ncoef, p)
  out <- matrix(NA_real_, k_use, 6)
  levs <- c(0.32, 0.10, 0.05, 0.01)
  cvs_mat <- if (cvs) matrix(NA_real_, k_use, 4) else NULL
  cvs_uncond <- if (cvs) vapply(levs, function(lv) .getcv(Omsfin, q, lv), 0.0) else NULL
  coef_names <- raw_coef_names[seq_len(k_use)]
  rownames(out) <- coef_names
  colnames(out) <- c("Coef", "Std_Err", "t", "P>|t|", "2.5 %", "97.5 %")

  for (j in seq_len(k_use)) {
    coef_j <- coef_vec[[j]]
    wj <- as.numeric(neff * bread_inv[j, ] %*% t(S)) + coef_j

    ## unconditional statistic -----------------------------------------------
    denom <- sqrt(sum((t(Wfin[, -1, drop = FALSE]) %*% wj)^2))
    if (!is.finite(denom) || denom <= sqrt(.Machine$double.eps)) {
      stop("SCPC variance projection is degenerate for coefficient: ", coef_names[[j]])
    }
    tau_u <- as.numeric(sqrt(q) * crossprod(Wfin[, 1], wj) / denom)
    SE <- as.numeric(denom / (sqrt(q) * sqrt(neff)))
    p_u <- .maxrp(Omsfin, q, abs(tau_u) / sqrt(q))$max

    ## conditional branch (skip when uncond = TRUE) --------------------------
    if (!uncond) {
      if (!is.null(cluster)) {
        xj_indiv <- as.numeric(neff * bread_inv[j, ] %*% t(model_mat_cond))
        if (is_fixest_iv) {
          Wx <- .orthogonalize_W_cluster_iv(
            Wfin = Wfin,
            cl_vec = cl_vec,
            xj_indiv = xj_indiv,
            residualize = iv_residualizer
          )
        } else {
          ## Clustered conditional: individual-level orthogonalization
          Wx <- .orthogonalize_W_cluster(
            Wfin, cl_vec, xj_indiv, model_mat_cond,
            include_intercept = cond_include_intercept
          )
        }
      } else if (is_fixest_iv) {
        xj <- as.numeric(neff * bread_inv[j, ] %*% t(model_mat_cond))
        xjs <- sign(xj)
        Wx <- .orthogonalize_W_iv(
          Wfin = Wfin,
          xj = xj,
          xjs = xjs,
          residualize = iv_residualizer
        )
      } else {
        ## Non-clustered conditional
        xj <- as.numeric(neff * bread_inv[j, ] %*% t(model_mat_cond))
        xjs <- sign(xj)
        Wx <- .orthogonalize_W(
          Wfin, xj, xjs, model_mat_cond,
          include_intercept = cond_include_intercept,
          fixef_id = cond_fixef_id
        )
      }

      if (!all(is.finite(Wx))) {
        stop("Conditional projection produced non-finite W values for coefficient: ", coef_names[[j]])
      }
      if (nrow(Wx) != nrow(D) || ncol(Wx) != ncol(Wfin)) {
        stop("Conditional projection produced W with incompatible dimensions for coefficient: ", coef_names[[j]])
      }
      if (
        ncol(Wx) > 1L &&
        all(colSums(abs(Wx[, -1, drop = FALSE])) <= sqrt(.Machine$double.eps))
      ) {
        stop("Conditional projection collapsed all non-constant W columns for coefficient: ", coef_names[[j]])
      }

      Omsx <- .getOms(D, spc$c0, spc$cmax, Wx, 1.2)
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
        cvs_cond <- vapply(levs, function(lv) .getcv(Omsx, q, lv), 0.0)
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
    call = match.call()
  ), class = "scpc")
}
