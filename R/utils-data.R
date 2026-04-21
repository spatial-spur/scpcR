# ---------------------------------------------------------------------------
# Data and model helpers
# ---------------------------------------------------------------------------

.get_obs_index <- function(model, data) {
  if (inherits(model, "fixest")) {
    if (!requireNamespace("fixest", quietly = TRUE)) {
      stop("Package 'fixest' is required for fixest model objects.")
    }
    obs_index <- fixest::obs(model)
  } else if (inherits(model, "lm")) {
    mf_rows <- rownames(stats::model.frame(model))
    if (!is.null(rownames(data))) {
      obs_index <- match(mf_rows, rownames(data))
      if (all(!is.na(obs_index))) {
        return(obs_index)
      }
    }

    obs_index <- suppressWarnings(as.integer(mf_rows))
    if (
      all(!is.na(obs_index)) &&
      all(obs_index >= 1L) &&
      all(obs_index <= nrow(data))
    ) {
      return(obs_index)
    }
    stop("Could not align model estimation sample to data rows.")
  } else {
    stop("Model class not supported. Provide an lm or fixest model.")
  }

  if (any(obs_index < 1L | obs_index > nrow(data))) {
    stop("Model observation indices are outside the row range of `data`.")
  }
  obs_index
}

.is_fixest_iv_second_stage <- function(model) {
  # fix: current fixest marks iv models with is_iv, so use that flag here.
  if (!inherits(model, "fixest") || !isTRUE(model$is_iv)) {
    return(FALSE)
  }
  stage <- tryCatch(as.integer(model$iv_stage[[1L]]), error = function(e) NA_integer_)
  identical(stage, 2L)
}

.align_scpc_model_matrix <- function(model_mat, coef_names, context) {
  model_mat <- as.matrix(model_mat)

  if (length(coef_names) == 0L) {
    stop(context, " has no coefficients to align.")
  }

  if (ncol(model_mat) != length(coef_names)) {
    if (
      !is.null(colnames(model_mat)) &&
      all(coef_names %in% colnames(model_mat))
    ) {
      model_mat <- model_mat[, coef_names, drop = FALSE]
    } else {
      stop(
        context, " column count cannot be aligned to the coefficient vector ",
        "(matrix columns = ", ncol(model_mat),
        ", coefficients = ", length(coef_names), ")."
      )
    }
  } else if (is.null(colnames(model_mat))) {
    stop(
      context, " has no column names, so its columns cannot be aligned ",
      "unambiguously to the coefficient vector."
    )
  } else if (!identical(colnames(model_mat), coef_names)) {
    if (all(coef_names %in% colnames(model_mat))) {
      model_mat <- model_mat[, coef_names, drop = FALSE]
    } else {
      stop(
        context, " column names do not align with model coefficients.\n",
        "Model matrix columns: ", paste(colnames(model_mat), collapse = ", "), "\n",
        "Coefficient names: ", paste(coef_names, collapse = ", ")
      )
    }
  }

  if (ncol(model_mat) != length(coef_names)) {
    stop(
      context, " remains misaligned after column matching ",
      "(matrix columns = ", ncol(model_mat),
      ", coefficients = ", length(coef_names), ")."
    )
  }

  model_mat
}

.get_fixest_matrix <- function(model, type, context) {
  mm <- tryCatch(
    stats::model.matrix(model, type = type),
    error = function(e) e
  )
  if (inherits(mm, "error")) {
    stop("Could not extract ", context, ": ", conditionMessage(mm))
  }
  as.matrix(mm)
}

.ensure_iv_intercept <- function(mat, n, include_intercept, context) {
  mat <- as.matrix(mat)
  cn <- colnames(mat)
  has_intercept <- !is.null(cn) && "(Intercept)" %in% cn

  if (isTRUE(include_intercept) && !has_intercept) {
    mat <- cbind("(Intercept)" = rep(1, n), mat)
  }
  if (!isTRUE(include_intercept) && has_intercept) {
    keep <- cn != "(Intercept)"
    mat <- mat[, keep, drop = FALSE]
  }

  cn <- colnames(mat)
  if (!is.null(cn) && sum(cn == "(Intercept)") > 1L) {
    stop(context, " contains multiple intercept columns.")
  }

  mat
}

.get_scpc_model_matrix <- function(model) {
  coef_names <- names(stats::coef(model))
  if (.is_fixest_iv_second_stage(model)) {
    mm <- .get_fixest_matrix(model, "iv.rhs2", "the fixest IV second-stage model matrix")
  } else {
    mm <- stats::model.matrix(model)
  }
  .align_scpc_model_matrix(mm, coef_names, "SCPC model matrix")
}

.has_fixest_fe <- function(model) {
  inherits(model, "fixest") &&
    !is.null(model$fixef_vars) &&
    length(model$fixef_vars) > 0L
}

.get_fixest_iv_design <- function(model) {
  if (!.is_fixest_iv_second_stage(model)) {
    stop("`.get_fixest_iv_design()` requires a fixest IV second-stage model.")
  }

  coef_names <- names(stats::coef(model))
  model_mat <- .align_scpc_model_matrix(
    .get_fixest_matrix(model, "iv.rhs2", "the fixest IV second-stage model matrix"),
    coef_names,
    "fixest IV second-stage model matrix"
  )

  exo <- .get_fixest_matrix(model, "iv.exo", "the fixest IV exogenous regressor matrix")
  endo <- .get_fixest_matrix(model, "iv.endo", "the fixest IV endogenous regressor matrix")
  inst <- .get_fixest_matrix(model, "iv.inst", "the fixest IV excluded-instrument matrix")

  n <- nrow(model_mat)
  include_intercept <- !is.null(colnames(model_mat)) &&
    "(Intercept)" %in% colnames(model_mat)

  exo <- .ensure_iv_intercept(
    exo,
    n = n,
    include_intercept = include_intercept,
    context = "The fixest IV exogenous regressor matrix"
  )

  if (!all(c(nrow(exo), nrow(endo), nrow(inst)) == n)) {
    stop(
      "Fixest IV design matrices do not share a common row count ",
      "(model_mat = ", n, ", exo = ", nrow(exo),
      ", endo = ", nrow(endo), ", inst = ", nrow(inst), ")."
    )
  }

  X <- cbind(exo, endo)
  Z <- cbind(exo, inst)

  if (nrow(X) != n || nrow(Z) != n) {
    stop("Internal error: fixest IV design matrices have incompatible row counts.")
  }

  if (any(!is.finite(model_mat)) || any(!is.finite(X)) || any(!is.finite(Z))) {
    stop("Fixest IV design extraction produced non-finite values.")
  }

  list(
    X = X,
    Z = Z,
    model_mat = model_mat,
    coef_names = coef_names,
    fixef_id = if (.has_fixest_fe(model)) model$fixef_id else NULL,
    has_fixef = .has_fixest_fe(model)
  )
}

.get_conditional_projection_setup <- function(model, model_mat, n, uncond) {
  setup <- list(model_mat = as.matrix(model_mat), include_intercept = TRUE, fixef_id = NULL)
  if (isTRUE(uncond) || !inherits(model, "fixest") || !.has_fixest_fe(model)) {
    return(setup)
  }

  mm <- .align_scpc_model_matrix(
    model_mat,
    names(stats::coef(model)),
    "Conditional projection model matrix"
  )

  if (nrow(mm) != n) {
    stop("Internal error: fixest conditional regressors have incompatible row count.")
  }

  mm <- tryCatch(
    as.matrix(fixest::demean(mm, f = model$fixef_id, nthreads = 1L)),
    error = function(e) e
  )
  if (inherits(mm, "error")) {
    stop("Could not FE-demean conditional regressors for fixest model: ", conditionMessage(mm))
  }
  list(model_mat = mm, include_intercept = FALSE, fixef_id = model$fixef_id)
}

.resolve_coords_input <- function(data, obs_index, lon, lat, coords_euclidean) {
  use_geodesic <- !is.null(lon) || !is.null(lat)
  use_euclidean <- !is.null(coords_euclidean)

  if (use_geodesic && use_euclidean) {
    stop("Specify either `lon`/`lat` or `coords_euclidean`, not both.")
  }
  if (!use_geodesic && !use_euclidean) {
    stop("Specify coordinates via `lon`/`lat` or `coords_euclidean`.")
  }

  if (use_geodesic) {
    if (is.null(lon) || is.null(lat)) {
      stop("For geodesic coordinates, provide both `lon` and `lat`.")
    }
    if (!is.character(lon) || length(lon) != 1L || !nzchar(lon) ||
        !is.character(lat) || length(lat) != 1L || !nzchar(lat)) {
      stop("`lon` and `lat` must each be a single column name.")
    }
    miss <- setdiff(c(lon, lat), names(data))
    if (length(miss) > 0) {
      stop("Coordinate variables not found in data: ", paste(miss, collapse = ", "))
    }
    coords <- data[obs_index, c(lon, lat), drop = FALSE]
    if (!all(vapply(coords, is.numeric, logical(1)))) {
      stop("`lon` and `lat` must reference numeric columns.")
    }
    if (any(!is.finite(as.matrix(coords)))) {
      stop("Geodesic coordinates must be finite.")
    }
    if (any(coords[[lon]] < -180 | coords[[lon]] > 180)) {
      stop("Longitude values must be in [-180, 180].")
    }
    if (any(coords[[lat]] < -90 | coords[[lat]] > 90)) {
      stop("Latitude values must be in [-90, 90].")
    }
    return(list(coords = as.matrix(coords), latlong = TRUE))
  }

  if (!is.character(coords_euclidean) || length(coords_euclidean) < 1L) {
    stop("`coords_euclidean` must be a character vector with at least one column name.")
  }
  miss <- setdiff(coords_euclidean, names(data))
  if (length(miss) > 0) {
    stop("Coordinate variables not found in data: ", paste(miss, collapse = ", "))
  }
  coords <- data[obs_index, coords_euclidean, drop = FALSE]
  if (!all(vapply(coords, is.numeric, logical(1)))) {
    stop("`coords_euclidean` columns must be numeric.")
  }
  if (any(!is.finite(as.matrix(coords)))) {
    stop("Euclidean coordinates must be finite.")
  }
  list(coords = as.matrix(coords), latlong = FALSE)
}
