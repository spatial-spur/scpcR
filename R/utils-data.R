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
    mf <- stats::model.frame(model)
    rn <- rownames(mf)
    suppressWarnings(obs_index <- as.integer(rn))
    if (anyNA(obs_index)) {
      obs_index <- match(rn, rownames(data))
      if (anyNA(obs_index)) {
        stop("Could not map lm model frame rows back to `data`.")
      }
    }
  } else {
    stop("Model class not supported. Provide an lm or fixest model.")
  }

  if (any(obs_index < 1L | obs_index > nrow(data))) {
    stop("Model observation indices are outside the row range of `data`.")
  }
  obs_index
}

.get_scpc_model_matrix <- function(model) {
  if (inherits(model, "fixest") && isTRUE(model$iv) && isTRUE(model$iv_stage == 2)) {
    mm <- tryCatch(
      utils::getFromNamespace("model.matrix.fixest", "fixest")(model, type = "iv.rhs2"),
      error = function(e) NULL
    )
    if (!is.null(mm)) {
      return(mm)
    }
  }
  stats::model.matrix(model)
}

.has_fixest_fe <- function(model) {
  inherits(model, "fixest") &&
    !is.null(model$fixef_vars) &&
    length(model$fixef_vars) > 0L
}

.get_conditional_projection_setup <- function(model, model_mat, n, uncond) {
  setup <- list(model_mat = as.matrix(model_mat), include_intercept = TRUE, fixef_id = NULL)
  if (isTRUE(uncond) || !inherits(model, "fixest") || !.has_fixest_fe(model)) {
    return(setup)
  }

  mm <- as.matrix(model_mat)
  coef_names <- names(stats::coef(model))
  if (!is.null(colnames(mm)) && all(coef_names %in% colnames(mm))) {
    mm <- mm[, coef_names, drop = FALSE]
  } else if (ncol(mm) == length(coef_names)) {
    colnames(mm) <- coef_names
  } else {
    stop(
      "Could not align demeaned fixest regressors with model coefficients ",
      "(coef count = ", length(coef_names), ", demeaned columns = ", ncol(mm), ")."
    )
  }

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
