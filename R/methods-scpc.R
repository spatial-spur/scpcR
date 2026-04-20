# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' @rdname scpc
#' @param x An object of class \code{"scpc"}.
#' @param \dots Further arguments (currently unused).
#' @export
print.scpc <- function(x, ...) {
  cat("\nSCPC Inference (ncoef =", nrow(x$scpcstats), ", q =", x$q, ")\n\n")
  stats::printCoefmat(x$scpcstats[, 1:4], P.values = TRUE, has.Pvalue = TRUE)
  if (!is.null(x$scpccvs)) {
    cat("\nTwo-sided critical values:\n")
    print(x$scpccvs)
  }
  invisible(x)
}

#' @rdname scpc
#' @param object An object of class \code{"scpc"}.
#' @export
summary.scpc <- function(object, ...) {
  cat("\nSCPC Inference (ncoef =", nrow(object$scpcstats), ", q =", object$q,
      ", avc =", object$avc, ")\n\n")
  stats::printCoefmat(object$scpcstats[, 1:4, drop = FALSE],
                      P.values = TRUE, has.Pvalue = TRUE,
                      signif.stars = TRUE)
  cat("\n95% Confidence Intervals:\n")
  print(object$scpcstats[, 5:6, drop = FALSE])
  if (!is.null(object$scpccvs)) {
    cat("\nTwo-sided critical values:\n")
    print(object$scpccvs)
  }
  invisible(object)
}

#' @rdname scpc
#' @export
coef.scpc <- function(object, ...) {
  out <- object$scpcstats[, "Coef"]
  names(out) <- rownames(object$scpcstats)
  out
}

#' @rdname scpc
#' @param parm Character vector of coefficient names, or numeric indices.
#'   Default is all coefficients.
#' @param level Confidence level.  Only 0.95 is available unless
#'   \code{cvs = TRUE} was used, in which case 0.68, 0.90, 0.95, and 0.99
#'   are supported.
#' @export
confint.scpc <- function(object, parm = NULL, level = 0.95, ...) {
  st <- object$scpcstats
  nms <- rownames(st)

  if (is.null(parm)) {
    idx <- seq_len(nrow(st))
  } else if (is.character(parm)) {
    idx <- match(parm, nms)
    if (anyNA(idx)) stop("Unknown coefficient(s): ", paste(parm[is.na(idx)], collapse = ", "))
  } else {
    idx <- parm
  }

  avail <- c(0.68, 0.90, 0.95, 0.99)
  avail_idx <- match(TRUE, abs(avail - level) < 1e-8)

  if (abs(level - 0.95) < 1e-8) {
    ci <- st[idx, c("2.5 %", "97.5 %"), drop = FALSE]
  } else if (!is.null(object$scpccvs) && !is.na(avail_idx)) {
    cv_vals <- object$scpccvs[idx, avail_idx]
    se_vals <- st[idx, "Std_Err"]
    co_vals <- st[idx, "Coef"]
    ci <- cbind(co_vals - cv_vals * se_vals, co_vals + cv_vals * se_vals)
  } else {
    avail_msg <- "0.95"
    if (!is.null(object$scpccvs)) avail_msg <- paste(c(avail_msg, "0.68, 0.90, 0.99"), collapse = ", ")
    stop("Confidence level ", level, " is not available. Available levels: ", avail_msg)
  }

  pct <- format(100 * c((1 - level) / 2, 1 - (1 - level) / 2), digits = 3)
  colnames(ci) <- paste0(pct, " %")
  rownames(ci) <- nms[idx]
  ci
}
