#' Set and get thread count for scpcR
#'
#' Configure the number of threads used for parallel computation in \pkg{scpcR}.
#'
#' @param num_threads Integer >= 1 specifying the number of threads, or \code{"auto"}
#'   / \code{NULL} to restore automatic thread detection.
#'
#' @return \code{scpc_set_num_threads} invisibly returns the previous package option (NULL for automatic detection).
#'   \code{scpc_get_num_threads} returns the currently configured or active thread count.
#'
#' @details
#' \pkg{scpcR} supports multithreading via \pkg{RcppParallel} for distance matrix calculations,
#' subsampling, and grid evaluations.
#'
#' Thread resolution follows this priority:
#' \enumerate{
#'   \item Explicit argument \code{nthreads} in \code{\link{scpc}}
#'   \item Package option \code{options(scpcR.threads = N)} (set directly or via \code{scpc_set_num_threads(N)})
#'   \item Environment variable \code{SCPCR_NUM_THREADS}
#'   \item General R parallel option \code{options(mc.cores = N)}
#'   \item \pkg{RcppParallel} environment variable \code{RCPP_PARALLEL_NUM_THREADS}
#'   \item Automatic hardware detection (all available logical cores)
#' }
#'
#' The CRAN \code{_R_CHECK_LIMIT_CORES_} limit caps every resolved setting at two threads.
#'
#' Setting \code{num_threads = 1} (or \code{options(scpcR.threads = 1)}) completely bypasses
#' background task submission, ensuring pure single-threaded execution. This is recommended
#' when parallelising outer loops with \code{parallel::mclapply}, \pkg{future}, or when running
#' on shared cluster nodes.
#'
#' @examples
#' \dontrun{
#' # Check currently active threads
#' scpc_get_num_threads()
#'
#' # Force single-threaded execution (e.g. inside mclapply or cluster jobs)
#' old_th <- scpc_set_num_threads(1)
#'
#' # Restore previous setting
#' scpc_set_num_threads(old_th)
#'
#' # Or restore default automatic threading
#' scpc_set_num_threads("auto")
#' }
#'
#' @export
scpc_set_num_threads <- function(num_threads = "auto") {
  old <- getOption("scpcR.threads")
  if (is.null(num_threads) || identical(num_threads, "auto")) {
    options(scpcR.threads = NULL)
  } else {
    if (!is.numeric(num_threads) || length(num_threads) != 1L ||
        !is.finite(num_threads) || num_threads < 1 ||
        num_threads > .Machine$integer.max || num_threads != floor(num_threads)) {
      stop("`num_threads` must be a single positive integer, NULL, or \"auto\".")
    }
    options(scpcR.threads = as.integer(num_threads))
  }
  invisible(old)
}

#' @rdname scpc_set_num_threads
#' @export
scpc_get_num_threads <- function() {
  .get_scpc_threads()
}

.get_scpc_threads <- function(nthreads = NULL) {
  valid_count <- function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value) &&
      value >= 1 && value <= .Machine$integer.max && value == floor(value)
  }
  if (!is.null(nthreads) && !valid_count(nthreads)) {
    stop("`nthreads` must be a single positive integer or NULL.")
  }
  count <- nthreads
  if (is.null(count)) count <- getOption("scpcR.threads")
  if (!is.null(count) && !valid_count(count)) {
    stop("`options(scpcR.threads)` must be a single positive integer.")
  }
  if (is.null(count)) {
    env <- Sys.getenv("SCPCR_NUM_THREADS", "")
    if (nzchar(env)) {
      count <- suppressWarnings(as.numeric(env))
      if (!valid_count(count)) stop("`SCPCR_NUM_THREADS` must be a positive integer.")
    }
  }
  if (is.null(count)) count <- getOption("mc.cores")
  if (!is.null(count) && !valid_count(count)) {
    stop("`options(mc.cores)` must be a single positive integer.")
  }
  if (is.null(count)) {
    env <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", "")
    if (nzchar(env)) {
      count <- suppressWarnings(as.numeric(env))
      if (!valid_count(count)) stop("`RCPP_PARALLEL_NUM_THREADS` must be a positive integer.")
    }
  }
  if (is.null(count)) count <- RcppParallel::defaultNumThreads()
  cran_limit <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  if (nzchar(cran_limit) && !cran_limit %in% c("false", "no", "0")) {
    count <- min(count, 2L)
  }
  as.integer(count)
}
