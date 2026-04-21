make_scpc_data <- function(n = 10, with_na = FALSE) {
  set.seed(42)
  x <- stats::rnorm(n)
  y <- 1 + 0.5 * x + stats::rnorm(n, sd = 0.2)
  lon <- seq(-125, -66, length.out = n) + stats::rnorm(n, sd = 0.01)
  lat <- seq(35, 45, length.out = n) + stats::rnorm(n, sd = 0.01)
  if (with_na) {
    x[3] <- NA_real_
  }
  data.frame(
    y = y,
    x = x,
    lon = lon,
    lat = lat,
    cluster = rep(seq_len(ceiling(n / 2)), each = 2)[seq_len(n)]
  )
}

make_python_scpc_data <- function() {
  data.frame(
    y = c(1.0, 1.8, 2.9, 3.7, 5.1),
    x = c(0.0, 1.0, 2.0, 3.0, 4.0),
    coord_x = c(0.0, 1.0, 0.5, 1.5, 2.0),
    coord_y = c(0.0, 0.0, 1.0, 1.0, 1.5)
  )
}

assert_columns_allclose_up_to_sign <- function(left, right, tolerance = ATOL) {
  expect_equal(dim(left), dim(right))
  expect_equal(left[, 1], right[, 1], tolerance = tolerance)
  if (ncol(left) > 1) {
    for (column in 2:ncol(left)) {
      expect_equal(abs(left[, column]), abs(right[, column]), tolerance = tolerance)
    }
  }
}

scpc_example_data_path <- function() {
  candidates <- c(
    system.file("extdata", "chetty_data_1.csv", package = "scpcR"),
    normalizePath(file.path(getwd(), "inst", "extdata", "chetty_data_1.csv"), mustWork = FALSE),
    normalizePath(file.path("..", "..", "inst", "extdata", "chetty_data_1.csv"), mustWork = FALSE),
    normalizePath(file.path(getwd(), "chetty_data_1.csv"), mustWork = FALSE),
    normalizePath(file.path("..", "..", "chetty_data_1.csv"), mustWork = FALSE)
  )
  hits <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(hits) == 0) {
    return(NA_character_)
  }
  hits[[1]]
}

with_fixest_single_thread <- function(code) {
  old_threads <- fixest::getFixest_nthreads()
  fixest::setFixest_nthreads(1)
  on.exit(fixest::setFixest_nthreads(old_threads), add = TRUE)
  eval.parent(substitute(code))
}
