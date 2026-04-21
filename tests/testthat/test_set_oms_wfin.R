test_that("setOmsWfin can force the exact path", {
  coords <- cbind(
    coord_x = c(0, 1, 0.5, 1.5, 2.0),
    coord_y = c(0, 0, 1.0, 1.0, 1.5)
  )

  out <- .setOmsWfin(
    coords,
    avc0 = 0.1,
    latlong = FALSE,
    method = "exact",
    large_n_seed = 1
  )

  expect_type(out, "list")
  expect_true(all(
    c("Wfin", "cvfin", "Omsfin", "c0", "cmax", "coords", "perm", "distmat", "method", "large_n", "random_state") %in%
      names(out)
  ))
  expect_identical(out$method, "exact")
  expect_false(out$large_n)
  expect_null(out$random_state)
  expect_equal(out$perm, seq_len(nrow(coords)))
  expect_equal(out$coords, coords)
  expect_equal(out$distmat, .getdistmat(coords, latlong = FALSE), tolerance = 1e-12)
  expect_equal(nrow(out$Wfin), nrow(coords))
  expect_true(ncol(out$Wfin) >= 2)
  expect_gt(out$cvfin, 0)
  expect_gt(out$c0, 0)
  expect_gt(out$cmax, out$c0)
  expect_true(length(out$Omsfin) >= 2)
})

test_that("setOmsWfin can force the large-n approximation path", {
  coords <- cbind(
    coord_x = seq(0, 1, length.out = 40),
    coord_y = seq(1, 0, length.out = 40)
  )

  out <- .setOmsWfin(
    coords,
    avc0 = 0.1,
    latlong = FALSE,
    method = "approx",
    large_n_seed = 1
  )

  expect_identical(out$method, "approx")
  expect_true(out$large_n)
  expect_null(out$distmat)
  expect_equal(sort(out$perm), seq_len(nrow(coords)))
  expect_equal(nrow(out$coords), nrow(coords))
  expect_equal(nrow(out$Wfin), nrow(coords))
  expect_true(all(is.finite(out$coords)))
  expect_true(all(is.finite(out$Wfin)))
  expect_true(all(vapply(out$Omsfin, function(Om) all(is.finite(Om)), logical(1))))
  expect_true(is.numeric(out$random_state))
  expect_length(out$random_state, 1)
  expect_true(is.finite(out$random_state))
})

test_that("setOmsWfin auto chooses exact on small problems", {
  coords <- cbind(
    coord_x = c(0, 1, 0.5, 1.5, 2.0),
    coord_y = c(0, 0, 1.0, 1.0, 1.5)
  )

  out <- .setOmsWfin(
    coords,
    avc0 = 0.1,
    latlong = FALSE,
    method = "auto",
    large_n_seed = 1
  )

  expect_identical(out$method, "exact")
  expect_false(out$large_n)
  expect_equal(out$perm, seq_len(nrow(coords)))
})

test_that("resolve_scpc_method routes auto at the large-n threshold", {
  expect_identical(.resolve_scpc_method(4499L, "auto"), "exact")
  expect_identical(.resolve_scpc_method(4500L, "auto"), "approx")
  expect_identical(.resolve_scpc_method(100L, "exact"), "exact")
  expect_identical(.resolve_scpc_method(100L, "approx"), "approx")
})

test_that("setOmsWfin validates method and large_n_seed", {
  coords <- cbind(
    coord_x = c(0, 1, 0.5, 1.5, 2.0),
    coord_y = c(0, 0, 1.0, 1.0, 1.5)
  )

  expect_error(
    .setOmsWfin(coords, avc0 = 0.1, latlong = FALSE, method = "bad", large_n_seed = 1),
    "`method` must be one of \"auto\", \"exact\", or \"approx\""
  )
  expect_error(
    .setOmsWfin(coords, avc0 = 0.1, latlong = FALSE, method = "exact", large_n_seed = 1.5),
    "`large_n_seed` must be a single integer-valued number in \\[0, 2\\^32\\)"
  )
})
