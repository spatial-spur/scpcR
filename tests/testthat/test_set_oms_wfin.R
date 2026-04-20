test_that("setOmsWfin returns a coherent spatial setup", {
  distmat <- matrix(
    c(0, 1, 2,
      1, 0, 1,
      2, 1, 0),
    nrow = 3,
    byrow = TRUE
  )

  out <- .setOmsWfin(distmat, 0.1)

  expect_type(out, "list")
  expect_true(all(c("Wfin", "cvfin", "Omsfin", "c0", "cmax") %in% names(out)))
  expect_equal(nrow(out$Wfin), nrow(distmat))
  expect_true(ncol(out$Wfin) >= 2)
  expect_gt(out$cvfin, 0)
  expect_gt(out$c0, 0)
  expect_gt(out$cmax, out$c0)
  expect_true(length(out$Omsfin) >= 2)
})
