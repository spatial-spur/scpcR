test_that("getnc matches the closed-form grid count", {
  expect_equal(.getnc(1, 4, 2), 2)
})
