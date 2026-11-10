test_that("the calibration plot driver covers the US only", {
  expect_error(calib(1, c(a = 0), "CA"), "^calib\\(\\) covers the US only; use calib_fit_for")
})
