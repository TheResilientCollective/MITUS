test_that("the calibration plot drivers cover the US only", {
  par <- c(a = 0)
  expect_error(calib(1, par, "CA"), "^calib\\(\\) covers the US only; use calib_fit_for")
  expect_error(calib_2020(1, par, "CA"), "^calib_2020\\(\\) covers the US only")
  expect_error(calib_noRB(1, par, "CA"), "^calib_noRB\\(\\) covers the US only")
  expect_error(calib_rblnc(1, par, "CA"), "^calib_rblnc\\(\\) covers the US only")
})
