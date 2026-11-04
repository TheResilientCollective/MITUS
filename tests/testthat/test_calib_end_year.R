# resolve_calib_end_year reads the global wts that model_load sets.
with_weights <- function(last_year, code) {
  g <- globalenv()
  had <- exists("wts", envir = g, inherits = FALSE)
  saved <- if (had) get("wts", envir = g)
  on.exit(if (had) assign("wts", saved, envir = g) else rm("wts", envir = g))
  assign("wts", MITUS:::make_impt_weights(last_year = last_year), envir = g)
  force(code)
}

test_that("calibration simulates through the last weighted year by default", {
  with_weights(2025, expect_identical(MITUS:::resolve_calib_end_year(NULL), 2025L))
})

test_that("a calib_end_year at or after the last weighted year is kept", {
  with_weights(2019, {
    expect_identical(MITUS:::resolve_calib_end_year(2019), 2019L)
    expect_identical(MITUS:::resolve_calib_end_year(2021), 2021L)
  })
})

test_that("a calib_end_year before the last weighted year is an error", {
  with_weights(2025, expect_error(MITUS:::resolve_calib_end_year(2021), "last weighted year"))
})
