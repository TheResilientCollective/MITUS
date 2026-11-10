test_that("model_load reads target globals for the US only", {
  local_model_load("US")
  expect_true(all(exists("CalibDat", envir = globalenv()), exists("CalibDatCases", envir = globalenv())))

  local_model_load("CA")
  expect_false(any(vapply(c("CalibDat", "CalibDatState", "CalibDatCases", "DeathByAge"),
                          exists, TRUE, envir = globalenv())))
})

test_that("model_load_demo covers the US only", {
  expect_error(model_load_demo("CA"), "^model_load_demo\\(\\) covers the US only$")
})
