test_that("the US, states and other locations resolve by type", {
  expect_identical(resolve_location("US")$loc_type, "national")
  ca <- resolve_location("CA")
  expect_identical(ca$loc_type, "state")
  expect_equal(unname(get_stateID()[ca$st, "USPS"]), "CA")
  other <- resolve_location("XX", data_dir = withr::local_tempdir())
  expect_identical(other[c("name", "code", "loc_type")],
                   list(name = "XX", code = "XX", loc_type = "county"))
  expect_null(other$state)
})

test_that("a location's Location file gives its name and state", {
  dir <- withr::local_tempdir()
  dir.create(file.path(dir, "XX"))
  saveRDS(list(name = "X County", state = "CA"), file.path(dir, "XX", "XX_Location.rds"))
  xx <- resolve_location("XX", data_dir = dir)
  expect_identical(xx[c("name", "loc_type", "state")],
                   list(name = "X County", loc_type = "county", state = "CA"))
})

test_that("loading a location leaves the state table unchanged", {
  before <- get_stateID()
  local_model_load("CA")
  expect_identical(get_stateID(), before)
  expect_identical(nrow(before), 51L)
})
