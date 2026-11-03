restore_stateID <- function() {
  env <- MITUS:::.mitus_state
  if (exists("stateID", envir = env, inherits = FALSE)) {
    saved <- get("stateID", envir = env)
    return(function() assign("stateID", saved, envir = env))
  }
  function() rm(list = intersect("stateID", ls(env)), envir = env)
}

test_that("a registered county still resolves as a county", {
  undo <- restore_stateID()
  on.exit(undo())
  first <- resolve_location("SanDiego")
  expect_equal(first$loc_type, "county")
  st <- register_location(first)

  again <- resolve_location("SanDiego")
  expect_equal(again$loc_type, "county")
  expect_equal(register_location(again), st)
})

test_that("states resolve as states", {
  ca <- resolve_location("CA")
  expect_equal(ca$loc_type, "state")
  expect_equal(unname(get_stateID()[ca$st, "USPS"]), "CA")
})
