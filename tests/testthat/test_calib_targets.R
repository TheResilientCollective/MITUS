states <- c("Alabama", "Alaska")

# All-state calibration data for two states, with each state's values marked
# by its position.
fake_all_states <- function() {
  per_state <- function(make) lapply(seq_along(states), make)
  list(
    cases_yr_st = per_state(function(i) cbind(year = 1993:1995, cases = 100 * i + 0:2)),
    cases_yr_ag_nat_st = per_state(function(i) array(i, c(2, 3, 2))),
    cases_yr_ag_nat_st_5yr = per_state(function(i) data.frame(usb = 0:1, X0.4 = i)),
    TLTBI_volume_state = per_state(function(i) c(i, i - 0.5, i + 0.5)),
    tbdeaths = per_state(function(i) data.frame(State = states[i], Year = 1999, Deaths = i)),
    homeless_pop = per_state(function(i) rep(i, 3)),
    pop_50_10 = per_state(function(i) matrix(i)),
    pop_00_19 = per_state(function(i) matrix(i)),
    cases_nat_st_5yr = data.frame(
      State = rep(states, each = 2), State.Code = rep(1:2, each = 2), usb = c(0, 1, 0, 1),
      X1996.2000 = 1:4
    ),
    hr_cases_sm = data.frame(st = factor(rep(states, each = 2)), fract_sm = 1:4 / 10),
    rt_fb_cases_sm = data.frame(st = factor(rep(states, each = 2)), p_recent = 1:4 / 10),
    rct_cases_sm = data.frame(st = states, fract_sm = c(0.2, 0.3)),
    cases_prop_change_53_94 = c(1.1, 0.9),
    tx_outcomes = "national tx_outcomes",
    TLTBI_dist = c(FB = 0.5, HR = 0.1),
    LTBI_prev_US_11_IGRA = "national LTBI_prev_US_11_IGRA",
    LTBI_prev_FB_11_IGRA = "national LTBI_prev_FB_11_IGRA",
    tbdeaths_age_yr = "national tbdeaths_age_yr",
    borgdorff_data = "borgdorff",
    ferebee_data = "ferebee",
    sutherland_data = "sutherland"
  )
}

# One county's calibration data, which provides two of the national targets.
fake_county <- function() {
  list(
    cases_yr_st = cbind(year = 1993, cases = 7),
    cases_yr_ag_nat_st = array(7, c(2, 3, 2)),
    cases_yr_ag_nat_st_5yr = data.frame(usb = 0:1, X0.4 = 7),
    TLTBI_volume_state = c(7, 6, 8),
    tbdeaths = data.frame(State = "County", Year = 1999, Deaths = 7),
    homeless_pop = rep(7, 3),
    pop_50_10 = matrix(7),
    pop_00_19 = matrix(7),
    cases_nat_st_5yr = data.frame(State = "County", State.Code = 52, usb = 0:1, X1996.2000 = 7:8),
    hr_cases_sm = data.frame(st = "County", fract_sm = 0.7),
    rt_fb_cases_sm = data.frame(st = "County", p_recent = 0.7),
    rct_cases_sm = data.frame(st = "County", fract_sm = 0.7),
    tx_outcomes = "county tx_outcomes",
    tbdeaths_age_yr = "county tbdeaths_age_yr"
  )
}

deaths_by_age <- list(matrix(1), matrix(2))

test_that("a state's targets are its own entries and rows", {
  targets <- calib_targets_state(fake_all_states(), deaths_by_age, 2, "Alaska")
  expect_identical(targets$cases_yr, cbind(year = 1993:1995, cases = 200 + 0:2))
  expect_identical(targets$TLTBI_volume, c(2, 1.5, 2.5))
  expect_identical(targets$deaths_by_age, matrix(2))
  expect_identical(targets$cases_nat_5yr$X1996.2000, 3:4)
  expect_identical(targets$hr_cases_sm$fract_sm, c(0.3, 0.4))
  expect_identical(targets$rt_fb_cases_sm$p_recent, c(0.3, 0.4))
  expect_identical(targets$rct_cases_sm$fract_sm, 0.3)
})

test_that("a state gets the national targets and study data", {
  targets <- calib_targets_state(fake_all_states(), deaths_by_age, 1, "Alabama")
  expect_identical(targets$tx_outcomes, "national tx_outcomes")
  expect_identical(targets$TLTBI_dist, c(FB = 0.5, HR = 0.1))
  expect_identical(targets$borgdorff_data, "borgdorff")
})

test_that("a state name missing from a table is an error", {
  expect_error(
    calib_targets_state(fake_all_states(), deaths_by_age, 1, "Arizona"),
    "no rows in cases_nat_5yr, hr_cases_sm, rt_fb_cases_sm, rct_cases_sm"
  )
})

test_that("a county's targets come from its own data", {
  targets <- suppressMessages(
    calib_targets_county(fake_county(), matrix(9), fake_all_states(), "Test")
  )
  expect_identical(targets$cases_yr, cbind(year = 1993, cases = 7))
  expect_identical(targets$cases_nat_5yr$X1996.2000, 7:8)
  expect_identical(targets$rct_cases_sm$fract_sm, 0.7)
  expect_identical(targets$deaths_by_age, matrix(9))
  expect_identical(targets$tx_outcomes, "county tx_outcomes")
  expect_identical(targets$tbdeaths_age_yr, "county tbdeaths_age_yr")
})

test_that("a county falls back to national targets it lacks, and says which", {
  expect_message(
    targets <- calib_targets_county(fake_county(), matrix(9), fake_all_states(), "Test"),
    paste(
      "Calibrating Test against national targets for cases_prop_change_53_94, TLTBI_dist,",
      "LTBI_prev_US_11_IGRA, LTBI_prev_FB_11_IGRA"
    )
  )
  expect_identical(targets$cases_prop_change_53_94, c(1.1, 0.9))
  expect_identical(targets$LTBI_prev_US_11_IGRA, "national LTBI_prev_US_11_IGRA")
  expect_identical(targets$sutherland_data, "sutherland")
})

test_that("a county with every national target gets no message", {
  county <- fake_county()
  national <- fake_all_states()
  county[national_target_fields] <- national[national_target_fields]
  expect_no_message(calib_targets_county(county, matrix(9), national, "Test"))
})

test_that("county tables must hold a single location", {
  county <- fake_county()
  county$hr_cases_sm <- data.frame(st = c("County", "Other"), fract_sm = 0.7)
  expect_error(
    calib_targets_county(county, matrix(9), fake_all_states(), "Test"),
    "more than one location in hr_cases_sm"
  )
})

test_that("a county missing one of its own targets is an error", {
  county <- fake_county()
  county$homeless_pop <- NULL
  expect_error(
    suppressMessages(calib_targets_county(county, matrix(9), fake_all_states(), "Test")),
    "Calibration targets for Test are missing homeless_pop"
  )
})

test_that("recent transmission needs exactly one row", {
  county <- fake_county()
  county$rct_cases_sm <- data.frame(st = "County", fract_sm = c(0.7, 0.8))
  expect_error(
    suppressMessages(calib_targets_county(county, matrix(9), fake_all_states(), "Test")),
    "need one row in rct_cases_sm"
  )
})

test_that("calib_targets takes only states and counties", {
  expect_error(calib_targets("US"), "takes a state or county code")
  expect_error(calib_targets("Atlantis"), "takes a state or county code")
})

test_that("every state's rows are the ones its position picks out", {
  calib_dat <- suppressWarnings(readRDS(find_loc_file("ST", "CalibDat")))
  dba <- suppressWarnings(readRDS(find_loc_file("ST", "deathbyAge")))
  state_ids <- get_stateID()
  for (st in seq_len(51)) {
    targets <- calib_targets_state(calib_dat, dba, st, state_ids[st, 1])
    cn <- calib_dat$cases_nat_st_5yr
    expect_identical(targets$cases_nat_5yr, cn[cn$State.Code == st, ], info = st)
    expect_identical(targets$rct_cases_sm, calib_dat$rct_cases_sm[st, ], info = st)
    expect_identical(targets$cases_yr, calib_dat$cases_yr_st[[st]], info = st)
    expect_identical(targets$deaths_by_age, dba[[st]], info = st)
  }
})
