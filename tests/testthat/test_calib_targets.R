states <- c("Alabama", "Alaska")
buckets <- c("X1996.2000", "X2001.2005", "X2006.2010", "X2011.2015", "X2016.2020")

# One location's own targets that the calibration data stores as lists, with
# every count set to `v` so different locations' data can be told apart.
location_values <- function(v) {
  list(
    cases_yr_st = cbind(year = 1993:1995, cases = v + 0:2),
    cases_yr_ag_nat_st = array(
      v, c(2, 3, 2), dimnames = list(NULL, c("year", "total", "0-4"), c("usb", "nusb"))
    ),
    cases_yr_ag_nat_st_5yr = data.frame(usb = rep(1:0, each = 5), X0.4 = v),
    TLTBI_volume_state = c(v, v - 0.5, v + 0.5),
    tbdeaths = data.frame(State = "x", Year = 1999, Deaths = v),
    homeless_pop = rep(v, 3),
    pop_50_10 = cbind(year = 1950, usb = 1, pop = v),
    pop_00_19 = cbind(age_group = 1, usb = 1, `2000` = v)
  )
}

# One location's own targets that the calibration data stores as tables.
location_tables <- function(name, code, v) {
  nativity <- data.frame(State = name, State.Code = code, usb = c(0, 1, 3), matrix(v, 3, 5))
  colnames(nativity)[-(1:3)] <- buckets
  list(
    cases_nat_st_5yr = nativity,
    hr_cases_sm = data.frame(st = name, ss = v, fract_sm = 0.1 * 1:5),
    rt_fb_cases_sm = data.frame(st = name, p_recent = 0.1 * 1:5),
    rct_cases_sm = data.frame(st = name, fract_sm = v / 1000)
  )
}

national_values <- list(
  cases_prop_change_53_94 = c(1.1, 0.9),
  tx_outcomes = data.frame(year = 1993, pct_disc = 0.1, pct_dead = 0.05, sample_size = 100),
  TLTBI_dist = c(FB = 0.5, HR = 0.1),
  LTBI_prev_US_11_IGRA = data.frame(SSp = 1, SSn = 9),
  LTBI_prev_FB_11_IGRA = data.frame(SSp = 2, SSn = 8),
  tbdeaths_age_yr = data.frame(year = 1999, X0_4 = 1),
  borgdorff_data = data.frame(Years = 0, mean = 1),
  ferebee_data = data.frame(Year = 1, At_Risk = 10, Count = 1),
  sutherland_data = data.frame(Year = 1, At_Risk = 10, Calc_Count = 0)
)

# All-state calibration data for two states, with each state's counts set to
# 100 times its position.
fake_all_states <- function() {
  values <- lapply(seq_along(states), function(i) location_values(100 * i))
  tables <- lapply(seq_along(states), function(i) location_tables(states[i], i, 100 * i))
  per_state <- lapply(names(values[[1]]), function(field) lapply(values, `[[`, field))
  names(per_state) <- names(values[[1]])
  stacked <- lapply(names(tables[[1]]), function(field) do.call(rbind, lapply(tables, `[[`, field)))
  names(stacked) <- names(tables[[1]])
  c(per_state, stacked, national_values)
}

# One county's calibration data, with its own treatment outcomes and TB deaths
# by age but none of the other national targets.
fake_county <- function() {
  own_national <- list(
    tx_outcomes = data.frame(year = 1993, pct_disc = 0.2, pct_dead = 0.1, sample_size = 7),
    tbdeaths_age_yr = data.frame(year = 1999, X0_4 = 7)
  )
  c(location_values(7), location_tables("County", 52, 7), own_national)
}

county_targets <- function(county) {
  suppressMessages(calib_targets_county(county, matrix(9), fake_all_states(), "Test"))
}

test_that("a state's targets are its own entries and rows", {
  targets <- calib_targets_state(fake_all_states(), list(matrix(1), matrix(2)), 2, "Alaska")
  expect_identical(targets$cases_yr, cbind(year = 1993:1995, cases = 200 + 0:2))
  expect_identical(targets$TLTBI_volume, c(200, 199.5, 200.5))
  expect_identical(targets$deaths_by_age, matrix(2))
  expect_identical(unique(targets$cases_nat_5yr$State), "Alaska")
  expect_identical(targets$cases_nat_5yr$X1996.2000, c(200, 200, 200))
  expect_identical(unique(targets$hr_cases_sm$st), "Alaska")
  expect_identical(unique(targets$rt_fb_cases_sm$st), "Alaska")
  expect_identical(targets$rct_cases_sm$fract_sm, 0.2)
})

test_that("a state gets the national targets and study data", {
  targets <- calib_targets_state(fake_all_states(), list(matrix(1), matrix(2)), 1, "Alabama")
  expect_identical(targets$tx_outcomes, national_values$tx_outcomes)
  expect_identical(targets$TLTBI_dist, c(FB = 0.5, HR = 0.1))
  expect_identical(targets$borgdorff_data, national_values$borgdorff_data)
})

test_that("a state name missing from the tables is an error", {
  expect_error(
    calib_targets_state(fake_all_states(), list(matrix(1), matrix(2)), 1, "Arizona"),
    "no rows in cases_nat_5yr, hr_cases_sm, rt_fb_cases_sm, rct_cases_sm"
  )
})

test_that("a county's targets come from its own data", {
  targets <- county_targets(fake_county())
  expect_identical(targets$cases_yr, cbind(year = 1993:1995, cases = 7 + 0:2))
  expect_identical(targets$cases_nat_5yr$X1996.2000, c(7, 7, 7))
  expect_identical(targets$rct_cases_sm$fract_sm, 0.007)
  expect_identical(targets$deaths_by_age, matrix(9))
  expect_identical(targets$tx_outcomes$sample_size, 7)
  expect_identical(targets$tbdeaths_age_yr$X0_4, 7)
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
  expect_identical(targets$LTBI_prev_US_11_IGRA, national_values$LTBI_prev_US_11_IGRA)
  expect_identical(targets$sutherland_data, national_values$sutherland_data)
})

test_that("a county that provides every national target gets no message", {
  county <- fake_county()
  county[national_target_fields] <- national_values[national_target_fields]
  expect_no_message(calib_targets_county(county, matrix(9), fake_all_states(), "Test"))
})

test_that("county tables must name a single location", {
  county <- fake_county()
  county$hr_cases_sm$st[2] <- "Other"
  expect_error(county_targets(county), "more than one location in hr_cases_sm")
  county <- fake_county()
  county$hr_cases_sm$st <- NULL
  expect_error(county_targets(county), "has no st column in hr_cases_sm")
})

test_that("a county missing one of its own targets is an error", {
  county <- fake_county()
  county$homeless_pop <- NULL
  expect_error(county_targets(county), "Calibration targets for Test are missing homeless_pop")
})

test_that("recent transmission needs exactly one row", {
  county <- fake_county()
  county$rct_cases_sm <- data.frame(st = "County", fract_sm = c(0.7, 0.8))
  expect_error(county_targets(county), "need one row in rct_cases_sm")
})

test_that("a target without a column the likelihood reads is an error", {
  county <- fake_county()
  county$tbdeaths <- data.frame(year = 1999, deaths = 7)
  expect_error(county_targets(county), "have no Year, Deaths in tbdeaths")
  county <- fake_county()
  county$tx_outcomes$sample_size <- NULL
  expect_error(county_targets(county), "have no sample_size in tx_outcomes")
})

test_that("case tables need one 5-year bucket for each simulated bucket", {
  county <- fake_county()
  county$cases_nat_st_5yr$X2016.2020 <- NULL
  expect_error(county_targets(county), "need 5 5-year buckets in cases_nat_5yr")
  county <- fake_county()
  county$hr_cases_sm <- county$hr_cases_sm[1:4, ]
  expect_error(county_targets(county), "need 5 5-year buckets in hr_cases_sm")
})

test_that("annual cases without 1993 are an error", {
  county <- fake_county()
  county$cases_yr_st <- cbind(year = 1994:1995, cases = 7)
  expect_error(county_targets(county), "have no 1993 row in cases_yr")
})

test_that("calib_targets takes only states and counties", {
  expect_error(calib_targets("US"), "takes a state or county code")
  expect_error(calib_targets("Atlantis"), "takes a state or county code")
})

test_that("calib_targets reads a state's targets from the all-state files", {
  calib_dat <- suppressWarnings(readRDS(find_loc_file("ST", "CalibDat")))
  deaths <- suppressWarnings(readRDS(find_loc_file("ST", "deathbyAge")))
  expect_identical(ca_targets(), calib_targets_state(calib_dat, deaths, 5, "California"))
})

test_that("calib_targets reads a county's targets from its own files", {
  data_dir <- withr::local_tempdir()
  dir.create(file.path(data_dir, "SanDiego"))
  saveRDS(fake_county(), file.path(data_dir, "SanDiego", "SanDiego_CalibDat.rds"))
  saveRDS(matrix(9), file.path(data_dir, "SanDiego", "SanDiego_deathbyAge.rds"))
  expect_message(
    targets <- suppressWarnings(calib_targets("SanDiego", data_dir = data_dir)),
    "Calibrating SanDiego against national targets"
  )
  expect_identical(targets$cases_yr, fake_county()$cases_yr_st)
  expect_identical(targets$deaths_by_age, matrix(9))
  national <- suppressWarnings(readRDS(find_loc_file("ST", "CalibDat")))
  expect_identical(targets$borgdorff_data, national$borgdorff_data)
})

test_that("every state's rows are the ones its position picks out", {
  calib_dat <- suppressWarnings(readRDS(find_loc_file("ST", "CalibDat")))
  deaths <- suppressWarnings(readRDS(find_loc_file("ST", "deathbyAge")))
  state_ids <- get_stateID()
  # Only the 51 states; a county registered by model_load() follows them
  for (st in seq_len(51)) {
    targets <- calib_targets_state(calib_dat, deaths, st, state_ids[st, 1])
    nativity <- calib_dat$cases_nat_st_5yr
    expect_identical(targets$cases_nat_5yr, nativity[nativity$State.Code == st, ], info = st)
    expect_identical(targets$rct_cases_sm, calib_dat$rct_cases_sm[st, ], info = st)
    expect_identical(targets$cases_yr, calib_dat$cases_yr_st[[st]], info = st)
    expect_identical(targets$deaths_by_age, deaths[[st]], info = st)
  }
})
