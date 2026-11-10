states <- c("Alabama", "Alaska")
age_cols <- c("X0.4", "X5.14", "X15.24", "X25.34", "X35.44", "X45.54", "X55.64", "X65.74",
              "X75.84", "X85p")
nativity_labels <- c("X1995.1999", "X2000.2004", "X2005.2009", "X2010.2014", "X2015.2019")
age_labels <- c("95-99", "00-04", "05-09", "10-14", "15-19")
homeless_labels <- c("95_99", "00_04", "05_09", "10_14", "15_19")
recent_labels <- c("94_98", "99_03", "04_08", "09_13", "14_18")

# One location's own targets that the calibration data stores as lists, with
# counts set from `v` so different locations' data can be told apart.
location_values <- function(v) {
  annual <- array(0, c(26, 2, 2), dimnames = list(NULL, c("year", "total"), c("usb", "nusb")))
  annual[, "year", ] <- 1993:2018
  annual[, "total", "nusb"] <- v
  by_age <- matrix(v, 10, 10, dimnames = list(NULL, age_cols))
  list(
    cases_yr_st = cbind(year = 1993:1995, cases = v + 0:2),
    cases_yr_ag_nat_st = annual,
    cases_yr_ag_nat_st_5yr = data.frame(usb = rep(1:0, each = 5), yr_txt = age_labels, by_age),
    TLTBI_volume_state = c(v, v - 0.5, v + 0.5),
    tbdeaths = data.frame(State = "x", Year = 1999, Deaths = v),
    homeless_pop = c(v, v - 1, v + 1),
    pop_50_10 = cbind(year = c(1950, 1950, 1960, 1960), usb = c(1, 0, 1, 0), pop = v),
    pop_00_19 = cbind(age_group = rep(1:11, 2), usb = rep(c(1, 0), each = 11), `2000` = v)
  )
}

# One location's own targets that the calibration data stores as tables.
location_tables <- function(name, code, v) {
  nativity <- data.frame(State = name, State.Code = code, usb = c(0, 1, 3), matrix(v, 3, 5))
  colnames(nativity)[-(1:3)] <- nativity_labels
  list(
    cases_nat_st_5yr = nativity,
    hr_cases_sm = data.frame(st = name, yr_txt = homeless_labels, ss = v, fract_sm = 0.1 * 1:5),
    rt_fb_cases_sm = data.frame(st = name, yr_txt = recent_labels, p_recent = 0.1 * 1:5),
    rct_cases_sm = data.frame(st = name, fract_sm = v / 1000)
  )
}

national_values <- list(
  cases_prop_change_53_94 = rep(c(1.02, 0.98), length.out = 41),
  tx_outcomes = data.frame(year = 1993, pct_disc = 0.1, pct_dead = 0.05, sample_size = 100),
  TLTBI_dist = c(FB = 0.5, HR = 0.1),
  LTBI_prev_US_11_IGRA = data.frame(age_group = "5_14", SSp = 1, SSn = 9),
  LTBI_prev_FB_11_IGRA = data.frame(age_group = "5_14", SSp = 2, SSn = 8),
  tbdeaths_age_yr = data.frame(year = 1999, matrix(1, 1, 10, dimnames = list(NULL, age_cols))),
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

# All-cause deaths by age: the year, ten age groups, then the total.
fake_deaths <- function(v) cbind(2019, matrix(v, 1, 10), 10 * v)

fake_state_targets <- function(st = 2) {
  calib_targets_state(fake_all_states(), list(fake_deaths(1), fake_deaths(2)), st, states[st])
}

# One county's calibration data, with its own treatment outcomes and TB deaths
# by age but none of the other national targets.
fake_county <- function() {
  own_national <- list(
    tx_outcomes = data.frame(year = 1993, pct_disc = 0.2, pct_dead = 0.1, sample_size = 7),
    tbdeaths_age_yr = data.frame(year = 1999, matrix(7, 1, 10, dimnames = list(NULL, age_cols)))
  )
  c(location_values(7), location_tables("County", 52, 7), own_national)
}

county_targets <- function(county) {
  suppressMessages(calib_targets_county(county, fake_deaths(9), fake_all_states(), "Test"))
}

test_that("period labels name their first and last years", {
  expect_identical(
    periods_from_labels(c("X1995.1999", "95-99", "00_04", "96-00", "17_21")),
    cbind(start = c(1995L, 1995L, 2000L, 1996L, 2017L), end = c(1999L, 1999L, 2004L, 2000L, 2021L))
  )
  expect_error(periods_from_labels("X2015.1019"), "end before they start: X2015.1019")
  expect_error(periods_from_labels("X14.May"), "must name two years: X14.May")
})

test_that("a state's targets are its own values", {
  targets <- fake_state_targets()
  expect_identical(target_obs(targets, "cases_total")$value, c(200, 201, 202))
  expect_identical(target_obs(targets, "cases_nativity", "usb")$value, rep(200, 5))
  expect_identical(target_obs(targets, "deaths_total")$value, 20)
  expect_identical(target_estimate(targets, "ltbi_tx_total")$estimate, 200)
  expect_identical(target_estimate(targets, "recent_transmission")$estimate, 0.2)
})

test_that("bucketed targets take their years from the labels", {
  targets <- fake_state_targets()
  starts <- c(1995L, 2000L, 2005L, 2010L, 2015L)
  expect_identical(target_periods(target_obs(targets, "cases_nativity", "nusb"))$start, starts)
  expect_identical(target_periods(target_obs(targets, "cases_age", "usb"))$end, starts + 4L)
  expect_identical(target_periods(target_obs(targets, "cases_homeless"))$start, starts)
  expect_identical(target_periods(target_obs(targets, "cases_recent_entry"))$start, starts - 1L)
})

test_that("derived counts are computed when the targets are built", {
  targets <- fake_state_targets()
  homeless <- target_obs(targets, "cases_homeless")
  expect_identical(homeless$value, 0.1 * (1:5) * 200)
  expect_identical(homeless$total, rep(200, 5))
  recent <- target_obs(targets, "cases_recent_entry")
  expect_identical(recent$total, rep(5 * 200, 5))
  expect_identical(recent$value, 0.1 * (1:5) * 1000)
  decline <- target_obs(targets, "cases_decline")
  expect_identical(decline$start, 1953:1993)
  expect_equal(decline$value[41], 200)
  outcomes <- target_obs(targets, "tx_outcomes")
  expect_equal(setNames(outcomes$value, outcomes$group),
               c(completed = 85, discontinued = 10, died = 5))
  ltbi <- target_obs(targets, "ltbi_prev", "nusb")
  expect_identical(c(ltbi$start, ltbi$value, ltbi$total), c(2011, 2, 10))
})

test_that("all-cause deaths are read as the year, ten age groups and the total", {
  deaths <- cbind(2019, matrix(1:10, 1, 10), 55)
  targets <- calib_targets_state(fake_all_states(), list(deaths), 1, "Alabama")
  expect_identical(target_obs(targets, "deaths_age")$value, as.numeric(1:10))
  expect_identical(target_obs(targets, "deaths_age")$group, data_age_groups)
  expect_identical(target_obs(targets, "deaths_total")$value, 55)
  expect_identical(target_obs(targets, "deaths_total")$start, 2019L)
})

test_that("population rows drop the empty 95+ age group", {
  population <- target_obs(fake_state_targets(), "population", "usb")
  expect_identical(population$group, data_age_groups)
  expect_identical(nrow(target_obs(fake_state_targets(), "population_census", "nusb")), 2L)
})

test_that("a state name missing from the tables is an error", {
  expect_error(
    calib_targets_state(fake_all_states(), list(fake_deaths(1)), 1, "Arizona"),
    "no rows in cases_nat_5yr, hr_cases_sm, rt_fb_cases_sm, rct_cases_sm"
  )
})

test_that("a county's targets come from its own data", {
  targets <- county_targets(fake_county())
  expect_identical(target_obs(targets, "cases_total")$value, c(7, 8, 9))
  expect_identical(target_obs(targets, "deaths_total")$value, 90)
  expect_identical(target_obs(targets, "tx_outcomes")$total, rep(7, 3))
  expect_identical(unique(target_obs(targets, "tb_deaths_age")$value), 7)
})

test_that("a county falls back to national targets it lacks, and says which", {
  expect_message(
    targets <- calib_targets_county(fake_county(), fake_deaths(9), fake_all_states(), "Test"),
    paste(
      "Calibrating Test against national targets for cases_prop_change_53_94, TLTBI_dist,",
      "LTBI_prev_US_11_IGRA, LTBI_prev_FB_11_IGRA"
    )
  )
  expect_identical(target_obs(targets, "ltbi_prev", "usb")$value, 1)
  expect_identical(targets$studies$sutherland_data, national_values$sutherland_data)
})

test_that("a county that provides every national target gets no message", {
  county <- fake_county()
  county[national_target_fields] <- national_values[national_target_fields]
  expect_no_message(calib_targets_county(county, fake_deaths(9), fake_all_states(), "Test"))
})

test_that("county tables must name a single location", {
  county <- fake_county()
  county$hr_cases_sm$st[2] <- "Other"
  expect_error(county_targets(county), "more than one location in hr_cases_sm")
  county <- fake_county()
  county$hr_cases_sm$st <- NULL
  expect_error(county_targets(county), "has no st column in hr_cases_sm")
})

test_that("missing targets, rows and columns are errors", {
  county <- fake_county()
  county$homeless_pop <- NULL
  expect_error(county_targets(county), "Calibration targets for Test are missing homeless_pop")
  county <- fake_county()
  county$rct_cases_sm <- data.frame(st = "County", fract_sm = c(0.7, 0.8))
  expect_error(county_targets(county), "need one row in rct_cases_sm")
  county <- fake_county()
  county$tbdeaths <- data.frame(year = 1999, deaths = 7)
  expect_error(county_targets(county), "have no Year, Deaths in tbdeaths")
  county <- fake_county()
  county$cases_yr_st <- cbind(year = 1994:1995, cases = 7)
  expect_error(county_targets(county), "have no 1993 row in cases_yr")
})

test_that("buckets that overlap, leave gaps or differ in length are errors", {
  county <- fake_county()
  county$hr_cases_sm$yr_txt[4] <- "05_09"
  expect_error(county_targets(county), "cases_homeless periods that overlap, leave gaps")
  county <- fake_county()
  county$rt_fb_cases_sm$yr_txt[2] <- "99_13"
  expect_error(county_targets(county), "cases_recent_entry periods")
  county <- fake_county()
  colnames(county$cases_nat_st_5yr)[8] <- "X2016.2020"
  expect_error(county_targets(county), "cases_nativity periods")
  county <- fake_county()
  county$hr_cases_sm$yr_txt[5] <- "15_20"
  expect_error(county_targets(county), "cases_homeless periods that overlap, leave gaps")
})

test_that("age buckets must match the nativity buckets", {
  county <- fake_county()
  county$cases_yr_ag_nat_st_5yr$yr_txt <- rep(c("96-00", "01-05", "06-10", "11-15", "16-20"), 2)
  expect_error(county_targets(county), "cases_age periods that differ from cases_nativity")
})

test_that("recent-arrival buckets need annual non-US-born cases for every year", {
  county <- fake_county()
  county$cases_yr_ag_nat_st <- county$cases_yr_ag_nat_st[1:25, , ]
  expect_error(county_targets(county), "lack annual non-US-born cases for 2018")
})

test_that("a bucket ending after the last weighted year is an error", {
  targets <- fake_state_targets()
  expect_silent(check_target_horizon(targets, flat_wts()))
  expect_error(check_target_horizon(targets, flat_wts()[as.character(1950:2018)]),
               "run past the last weighted year, 2018, in cases_nativity, cases_age")
})

test_that("single years after the last weighted year are allowed", {
  targets <- fake_state_targets()
  targets$obs <- rbind(targets$obs, obs_rows("cases_total", 2020, value = 300))
  expect_silent(check_target_horizon(targets, flat_wts()))
})

test_that("calib_targets takes only states and counties", {
  expect_error(calib_targets("US"), "takes a state or county code")
  expect_error(calib_targets("Atlantis"), "takes a state or county code")
})

test_that("calib_targets reads the calibration data file it is given", {
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(state_calib_data(), path)
  expect_identical(suppressWarnings(calib_targets("CA", calib_data_file = path)), ca_targets())
})

test_that("calib_targets builds a state's targets from the default data", {
  expect_identical(suppressWarnings(calib_targets("CA")), ca_targets())
})

test_that("the 04-20-22 state calibration data has buckets that are not contiguous", {
  path <- system.file("ST", "ST_CalibDat_04-20-22.rds", package = "MITUS")
  expect_error(suppressWarnings(calib_targets("CA", calib_data_file = path)),
               "cases_homeless periods that overlap")
})

test_that("calib_targets reads a county's targets from its own files", {
  data_dir <- withr::local_tempdir()
  dir.create(file.path(data_dir, "SanDiego"))
  saveRDS(fake_county(), file.path(data_dir, "SanDiego", "SanDiego_CalibDat.rds"))
  saveRDS(fake_deaths(9), file.path(data_dir, "SanDiego", "SanDiego_deathbyAge.rds"))
  expect_message(
    targets <- suppressWarnings(calib_targets("SanDiego", data_dir = data_dir)),
    "Calibrating SanDiego against national targets"
  )
  expect_identical(target_obs(targets, "cases_total")$value, c(7, 8, 9))
  national <- suppressWarnings(readRDS(find_loc_file("ST", "CalibDat")))
  expect_identical(targets$studies$borgdorff_data, national$borgdorff_data)
})

test_that("every state's targets are the values its position picks out", {
  calib_dat <- state_calib_data()
  deaths <- state_deaths_by_age()
  state_ids <- get_stateID()
  # Only the 51 states; a county registered by model_load() follows them
  for (st in seq_len(51)) {
    targets <- calib_targets_state(calib_dat, deaths, st, state_ids[st, 1])
    nativity <- calib_dat$cases_nat_st_5yr
    usb <- nativity[nativity$State.Code == st & nativity$usb == 1, -(1:3)]
    expect_identical(target_obs(targets, "cases_nativity", "usb")$value,
                     as.numeric(unlist(usb)), info = st)
    expect_identical(target_estimate(targets, "recent_transmission")$estimate,
                     as.numeric(calib_dat$rct_cases_sm[st, "fract_sm"]), info = st)
    expect_identical(target_obs(targets, "cases_total")$value,
                     as.numeric(calib_dat$cases_yr_st[[st]][, "cases"]), info = st)
    expect_equal(target_obs(targets, "deaths_total")$value, deaths[[st]][, 12], info = st)
  }
})
