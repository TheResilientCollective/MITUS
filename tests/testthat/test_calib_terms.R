years <- 1950:2019
flat_wts <- setNames(rep(1, length(years)), years)
row_of <- function(year) year_to_idx(year)

# An annual series in millions, zero except in the given years.
annual <- function(values_by_year) {
  sim <- numeric(length(years))
  sim[row_of(as.integer(names(values_by_year)))] <- values_by_year / 1e6
  sim
}

bucket_cols <- c("X1995.1999", "X2000.2004", "X2005.2009", "X2010.2014", "X2015.2019")

# Five-year case totals by nativity: non-US-born, US-born, all.
nat_5yr <- function(nusb, usb) {
  d <- data.frame(State = "Test", usb = c(0, 1, 3), rbind(nusb, usb, nusb + usb))
  colnames(d)[3:7] <- bucket_cols
  d
}

# An annual series whose sums over each bucket equal the bucket totals.
spread_over_buckets <- function(totals, buckets = case_bucket_years) {
  sim <- numeric(length(years))
  for (b in seq_along(buckets)) {
    sim[row_of(buckets[[b]])] <- totals[b] / length(buckets[[b]]) / 1e6
  }
  sim
}

test_that("llik_normal is zero at a perfect fit and scales with the weights", {
  expect_identical(llik_normal(c(5, 9), c(5, 9), c(1, 2), c(1, 3)), 0)
  expect_equal(llik_normal(10, 12, 1, 2), -4)
})

test_that("llik_dirmult is zero when the simulation matches the padded counts", {
  obs <- rbind(c(10, 0, 5), c(3, 3, 3))
  expect_equal(llik_dirmult(obs, obs + 0.1, rho = 0.01, w = c(1, 1), pad = 0.1), 0)
  expect_lt(llik_dirmult(obs, obs[, 3:1] + 0.1, rho = 0.01, w = c(1, 1), pad = 0.1), 0)
})

test_that("sum_buckets sums vectors and matrix rows by year", {
  sim <- as.numeric(seq_along(years))
  expect_identical(sum_buckets(sim, list(1950:1951, 1952)), c(3, 3))
  m <- matrix(c(sim, 2 * sim), ncol = 2)
  expect_identical(sum_buckets(m, list(1950:1951, 1952)), rbind(c(3, 6), c(3, 6)))
})

test_that("merge_85_plus adds the last two age groups", {
  m <- matrix(1:22, 2, 11)
  merged <- merge_85_plus(m)
  expect_identical(dim(merged), c(2L, 10L))
  expect_identical(merged[, 10], m[, 10] + m[, 11])
  expect_identical(merged[, 1:9], m[, 1:9])
})

test_that("total cases skip zero and missing years and double the latest", {
  targets <- list(cases_yr = cbind(
    year = 2016:2020, cases = c(100, 0, NA, 200, 50)
  ))
  sim <- annual(c(`2016` = 110, `2019` = 200))
  sd_2016 <- 100 * 0.1 / 1.96
  expected <- stats::dnorm(100, 110, sd_2016, log = TRUE) -
    stats::dnorm(100, 100, sd_2016, log = TRUE)
  expect_equal(llik_cases_total(sim, targets, flat_wts), expected)

  weighted <- flat_wts
  weighted["2019"] <- 3
  off_2019 <- annual(c(`2016` = 100, `2019` = 210))
  sd_2019 <- 200 * 0.1 / 1.96
  expected_2019 <- 2 * 3 * (
    stats::dnorm(200, 210, sd_2019, log = TRUE) - stats::dnorm(200, 200, sd_2019, log = TRUE)
  )
  expect_equal(llik_cases_total(off_2019, targets, weighted), expected_2019)
})

test_that("the case decline ends at the location's 1993 cases", {
  change <- rep(0.95, 41)
  targets <- list(cases_yr = cbind(year = 1993:1994, cases = c(400, 380)),
                  cases_prop_change_53_94 = change)
  obs <- cumprod(change) / prod(change) * 400
  expect_equal(obs[41], 400)
  sim <- numeric(length(years))
  sim[row_of(1953:1993)] <- obs / 1e6
  expect_equal(llik_cases_decline(sim, targets, flat_wts), 0)
  sim[row_of(1993)] <- 420 / 1e6
  expect_lt(llik_cases_decline(sim, targets, flat_wts), 0)
})

test_that("cases by nativity compare bucket sums, weighted by each bucket's year", {
  nusb <- c(500, 450, 400, 350, 300)
  usb <- c(200, 180, 160, 140, 120)
  targets <- list(cases_nat_5yr = nat_5yr(nusb, usb))
  expect_equal(llik_cases_nat_5yr(spread_over_buckets(nusb), 0, targets, flat_wts), 0)
  expect_equal(llik_cases_nat_5yr(spread_over_buckets(usb), 1, targets, flat_wts), 0)

  weighted <- flat_wts
  weighted[as.character(case_bucket_wt_years)] <- 1:5
  off <- spread_over_buckets(replace(nusb, 2, 460))
  sd <- 450 * 0.05 / 1.96
  expected <- 2 * (stats::dnorm(450, 460, sd, log = TRUE) - stats::dnorm(450, 450, sd, log = TRUE))
  expect_equal(llik_cases_nat_5yr(off, 0, targets, weighted), expected)
})

test_that("the nativity trend compares the last two buckets", {
  nusb <- c(500, 450, 400, 100, 80)
  usb <- c(200, 180, 160, 50, 40)
  targets <- list(cases_nat_5yr = nat_5yr(nusb, usb))
  sim_nusb <- spread_over_buckets(nusb)
  sim_usb <- spread_over_buckets(usb)
  expect_equal(llik_cases_nat_trend(sim_nusb, sim_usb, targets), 0)

  fixed <- llik_cases_nat_trend(sim_nusb, sim_usb, targets, usb_trend = -0.1)
  expected <- stats::dnorm(-0.1, log(0.8), 0.005, log = TRUE) -
    stats::dnorm(-0.1, -0.1, 0.005, log = TRUE)
  expect_equal(fixed, expected)
})

# Five-year age table and matching totals for one nativity group.
age_targets <- function(counts, usb = 1) {
  ages <- c("X0.4", "X5.14", "X15.24", "X25.34", "X35.44", "X45.54", "X55.64", "X65.74",
            "X75.84", "X85p")
  shown <- data.frame(state = "Test", usb = usb, counts)
  colnames(shown)[-(1:2)] <- ages
  list(cases_yr_ag_nat_5yr = shown,
       cases_nat_5yr = nat_5yr(rowSums(counts), rowSums(counts)))
}

# Simulated cases by model age group whose bucket sums are proportional to
# `by_bucket`, with the 85+ group split over the model's last two groups.
age_sim <- function(by_bucket) {
  sim <- matrix(0, length(years), 11)
  for (b in seq_along(case_bucket_years)) {
    first <- row_of(case_bucket_years[[b]][1])
    sim[first, ] <- c(by_bucket[b, 1:9], by_bucket[b, 10] / 2, by_bucket[b, 10] / 2) / 1e6
  }
  sim
}

test_that("the case age distribution is zero when the simulation matches", {
  set.seed(1)
  counts <- matrix(sample(20:200, 50, replace = TRUE), 5, 10)
  targets <- age_targets(counts)
  expect_equal(llik_cases_age_5yr(age_sim(counts + 0.01), 1, targets, flat_wts), 0,
               tolerance = 1e-8)
  expect_lt(llik_cases_age_5yr(age_sim(counts[, 10:1] + 0.01), 1, targets, flat_wts), -1)
})

test_that("suppressed age groups are pooled into the bucket's remaining cases", {
  set.seed(2)
  counts <- matrix(sample(20:200, 50, replace = TRUE), 5, 10)
  targets <- age_targets(counts)
  shown <- counts
  shown[2, c(3, 7)] <- NA
  targets$cases_yr_ag_nat_5yr[, -(1:2)] <- shown

  # The simulation matches the shown groups and the pooled remainder
  sim_counts <- counts + 0.01
  sim_counts[2, c(3, 7)] <- c(sum(counts[2, c(3, 7)]) + 0.01, 0)
  matched <- llik_cases_age_5yr(age_sim(sim_counts), 1, targets, flat_wts)
  expect_equal(matched, 0, tolerance = 1e-8)

  # How the pooled cases split among the suppressed groups does not matter
  sim_counts[2, c(3, 7)] <- c(0, sum(counts[2, c(3, 7)]) + 0.01)
  expect_equal(llik_cases_age_5yr(age_sim(sim_counts), 1, targets, flat_wts), matched)
})

test_that("homeless cases compare bucket sums with the smoothed counts", {
  hr <- data.frame(st = "Test", ss = c(1000, 900, 800, 700, 600), fract_sm = 0.1)
  targets <- list(hr_cases_sm = hr)
  expect_equal(llik_cases_homeless(spread_over_buckets(hr$ss * 0.1), targets, flat_wts), 0)
  weighted <- flat_wts
  weighted[as.character(smoothed_bucket_wt_years)] <- 1:5
  off <- spread_over_buckets(c(100, 90, 80, 70, 66))
  sd <- 60 * 0.1 / 1.96
  expected <- 5 * (stats::dnorm(60, 66, sd, log = TRUE) - stats::dnorm(60, 60, sd, log = TRUE))
  expect_equal(llik_cases_homeless(off, targets, weighted), expected)
})

test_that("recent transmission is zero at the observed fraction", {
  targets <- list(rct_cases_sm = data.frame(st = "Test", fract_sm = 0.25))
  expect_identical(llik_recent_transmission(0.25, targets), 0)
  expect_lt(llik_recent_transmission(0.3, targets), 0)
})

test_that("recent-entry cases are the recent fraction of non-US-born cases", {
  annual_cases <- array(
    0, c(26, 2, 2), dimnames = list(NULL, c("year", "total"), c("usb", "nusb"))
  )
  annual_cases[, "year", ] <- 1993:2018
  annual_cases[, "total", "nusb"] <- 10
  targets <- list(cases_yr_ag_nat = annual_cases,
                  rt_fb_cases_sm = data.frame(st = "Test", p_recent = c(0.5, 0.4, 0.3, 0.2, 0.1)))
  obs <- c(0.5, 0.4, 0.3, 0.2, 0.1) * 50
  sim <- spread_over_buckets(obs, recent_entry_bucket_years)
  expect_equal(llik_cases_recent_entry(sim, targets, flat_wts), 0)
  expect_lt(llik_cases_recent_entry(sim * 1.2, targets, flat_wts), 0)
})

test_that("treatment outcomes use only weighted years", {
  outcomes <- data.frame(year = c(2018, 2019, 2020), pct_disc = 0.1, pct_dead = 0.05,
                         sample_size = 1000)
  targets <- list(tx_outcomes = outcomes)
  obs <- c(850, 100, 50)
  sim <- matrix(0, length(years), 3)
  sim[row_of(2018:2019), ] <- rep(obs + 0.01, each = 2)
  expect_equal(llik_tx_outcomes(sim, targets, flat_wts), 0, tolerance = 1e-10)
  sim[row_of(2019), ] <- c(800, 150, 50)
  expect_lt(llik_tx_outcomes(sim, targets, flat_wts), 0)
})

test_that("LTBI treatment starts and their distribution are zero at the estimates", {
  targets <- list(TLTBI_volume = c(1000, 900, 1100), TLTBI_dist = c(FB = 0.5, HR = 0.1, HV = 0))
  expect_equal(llik_ltbi_tx_total(1000 / 1e6, targets), 0)
  sd <- 200 / 1.96
  expect_equal(
    llik_ltbi_tx_total(1100 / 1e6, targets),
    stats::dnorm(1000, 1100, sd, log = TRUE) - stats::dnorm(1000, 1000, sd, log = TRUE)
  )
  expect_identical(llik_ltbi_tx_dist(c(0.5, 0.1), targets), 0)
  expect_lt(llik_ltbi_tx_dist(c(0.4, 0.1), targets), 0)
})

test_that("igra_results applies sensitivity and specificity", {
  expect_equal(igra_results(10, 90, sens = 0.8, spec = 0.9), cbind(17, 83))
})

# Survey counts and expected test results that match them in ages 5 to 75+.
ltbi_survey <- data.frame(age_group = 1:8, SSp = c(4, 8, 12, 20, 30, 40, 50, 60),
                          SSn = c(500, 380, 300, 250, 200, 150, 120, 100))
matching_tests <- function() {
  tested <- matrix(1, 11, 2)
  tested[2:8, ] <- as.matrix(ltbi_survey[1:7, c("SSp", "SSn")])
  tested[9:11, ] <- rbind(c(30, 50), c(20, 30), c(10, 20))
  tested
}

test_that("LTBI prevalence is zero when the tested fractions match the survey", {
  expect_equal(llik_ltbi_prev(matching_tests(), ltbi_survey), 0)
  expect_equal(llik_ltbi_prev(matching_tests() * 3, ltbi_survey), 0)
})

test_that("LTBI prevalence counts four times and ignores children under 5", {
  tested <- matching_tests()
  tested[1, ] <- c(1000, 1)
  expect_equal(llik_ltbi_prev(tested, ltbi_survey), 0)
  tested[3, ] <- c(10, 380)
  expected <- 4 * (
    stats::dbeta(10 / 390, 8, 380, log = TRUE) - stats::dbeta(8 / 388, 8, 380, log = TRUE)
  )
  expect_equal(llik_ltbi_prev(tested, ltbi_survey), expected)
})

test_that("LTBI prevalence applies each group's IGRA accuracy", {
  ltbi <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
  pop <- 100 * ltbi
  targets <- list(LTBI_prev_US_11_IGRA = ltbi_survey, LTBI_prev_FB_11_IGRA = ltbi_survey)
  expect_equal(
    llik_ltbi_prev_usb(ltbi, pop, targets),
    llik_ltbi_prev(igra_results(ltbi, pop - ltbi, 0.780, 0.979), ltbi_survey)
  )
  tested_nusb <- igra_results(ltbi, pop - ltbi, 0.789, 0.985)
  tested_nusb[1, ] <- igra_results(ltbi[1], pop[1] - ltbi[1], 0.712, 0.989)
  expect_equal(llik_ltbi_prev_nusb(ltbi, pop, targets), llik_ltbi_prev(tested_nusb, ltbi_survey))
})

test_that("TB deaths skip missing years and years outside the weights", {
  targets <- list(tbdeaths = data.frame(
    State = "Test", Year = c("2017", "2018", "2019", "2020"), Deaths = c("12", NA, "20", "5")
  ))
  sim <- matrix(0, length(years), 11)
  sim[row_of(2017), 1] <- 12 / 1e6
  sim[row_of(2019), 2] <- 20 / 1e6
  expect_equal(llik_tb_deaths_total(sim, targets, flat_wts), 0)
  sim[row_of(2019), 2] <- 22 / 1e6
  sd <- 20 * 0.1 / 1.96
  expected <- stats::dnorm(20, 22, sd, log = TRUE) - stats::dnorm(20, 20, sd, log = TRUE)
  expect_equal(llik_tb_deaths_total(sim, targets, flat_wts), expected)
})

test_that("the TB death age distribution pools the model's 85+ groups", {
  counts <- rbind(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2))
  targets <- list(tbdeaths_age_yr = data.frame(year = 2018:2019, counts))
  sim <- matrix(0, length(years), 11)
  padded <- counts + 0.005
  sim[row_of(2018:2019), ] <- cbind(padded[, 1:9], padded[, 10] / 4, padded[, 10] * 3 / 4)
  expect_equal(llik_tb_deaths_age(sim, targets, flat_wts), 0, tolerance = 1e-10)
})

progression <- list(
  Mpfast = matrix(0.1, 4, 4),
  Mrslow = matrix(c(0.0001, 0.0002, 0.0003, 0.0004), 4, 4, byrow = TRUE),
  rfast = 0.1,
  rRecov = 0.002
)

test_that("everyone is TB-free at infection, and fewer are later", {
  p <- tb_free_probability(progression, 0:10)
  expect_equal(p[1], sum(dist_gen))
  expect_true(all(diff(p) < 0))
})

test_that("study_llik replaces failed or impossible fits with -10^4", {
  expect_identical(study_llik(-3), -3)
  expect_identical(study_llik(-Inf), -1e4)
  expect_identical(study_llik(NaN), -1e4)
  expect_identical(study_llik(stop("no fit")), -1e4)
})

test_that("the published studies give finite likelihoods", {
  calib_dat <- suppressWarnings(readRDS(find_loc_file("ST", "CalibDat")))
  targets <- calib_dat[shared_target_fields]
  expect_true(any(targets$sutherland_data[, 3] == 0))
  for (llik in list(llik_borgdorff, llik_ferebee, llik_sutherland)) {
    value <- llik(progression, targets)
    expect_true(is.finite(value) && value != -1e4)
  }
})

test_that("Tiemersma is zero at a 3-year duration with 90% fatality", {
  mortality <- 0.3 / 12
  self_cure <- (1 / 3 - 0.3) / 12
  expect_equal(llik_tiemersma(self_cure, mortality), 0)
  expect_lt(llik_tiemersma(self_cure * 2, mortality), 0)
})

# Census populations by decade, and annual populations by age group.
pop_targets <- function() {
  census <- rbind(cbind(year = seq(1950, 2010, 10), usb = 1, pop = seq(1e6, 7e6, 1e6)),
                  cbind(year = seq(1950, 2010, 10), usb = 0, pop = seq(1e5, 7e5, 1e5)))
  annual_pop <- cbind(age_group = rep(1:11, 2), usb = rep(c(1, 0), each = 11),
                      matrix(1000 * rep(1:22, 20), 22, 20, dimnames = list(NULL, 2000:2019)))
  list(pop_50_10 = census, pop_00_19 = annual_pop)
}

test_that("population targets are the censuses after 1950 and the latest year", {
  tg <- pop_total_obs(pop_targets(), usb = 1, flat_wts)
  expect_identical(tg$years, c(seq(1960L, 2010L, 10L), 2019L))
  expect_identical(tg$pop, c(seq(2e6, 7e6, 1e6), 1000 * sum(1:10)))
  expect_identical(tg$scale, 7e6)
  short_wts <- flat_wts[as.character(1950:2015)]
  expect_identical(pop_total_obs(pop_targets(), 1, short_wts)$years[7], 2015L)
})

test_that("total population is zero at the targets", {
  targets <- pop_targets()
  tg <- pop_total_obs(targets, usb = 0, flat_wts)
  sim <- numeric(length(years))
  sim[row_of(tg$years)] <- tg$pop / 1e6
  expect_equal(llik_pop_total(sim, 0, targets, flat_wts), 0)
  sim[row_of(2019)] <- 1.1 * sim[row_of(2019)]
  expect_lt(llik_pop_total(sim, 0, targets, flat_wts), 0)
})

test_that("the population age distribution is zero when proportions match", {
  # The annual table's age groups 1-10 hold 1-10 (US-born) and 12-21
  # (non-US-born) thousand people; the model splits group 10 into two.
  targets <- pop_targets()
  sim_usb <- sim_nusb <- matrix(0, length(years), 11)
  sim_usb[row_of(2017:2019), ] <- rep(c(1:9, 4, 6), each = 3)
  sim_nusb[row_of(2017:2019), ] <- rep(c(12:20, 10, 11), each = 3)
  expect_equal(llik_pop_age(sim_usb, sim_nusb, targets, flat_wts), 0)
  sim_usb[row_of(2019), 1] <- 5
  expect_lt(llik_pop_age(sim_usb, sim_nusb, targets, flat_wts), 0)
})

test_that("the homeless population is zero at the estimate", {
  targets <- list(homeless_pop = c(0.02, 0.01, 0.03))
  expect_identical(llik_homeless_pop(0.02, targets), 0)
  sd <- 0.02 / 2 / 1.96
  expect_equal(
    llik_homeless_pop(0.025, targets),
    stats::dnorm(0.02, 0.025, sd, log = TRUE) - stats::dnorm(0.02, 0.02, sd, log = TRUE)
  )
})

test_that("all-cause deaths use the latest weighted years", {
  by_age_counts <- matrix(1:10 * 100, 4, 10, byrow = TRUE)
  dba <- cbind(2017:2020, by_age_counts, rowSums(by_age_counts))
  targets <- list(deaths_by_age = dba)
  sim <- matrix(0, length(years), 11)
  sim[row_of(2018:2019), ] <- cbind(by_age_counts[1:2, 1:9], 500, 500) / 1e6
  expect_equal(llik_deaths_total(sim, targets, flat_wts), 0)
  expect_equal(llik_deaths_age(sim, targets, flat_wts), llik_dirmult(
    by_age_counts[1:2, ], merge_85_plus(sim[row_of(2018:2019), ]), 0.01, c(1, 1), pad = 0.1
  ))
  sim[row_of(2019), 1] <- 2 * sim[row_of(2019), 1]
  expect_lt(llik_deaths_total(sim, targets, flat_wts), 0)
})

test_that("the mortality risk distribution is zero at the model's assumption", {
  share <- rowSums(dist_gen) + 0.01
  sim <- matrix(0, length(years), 44)
  sim[row_of(2015:2016), ] <- rep(rep(share, 11), each = 2) / 1e6
  expect_equal(llik_mort_dist(sim, flat_wts), 0, tolerance = 1e-10)
  sim[row_of(2015), 1] <- 2 * sim[row_of(2015), 1]
  expect_lt(llik_mort_dist(sim, flat_wts), 0)
})
