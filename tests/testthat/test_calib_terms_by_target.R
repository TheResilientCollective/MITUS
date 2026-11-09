years <- 1950:2019

# An annual series in millions, zero except in the given years.
series_at <- function(values_by_year) {
  sim <- numeric(length(years))
  sim[year_to_idx(as.integer(names(values_by_year)))] <- values_by_year / 1e6
  sim
}

# Bucket labels as in the calibration data, which label each bucket a year later
# than the years the likelihood sums.
data_bucket_labels <- c("X1996.2000", "X2001.2005", "X2006.2010", "X2011.2015", "X2016.2020")

# Five-year case totals by nativity: non-US-born, US-born, all.
nativity_targets <- function(nusb, usb) {
  d <- data.frame(State = "Test", usb = c(0, 1, 3), rbind(nusb, usb, nusb + usb))
  colnames(d)[-(1:2)] <- data_bucket_labels
  list(cases_nat_5yr = d)
}

# An annual series whose sums over each bucket of years equal the bucket totals.
spread_over_buckets <- function(totals, buckets = case_bucket_years) {
  sim <- numeric(length(years))
  for (b in seq_along(buckets)) {
    sim[year_to_idx(buckets[[b]])] <- totals[b] / length(buckets[[b]]) / 1e6
  }
  sim
}

# Weights of 1 except at the given years.
wts_with <- function(values_by_year) {
  w <- flat_wts()
  w[names(values_by_year)] <- values_by_year
  w
}

normal_diff <- function(obs, sim, sd) {
  stats::dnorm(obs, sim, sd, log = TRUE) - stats::dnorm(obs, obs, sd, log = TRUE)
}

test_that("llik_normal is zero when the simulation matches and scales with the weights", {
  expect_identical(llik_normal(c(5, 9), c(5, 9), c(1, 2), c(1, 3)), 0)
  expect_equal(llik_normal(10, 12, 1, 2), -4)
})

test_that("sd_from_interval gives the sd of a 95% interval of +/- a share of the value", {
  expect_equal(sd_from_interval(100, 0.1), 10 / 1.96)
})

test_that("llik_dirmult is zero when the simulation matches the padded counts", {
  obs <- rbind(c(10, 0, 5), c(3, 3, 3))
  expect_equal(llik_dirmult(obs, obs + 0.1, rho = 0.01, w = c(1, 1), pad = 0.1), 0)
  expect_lt(llik_dirmult(obs, obs[, 3:1] + 0.1, rho = 0.01, w = c(1, 1), pad = 0.1), 0)
})

test_that("years stored as numbers, text or factors are read alike", {
  expect_identical(as_year(factor(c("2019", "1993"))), c(2019L, 1993L))
  expect_identical(weighted_years(factor(c("2019", "2020")), flat_wts()), c(TRUE, FALSE))
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
  targets <- list(cases_yr = cbind(year = 2016:2020, cases = c(100, 0, NA, 200, 50)))
  sim <- series_at(c(`2016` = 110, `2019` = 200))
  expect_equal(llik_cases_total(sim, targets, flat_wts()), normal_diff(100, 110, 10 / 1.96))

  off_2019 <- series_at(c(`2016` = 100, `2019` = 210))
  expected <- 2 * 3 * normal_diff(200, 210, 20 / 1.96)
  expect_equal(llik_cases_total(off_2019, targets, wts_with(c(`2019` = 3))), expected)
})

test_that("total cases read factor years as the years they name", {
  as_numbers <- list(cases_yr = data.frame(year = 2017:2019, cases = c(100, 120, 140)))
  as_factors <- as_numbers
  as_factors$cases_yr$year <- factor(as_factors$cases_yr$year)
  sim <- series_at(c(`2017` = 90, `2018` = 125, `2019` = 150))
  expect_identical(
    llik_cases_total(sim, as_factors, flat_wts()),
    llik_cases_total(sim, as_numbers, flat_wts())
  )
})

test_that("the case decline ends at the location's 1993 cases", {
  change <- rep(0.95, 41)
  targets <- list(cases_yr = cbind(year = 1993:1994, cases = c(400, 380)),
                  cases_prop_change_53_94 = change)
  obs <- cumprod(change) / prod(change) * 400
  expect_equal(obs[41], 400)
  sim <- numeric(length(years))
  sim[year_to_idx(1953:1993)] <- obs / 1e6
  expect_equal(llik_cases_decline(sim, targets, flat_wts()), 0)
  sim[year_to_idx(1993)] <- 420 / 1e6
  expect_equal(llik_cases_decline(sim, targets, flat_wts()), normal_diff(400, 420, 80 / 1.96))
})

test_that("case buckets are summed over 1995-1999 to 2015-2019, whatever their labels", {
  nusb <- c(500, 450, 400, 350, 300)
  usb <- c(200, 180, 160, 140, 120)
  targets <- nativity_targets(nusb, usb)
  expect_equal(llik_cases_nativity_5yr(spread_over_buckets(nusb), 0, targets, flat_wts()), 0)
  expect_equal(llik_cases_nativity_5yr(spread_over_buckets(usb), 1, targets, flat_wts()), 0)
  shifted <- spread_over_buckets(nusb, lapply(case_bucket_years, `+`, 1))
  expect_lt(llik_cases_nativity_5yr(shifted, 0, targets, flat_wts()), -1)
})

test_that("each nativity bucket carries the weight of 1998, 2003, 2008, 2013 or 2019", {
  nusb <- c(500, 450, 400, 350, 300)
  targets <- nativity_targets(nusb, nusb)
  bucket_wts <- wts_with(c(`1998` = 1, `2003` = 2, `2008` = 3, `2013` = 4, `2019` = 5))
  for (b in 1:5) {
    off <- spread_over_buckets(replace(nusb, b, nusb[b] + 10))
    expected <- b * normal_diff(nusb[b], nusb[b] + 10, nusb[b] * 0.05 / 1.96)
    expect_equal(llik_cases_nativity_5yr(off, 0, targets, bucket_wts), expected, info = b)
  }
})

test_that("the nativity trend compares the last two buckets of each group", {
  nusb <- c(500, 450, 400, 100, 80)
  usb <- c(200, 180, 160, 50, 45)
  targets <- nativity_targets(nusb, usb)
  sim_nusb <- spread_over_buckets(nusb)
  sim_usb <- spread_over_buckets(usb)
  expect_equal(llik_cases_nativity_trend(sim_usb, sim_nusb, targets), 0)

  flat_usb <- spread_over_buckets(c(200, 180, 160, 50, 50))
  expect_equal(
    llik_cases_nativity_trend(flat_usb, sim_nusb, targets),
    normal_diff(log(0.9), 0, 0.005)
  )
})

test_that("a fixed US-born trend replaces only the US-born target", {
  nusb <- c(500, 450, 400, 100, 80)
  usb <- c(200, 180, 160, 50, 45)
  targets <- nativity_targets(nusb, usb)
  fixed <- llik_cases_nativity_trend(
    spread_over_buckets(usb), spread_over_buckets(nusb), targets, usb_trend = -0.1
  )
  expect_equal(fixed, normal_diff(-0.1, log(0.9), 0.005))
})

# Five-year age table and matching totals for one nativity group.
age_targets <- function(counts, usb = 1) {
  ages <- c("X0.4", "X5.14", "X15.24", "X25.34", "X35.44", "X45.54", "X55.64", "X65.74",
            "X75.84", "X85p")
  shown <- data.frame(state = "Test", usb = usb, counts)
  colnames(shown)[-(1:2)] <- ages
  c(list(cases_yr_ag_nat_5yr = shown), nativity_targets(rowSums(counts), rowSums(counts)))
}

# Simulated cases by model age group, in millions, whose sums over each case
# bucket equal `by_bucket`. Each bucket's cases fall in its first year, and the
# 85+ group is split over the model's last two groups.
age_sim <- function(by_bucket) {
  sim <- matrix(0, length(years), 11)
  for (b in seq_along(case_bucket_years)) {
    first <- year_to_idx(case_bucket_years[[b]][1])
    sim[first, ] <- c(by_bucket[b, 1:9], by_bucket[b, 10] / 2, by_bucket[b, 10] / 2) / 1e6
  }
  sim
}

test_that("the case age distribution is zero when the simulation matches", {
  set.seed(1)
  counts <- matrix(sample(20:200, 50, replace = TRUE), 5, 10)
  targets <- age_targets(counts)
  expect_equal(llik_cases_age_5yr(age_sim(counts + 0.01), 1, targets, flat_wts()), 0,
               tolerance = 1e-8)
  expect_lt(llik_cases_age_5yr(age_sim(counts[, 10:1] + 0.01), 1, targets, flat_wts()), -1)
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
  matched <- llik_cases_age_5yr(age_sim(sim_counts), 1, targets, flat_wts())
  expect_equal(matched, 0, tolerance = 1e-8)

  # How the pooled cases split among the suppressed groups does not matter
  sim_counts[2, c(3, 7)] <- c(0, sum(counts[2, c(3, 7)]) + 0.01)
  expect_equal(llik_cases_age_5yr(age_sim(sim_counts), 1, targets, flat_wts()), matched)
})

test_that("the case age distribution uses rho 0.1 and each bucket's weight", {
  set.seed(4)
  counts <- matrix(sample(20:200, 50, replace = TRUE), 5, 10)
  reversed <- counts[, 10:1] + 0.01
  bucket_wts <- wts_with(c(`1998` = 1, `2003` = 2, `2008` = 3, `2013` = 4, `2019` = 5))
  expected <- sum(vapply(1:5, function(b) {
    llik_dirmult(counts[b, ], reversed[b, ], rho = 0.1, w = b, pad = 0.01)
  }, numeric(1)))
  expect_equal(
    llik_cases_age_5yr(age_sim(reversed), 1, age_targets(counts), bucket_wts), expected
  )
})

test_that("homeless cases are weighted by the year before each bucket", {
  hr <- data.frame(st = "Test", ss = c(1000, 900, 800, 700, 600), fract_sm = 0.1)
  targets <- list(hr_cases_sm = hr)
  obs <- hr$ss * 0.1
  expect_equal(llik_cases_homeless(spread_over_buckets(obs), targets, flat_wts()), 0)
  bucket_wts <- wts_with(c(`1994` = 1, `1999` = 2, `2004` = 3, `2009` = 4, `2014` = 5))
  for (b in 1:5) {
    off <- spread_over_buckets(replace(obs, b, obs[b] + 6))
    expected <- b * normal_diff(obs[b], obs[b] + 6, obs[b] * 0.1 / 1.96)
    expect_equal(llik_cases_homeless(off, targets, bucket_wts), expected, info = b)
  }
})

test_that("recent transmission is zero at the observed fraction", {
  targets <- list(rct_cases_sm = data.frame(st = "Test", fract_sm = 0.25))
  expect_identical(llik_recent_transmission(0.25, targets), 0)
  expect_lt(llik_recent_transmission(0.3, targets), 0)
})

test_that("recent-entry cases are the recent fraction of non-US-born cases, 1994-2018", {
  annual_cases <- array(
    0, c(26, 2, 2), dimnames = list(NULL, c("year", "total"), c("usb", "nusb"))
  )
  annual_cases[, "year", ] <- 1993:2018
  annual_cases[, "total", "nusb"] <- (1993:2018) - 1990
  p_recent <- c(0.5, 0.4, 0.3, 0.2, 0.1)
  targets <- list(cases_yr_ag_nat = annual_cases,
                  rt_fb_cases_sm = data.frame(st = "Test", p_recent = p_recent))
  nusb <- vapply(recent_entry_bucket_years, function(y) sum(y - 1990), numeric(1))
  obs <- p_recent * nusb
  sim <- spread_over_buckets(obs, recent_entry_bucket_years)
  expect_equal(llik_cases_recent_entry(sim, targets, flat_wts()), 0)

  bucket_wts <- wts_with(c(`1994` = 1, `1999` = 2, `2004` = 3, `2009` = 4, `2014` = 5))
  off <- spread_over_buckets(replace(obs, 3, obs[3] * 1.1), recent_entry_bucket_years)
  expected <- 3 * normal_diff(obs[3], obs[3] * 1.1, obs[3] * 0.1 / 1.96)
  expect_equal(llik_cases_recent_entry(off, targets, bucket_wts), expected)
})

test_that("treatment outcomes use only weighted years", {
  outcomes <- data.frame(year = c(2018, 2019, 2020), pct_disc = 0.1, pct_dead = 0.05,
                         sample_size = 1000)
  targets <- list(tx_outcomes = outcomes)
  obs <- c(850, 100, 50)
  sim <- matrix(0, length(years), 3)
  sim[year_to_idx(2018:2019), ] <- rep(obs + 0.01, each = 2)
  expect_equal(llik_tx_outcomes(sim, targets, flat_wts()), 0, tolerance = 1e-10)
  sim[year_to_idx(2019), ] <- c(800, 150, 50)
  expected <- llik_dirmult(rbind(obs, obs), sim[year_to_idx(2018:2019), ], rho = 0.01,
                           w = c(1, 1), pad = 0.01)
  expect_equal(llik_tx_outcomes(sim, targets, flat_wts()), expected)
  expect_lt(expected, 0)
})

test_that("LTBI treatment starts and their distribution are zero at the estimates", {
  targets <- list(TLTBI_volume = c(1000, 900, 1100), TLTBI_dist = c(FB = 0.5, HR = 0.1, HV = 0))
  expect_equal(llik_ltbi_tx_total(1000 / 1e6, targets), 0)
  expect_equal(llik_ltbi_tx_total(1100 / 1e6, targets), normal_diff(1000, 1100, 200 / 1.96))
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

test_that("LTBI prevalence uses the US-born and non-US-born IGRA accuracy", {
  ltbi <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
  pop <- 100 * ltbi
  targets <- list(LTBI_prev_US_11_IGRA = ltbi_survey, LTBI_prev_FB_11_IGRA = ltbi_survey)
  expect_equal(
    llik_ltbi_prev_usb(ltbi, pop, targets),
    llik_ltbi_prev(igra_results(ltbi, pop - ltbi, 0.780, 0.979), ltbi_survey)
  )
  expect_equal(
    llik_ltbi_prev_nusb(ltbi, pop, targets),
    llik_ltbi_prev(igra_results(ltbi, pop - ltbi, 0.789, 0.985), ltbi_survey)
  )
})

tb_deaths_targets <- function() {
  list(tbdeaths = data.frame(
    State = "Test", Year = c("2017", "2018", "2019", "2020"), Deaths = c("12", NA, "20", "5")
  ))
}

test_that("TB deaths skip missing years and years outside the weights", {
  sim <- matrix(0, length(years), 11)
  sim[year_to_idx(2017), 1] <- 12 / 1e6
  sim[year_to_idx(2019), 2] <- 20 / 1e6
  expect_equal(llik_tb_deaths_total(sim, tb_deaths_targets(), flat_wts()), 0)
  sim[year_to_idx(2019), 2] <- 22 / 1e6
  expect_equal(
    llik_tb_deaths_total(sim, tb_deaths_targets(), flat_wts()), normal_diff(20, 22, 2 / 1.96)
  )
})

test_that("TB deaths stored as a matrix are read alike", {
  sim <- matrix(0, length(years), 11)
  sim[year_to_idx(2017:2019), 1] <- c(13, 0, 19) / 1e6
  as_matrix <- list(tbdeaths = as.matrix(tb_deaths_targets()$tbdeaths))
  expect_identical(
    llik_tb_deaths_total(sim, as_matrix, flat_wts()),
    llik_tb_deaths_total(sim, tb_deaths_targets(), flat_wts())
  )
})

test_that("the TB death age distribution pools 85+ and skips years outside the weights", {
  counts <- rbind(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
                  c(9, 9, 9, 9, 9, 9, 9, 9, 9, 9))
  targets <- list(tbdeaths_age_yr = data.frame(year = 2018:2020, counts))
  sim <- matrix(0, length(years), 11)
  padded <- counts[1:2, ] + 0.005
  sim[year_to_idx(2018:2019), ] <- cbind(padded[, 1:9], padded[, 10] / 4, padded[, 10] * 3 / 4)
  expect_equal(llik_tb_deaths_age(sim, targets, flat_wts()), 0, tolerance = 1e-10)
  sim[year_to_idx(2019), 1] <- 5
  expected <- llik_dirmult(counts[1:2, ], merge_85_plus(sim[year_to_idx(2018:2019), ]),
                           rho = 0.005, w = c(1, 1), pad = 0.005)
  expect_equal(llik_tb_deaths_age(sim, targets, flat_wts()), expected)
})

test_that("the TB-free probability uses row 3 and annual rates", {
  progression <- fake_progression()
  times <- c(0, 1, 5)
  pfast <- 0.1
  rslow <- c(0.0001, 0.0002, 0.0003, 0.0004) * 12
  rate <- (0.1 + 0.002) * 12
  share <- colSums(dist_gen)
  expected <- vapply(times, function(t) {
    fast <- pfast * (1 - (1 - exp(-rate * t)) * 0.1 / 0.102)
    slow <- (1 - pfast) * (1 - (1 - exp(-(rslow + 0.002 * 12) * t)) * rslow / (rslow + 0.002 * 12))
    sum((fast + slow) * share)
  }, numeric(1))
  expect_equal(tb_free_probability(progression, times), expected)
})

test_that("study_llik replaces failed or impossible fits with -10^4", {
  expect_identical(study_llik(-3), -3)
  expect_identical(study_llik(-Inf), -1e4)
  expect_identical(study_llik(NaN), -1e4)
  expect_identical(study_llik(stop("no fit")), -1e4)
})

test_that("Borgdorff is binomial on each interval's share of cases", {
  study <- list(borgdorff_data = data.frame(Years = c(0, 1, 2, 4), mean = c(1, 0.6, 0.3, 0)))
  progression <- fake_progression()
  ss <- 854.5463 / 4
  obs <- c(0.4, 0.3, 0.3)
  p_free <- tb_free_probability(progression, c(0, 1, 2, 4))
  sim <- diff(-(1 - (1 - p_free) / (1 - p_free)[4]))
  expected <- sum(obs * ss * log(sim) + (1 - obs) * ss * log(1 - sim)) -
    sum(obs * ss * log(obs) + (1 - obs) * ss * log(1 - obs))
  expect_equal(llik_borgdorff(progression, study), expected)
})

test_that("cohort studies are binomial on each year's cases, divided by 4", {
  cohort <- data.frame(Year = 1:2, At_Risk = c(1000, 900), Count = c(20, 0))
  progression <- fake_progression()
  p_free <- tb_free_probability(progression, 0:2)
  rate <- -log(1 - diff(-p_free) / p_free[1:2])
  at_risk <- cohort$At_Risk
  fit <- c(20, 0.01)
  adj <- sum(fit * log(fit / at_risk) + (at_risk - fit) * log(1 - fit / at_risk))
  lik <- sum(cohort$Count * log(rate) + (at_risk - cohort$Count) * log(1 - rate))
  expect_equal(llik_ferebee(progression, list(ferebee_data = cohort)), lik / 4 - adj / 4)
})

test_that("the published studies give finite likelihoods", {
  targets <- ca_targets()
  expect_true(any(targets$sutherland_data[, 3] == 0))
  for (llik in list(llik_borgdorff, llik_ferebee, llik_sutherland)) {
    value <- llik(fake_progression(), targets)
    expect_true(is.finite(value) && value != -1e4)
  }
})

test_that("the Tiemersma term is zero at a 3-year duration with 90% fatality", {
  mortality <- 0.3 / 12
  self_cure <- (1 / 3 - 0.3) / 12
  expect_equal(llik_tiemersma(self_cure, mortality), 0)
  expect_lt(llik_tiemersma(self_cure * 2, mortality), 0)
})

# Census populations by decade, and annual populations by age group, in which
# age group 11 is empty as in the calibration data.
pop_targets <- function() {
  census <- rbind(cbind(year = seq(1950, 2010, 10), usb = 1, pop = seq(1e6, 7e6, 1e6)),
                  cbind(year = seq(1950, 2010, 10), usb = 0, pop = seq(1e5, 7e5, 1e5)))
  by_age <- 1000 * c(1:10, 0, 12:21, 0)
  annual_pop <- cbind(age_group = rep(1:11, 2), usb = rep(c(1, 0), each = 11),
                      matrix(by_age, 22, 20, dimnames = list(NULL, 2000:2019)))
  list(pop_50_10 = census, pop_00_19 = annual_pop)
}

test_that("population targets are the censuses after 1950 and the latest year", {
  tg <- pop_total_obs(pop_targets(), usb = 1, flat_wts())
  expect_identical(tg$years, c(seq(1960L, 2010L, 10L), 2019L))
  expect_identical(tg$pop, c(seq(2e6, 7e6, 1e6), 1000 * sum(1:10)))
  expect_identical(tg$last_census_pop, 7e6)
  short_wts <- flat_wts()[as.character(1950:2015)]
  expect_identical(pop_total_obs(pop_targets(), 1, short_wts)$years[7], 2015L)
})

test_that("the total population term is zero at the targets", {
  targets <- pop_targets()
  tg <- pop_total_obs(targets, usb = 0, flat_wts())
  sim <- numeric(length(years))
  sim[year_to_idx(tg$years)] <- tg$pop / 1e6
  expect_equal(llik_pop_total(sim, 0, targets, flat_wts()), 0)
  sim[year_to_idx(2019)] <- 1.1 * sim[year_to_idx(2019)]
  expect_equal(
    llik_pop_total(sim, 0, targets, flat_wts()),
    normal_diff(tg$pop[7], 1.1 * tg$pop[7], 7e5 * 0.05 / 1.96)
  )
})

test_that("the population age term is zero when proportions match", {
  # Age groups 1-10 hold 1-10 (US-born) and 12-21 (non-US-born) thousand
  # people; the model splits group 10 into two.
  targets <- pop_targets()
  sim_usb <- sim_nusb <- matrix(0, length(years), 11)
  sim_usb[year_to_idx(2017:2019), ] <- rep(c(1:9, 4, 6), each = 3)
  sim_nusb[year_to_idx(2017:2019), ] <- rep(c(12:20, 10, 11), each = 3)
  expect_equal(llik_pop_age(sim_usb, sim_nusb, targets, flat_wts()), 0)
})

test_that("the population age term pools the latest 3 years with sample size 500", {
  targets <- pop_targets()
  sim_usb <- sim_nusb <- matrix(0, length(years), 11)
  sim_usb[year_to_idx(2016:2019), ] <- rep(c(1:9, 4, 6), each = 4)
  sim_nusb[year_to_idx(2016:2019), ] <- rep(c(12:20, 10, 11), each = 4)
  sim_usb[year_to_idx(2019), 1] <- 4
  pooled_usb <- c(1 + 1 + 4, 2:9 * 3, 10 * 3)
  obs_usb <- (1:10) / sum(1:10)
  expected <- 500 * (sum(log(pooled_usb / sum(pooled_usb)) * obs_usb) - sum(log(obs_usb) * obs_usb))
  expect_equal(llik_pop_age(sim_usb, sim_nusb, targets, flat_wts()), expected)
  sim_usb[year_to_idx(2016), ] <- 99
  expect_equal(llik_pop_age(sim_usb, sim_nusb, targets, flat_wts()), expected)
})

test_that("the homeless population term is zero at the estimate", {
  targets <- list(homeless_pop = c(0.02, 0.01, 0.03))
  expect_identical(llik_homeless_pop(0.02, targets), 0)
  expect_equal(llik_homeless_pop(0.025, targets), normal_diff(0.02, 0.025, 0.02 / 2 / 1.96))
})

test_that("all-cause deaths use the latest weighted years", {
  by_age_counts <- matrix(1:10 * 100, 4, 10, byrow = TRUE)
  by_age_counts[4, ] <- 1
  deaths <- cbind(2017:2020, by_age_counts, rowSums(by_age_counts))
  targets <- list(deaths_by_age = deaths)
  sim <- matrix(0, length(years), 11)
  sim[year_to_idx(2018:2019), ] <- cbind(by_age_counts[2:3, 1:9], 500, 500) / 1e6
  expect_equal(llik_deaths_total(sim, targets, flat_wts()), 0)
  expect_equal(llik_deaths_age(sim, targets, flat_wts()), llik_dirmult(
    by_age_counts[2:3, ], merge_85_plus(sim[year_to_idx(2018:2019), ]), 0.01, c(1, 1), pad = 0.1
  ))
  sim[year_to_idx(2019), 1] <- 2 * sim[year_to_idx(2019), 1]
  expect_equal(llik_deaths_total(sim, targets, flat_wts()), normal_diff(5500, 5600, 550 / 1.96))
})

test_that("the mortality risk term is zero at the model's assumed distribution", {
  share <- rowSums(dist_gen) + 0.01
  sim <- matrix(0, length(years), 44)
  sim[year_to_idx(2015:2016), ] <- rep(rep(share, 11), each = 2) / 1e6
  expect_equal(llik_mort_dist(sim, flat_wts()), 0, tolerance = 1e-10)
  sim[year_to_idx(2015), 1] <- 2 * sim[year_to_idx(2015), 1]
  expect_lt(llik_mort_dist(sim, flat_wts()), 0)
})
