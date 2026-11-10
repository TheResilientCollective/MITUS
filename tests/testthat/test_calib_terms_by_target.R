bucket_starts <- c(1995L, 2000L, 2005L, 2010L, 2015L)

# Fit rows of one target.
fit_of <- function(target, observed, simulated, start = NA_integer_, end = start,
                   nativity = "all", group = "all", total = NA_real_) {
  data.frame(target = target, start = as.integer(start), end = as.integer(end),
             nativity = nativity, group = group, observed = observed, simulated = simulated,
             total = total, stringsAsFactors = FALSE)
}

# Fit rows of 5-year buckets.
bucket_fit <- function(target, observed, simulated, starts = bucket_starts, nativity = "all") {
  fit_of(target, observed, simulated, starts, starts + 4L, nativity)
}

# Weights of 1 except at the given years.
wts_with <- function(values_by_year) {
  w <- flat_wts()
  w[names(values_by_year)] <- values_by_year
  w
}

# Weights that vary within each bucket and average to the bucket's position.
bucket_mean_wts <- function(starts) {
  w <- flat_wts()
  for (b in seq_along(starts)) {
    w[as.character(starts[b] + 0:4)] <- b + c(-0.4, -0.2, 0, 0.2, 0.4)
  }
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

test_that("a period's weight is the mean of its years' weights", {
  w <- wts_with(c(`2000` = 1, `2001` = 2, `2002` = 6))
  expect_equal(period_wts(w, c(2000, 2002), c(2002, 2002)), c(3, 6))
})

test_that("sum_periods sums vectors and matrix rows over each period's years", {
  sim <- as.numeric(seq_len(70))
  expect_identical(sum_periods(sim, c(1950, 1952), c(1951, 1952)), c(3, 3))
  m <- matrix(c(sim, 2 * sim), ncol = 2)
  expect_identical(sum_periods(m, c(1950, 1952), c(1951, 1952)), rbind(c(3, 6), c(3, 6)))
})

test_that("merge_85_plus adds the last two age groups", {
  m <- matrix(1:22, 2, 11)
  merged <- merge_85_plus(m)
  expect_identical(dim(merged), c(2L, 10L))
  expect_identical(merged[, 10], m[, 10] + m[, 11])
  expect_identical(merged[, 1:9], m[, 1:9])
})

test_that("total cases skip zero, missing and unweighted years and double the latest", {
  f <- fit_of("cases_total", c(100, 0, NA, 200, 50), c(110, 5, 5, 200, 9), 2016:2020)
  expect_equal(llik_cases_total(f, flat_wts()), normal_diff(100, 110, 10 / 1.96))
  f$simulated <- c(100, 5, 5, 210, 9)
  expect_equal(llik_cases_total(f, wts_with(c(`2019` = 3))),
               2 * 3 * normal_diff(200, 210, 20 / 1.96))
})

test_that("the case decline weights each year by its importance weight", {
  f <- fit_of("cases_decline", c(400, 380), c(400, 400), 1992:1993)
  expect_equal(llik_cases_decline(f, wts_with(c(`1993` = 3))), 3 * normal_diff(380, 400, 76 / 1.96))
})

test_that("each nativity bucket carries the mean weight of its years", {
  obs <- c(500, 450, 400, 350, 300)
  for (b in 1:5) {
    f <- bucket_fit("cases_nativity", obs, replace(obs, b, obs[b] + 10))
    expected <- b * normal_diff(obs[b], obs[b] + 10, obs[b] * 0.05 / 1.96)
    expect_equal(llik_cases_nativity_5yr(f, bucket_mean_wts(bucket_starts)), expected, info = b)
  }
})

test_that("the nativity trend compares the last two buckets of each group", {
  usb <- bucket_fit("cases_nativity", c(200, 180, 160, 50, 45), c(200, 180, 160, 50, 50),
                    nativity = "usb")
  nusb <- bucket_fit("cases_nativity", c(500, 450, 400, 100, 80), c(1, 2, 3, 100, 80),
                     nativity = "nusb")
  expect_equal(llik_cases_nativity_trend(usb, nusb), normal_diff(log(0.9), 0, 0.005))
  expect_equal(llik_cases_nativity_trend(usb[5:1, ], nusb), normal_diff(log(0.9), 0, 0.005))
})

test_that("a fixed US-born trend replaces only the US-born target", {
  usb <- bucket_fit("cases_nativity", c(200, 180, 160, 50, 45), c(200, 180, 160, 50, 45),
                    nativity = "usb")
  nusb <- bucket_fit("cases_nativity", c(500, 450, 400, 100, 80), c(500, 450, 400, 100, 80),
                     nativity = "nusb")
  expect_equal(llik_cases_nativity_trend(usb, nusb, usb_trend = -0.1),
               normal_diff(-0.1, log(0.9), 0.005))
})

# Fit rows of cases by age for one nativity group, with the group's bucket totals.
age_fit <- function(shown, simulated, totals = rowSums(shown)) {
  n_ages <- length(data_age_groups)
  list(
    age = fit_of("cases_age", as.vector(t(shown)), as.vector(t(simulated)),
                 rep(bucket_starts, each = n_ages), rep(bucket_starts + 4L, each = n_ages),
                 nativity = "usb", group = rep(data_age_groups, length(bucket_starts))),
    totals = bucket_fit("cases_nativity", totals, totals, nativity = "usb")
  )
}

test_that("the case age distribution is zero when the simulation matches", {
  set.seed(1)
  counts <- matrix(sample(20:200, 50, replace = TRUE), 5, 10)
  f <- age_fit(counts, counts + 0.01)
  expect_equal(llik_cases_age_5yr(f$age, f$totals, flat_wts()), 0, tolerance = 1e-10)
  f <- age_fit(counts, counts[, 10:1] + 0.01)
  expect_lt(llik_cases_age_5yr(f$age, f$totals, flat_wts()), -1)
})

test_that("the case age distribution uses rho 0.1 and each bucket's mean weight", {
  set.seed(4)
  counts <- matrix(sample(20:200, 50, replace = TRUE), 5, 10)
  reversed <- counts[, 10:1] + 0.01
  f <- age_fit(counts, reversed)
  expected <- sum(vapply(1:5, function(b) {
    llik_dirmult(counts[b, ], reversed[b, ], rho = 0.1, w = b, pad = 0.01)
  }, numeric(1)))
  expect_equal(llik_cases_age_5yr(f$age, f$totals, bucket_mean_wts(bucket_starts)), expected)
})

test_that("suppressed age groups are pooled into the bucket's remaining cases", {
  set.seed(2)
  counts <- matrix(sample(20:200, 50, replace = TRUE), 5, 10)
  shown <- counts
  shown[2, c(3, 7)] <- NA

  # The simulation matches the shown groups and the pooled remainder
  sim <- counts + 0.01
  sim[2, c(3, 7)] <- c(sum(counts[2, c(3, 7)]) + 0.01, 0)
  f <- age_fit(shown, sim, totals = rowSums(counts))
  matched <- llik_cases_age_5yr(f$age, f$totals, flat_wts())
  expect_equal(matched, 0, tolerance = 1e-10)

  # How the pooled cases split among the suppressed groups does not matter
  sim[2, c(3, 7)] <- c(0, sum(counts[2, c(3, 7)]) + 0.01)
  f <- age_fit(shown, sim, totals = rowSums(counts))
  expect_equal(llik_cases_age_5yr(f$age, f$totals, flat_wts()), matched)
})

test_that("homeless and recent-arrival buckets carry the mean weight of their years", {
  obs <- c(100, 90, 80, 70, 60)
  starts <- bucket_starts - 2L
  for (b in 1:5) {
    sim <- replace(obs, b, obs[b] + 6)
    expected <- b * normal_diff(obs[b], obs[b] + 6, obs[b] * 0.1 / 1.96)
    w <- bucket_mean_wts(starts)
    expect_equal(llik_cases_homeless(bucket_fit("cases_homeless", obs, sim, starts), w),
                 expected, info = b)
    expect_equal(llik_cases_recent_entry(bucket_fit("cases_recent_entry", obs, sim, starts), w),
                 expected, info = b)
  }
})

test_that("recent transmission is zero at the estimated fraction", {
  expect_identical(llik_recent_transmission(fit_of("recent_transmission", 0.25, 0.25)), 0)
  expect_lt(llik_recent_transmission(fit_of("recent_transmission", 0.25, 0.3)), 0)
})

test_that("treatment outcomes use only weighted years and rho 0.01", {
  obs <- c(850, 100, 50)
  sim <- c(800, 150, 50)
  f <- fit_of("tx_outcomes", rep(obs, 3), c(obs + 0.01, sim, sim), rep(2018:2020, each = 3),
              group = tx_outcome_groups)
  expected <- llik_dirmult(rbind(obs, obs), rbind(obs + 0.01, sim), rho = 0.01, w = c(1, 1),
                           pad = 0.01)
  expect_equal(llik_tx_outcomes(f, flat_wts()), expected)
  expect_equal(llik_tx_outcomes(f[c(4:6, 1:3), ], flat_wts()), expected)
})

test_that("LTBI treatment starts and their shares are zero at the estimates", {
  targets <- list(estimates = estimate_rows("ltbi_tx_total", 2002, estimate = 1000, lower = 900,
                                            upper = 1100))
  expect_identical(llik_ltbi_tx_total(fit_of("ltbi_tx_total", 1000, 1000), targets), 0)
  expect_equal(llik_ltbi_tx_total(fit_of("ltbi_tx_total", 1000, 1100), targets),
               normal_diff(1000, 1100, 200 / 1.96))
  shares <- fit_of("ltbi_tx_share", c(0.1, 0.5), c(0.1, 0.5), group = c("homeless", "nusb"))
  expect_identical(llik_ltbi_tx_dist(shares), 0)
  shares$simulated <- c(0.5, 0.1)
  expect_lt(llik_ltbi_tx_dist(shares), 0)
})

test_that("igra_results applies sensitivity and specificity", {
  expect_equal(igra_results(10, 90, sens = 0.8, spec = 0.9), cbind(17, 83))
})

test_that("the IGRA-positive fraction pools 75+ and skips children under 5", {
  tested <- matrix(1, 11, 2)
  tested[1, ] <- c(1000, 1)
  tested[2:8, ] <- cbind(1:7, 9:3)
  tested[9:11, ] <- rbind(c(3, 7), c(2, 3), c(1, 4))
  expect_equal(igra_positive_fraction(tested), c((1:7) / 10, 6 / 20))
})

test_that("LTBI prevalence is beta on the survey counts and counts four times", {
  positive <- c(8, 12)
  tested <- c(388, 312)
  f <- fit_of("ltbi_prev", positive, positive, group = c("15_24", "25_34"), total = tested)
  expect_equal(llik_ltbi_prev(f), 0)
  f$simulated <- c(10 / 390 * 388, 12)
  expected <- 4 * (
    stats::dbeta(10 / 390, 8, 380, log = TRUE) - stats::dbeta(8 / 388, 8, 380, log = TRUE)
  )
  expect_equal(llik_ltbi_prev(f), expected)
})

test_that("TB deaths skip missing years and years outside the weights", {
  f <- fit_of("tb_deaths_total", c(12, NA, 20, 5), c(12, 99, 22, 99), 2017:2020)
  expect_equal(llik_tb_deaths_total(f, flat_wts()), normal_diff(20, 22, 2 / 1.96))
})

test_that("the TB death age distribution skips unweighted years and uses rho 0.005", {
  counts <- rbind(1:10, rep(2, 10), rep(9, 10))
  sim <- rbind(1:10 + 0.005, c(5, rep(2.005, 9)), rep(1, 10))
  f <- fit_of("tb_deaths_age", as.vector(t(counts)), as.vector(t(sim)), rep(2018:2020, each = 10),
              group = data_age_groups)
  expected <- llik_dirmult(counts[1:2, ], sim[1:2, ], rho = 0.005, w = c(1, 1), pad = 0.005)
  expect_equal(llik_tb_deaths_age(f, flat_wts()), expected)
})

test_that("the TB-free probability uses row 3 and annual rates", {
  times <- c(0, 1, 5)
  rslow <- c(0.0001, 0.0002, 0.0003, 0.0004) * 12
  recov <- 0.002 * 12
  expected <- vapply(times, function(t) {
    fast <- 0.1 * (1 - (1 - exp(-(1.2 + recov) * t)) * 1.2 / (1.2 + recov))
    slow <- 0.9 * (1 - (1 - exp(-(rslow + recov) * t)) * rslow / (rslow + recov))
    sum((fast + slow) * colSums(dist_gen))
  }, numeric(1))
  expect_equal(tb_free_probability(fake_progression(), times), expected)
})

test_that("study_llik replaces failed or impossible fits with -10^4", {
  expect_identical(study_llik(-3), -3)
  expect_identical(study_llik(-Inf), -1e4)
  expect_identical(study_llik(NaN), -1e4)
  expect_identical(study_llik(stop("no fit")), -1e4)
})

test_that("Borgdorff is binomial on each interval's share of cases", {
  obs <- c(0.4, 0.3, 0.3)
  sim <- c(0.5, 0.3, 0.2)
  ss <- 854.5463 / 4
  expected <- sum(obs * ss * log(sim) + (1 - obs) * ss * log(1 - sim)) -
    sum(obs * ss * log(obs) + (1 - obs) * ss * log(1 - obs))
  expect_equal(llik_borgdorff(fit_of("borgdorff", obs, sim)), expected)
  expect_identical(llik_borgdorff(fit_of("borgdorff", obs, c(NaN, 0.3, 0.2))), -1e4)
})

test_that("cohort studies are binomial on each year's cases, divided by 4", {
  at_risk <- c(1000, 900)
  cases <- c(20, 0)
  rate <- c(0.03, 0.001)
  fit <- c(20, 0.01)
  adj <- sum(fit * log(fit / at_risk) + (at_risk - fit) * log(1 - fit / at_risk))
  lik <- sum(cases * log(rate) + (at_risk - cases) * log(1 - rate))
  f <- fit_of("ferebee", cases, rate, group = c("1", "2"), total = at_risk)
  expect_equal(llik_cohort(f), lik / 4 - adj / 4)
})

test_that("the Tiemersma term is zero at a 3-year duration with 90% fatality", {
  f <- fit_of("tiemersma", c(3, 0.9), c(3, 0.9), group = c("duration", "fatality"))
  expect_equal(llik_tiemersma(f), 0)
  f$simulated <- c(2, 0.9)
  expect_equal(llik_tiemersma(f), normal_diff(3, 2, 0.5 / 1.96))
})

test_that("the total population compares the censuses after the first and the latest year", {
  census <- fit_of("population_census", c(1e6, 2e6, 3e6), c(9e6, 2e6, 3.1e6), c(1990, 2000, 2010))
  total <- fit_of("population_total", c(3.2e6, 3.3e6, 3.4e6), c(3.2e6, 3.4e6, 1), 2018:2020)
  sd <- 3e6 * 0.05 / 1.96
  expected <- normal_diff(3e6, 3.1e6, sd) + normal_diff(3.3e6, 3.4e6, sd)
  expect_equal(llik_pop_total(census, total, flat_wts()), expected)
})

# Fit rows of the annual population by age for one nativity group, for
# 2016-2019.
pop_age_fit <- function(observed, simulated_by_year, nativity) {
  fit_of("population", rep(observed, 4), as.vector(t(simulated_by_year)),
         rep(2016:2019, each = 10), nativity = nativity, group = rep(data_age_groups, 4))
}

test_that("the population age term pools the latest 3 years with sample size 500", {
  usb_sim <- matrix(rep(1:10, 4), 4, 10, byrow = TRUE)
  usb_sim[1, ] <- 99
  usb_sim[4, 1] <- 4
  usb <- pop_age_fit(1:10, usb_sim, "usb")
  nusb <- pop_age_fit(12:21, matrix(rep(12:21, 4), 4, 10, byrow = TRUE), "nusb")
  pooled <- c(1 + 1 + 4, 2:10 * 3)
  obs <- (1:10) / sum(1:10)
  expected <- 500 * (sum(log(pooled / sum(pooled)) * obs) - sum(log(obs) * obs))
  expect_equal(llik_pop_age(usb, nusb, flat_wts()), expected)
})

test_that("the homeless population term is zero at the estimate", {
  targets <- list(estimates = estimate_rows("homeless_pop", 2010, estimate = 0.02, lower = 0.01,
                                            upper = 0.03))
  expect_identical(llik_homeless_pop(fit_of("homeless_pop", 0.02, 0.02), targets), 0)
  expect_equal(llik_homeless_pop(fit_of("homeless_pop", 0.02, 0.025), targets),
               normal_diff(0.02, 0.025, 0.02 / 2 / 1.96))
})

test_that("all-cause deaths use the latest weighted years", {
  totals <- fit_of("deaths_total", c(1, 5500, 5500, 1), c(9, 9, 5600, 9), 2017:2020)
  expect_equal(llik_deaths_total(totals, flat_wts()), normal_diff(5500, 5600, 550 / 1.96))
  by_age <- rbind(rep(1, 10), 1:10 * 100, 1:10 * 100, rep(1, 10))
  sim <- rbind(10:1, 1:10 * 100, c(200, 2:10 * 100), rep(5, 10))
  ages <- fit_of("deaths_age", as.vector(t(by_age)), as.vector(t(sim)), rep(2017:2020, each = 10),
                 group = data_age_groups)
  expect_equal(llik_deaths_age(ages, flat_wts()),
               llik_dirmult(by_age[2:3, ], sim[2:3, ], 0.01, c(1, 1), pad = 0.1))
})

test_that("the mortality risk term is zero at the model's assumed distribution", {
  share <- rowSums(dist_gen)
  groups <- paste0(rep(model_age_groups, each = 4), "_M", 1:4)
  f <- fit_of("mort_dist", rep(rep(share, 11), each = 2), rep(rep(share + 0.01, 11), each = 2),
              start = rep(2015:2016, 44), group = rep(groups, each = 2))
  expect_equal(llik_mort_dist(f, flat_wts()), 0, tolerance = 1e-10)
  f$simulated[1] <- 2 * f$simulated[1]
  expect_lt(llik_mort_dist(f, flat_wts()), 0)
})
