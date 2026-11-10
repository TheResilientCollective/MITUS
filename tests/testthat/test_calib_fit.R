rows_of <- function(fit, target, nativity = NULL) {
  f <- fit[fit$target == target, , drop = FALSE]
  if (!is.null(nativity)) f <- f[f$nativity == nativity, , drop = FALSE]
  f
}

ca_fit <- function(out = fake_outputs(), targets = ca_targets()) {
  calib_fit(out, fake_progression(), c(rSlfCur = 0.01, muIp = 0.02), targets)
}

test_that("fit rows past the simulated years are left out", {
  d <- obs_rows("cases_total", 2018:2021, value = 1:4)
  expect_identical(simulated_rows(d, 70)$start, 2018:2019)
})

test_that("period_fit sums the simulated series over each row's years", {
  targets <- list(obs = obs_rows("cases_homeless", c(2010, 2015), c(2014, 2019), c(5, 6)))
  sim <- seq_len(70) / 1e6
  f <- period_fit(targets, "cases_homeless", sim)
  expect_equal(f$simulated, c(sum(61:65), sum(66:70)))
  expect_identical(f$observed, c(5, 6))
})

test_that("group_fit sums each row's group over its years", {
  targets <- list(obs = obs_rows("deaths_age", c(2018, 2019), value = 1:2, group = c("b", "a")))
  sim <- cbind(a = seq_len(70), b = 100 * seq_len(70))
  expect_identical(group_fit(targets, "deaths_age", sim)$simulated, c(6900, 70))
})

test_that("by_data_age pools the model's 85+ groups and names the data's groups", {
  m <- matrix(1, 2, 11)
  merged <- by_data_age(m)
  expect_identical(colnames(merged), data_age_groups)
  expect_identical(unname(merged[1, "85p"]), 2)
})

test_that("the fit has rows for every target the likelihood reads", {
  fit <- ca_fit()
  expect_setequal(unique(fit$target), c(
    "cases_total", "cases_decline", "cases_nativity", "cases_age", "cases_homeless",
    "cases_recent_entry", "recent_transmission", "tx_outcomes", "ltbi_tx_total",
    "ltbi_tx_share", "ltbi_prev", "tb_deaths_total", "tb_deaths_age", "borgdorff", "ferebee",
    "sutherland", "tiemersma", "population_census", "population_total", "population",
    "homeless_pop", "deaths_total", "deaths_age", "mort_dist"
  ))
  expect_true(all(is.finite(fit$simulated)))
})

test_that("case targets get the simulated cases they compare", {
  out <- fake_outputs()
  fit <- ca_fit(out)
  usb <- rows_of(fit, "cases_nativity", "usb")
  usb_cases <- out[, "NOTIF_US"] + out[, "NOTIF_MORT_US"]
  expect_equal(usb$simulated, sum_periods(usb_cases, usb$start, usb$end) * 1e6)
  age <- rows_of(fit, "cases_age", "nusb")
  nusb_age <- by_age(out, "NOTIF_") + by_age(out, "NOTIF_MORT_") -
    (by_age(out, "NOTIF_US_") + by_age(out, "NOTIF_US_MORT_"))
  first <- age[1, ]
  expect_equal(first$simulated,
               sum(by_data_age(nusb_age)[year_to_idx(first$start:first$end), first$group]) * 1e6)
})

test_that("single-year and fixed-window targets take their years from the targets", {
  out <- fake_outputs()
  targets <- ca_targets()
  est <- targets$estimates
  est[est$target == "recent_transmission", c("start", "end")] <- c(2012L, 2016L)
  est[est$target %in% c("ltbi_tx_total", "ltbi_tx_share"), c("start", "end")] <- 2005L
  est[est$target == "homeless_pop", c("start", "end")] <- 2013L
  targets$estimates <- est
  targets$obs[targets$obs$target == "ltbi_prev", c("start", "end")] <- 2014L
  fit <- ca_fit(out, targets)

  window <- year_to_idx(2012:2016)
  recent <- sum(out[window, c("INCID_REC_US", "INCID_REC_FB")]) /
    sum(out[window, c("INCID_ALL_US", "INCID_ALL_FB")])
  expect_identical(rows_of(fit, "recent_transmission")$simulated, recent)
  expect_identical(rows_of(fit, "ltbi_tx_total")$simulated,
                   unname(out[year_to_idx(2005), "TLTBI_INITS"] * 1e6))
  expect_identical(rows_of(fit, "homeless_pop")$simulated, unname(out[year_to_idx(2013), "N_HR"]))
  y2014 <- year_to_idx(2014)
  ltbi <- by_age(out, "N_US_LTBI_")[y2014, ]
  tested <- igra_results(ltbi, by_age(out, "N_US_")[y2014, ] - ltbi, 0.780, 0.979)
  survey <- rows_of(fit, "ltbi_prev", "usb")
  expect_identical(survey$simulated, unname(igra_positive_fraction(tested) * survey$total))
})

test_that("non-US-born children under 5 take their own IGRA accuracy", {
  out <- fake_outputs()
  y2011 <- year_to_idx(2011)
  ltbi <- by_age(out, "N_FB_LTBI_")[y2011, ]
  pop <- by_age(out, "N_FB_")[y2011, ]
  tested <- igra_results(ltbi, pop - ltbi, 0.789, 0.985)
  tested[1, ] <- igra_results(ltbi[1], pop[1] - ltbi[1], 0.712, 0.989)
  survey <- rows_of(ca_fit(out), "ltbi_prev", "nusb")
  expect_equal(survey$simulated, unname(igra_positive_fraction(tested) * survey$total))
})

test_that("the population totals sum the annual population over age groups", {
  out <- fake_outputs()
  targets <- ca_targets()
  totals <- rows_of(ca_fit(out, targets), "population_total", "nusb")
  annual <- target_obs(targets, "population", "nusb")
  expect_identical(totals$observed[1], sum(annual$value[annual$start == totals$start[1]]))
  expect_identical(totals$simulated[1],
                   (out[, "N_F1"] + out[, "N_F2"])[year_to_idx(totals$start[1])] * 1e6)
})

test_that("published studies compare the simulated progression with each study", {
  study <- data.frame(Years = c(0, 1, 2), mean = c(1, 0.6, 0.3))
  f <- borgdorff_fit(fake_progression(), study)
  p_free <- tb_free_probability(fake_progression(), c(0, 1, 2))
  expect_identical(f$simulated, diff(-(1 - (1 - p_free) / (1 - p_free)[3])))
  expect_identical(f$observed, c(0.4, 0.3))
  bad <- fake_progression()
  bad$Mpfast <- bad$Mpfast[, 1:2]
  expect_identical(borgdorff_fit(bad, study)$simulated, c(NaN, NaN))
  cohort <- data.frame(Year = 1:2, At_Risk = c(1000, 900), Count = c(20, 0))
  f <- cohort_fit("ferebee", fake_progression(), cohort)
  expect_identical(f$total, c(1000, 900))
  expect_identical(f$observed, c(20, 0))
})

test_that("Tiemersma compares the simulated duration and fatality with 3 years and 0.9", {
  f <- tiemersma_fit(self_cure = 0.01, mortality = 0.02)
  expect_identical(f$observed, c(3, 0.9))
  expect_equal(f$simulated, c(1 / 0.36, 2 / 3))
})

test_that("the mortality risk rows pair each age and risk group's simulated count", {
  sim <- matrix(seq_len(70 * 44), 70, 44)
  f <- mort_dist_fit(sim)
  row <- f[f$group == "5_14_M3" & f$start == 2016, ]
  expect_identical(row$simulated, sim[year_to_idx(2016), 7] * 1e6)
  expect_identical(row$observed, rowSums(dist_gen)[[3]])
})
