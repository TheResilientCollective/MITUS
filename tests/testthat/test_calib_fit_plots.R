# Fit rows of one target, for the plot tests.
plot_rows <- function(target, observed, simulated, start = NA_integer_, end = start,
                      nativity = "all", group = "all", total = NA_real_) {
  data.frame(target = target, start = as.integer(start), end = as.integer(end),
             nativity = nativity, group = group, observed = observed, simulated = simulated,
             total = total, stringsAsFactors = FALSE)
}

test_that("distributions are drawn as each period's shares", {
  fit <- plot_rows("deaths_age", c(1, 3, 2, 2), c(10, 30, 50, 50), rep(2018:2019, each = 2),
                   group = rep(c("0_4", "5_14"), 2))
  v <- fit_plot_values(fit)
  expect_equal(v$obs_y, c(0.25, 0.75, 0.5, 0.5))
  expect_equal(v$sim_y, c(0.25, 0.75, 0.5, 0.5))
})

test_that("complete periods take observed shares of their sum, whatever their total", {
  fit <- plot_rows("tx_outcomes", c(60, 20), c(3, 1), 2015, group = c("completed", "died"),
                   total = 100)
  expect_equal(fit_plot_values(fit)$obs_y, c(0.75, 0.25))
})

test_that("suppressed groups leave observed shares of the period's total", {
  fit <- plot_rows("cases_age", c(2, NA, 6), c(1, 1, 2), 2015, 2019, group = c("a", "b", "c"),
                   total = 10)
  v <- fit_plot_values(fit)
  expect_equal(v$obs_y, c(0.2, NA, 0.6))
  expect_equal(v$sim_y, c(0.25, 0.25, 0.5))
})

test_that("mortality risk shares are within each age group", {
  fit <- plot_rows("mort_dist", c(0.8, 0.2, 0.8, 0.2), c(3, 1, 1, 1), 2015,
                   group = c("0_4_M1", "0_4_M2", "5_14_M1", "5_14_M2"))
  expect_equal(fit_plot_values(fit)$sim_y, c(0.75, 0.25, 0.5, 0.5))
})

test_that("rate targets are drawn per total and counts as they are", {
  fit <- rbind(
    plot_rows("ltbi_prev", 8, 10, 2011, group = "15_24", total = 400),
    plot_rows("ferebee", 20, 30, group = "1", total = 1000),
    plot_rows("sutherland", 5, 4, group = "1", total = 100),
    plot_rows("cases_homeless", 500, 400, 2015, 2019, total = 9000)
  )
  v <- fit_plot_values(fit)
  expect_equal(v$obs_y, c(0.02, 0.02, 0.05, 500))
  expect_equal(v$sim_y, c(0.025, 0.03, 0.04, 400))
})

test_that("labels read as ages, periods and compact numbers", {
  expect_identical(group_label(c("0_4", "85p", "95p", "15_24_M3", "completed")),
                   c("0-4", "85+", "95+", "15-24 M3", "completed"))
  expect_identical(period_label(c(2015, 2015), c(2015, 2018)), c("2015", "2015-2018"))
  expect_identical(compact_labels(c(0, 2e6, 4e6)), c("0", "2M", "4M"))
  expect_identical(compact_labels(c(0, 2e4)), c("0", "20k"))
  expect_identical(compact_labels(c(0, 0.5)), c("0", "0.5"))
})

test_that("figures with enough periods are drawn over time, one panel per group", {
  years <- rep(2010:2012, each = 2)
  d <- fit_plot_values(plot_rows("tb_deaths_age", 1:6, 1:6, years, group = c("0_4", "5_14")))
  panels <- figure_panels(d, "tb_deaths_age")
  expect_identical(vapply(panels, `[[`, "", "kind"), c("time", "time"))
  expect_identical(vapply(panels, `[[`, "", "title"), c("0-4", "5-14"))
  single <- figure_panels(fit_plot_values(plot_rows("cases_total", 1:3, 1:3, 2010:2012)),
                          "cases_total")
  expect_length(single, 1)
  expect_identical(single[[1]]$title, "")
})

test_that("figures with few periods are drawn by group, one panel per period", {
  d <- fit_plot_values(plot_rows("mort_dist", rep(0.5, 4), 1:4, rep(2015:2016, each = 2),
                                 group = c("0_4_M1", "0_4_M2")))
  panels <- figure_panels(d, "mort_dist")
  expect_identical(vapply(panels, `[[`, "", "kind"), c("group", "group"))
  expect_identical(vapply(panels, `[[`, "", "title"), c("2015", "2016"))
  studies <- figure_panels(fit_plot_values(plot_rows("borgdorff", 1:3, 1:3, group = 1:3)),
                           "borgdorff")
  expect_length(studies, 1)
  expect_identical(nrow(studies[[1]]$data), 3L)
})

test_that("each group of a mixed-unit target gets its own panel", {
  d <- fit_plot_values(plot_rows("tiemersma", c(3, 0.9), c(2, 0.8),
                                 group = c("duration", "fatality")))
  panels <- figure_panels(d, "tiemersma")
  expect_identical(vapply(panels, `[[`, "", "title"), c("duration", "fatality"))
})

test_that("figures are titled by target and nativity group", {
  fit <- rbind(
    plot_rows("cases_nativity", 1:3, 1:3, c(2005, 2010, 2015), c(2009, 2014, 2019),
              nativity = "usb"),
    plot_rows("cases_total", 1:3, 1:3, 2017:2019)
  )
  expect_identical(vapply(fit_figures(fit), `[[`, "", "title"),
                   c("TB cases, US-born", "TB cases"))
})

test_that("single-panel figures share pages ahead of multi-panel figures", {
  one <- function(title) list(title = title, panels = list(list(title = "2010", kind = "group")))
  many <- list(title = "Ages", panels = rep(list(list(title = "", kind = "time")), 13))
  figures <- c(list(one("A")), list(many), rep(list(one("B")), 12))
  pages <- fit_pages(figures, "Loc")
  expect_identical(vapply(pages, `[[`, "", "title"), c("Loc", "Loc", "Loc: Ages", "Loc: Ages"))
  expect_identical(vapply(pages, function(p) length(p$panels), 0L), c(12L, 1L, 12L, 1L))
  expect_identical(pages[[1]]$panels[[1]]$title, "A, 2010")
})

test_that("plot_calib_fit draws every page to the file", {
  fit <- calib_fit(fake_outputs(), fake_progression(), c(rSlfCur = 0.01, muIp = 0.02),
                   ca_targets())
  file <- withr::local_tempfile(fileext = ".pdf")
  expect_identical(plot_calib_fit(fit, file, "CA"), file)
  expect_gt(file.size(file), 0)

  drawn <- character(0)
  local_mocked_bindings(draw_fit_page = function(page) {
    drawn <<- c(drawn, page$title)
    graphics::plot.new()
  })
  plot_calib_fit(fit, file, "CA")
  expect_identical(drawn, vapply(fit_pages(fit_figures(fit), "CA"), `[[`, "", "title"))
})
