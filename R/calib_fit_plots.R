# Plots of a calibration fit: each target's observed values beside the
# simulated values they are compared with, drawn from the fit table of
# calib_fit(). Observed values are black points, with a line across the years
# of multi-year buckets; simulated values are blue.

fit_plot_colors <- c(
  observed = "#0b0b0b", simulated = "#2a78d6", grid = "#e4e3df", text = "#52514e"
)

# Targets compared as each period's distribution over groups. Their simulated
# values are on the model's own scale, so the plots show each group's share.
distribution_targets <- c(
  "cases_age", "tx_outcomes", "tb_deaths_age", "deaths_age", "population", "mort_dist"
)

# Targets compared as values per observed total: positive tests per person
# tested, and cases per person at risk.
rate_targets <- c("ltbi_prev", "ferebee", "sutherland")

# Targets whose groups are different quantities, each drawn in its own panel.
mixed_unit_targets <- "tiemersma"

# Figures with at least this many periods are drawn over time; others are drawn
# by group, one panel per period.
min_periods_over_time <- 3

fit_plot_titles <- c(
  cases_total = "TB cases",
  cases_decline = "TB cases in the decline from 1953",
  cases_nativity = "TB cases",
  cases_age = "Share of TB cases by age",
  cases_homeless = "TB cases among people experiencing homelessness",
  cases_recent_entry = "TB cases among recent arrivals",
  recent_transmission = "Share of TB cases from recent transmission",
  tx_outcomes = "Share of TB treatment outcomes",
  ltbi_tx_total = "LTBI treatment starts",
  ltbi_tx_share = "Shares of LTBI treatment starts",
  ltbi_prev = "IGRA-positive fraction by age",
  tb_deaths_total = "Deaths with TB",
  tb_deaths_age = "Share of deaths with TB by age",
  borgdorff = "Borgdorff: share of TB cases by years since infection",
  ferebee = "Ferebee: annual TB rate by year since infection",
  sutherland = "Sutherland: annual TB rate by year since infection",
  tiemersma = "Tiemersma: untreated TB",
  population_census = "Census population",
  population_total = "Population",
  population = "Share of population by age",
  homeless_pop = "Population experiencing homelessness",
  deaths_total = "All-cause deaths",
  deaths_age = "Share of all-cause deaths by age",
  mort_dist = "Share of each age group by mortality risk group"
)

nativity_labels <- c(usb = "US-born", nusb = "non-US-born")

#' The values a fit plot draws
#'
#' Distribution targets show each group's share of its period, and the
#' mortality risk distribution shows shares within each age group. Observed
#' shares are of the period's observed sum, or of its total when some groups
#' are suppressed. Rate targets show values per total, and the rest show the
#' values themselves.
#'
#' @param fit fit table from \code{calib_fit}
#' @return \code{fit} with the drawn values in \code{obs_y} and \code{sim_y}
#' @noRd
fit_plot_values <- function(fit) {
  fit$obs_y <- fit$observed
  fit$sim_y <- fit$simulated
  rates <- fit$target %in% rate_targets
  fit$obs_y[rates] <- fit$observed[rates] / fit$total[rates]
  fit$sim_y[rates] <- fit$simulated[rates] / fit$total[rates]
  shares <- fit$target %in% distribution_targets
  d <- fit[shares, , drop = FALSE]
  age_set <- ifelse(d$target == "mort_dist", sub("_M[0-9]+$", "", d$group), "")
  key <- paste(d$target, d$nativity, d$start, age_set)
  obs_sum <- stats::ave(d$observed, key, FUN = sum)
  obs_total <- ifelse(is.na(obs_sum), d$total, obs_sum)
  fit$obs_y[shares] <- d$observed / obs_total
  fit$sim_y[shares] <- d$simulated / stats::ave(d$simulated, key, FUN = sum)
  fit
}

#' Label of a group in plots: age groups as "15-24" and "85+"
#' @noRd
group_label <- function(group) {
  label <- sub("^([0-9]+)p", "\\1+", group)
  label <- sub("^([0-9]+)_([0-9]+)", "\\1-\\2", label)
  sub("_M([0-9]+)$", " M\\1", label)
}

#' Label of a period in plots: "2015" or "2015-2018"
#' @noRd
period_label <- function(start, end) {
  ifelse(start == end, as.character(start), paste0(start, "-", end))
}

#' The panels of a fit plot
#'
#' Each target and nativity group is one figure. A figure is drawn over time
#' when it covers enough periods, with one panel per group, and otherwise by
#' group, with one panel per period. The groups of mixed-unit targets each get
#' their own panel.
#'
#' @param fit fit table from \code{calib_fit}
#' @return list of figures, each a list of a \code{title} and its \code{panels};
#'   a panel is a list of a \code{title}, a \code{kind} of \code{"time"} or
#'   \code{"group"}, and its rows of the fit with the drawn values
#' @noRd
fit_figures <- function(fit) {
  fit <- fit_plot_values(fit)
  figure_keys <- unique(fit[, c("target", "nativity")])
  lapply(seq_len(nrow(figure_keys)), function(i) {
    target <- figure_keys$target[i]
    nativity <- figure_keys$nativity[i]
    d <- fit[fit$target == target & fit$nativity == nativity, , drop = FALSE]
    title <- fit_plot_titles[[target]]
    if (nativity %in% names(nativity_labels)) {
      title <- paste0(title, ", ", nativity_labels[[nativity]])
    }
    list(title = title, panels = figure_panels(d, target))
  })
}

#' The panels of one figure
#' @noRd
figure_panels <- function(d, target) {
  panel <- function(title, kind, rows) list(title = title, kind = kind, data = rows)
  groups <- unique(d$group)
  if (target %in% mixed_unit_targets) {
    return(lapply(groups, function(g) panel(group_label(g), "group", d[d$group == g, ])))
  }
  periods <- unique(d$start)
  if (length(periods) >= min_periods_over_time) {
    if (identical(groups, "all")) {
      return(list(panel("", "time", d)))
    }
    return(lapply(groups, function(g) panel(group_label(g), "time", d[d$group == g, ])))
  }
  lapply(periods, function(start) {
    rows <- d[d$start %in% start, , drop = FALSE]
    title <- if (is.na(start)) "" else period_label(start, rows$end[1])
    panel(title, "group", rows)
  })
}

#' Axis labels in thousands or millions for large values
#' @noRd
compact_labels <- function(ticks) {
  top <- max(abs(ticks))
  labels <- if (top >= 1e6) {
    paste0(format(ticks / 1e6, trim = TRUE), "M")
  } else if (top >= 1e4) {
    paste0(format(ticks / 1e3, trim = TRUE), "k")
  } else {
    format(ticks, trim = TRUE)
  }
  replace(labels, ticks == 0, "0")
}

#' Set up a panel's plotting region, grid and value axis
#'
#' @param bottom lines of margin below the panel, for its labels
#' @noRd
panel_frame <- function(xlim, y, title, bottom = 2.6) {
  top <- suppressWarnings(max(y, na.rm = TRUE))
  if (!is.finite(top) || top <= 0) top <- 1
  graphics::par(mar = c(bottom, 3, 2, 0.8))
  graphics::plot.new()
  graphics::plot.window(xlim, c(0, top * 1.08))
  ticks <- graphics::axTicks(2)
  graphics::abline(h = ticks, col = fit_plot_colors[["grid"]], lwd = 0.5)
  graphics::axis(2, at = ticks, labels = compact_labels(ticks), las = 1, lwd = 0,
                 col.axis = fit_plot_colors[["text"]], cex.axis = 0.7)
  graphics::title(main = title, cex.main = 0.8, font.main = 1, line = 0.6)
}

#' Draw a panel over time
#' @noRd
draw_time_panel <- function(p) {
  d <- p$data[order(p$data$start), , drop = FALSE]
  mid <- (d$start + d$end) / 2
  panel_frame(range(c(d$start, d$end)), c(d$obs_y, d$sim_y), p$title)
  graphics::axis(1, lwd = 0, lwd.ticks = 0.5, col.axis = fit_plot_colors[["text"]],
                 cex.axis = 0.7, padj = -1)
  graphics::lines(mid, d$sim_y, col = fit_plot_colors[["simulated"]], lwd = 1.5)
  graphics::points(mid, d$sim_y, col = fit_plot_colors[["simulated"]], pch = 19, cex = 0.35)
  span <- d$end > d$start
  graphics::segments(d$start[span], d$obs_y[span], d$end[span], d$obs_y[span],
                     col = fit_plot_colors[["observed"]], lwd = 1)
  graphics::points(mid, d$obs_y, col = fit_plot_colors[["observed"]], pch = 19, cex = 0.55)
}

#' Draw a panel by group
#' @noRd
draw_group_panel <- function(p) {
  d <- p$data
  x <- seq_len(nrow(d))
  labels <- if (identical(unique(d$group), "all")) "" else group_label(d$group)
  rotated <- nrow(d) > 6
  panel_frame(c(0.5, nrow(d) + 0.5), c(d$obs_y, d$sim_y), p$title,
              bottom = if (rotated) 4.5 else 2.6)
  graphics::axis(1, at = x, labels = labels, las = if (rotated) 2 else 1, lwd = 0,
                 col.axis = fit_plot_colors[["text"]], cex.axis = 0.65,
                 padj = if (rotated) 0.5 else -1)
  graphics::points(x - 0.12, d$obs_y, col = fit_plot_colors[["observed"]], pch = 19, cex = 0.7)
  graphics::points(x + 0.12, d$sim_y, col = fit_plot_colors[["simulated"]], pch = 19, cex = 0.7)
}

#' Rows and columns of a page holding \code{n} panels
#' @noRd
page_grid <- function(n) {
  if (n <= 2) c(n, 1) else if (n <= 4) c(2, 2) else if (n <= 6) c(2, 3) else if (n <= 9) c(3, 3)
  else c(3, 4)
}

#' Pages of a fit plot: the figures with one panel together, up to twelve to
#' a page, then each figure with several panels on its own pages
#'
#' @param figures figures from \code{fit_figures}
#' @param title title of the pages
#' @return list of pages, each a list of a \code{title} and its \code{panels}
#' @noRd
fit_pages <- function(figures, title) {
  per_page <- 12
  pages_of <- function(page_title, panels) {
    lapply(split(panels, ceiling(seq_along(panels) / per_page)), function(chunk) {
      list(title = page_title, panels = unname(chunk))
    })
  }
  single <- vapply(figures, function(figure) length(figure$panels) == 1, logical(1))
  single_panels <- lapply(figures[single], function(figure) {
    panel <- figure$panels[[1]]
    panel$title <- paste(c(figure$title, panel$title[nzchar(panel$title)]), collapse = ", ")
    panel
  })
  multi_pages <- lapply(figures[!single], function(figure) {
    pages_of(paste0(title, ": ", figure$title), figure$panels)
  })
  unname(c(pages_of(title, single_panels), unlist(multi_pages, recursive = FALSE)))
}

#' Draw one page of a fit plot, with its title and a legend
#' @noRd
draw_fit_page <- function(page) {
  grid <- page_grid(length(page$panels))
  graphics::par(mfrow = grid, oma = c(2, 0, 2.5, 0))
  for (p in page$panels) {
    if (p$kind == "time") draw_time_panel(p) else draw_group_panel(p)
  }
  graphics::mtext(page$title, side = 3, outer = TRUE, line = 0.8, cex = 1, font = 2)
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot.new()
  graphics::legend("bottom", legend = c("Observed", "Simulated"), horiz = TRUE, bty = "n",
                   col = fit_plot_colors[c("observed", "simulated")], pch = 19,
                   lty = c(0, 1), lwd = c(0, 1.5), pt.cex = c(0.8, 0.6), cex = 0.8,
                   text.col = fit_plot_colors[["text"]])
}

#' Plot a calibration fit
#'
#' Draws each calibration target's observed values beside the simulated values
#' they are compared with, one panel per target, nativity group and, for
#' distributions over time, group.
#'
#' @param fit fit table, from \code{calib_fit_for}
#' @param file PDF file to write
#' @param title title of each page, such as the location
#' @return \code{file}, invisibly
#' @export
plot_calib_fit <- function(fit, file, title = "Calibration fit") {
  pages <- fit_pages(fit_figures(fit), title)
  grDevices::pdf(file, width = 11, height = 8.5)
  on.exit(grDevices::dev.off())
  for (page in pages) {
    draw_fit_page(page)
  }
  invisible(file)
}
