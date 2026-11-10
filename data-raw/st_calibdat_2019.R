# Builds inst/ST/ST_CalibDat_2026-11-09.rds, the state calibration targets
# through 2019. Run from the package root: Rscript data-raw/st_calibdat_2019.R
#
# Starts from ST_CalibDat_10-29-21.rds, whose 5-year case buckets end in 2019,
# and
# - adds 2019 to tbdeaths and tbdeaths_age_yr from ST_CalibDat_04-20-22.rds,
#   after checking that the two files agree through 2018
# - corrects the labels of the 2015-2019 nativity bucket, the 1999-2003
#   recent-arrival bucket and the 5-14 age group

old <- readRDS("inst/ST/ST_CalibDat_10-29-21.rds")
new <- readRDS("inst/ST/ST_CalibDat_04-20-22.rds")
out <- old

through_2018 <- function(d, year_col) d[d[, year_col] <= 2018, , drop = FALSE]
for (st in seq_along(old$tbdeaths)) {
  stopifnot(isTRUE(all.equal(
    through_2018(old$tbdeaths[[st]], "Year"), through_2018(new$tbdeaths[[st]], "Year"),
    check.attributes = FALSE
  )))
}
stopifnot(isTRUE(all.equal(
  old$tbdeaths_age_yr, through_2018(new$tbdeaths_age_yr, "year"),
  check.attributes = FALSE
)))
out$tbdeaths <- new$tbdeaths
out$tbdeaths_age_yr <- new$tbdeaths_age_yr

relabel <- function(x, from, to) {
  stopifnot(sum(x == from) > 0)
  replace(x, x == from, to)
}
colnames(out$cases_nat_st_5yr) <-
  relabel(colnames(out$cases_nat_st_5yr), "X2015.1019", "X2015.2019")
out$rt_fb_cases_sm$yr_txt <- relabel(out$rt_fb_cases_sm$yr_txt, "99_13", "99_03")
out$cases_yr_ag_nat_st_5yr <- lapply(out$cases_yr_ag_nat_st_5yr, function(d) {
  colnames(d) <- relabel(colnames(d), "X14.May", "X5.14")
  d
})

saveRDS(out, "inst/ST/ST_CalibDat_2026-11-09.rds")
