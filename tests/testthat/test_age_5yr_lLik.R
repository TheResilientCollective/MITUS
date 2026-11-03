# The 5-year age likelihoods read CalibDatCases and wts from the global
# environment, where model_load() puts them.
with_age_data <- function(counts_us, counts_nus, shown_us, shown_nus, code) {
  g <- globalenv()
  saved <- mget(c("CalibDatCases", "wts"), envir = g, ifnotfound = list(NULL))
  on.exit(for (nm in names(saved)) {
    if (is.null(saved[[nm]])) rm(list = intersect(nm, ls(g)), envir = g)
    else assign(nm, saved[[nm]], envir = g)
  })
  ages <- paste0("X", 1:10)
  ag5 <- rbind(data.frame(usb = 1, setNames(as.data.frame(shown_us), ages)),
               data.frame(usb = 0, setNames(as.data.frame(shown_nus), ages)))
  buckets <- paste0("X", 1995 + 5 * (0:4), ".", 1999 + 5 * (0:4))
  cn <- data.frame(State.Code = 1, usb = c(1, 0),
                   rbind(rowSums(counts_us), rowSums(counts_nus)))
  colnames(cn)[3:7] <- buckets
  assign("CalibDatCases", list(cases_yr_ag_nat_st_5yr = list(ag5),
                               cases_nat_st_5yr = cn), envir = g)
  assign("wts", rep(1, 70), envir = g)
  code
}

# Model output whose 5-year sums match the counts, at an arbitrary scale.
matching_output <- function(counts) {
  V <- matrix(0, 25, 11)
  for (b in 1:5) V[(b - 1) * 5 + 1, 1:10] <- 3 * counts[b, ]
  V
}

set.seed(1)
us  <- matrix(sample(20:200, 50, replace = TRUE), 5, 10)
nus <- matrix(sample(20:200, 50, replace = TRUE), 5, 10)
us_na <- us;   us_na[2, 9:10] <- NA
nus_na <- nus; nus_na[4, 1] <- NA

test_that("a matching model scores near zero with complete data", {
  with_age_data(us, nus, us, nus, {
    expect_lt(abs(notif_age_us_5yr_lLik_st(matching_output(us), st = 1)), 0.01)
    expect_lt(abs(notif_age_nus_5yr_lLik_st(matching_output(nus), st = 1)), 0.01)
  })
})

test_that("a matching model scores near zero with suppressed age cells", {
  with_age_data(us, nus, us_na, nus_na, {
    expect_lt(abs(notif_age_us_5yr_lLik_st(matching_output(us), st = 1)), 0.01)
    expect_lt(abs(notif_age_nus_5yr_lLik_st(matching_output(nus), st = 1)), 0.01)
  })
})
