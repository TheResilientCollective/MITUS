# Single years of age, 0 to 110, in each of the model's 11 age groups.
age_group_rows <- list(1:5, 6:15, 16:25, 26:35, 36:45, 46:55, 56:65, 66:75, 76:85, 86:95, 96:111)

#' Background mortality rates by age group for a location
#'
#' For the US, the national death counts over the national population. For a
#' state or county, the national rates scaled by the location's life table,
#' over the years the life table covers, weighted by the national population.
#' For the 5-14 age group the life table's ratios for ages 4-14 are matched to
#' the national rates for ages 5-14, with age 5 used twice.
#'
#' @name weight_mort
#' @param loc location code
#' @return matrix of yearly mortality rates, one row per year, named by year,
#'   and one column per age group
#' @export
weight_mort <- function(loc) {
  death_age <- national_by_age("US/US_MortalityCountsByAge.rds")
  popdist <- national_by_age("US/US_PopCountsByAge.rds")
  years <- intersect(colnames(death_age), colnames(popdist))
  death_age <- death_age[, years]
  popdist <- popdist[, years]
  if (loc == "US") {
    return(sapply(age_group_rows, function(rows) {
      colSums(death_age[rows, ]) / colSums(popdist[rows, ])
    }))
  }
  life_table <- location_life_table(loc)
  life_years <- colnames(life_table)
  if (!all(life_years %in% years)) {
    stop("The life table for ", loc, " covers years without national mortality data",
         call. = FALSE)
  }
  nat_mort <- (death_age / popdist)[, life_years]
  popdist <- popdist[, life_years]
  ratio <- life_table / nat_mort
  ratio_rows <- age_group_rows
  ratio_rows[[2]] <- 5:15
  sapply(seq_along(age_group_rows), function(k) {
    rows <- age_group_rows[[k]]
    colSums(nat_mort[rows, ] * popdist[rows, ] * ratio[ratio_rows[[k]], ]) /
      colSums(popdist[rows, ])
  })
}

#' Mortality rates by single year of age and year for a state or county
#'
#' Read from the location's own life table file if it has one, otherwise from
#' the all-state life tables.
#'
#' @return data frame with one column per year
#' @noRd
location_life_table <- function(loc) {
  loc_file <- find_loc_file(loc, "SingleYearLifeTable", data_dir = get_data_dir(),
                            required = FALSE)
  if (!is.null(loc_file)) {
    return(readRDS(loc_file)[, -1])
  }
  index <- which(as.data.frame(get_stateID())$USPS == loc)
  readRDS(system.file("ST/ST_SingleYearLifeTable.rds", package = "MITUS"))[[index]][, -1]
}

#' National counts by single year of age, one column per year, from a package
#' file with an age column followed by yearly columns
#'
#' Year columns are read by their leading four digits: the national population
#' counts label 1959 as "1959+".
#' @noRd
national_by_age <- function(file) {
  table <- readRDS(system.file(file, package = "MITUS"))
  year_cols <- grep("^[0-9]{4}", colnames(table))
  years <- substr(colnames(table)[year_cols], 1, 4)
  if (any(diff(as.integer(years)) != 1)) {
    stop(file, " must have one column per year", call. = FALSE)
  }
  counts <- as.matrix(table[, year_cols])
  dimnames(counts) <- list(table[, 1], years)
  counts
}
