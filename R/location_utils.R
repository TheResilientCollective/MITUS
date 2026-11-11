#' Location and Year Index Utilities
#'
#' Helper functions for converting between calendar years and model indices,
#' and for discovering location-specific data files within the MITUS package.
#' The MITUS model uses 1950 as the base year (index 1).

#' Convert a calendar year to a model index
#'@name year_to_idx
#'@param year Calendar year (numeric)
#'@param base_year Base year of the model (default: 1950)
#'@return Integer index (1-based, where base_year = 1)
#'@export
year_to_idx<-function(year,base_year=1950){
  year-base_year+1
}

#' Convert a calendar year and month to a monthly model index
#'@name year_month_to_idx
#'@param year Calendar year (numeric)
#'@param month Month (1-12)
#'@param base_year Base year of the model (default: 1950)
#'@return Integer monthly index (1-based, where base_year January = 1)
#'@export
year_month_to_idx<-function(year,month,base_year=1950){
  (year-base_year)*12+month
}

#' Importance weights used to weight calibration targets by year.
#'
#' Returns a vector of per-year weights spanning \code{base_year:last_year},
#' named by calendar year. Weights follow a logistic curve rising from
#' \code{floor_val} at \code{base_year} to \code{ramp_end_val + floor_val} at
#' \code{ramp_end}, so recent years carry more weight in the objective.
#'
#'@param last_year Last calendar year covered by the weights
#'@param base_year First calendar year covered; the model base year
#'@param ramp_start,ramp_end Years over which the logistic curve rises
#'@param ramp_end_val Curve value attained at \code{ramp_end}
#'@param floor_val Constant added to every weight
#'@return Named numeric vector of length \code{last_year - base_year + 1}
#'@noRd
make_impt_weights <- function(last_year    = 2019,
                              base_year    = 1950,
                              ramp_start   = 2000,
                              ramp_end     = 2019,
                              ramp_end_val = 0.95,
                              floor_val    = 0.05) {
  if (last_year < base_year) {
    stop("last_year (", last_year, ") must be >= base_year (", base_year, ")")
  }
  z    <- log(1/0.005 - 1)
  step <- (2*z)/(ramp_end - ramp_start)
  from <- -z*(1 + 2*(ramp_start - base_year)/(ramp_end - ramp_start))
  yrs  <- base_year:last_year
  wts  <- as.numeric(ramp_end_val)/(1 + exp(-(from + (yrs - base_year)*step))) + floor_val
  names(wts) <- yrs
  wts
}

#' Last simulated year for a calibration run.
#'
#' The likelihood reads simulator output for every year in \code{names(wts)},
#' so the simulation must run at least through the last weighted year.
#'
#'@param calib_end_year Last simulated year, or NULL for the last year of the
#'  importance weights \code{wts} set by \code{model_load}
#'@return Integer calendar year
#'@noRd
resolve_calib_end_year <- function(calib_end_year = NULL) {
  last_weight_year <- max(as.integer(names(wts)))
  if (is.null(calib_end_year)) {
    return(last_weight_year)
  }
  if (calib_end_year < last_weight_year) {
    stop("calib_end_year (", calib_end_year, ") is before the last weighted year (",
         last_weight_year, "); simulate through at least that year or lower ",
         "last_weight_year in model_load()")
  }
  as.integer(calib_end_year)
}

#' Find a location-specific data file
#'
#' Searches for RDS files matching {loc}_{filetype}[_.]*.rds. When data_dir
#' is provided, searches there first (e.g., an external input directory).
#' Falls back to the installed package inst/ directory via system.file().
#' Supports dated filenames (ST_CalibDat_04-20-22.rds) and plain filenames
#' (SanDiego_CalibDat.rds). If multiple matches, returns most recently modified.
#'@name find_loc_file
#'@param loc Location code (e.g., "CA", "SanDiego", "US", "ST")
#'@param filetype File type pattern (e.g., "ModelInputs", "CalibDat")
#'@param fallback_loc Optional fallback location to try if no file found
#'@param data_dir Optional external directory to search before package inst/
#'@param required If TRUE (default), stop with error when no file found
#'@return Full file path to the matched RDS file
#'@export
find_loc_file<-function(loc,filetype,fallback_loc=NULL,data_dir=NULL,required=TRUE){
  pattern<-paste0("^",loc,"_?",filetype,"([_.].*)?\\.rds$")
  fb_pattern<-if(!is.null(fallback_loc)) paste0("^",fallback_loc,"_?",filetype,"([_.].*)?\\.rds$")

  # Parse calibration-vintage date from a filename, trying common formats:
  # YYYY-MM-DD (e.g. 2022-07-08), MM-DD-YY (e.g. 04-20-22), MMDDYY (e.g. 111620).
  # Returns Date or NA. We prefer date-in-filename over filesystem mtime because
  # mtime is non-deterministic across checkouts (a freshly-cloned tree gives all
  # files identical mtimes), while filename dates are stable and meaningful.
  # Future-dated matches (e.g. ST_CalibDat_09-22-28 typo for ...18) are treated
  # as NA — a calibration file dated past today is almost certainly a typo.
  today<-Sys.Date()
  parse_filename_date<-function(filename){
    bn<-sub("\\.rds$","",basename(filename))
    check<-function(d){
      if(is.na(d)) return(as.Date(NA))
      if(d>today) return(as.Date(NA))
      d
    }
    m<-regmatches(bn,regexpr("\\d{4}-\\d{2}-\\d{2}",bn))
    if(length(m)>0 && nchar(m)>0){
      d<-check(suppressWarnings(as.Date(m,"%Y-%m-%d")))
      if(!is.na(d)) return(d)
    }
    m<-regmatches(bn,regexpr("\\d{2}-\\d{2}-\\d{2}",bn))
    if(length(m)>0 && nchar(m)>0){
      d<-check(suppressWarnings(as.Date(m,"%m-%d-%y")))
      if(!is.na(d)) return(d)
    }
    m<-regmatches(bn,regexpr("(?<![0-9])\\d{6}(?![0-9])",bn,perl=TRUE))
    if(length(m)>0 && nchar(m)>0){
      d<-check(suppressWarnings(as.Date(m,"%m%d%y")))
      if(!is.na(d)) return(d)
    }
    as.Date(NA)
  }
  # Pick the file with the latest filename-date. Files with no parseable date
  # rank below any dated file; ties (or no dated file at all) fall back to mtime.
  pick_canonical<-function(files){
    if(length(files)==1) return(files)
    dates<-as.Date(sapply(files,parse_filename_date),origin="1970-01-01")
    has_date<-!is.na(dates)
    if(any(has_date)){
      candidates<-files[has_date]
      cand_dates<-dates[has_date]
      max_date<-max(cand_dates)
      finalists<-candidates[cand_dates==max_date]
      if(length(finalists)==1) return(finalists)
      return(finalists[which.max(file.mtime(finalists))])
    }
    files[which.max(file.mtime(files))]
  }

  # helper: search a directory for matching files
  search_dir<-function(base_dir,prefix,pat){
    d<-file.path(base_dir,prefix)
    if(!dir.exists(d)) return(NULL)
    files<-list.files(d,pattern=pat,full.names=TRUE,ignore.case=TRUE)
    if(length(files)==0) return(NULL)
    picked<-pick_canonical(files)
    if(length(files)>1){
      warning(paste0("Multiple files matched for ",prefix,"/",filetype,
                     "; using ",basename(picked)))
    }
    picked
  }

  # search external data_dir first (e.g., SIMULATOR_INPUT_DIR)
  if(!is.null(data_dir)){
    f<-search_dir(data_dir,loc,pattern)
    if(!is.null(f)) return(f)
    if(!is.null(fallback_loc)){
      f<-search_dir(data_dir,fallback_loc,fb_pattern)
      if(!is.null(f)) return(f)
    }
  }
  # fall back to installed package inst/ directory
  pkg_dir<-system.file(loc,package="MITUS")
  if(nchar(pkg_dir)>0){
    files<-list.files(pkg_dir,pattern=pattern,full.names=TRUE,ignore.case=TRUE)
    if(length(files)>0){
      picked<-pick_canonical(files)
      if(length(files)>1){
        warning(paste0("Multiple files matched for ",loc,"/",filetype,
                       "; using ",basename(picked)))
      }
      return(picked)
    }
  }
  if(!is.null(fallback_loc)){
    fb_pkg_dir<-system.file(fallback_loc,package="MITUS")
    if(nchar(fb_pkg_dir)>0){
      files<-list.files(fb_pkg_dir,pattern=fb_pattern,full.names=TRUE,ignore.case=TRUE)
      if(length(files)>0){
        picked<-pick_canonical(files)
        if(length(files)>1){
          warning(paste0("Multiple files matched for ",fallback_loc,"/",filetype,
                         "; using ",basename(picked)))
        }
        return(picked)
      }
    }
  }
  if(required){
    searched<-loc
    if(!is.null(fallback_loc)) searched<-paste0(loc," or ",fallback_loc)
    stop(paste0("Required data file not found: ",filetype," for location '",
                searched,"'"))
  }
  return(NULL)
}

# Package state: values model_load sets for the loaded location, kept outside
# the user's workspace.
.mitus_state <- new.env(parent = emptyenv())

#' Get the data directory passed to model_load (or NULL)
#'
#' Returns the external data directory \code{model_load} was called with, or
#' NULL if it was called without one or has not been called.
#'@name get_data_dir
#'@return Directory path (string) or NULL
#'@export
get_data_dir <- function() {
  if (exists("data_dir", envir = .mitus_state, inherits = FALSE)) {
    return(get("data_dir", envir = .mitus_state, inherits = FALSE))
  }
  NULL
}

#' Get the state table
#'
#' Returns the package's \code{stateID} table of the 50 states and DC: one row
#' per state, in the order of the all-state data files, with its name, FIPS
#' code and USPS code.
#'
#'@name get_stateID
#'@return stateID matrix
#'@export
get_stateID <- function() {
  e <- new.env()
  utils::data("stateID", package = "MITUS", envir = e)
  e$stateID
}

#' Look up a location code
#'
#' Returns the US as national, a USPS code in the state table as a state, and
#' any other code as a county. A county's name and state come from its
#' optional \code{Location} file, a list with \code{name} and \code{state}, the
#' USPS code of its state; without one, its name is its code.
#'
#'@name resolve_location
#'@param loc location code
#'@param data_dir directory searched for a county's \code{Location} file before
#'  the package's \code{inst/}
#'@return list with \code{name}, \code{code} and \code{loc_type}
#'  (\code{"national"}, \code{"state"} or \code{"county"}); for a state also
#'  \code{fips} and \code{st}, its row in the state table; for a county also
#'  \code{state}, NULL when not given
#'@export
resolve_location <- function(loc, data_dir = get_data_dir()) {
  if (loc == "US") {
    return(list(name = "United States", code = "US", loc_type = "national"))
  }
  states <- as.data.frame(get_stateID(), stringsAsFactors = FALSE)
  st <- which(states$USPS == loc)
  if (length(st) > 0) {
    return(list(name = states$Name[st], code = loc, loc_type = "state",
                fips = states$FIPS[st], st = st))
  }
  file <- find_loc_file(loc, "Location", data_dir = data_dir, required = FALSE)
  about <- if (is.null(file)) list() else readRDS(file)
  list(name = if (is.null(about$name)) loc else about$name, code = loc,
       loc_type = "county", state = about$state)
}
