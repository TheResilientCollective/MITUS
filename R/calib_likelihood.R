# Calibration log-likelihood for one state or county: simulate the location
# with a parameter set, then compare the simulation with the location's
# calibration targets.

# Locations whose calibration leaves out a target.
locs_without_nativity_trend <- "MT"
locs_without_homeless_cases <- c("RI", "ID", "NH", "VT", "WY")

# Locations whose US-born case trend target is set in place of the data's.
fixed_usb_trend <- c(ND = -0.1)

# Adjustments for the 2020 pandemic: the fractional reductions in diagnosis and
# in transmission from March 2020 through 2021, and the relative risk of TB
# death over that period.
pandemic_2020 <- c(Dxt = 0.3707595, Trans = 0.1984619, CaseFat = 1.1158255)

# The latest year the TB and demographic terms compare. The importance weights
# must reach it.
last_compared_year <- 2019

# The likelihood is set to this value when a parameter set gives an invalid
# simulation, and to this value less a distance penalty when the likelihood
# cannot be evaluated.
invalid_llik <- -1e12

#' Calibration log-likelihood for parameter sets
#'
#' Simulates the location with each parameter set and compares the simulation
#' with the location's calibration targets.
#'
#' @param par_matrix parameter sets on the standard normal scale of the
#'   calibration, one per row; a single set may be a vector or a one-column
#'   matrix
#' @param loc state USPS code, or county code from \code{county_ID.csv}
#' @param targets calibration targets for \code{loc}, from \code{calib_targets}
#' @param n_cores number of parameter sets evaluated in parallel
#' @param tb whether to include the TB targets, or only the demographic ones
#' @param calib_end_year last simulated year; defaults to the last weighted year
#' @return log-likelihood of each parameter set
#' @noRd
calib_llikelihood <- function(par_matrix, loc, targets, n_cores = 1, tb = TRUE,
                              calib_end_year = NULL) {
  calib_end_year <- resolve_calib_end_year(calib_end_year)
  setup <- calib_setup()
  check_weight_years(setup$wts)
  par_matrix <- as.matrix(par_matrix)
  if (ncol(par_matrix) == 1) {
    par_matrix <- t(par_matrix)
  }
  if (ncol(par_matrix) != nrow(setup$ParamInitZ)) {
    stop("Each parameter set needs ", nrow(setup$ParamInitZ), " values, not ",
         ncol(par_matrix), call. = FALSE)
  }
  llik <- parallel::mclapply(
    seq_len(nrow(par_matrix)),
    function(i) {
      calib_llikelihood_one(par_matrix[i, ], loc, targets, setup, tb, calib_end_year)
    },
    mc.cores = n_cores
  )
  failed <- !vapply(llik, function(x) is.numeric(x) && length(x) == 1, logical(1))
  if (any(failed)) {
    stop("Calibration likelihood failed for parameter set ", which(failed)[1], ": ",
         format(llik[[which(failed)[1]]]), call. = FALSE)
  }
  unlist(llik)
}

#' Calibration log-likelihood for one parameter set
#'
#' @param par parameter set on the standard normal scale of the calibration
#' @param setup calibration setup from \code{calib_setup}
#' @inheritParams calib_llikelihood
#' @return log-likelihood; \code{invalid_llik} when the simulation is invalid,
#'   and \code{invalid_llik} less the squared distance of \code{par} from the
#'   calibration's starting values when the likelihood cannot be evaluated
#' @noRd
calib_llikelihood_one <- function(par, loc, targets, setup, tb, calib_end_year) {
  p <- calib_params(par, setup)
  llik <- tryCatch(
    {
      terms <- calib_llik_terms(p, loc, targets, setup$wts, tb, calib_end_year)
      if (is.null(terms)) invalid_llik else sum(terms)
    },
    error = function(e) NA_real_
  )
  if (is.na(llik) || is.infinite(llik)) {
    message("Calibration likelihood could not be evaluated; using the distance penalty")
    llik <- invalid_llik - sum((setup$ParamInitZ[, 8] - par)^2)
  }
  llik
}

#' Calibration log-likelihood terms for one set of model parameters
#'
#' Called by \code{calib_llikelihood_one}, and useful on its own to see each
#' term.
#'
#' @param p model parameter vector, from \code{calib_params}
#' @param wts importance weights by year
#' @inheritParams calib_llikelihood
#' @return named vector of log-likelihood terms, or NULL when the simulation is
#'   invalid
#' @noRd
calib_llik_terms <- function(p, loc, targets, wts, tb = TRUE, calib_end_year = NULL) {
  check_weight_years(wts)
  calib_end_year <- resolve_calib_end_year(calib_end_year)
  sim <- simulate_calibration(p, loc, calib_end_year)
  if (is.null(sim)) {
    return(NULL)
  }
  out <- sim$out
  terms <- if (tb) tb_llik_terms(out, sim$prms, p, loc, targets, wts) else numeric(0)
  c(terms, demographic_llik_terms(out, targets, wts))
}

#' Stop unless the importance weights reach the latest year the terms compare
#' @noRd
check_weight_years <- function(wts) {
  last_weight_year <- max(as.integer(names(wts)))
  if (last_weight_year < last_compared_year) {
    stop("The calibration likelihood compares years through ", last_compared_year,
         ", but the importance weights end in ", last_weight_year,
         "; call model_load() with last_weight_year of at least ", last_compared_year,
         call. = FALSE)
  }
}

#' Model globals the calibration reads, as set by \code{model_load}
#'
#' @return list of the background parameter values \code{P}; the logical index
#'   \code{ii} of the calibrated parameters within \code{P}; their priors
#'   \code{ParamInitZ}, one row per calibrated parameter, with the prior's two
#'   shape parameters in columns 6 and 7 and the starting value on the normal
#'   scale in column 8; the logical indexes \code{idZ0}, \code{idZ1} and
#'   \code{idZ2} of the parameters with beta, gamma and normal priors; and the
#'   importance weights \code{wts}
#' @noRd
calib_setup <- function() {
  needed <- c("P", "ii", "ParamInitZ", "idZ0", "idZ1", "idZ2", "wts")
  found <- vapply(needed, exists, logical(1), envir = globalenv(), inherits = FALSE)
  if (!all(found)) {
    stop("Run model_load() before calibrating; missing ", toString(needed[!found]),
         call. = FALSE)
  }
  mget(needed, envir = globalenv())
}

#' Model parameters for a parameter set on the calibration's normal scale
#'
#' Maps each calibrated parameter through the standard normal distribution to
#' a quantile of its prior.
#'
#' @param par calibrated parameters on the standard normal scale
#' @param setup calibration setup from \code{calib_setup}
#' @return the full model parameter vector
#' @noRd
calib_params <- function(par, setup) {
  u <- stats::pnorm(par, 0, 1)
  x <- u
  prior <- setup$ParamInitZ
  beta_prior <- setup$idZ0
  gamma_prior <- setup$idZ1
  normal_prior <- setup$idZ2
  x[beta_prior] <- stats::qbeta(
    u[beta_prior], shape1 = prior[beta_prior, 6], shape2 = prior[beta_prior, 7]
  )
  x[gamma_prior] <- stats::qgamma(
    u[gamma_prior], shape = prior[gamma_prior, 6], rate = prior[gamma_prior, 7]
  )
  x[normal_prior] <- stats::qnorm(
    u[normal_prior], mean = prior[normal_prior, 6], sd = prior[normal_prior, 7]
  )
  p <- setup$P
  p[setup$ii] <- x
  p
}

#' Apply the 2020 pandemic adjustments to the model parameters
#'
#' Diagnosis and transmission are reduced from March 2020 through 2021, then
#' return linearly to their February 2020 values over 2022-2023. The relative
#' risk of TB death is \code{CaseFat} from March 2020 through 2021, then falls
#' linearly to 1 over 2022-2023.
#'
#' @param prms model parameters from \code{param_init}, with monthly series
#' @return \code{prms} with adjusted \code{rDxt} and \code{NixTrans}, and
#'   \code{RRmuTBPand} added
#' @noRd
apply_pandemic_2020 <- function(prms) {
  lockdown <- year_month_to_idx(2020, 3):year_month_to_idx(2021, 12)
  recovery <- year_month_to_idx(2022, 1):year_month_to_idx(2023, 12)
  pre_pandemic <- year_month_to_idx(2020, 2)
  last_lockdown <- year_month_to_idx(2021, 12)
  prms$rDxt[lockdown, ] <- prms$rDxt[lockdown, ] -
    prms$rDxt[lockdown, ] * pandemic_2020[["Dxt"]]
  prms$NixTrans[lockdown] <- 1 - pandemic_2020[["Trans"]]
  for (group in seq_len(ncol(prms$rDxt))) {
    prms$rDxt[recovery, group] <- seq(
      prms$rDxt[last_lockdown, group], prms$rDxt[pre_pandemic, group],
      length.out = length(recovery)
    )
  }
  prms$NixTrans[recovery] <- seq(
    prms$NixTrans[last_lockdown], prms$NixTrans[pre_pandemic],
    length.out = length(recovery)
  )
  rr_mu <- rep(1, year_month_to_idx(2100, 12))
  rr_mu[c(lockdown, recovery)] <- c(
    rep(pandemic_2020[["CaseFat"]], length(lockdown)),
    seq(pandemic_2020[["CaseFat"]], 1, length.out = length(recovery))
  )
  prms$RRmuTBPand <- rr_mu
  prms
}

#' Simulate a location for calibration
#'
#' Besides the model parameters, the simulation uses the package constants
#' \code{can_go}, \code{HRdist} and \code{dist_gen_v} to rebalance the
#' population across risk groups. A message reports any rebalancing transition
#' probability above 1.
#'
#' @param p model parameter vector
#' @param loc location code
#' @param calib_end_year last simulated year
#' @return list of the annual outputs \code{out}, with named columns, and the
#'   model parameters \code{prms}; NULL when the simulation is invalid, meaning
#'   its 2017 outputs have missing or negative values or its state array
#'   \code{V1} has negative values
#' @noRd
simulate_calibration <- function(p, loc, calib_end_year) {
  prms <- param_init(p, loc, prg_chng = def_prgchng(p), ttt_list = def_ttt())
  prms <- apply_pandemic_2020(prms)
  trans_mat <- reblncd(
    mubt = prms[["mubt"]], can_go = can_go, RRmuHR = prms[["RRmuHR"]][2],
    RRmuRF = prms[["RRmuRF"]], HRdist = HRdist, dist_gen_v = dist_gen_v,
    adj_fact = prms[["adj_fact"]]
  )
  if (any(trans_mat > 1)) {
    message("Transition probabilities are too high")
  }
  cascade <- def_care_cascade()
  run <- cSim(
    setup_pars = c(year_to_idx(calib_end_year), length(func_ResNam()), length(model_age_groups)),
    rDxt = prms[["rDxt"]],
    TxQualt = prms[["TxQualt"]],
    InitPop = prms[["InitPop"]],
    Mpfast = prms[["Mpfast"]],
    ExogInf = prms[["ExogInf"]],
    MpfastPI = prms[["MpfastPI"]],
    Mrslow = prms[["Mrslow"]],
    rrSlowFB = prms[["rrSlowFB"]],
    rfast = prms[["rfast"]],
    RRcurDef = prms[["RRcurDef"]],
    rSlfCur = prms[["rSlfCur"]],
    p_HR = prms[["p_HR"]],
    vTMort = prms[["vTMort"]],
    RRmuRF = prms[["RRmuRF"]],
    RRmuHR = prms[["RRmuHR"]],
    Birthst = prms[["Birthst"]],
    HrEntEx = prms[["HrEntEx"]],
    ImmNon = prms[["ImmNon"]],
    ImmLat = prms[["ImmLat"]],
    ImmAct = prms[["ImmAct"]],
    ImmFst = prms[["ImmFst"]],
    Int1Test = prms[["Int1Test"]],
    Int1Init = prms[["Int1Init"]],
    Int1Tx = prms[["Int1Tx"]],
    net_mig_usb = prms[["net_mig_usb"]],
    net_mig_nusb = prms[["net_mig_nusb"]],
    RRmuTBPand = prms[["RRmuTBPand"]],
    mubt = prms[["mubt"]],
    RelInf = prms[["RelInf"]],
    RelInfRg = prms[["RelInfRg"]],
    RRcrAG = prms[["RRcrAG"]],
    Vmix = prms[["Vmix"]],
    rEmmigFB = prms[["rEmmigFB"]],
    TxVec = prms[["TxVec"]],
    TunTxMort = prms[["TunTxMort"]],
    rDeft = prms[["rDeft"]],
    ttt_samp_dist = prms[["ttt_sampling_dist"]],
    ttt_ag = prms[["ttt_ag"]],
    ttt_na = prms[["ttt_na"]],
    ttt_month = prms[["ttt_month"]],
    ttt_pop_scrn = prms[["ttt_pop_scrn"]],
    ttt_ltbi = prms[["ttt_ltbi"]],
    LtTxPar = prms[["LtTxPar"]],
    LtDxPar_lt = prms[["LtDxPar_lt"]],
    LtDxPar_nolt = prms[["LtDxPar_nolt"]],
    rrTestLrNoTb = prms[["rrTestLrNoTb"]],
    rrTestHr = prms[["rrTestHr"]],
    rLtScrt = prms[["rLtScrt"]],
    RRdxAge = prms[["RRdxAge"]],
    ttt_ltbi_init = cascade[1],
    ttt_ltbi_comp = cascade[2],
    ttt_ltbi_eff = cascade[3],
    ttt_ltbi_sens = cascade[4],
    ttt_ltbi_spec = cascade[5],
    ttt_ltbi_accept = cascade[6],
    rRecov = prms[["rRecov"]],
    pImmScen = prms[["pImmScen"]],
    EarlyTrend = prms[["EarlyTrend"]],
    pReTx = prms[["pReTx"]],
    ag_den = prms[["aging_denom"]],
    NixTrans = prms[["NixTrans"]],
    dist_gen = prms[["dist_gen"]],
    trans_mat_tot_ages = trans_mat
  )
  check_year <- run$Outputs[year_to_idx(2017), ]
  if (any(is.na(check_year)) || min(check_year) < 0 || min(run$V1) < 0) {
    message("Simulation has missing or negative values")
    return(NULL)
  }
  out <- run$Outputs
  colnames(out) <- prms[["ResNam"]]
  list(out = out, prms = prms)
}

#' Simulated outputs by model age group
#'
#' @param out annual simulation outputs with named columns
#' @param prefix output name before the age group, e.g. \code{"NOTIF_"}
#' @noRd
by_age <- function(out, prefix) out[, paste0(prefix, model_age_groups), drop = FALSE]

#' TB log-likelihood terms
#'
#' @param out annual simulation outputs with named columns
#' @param prms model parameters from \code{param_init}
#' @param p model parameter vector
#' @inheritParams calib_llikelihood
#' @param wts importance weights by year
#' @return named vector of log-likelihood terms
#' @noRd
tb_llik_terms <- function(out, prms, p, loc, targets, wts) {
  y2002 <- year_to_idx(2002)
  y2011 <- year_to_idx(2011)
  cases <- out[, "NOTIF_ALL"] + out[, "NOTIF_MORT_ALL"]
  cases_usb_age <- by_age(out, "NOTIF_US_") + by_age(out, "NOTIF_US_MORT_")
  cases_nusb_age <- (by_age(out, "NOTIF_") + by_age(out, "NOTIF_MORT_")) - cases_usb_age
  cases_usb <- out[, "NOTIF_US"] + out[, "NOTIF_MORT_US"]
  cases_nusb <- out[, "NOTIF_F1"] + out[, "NOTIF_F2"] +
    (out[, "NOTIF_MORT_F1"] + out[, "NOTIF_MORT_F2"])
  recent_years <- year_to_idx(2015:2018)
  recent_transmission <- sum(out[recent_years, c("INCID_REC_US", "INCID_REC_FB")]) /
    sum(out[recent_years, c("INCID_ALL_US", "INCID_ALL_FB")])
  tb_deaths <- by_age(out, "TOTMORT_W_TB_")
  progression <- prms[c("Mpfast", "Mrslow", "rfast", "rRecov")]

  terms <- c(
    cases_total = llik_cases_total(cases, targets, wts),
    cases_decline = llik_cases_decline(cases, targets, wts),
    cases_age_usb = llik_cases_age_5yr(cases_usb_age, 1, targets, wts),
    cases_age_nusb = llik_cases_age_5yr(cases_nusb_age, 0, targets, wts),
    cases_nusb = llik_cases_nativity_5yr(cases_nusb, 0, targets, wts),
    cases_usb = llik_cases_nativity_5yr(cases_usb, 1, targets, wts)
  )
  if (!loc %in% locs_without_nativity_trend) {
    usb_trend <- if (loc %in% names(fixed_usb_trend)) fixed_usb_trend[[loc]] else NA
    terms["cases_nativity_trend"] <- llik_cases_nativity_trend(
      cases_usb, cases_nusb, targets, usb_trend
    )
  }
  if (!loc %in% locs_without_homeless_cases) {
    cases_homeless <- out[, "NOTIF_HR"] + out[, "NOTIF_MORT_HR"]
    terms["cases_homeless"] <- llik_cases_homeless(cases_homeless, targets, wts)
  }
  c(
    terms,
    recent_transmission = llik_recent_transmission(recent_transmission, targets),
    cases_recent_entry = llik_cases_recent_entry(
      out[, "NOTIF_F1"] + out[, "NOTIF_MORT_F1"], targets, wts
    ),
    tx_outcomes = llik_tx_outcomes(
      out[, c("TBTX_COMPLT", "TBTX_DISCONT", "TBTX_DIED")], targets, wts
    ),
    ltbi_tx_total = llik_ltbi_tx_total(out[y2002, "TLTBI_INITS"], targets),
    ltbi_tx_dist = llik_ltbi_tx_dist(
      out[y2002, c("TLTBI_INITS_FB", "TLTBI_INITS_HR")] / out[y2002, "TLTBI_INITS"],
      targets
    ),
    ltbi_prev_usb = llik_ltbi_prev_usb(
      by_age(out, "N_US_LTBI_")[y2011, ], by_age(out, "N_US_")[y2011, ], targets
    ),
    ltbi_prev_nusb = llik_ltbi_prev_nusb(
      by_age(out, "N_FB_LTBI_")[y2011, ], by_age(out, "N_FB_")[y2011, ], targets
    ),
    tb_deaths_total = llik_tb_deaths_total(tb_deaths, targets, wts),
    tb_deaths_age = llik_tb_deaths_age(tb_deaths, targets, wts),
    borgdorff = llik_borgdorff(progression, targets),
    ferebee = llik_ferebee(progression, targets),
    sutherland = llik_sutherland(progression, targets),
    tiemersma = llik_tiemersma(p[["rSlfCur"]], p[["muIp"]])
  )
}

#' Demographic log-likelihood terms
#'
#' @inheritParams tb_llik_terms
#' @return named vector of log-likelihood terms
#' @noRd
demographic_llik_terms <- function(out, targets, wts) {
  deaths <- by_age(out, "TOTMORT_")
  mort_groups <- out[, paste0("N_ag_", rep(seq_along(model_age_groups), each = 4), "_M", 1:4)]
  c(
    pop_nusb = llik_pop_total(out[, "N_F1"] + out[, "N_F2"], 0, targets, wts),
    pop_usb = llik_pop_total(out[, "N_US"], 1, targets, wts),
    pop_age = llik_pop_age(by_age(out, "N_US_"), by_age(out, "N_FB_"), targets, wts),
    homeless_pop = llik_homeless_pop(out[year_to_idx(2010), "N_HR"], targets),
    deaths_total = llik_deaths_total(deaths, targets, wts),
    deaths_age = llik_deaths_age(deaths, targets, wts),
    mort_dist = llik_mort_dist(mort_groups, wts)
  )
}
