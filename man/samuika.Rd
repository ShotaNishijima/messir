% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samuika.R
\encoding{UTF-8}
\name{samuika}
\alias{samuika}
\title{State-Space Assessment model Used for surume-IKA (samuika) ----}
\usage{
samuika(
  catch_data,
  weight_data,
  index_data,
  SR = "BH",
  regime_year = NULL,
  regime_par = c("a", "b", "sd")[c(1, 2, 3)],
  regime_key = 0:length(regime_year),
  stock_shared_par = c("a", "b", "sd", "SDlogF", "SDlogC")[c(1:5)],
  SDlogCPUE_key = rep(0, length(unique(index_data$Index_ID))),
  beta_key = rep(0, length(unique(index_data$Index_ID))),
  beta_fix = 1,
  logZ_mean = NULL,
  logZ_sd = NULL,
  logZ_weight = NULL,
  restrict_mean = FALSE,
  add_cpue = NULL,
  add_cpue_info = NULL,
  add_cpue_tday = NULL,
  fish_days = 180,
  add_cpue_covariate = NULL,
  add_cpue_SDkey = NULL,
  add_cpue_error_type = 0,
  add_cpue_all = 1,
  Fprocess_remove_year = NULL,
  p0_list = NULL,
  M = 0.6,
  Pope = FALSE,
  scale_num_to_mass = 0.1,
  bias_correct = TRUE,
  bias_correct_control = list(sd = FALSE),
  fixed_par = c("a", "b", "sd", "rec_rho", "SDlogF", "rho_SDlogF", "SDlogC", "q",
    "SDlogCPUE")[6],
  map_add = NULL,
  logF_diff = 0,
  SDlogF_init = 0.2,
  rho_SDlogF_init = 0,
  reca_init = c(3.5),
  recb_init = ifelse(SR == "HS", 14, 0.04),
  recSD_init = 0.4,
  rec_rho_init = 0,
  SDlogC_init = 0.2,
  q_init = c(0.22, 0.5),
  SDcpue_init = 0.2,
  N_init = 23.7,
  F_init = 0.38,
  silent = FALSE,
  maxrep = 100,
  HS_restrict = TRUE,
  nlminb_control = list(eval.max = 2000, iter.max = 2000)
)
}
\arguments{
\item{catch_data}{time-series of catch}

\item{weight_data}{time-series of weight}

\item{index_data}{time-series of indices}

\item{SR}{stock-recruit function ("BH", "RI", or "HS")}

\item{regime_year}{years when a regime shift occurred}

\item{regime_par}{parameters that changed by the regime shift (\code{c("a","b","sd")})}

\item{regime_key}{KEY representing regime (e.g., \code{c(0,1,0)} indicates the first and the third regimes were identical)}

\item{stock_shared_par}{parameter(s) that were shared between stocks (\code{c("a","b","sd","SDlogF","SDlogC")})}

\item{logZ_mean}{prior mean of logF (NOTICE: Not logZ!!)}

\item{logZ_sd}{prior SD of logF (NOTICE: Not logZ!!)}

\item{logZ_weight}{matrix for which years are included when using prior of \code{logZ_mean} and \code{logZ_sd} or \code{add_cpue}}

\item{add_cpue}{VECTOR of additional cpue time-series}

\item{add_cpue_info}{column 1: CPUE_ID, 2: Stock_ID, 3: Year_ID}

\item{add_cpue_tday}{VECTOR of days since the beginning of fishing season}

\item{fish_days}{duration of fishing season when using \code{add_cpue} (default: 180)}

\item{add_cpue_covariate}{array used as covariate(s) for additional CPUE}

\item{add_cpue_error_type}{error distribution for fitting additional CPUE; 0: lonormal, 1:log-Laplace, 2: gamma, 3: normal, 4: Laplace}

\item{add_cpue_all}{fitting additional CPUE not separately per-year (0: FALSE, 1 (default): using mean(F+Z), 2: using mean(log(F+Z))}

\item{M}{natural mortality coefficient (default: 0.6)}

\item{Pope}{whether the Pope approximation is used (TRUE) or not (FALSE: default)}

\item{scale_num_to_mass}{scaling multiplier in conversion of number to mass}

\item{bias_correct}{bias correct option in \code{sdreport}}

\item{fixed_par}{which parameters among \code{c("a","b","sd","rec_rho","SDlogF","rho_SDlogF","SDlogC","q","SDlogCPUE")} are fixed at their initial values}

\item{map_add}{added factor to map other than fixed_par}
}
\description{
State-Space Assessment model Used for surume-IKA (samuika) ----
}
