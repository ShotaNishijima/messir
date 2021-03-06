% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/messir.R
\encoding{UTF-8}
\name{make_tmb_data}
\alias{make_tmb_data}
\title{Function for making input data for est_MSY with TMB}
\usage{
make_tmb_data(future_sim_res, obj_catch, objective, objective_value)
}
\arguments{
\item{future_sim_res}{\code{future_sim} object}

\item{obj_catch}{objective of maximum catch: mean(0) or geomean(1)}

\item{objective}{# 0: MSY, 1: PGY, 2: percentB0}

\item{objective_value}{# for PGY and percentB0}
}
\description{
Function for making input data for est_MSY with TMB
}
