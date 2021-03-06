% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utlitika.R
\encoding{UTF-8}
\name{plot_catch}
\alias{plot_catch}
\title{Function for plotting time series of catch biomass}
\usage{
plot_catch(
  model,
  CI = 0.8,
  plot_CI = TRUE,
  Stock_name = NULL,
  plot_obs = TRUE,
  base_size = 18,
  path_size = 1.5,
  colour = "blue",
  point_size = 2,
  alpha = 0.4,
  scales = "fix",
  add_MSY = TRUE,
  MSY_path_size = 1,
  MSY_path_colour = "orange",
  MSY_path_linetype = "dashed"
)
}
\arguments{
\item{model}{\code{samuika} object}
}
\description{
Function for plotting time series of catch biomass
}
