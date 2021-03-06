% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{COVIDmodel_run_and_mutate}
\alias{COVIDmodel_run_and_mutate}
\title{Runs the conisi COVID-19 model and applies mutate_model to the result}
\usage{
COVIDmodel_run_and_mutate(parm_table, pop_size, num_days, start_date = NULL)
}
\arguments{
\item{parm_table}{A data frame containing the time-varying parameters (in long format).}

\item{pop_size}{Integer The size of the population being modelled.}

\item{num_days}{Integer How many days to run the simulation for.}
}
\value{
A data frame with the values of the various compartments, parameters and other statistics over time
}
\description{
Runs the conisi COVID-19 model and applies mutate_model to the result
}
