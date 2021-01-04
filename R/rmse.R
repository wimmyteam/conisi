
se <- function(actual, predicted){
  return((actual - predicted) ^ 2)
}

mse <- function(actual, predicted, weights = rep(1, length(actual))){
  return(weighted.mean(x = se(actual, predicted), w = weights))
}

rmse <- function(actual, predicted, weights){
  return(sqrt(mse(actual, predicted, weights)))
}

myrmse <- function(model.features, target.features, weights){
  return(rmse(target.features, model.features, weights))
}

#' Calculates model RMSE
#'
#' @param modelOutput results from the model
#' @param start_date date The date which the first model time-point is being compared to.
#' @param target_data Data frame The observed data that the model is targeted to match.
#' @param weights Vector Vector of length 2 containing weights indicating the relative weighting to give cases and deaths.
#'   For example  c(1.0, 0.5)
#' @param under_report_factor fraction Estimated fraction of deaths that get reported as covid-19 deaths. This is deprecated and no longer used.
#'   The estimated_deaths in the variable "cum_est_deaths" allows for adjusting the death reporting rates.
#'
#' @export
#'
modelrmse <- function(modelOutput,
                      start_date,
                      target_data,
                      weights = c(1.0, 0.5),
                      under_report_factor = 1)
{

  # The last day of observed data to be used
  end_date <- start_date + max(modelOutput$time)

  last_observed_data <- max(target_data$date)

  # Filter the observed data so only looking at the part that overlaps with with model prediction
  target_data <- target_data %>%
  dplyr::filter(date >= start_date & date <= end_date)

  target_features <- c(target_data$cum_cases / tail(target_data$cum_cases, 1),
                       target_data$cum_est_deaths / tail(target_data$cum_est_deaths, 1))

  target_features[is.na(target_features)] <- 0

  model.df_current <- modelOutput %>%
    dplyr::mutate(AllDeaths = D_s + D_h + D_c,
#           local_epi_start_date = local_epi_start_date_0,
           date = start_date + time)

  model.df_current <- model.df_current %>%
    dplyr::filter(date >= start_date & date <= max(target_data$date))

  model.df_current <- model.df_current %>%
    dplyr::mutate(ConfirmedCasesRescaled = ConfirmedCases / tail(target_data$cum_cases, 1),
           AllDeathsRescaled = AllDeaths / tail(target_data$cum_est_deaths, 1))


  if(under_report_factor != 1) {
    model.df_current <- model.df_current %>%
      mutate(AllDeathsRescaled = AllDeathsRescaled * under_report_factor)
  }

  if("experiment" %in% model.df_current)
  {
    model.df_current <- model.df_current %>%
      dplyr::arrange(experiment, time) %>%
      dplyr::group_by(experiment) %>%
      dplyr::select(experiment, time, ConfirmedCasesRescaled, AllDeathsRescaled) %>%
      dplyr::ungroup()
  } else {
    model.df_current <- model.df_current %>%
      dplyr::arrange(time) %>%
      dplyr::select(time, ConfirmedCasesRescaled, AllDeathsRescaled)
  }

  model_current_wider.df <- model.df_current %>%
    tidyr::pivot_wider(names_from = "time",
                values_from = c("ConfirmedCasesRescaled", "AllDeathsRescaled"))

  weight_vector <- c(rep(weights[1], (length(target_features) / 2)),
                     rep(weights[2], (length(target_features) / 2)))

  if("experiment" %in% model.df_current)
  {
    model_rmse <- apply(dplyr::select(model_current_wider.df, -experiment),
                        MARGIN = 1,
                        FUN = myrmse,
                        target_features,
                        weight_vector)
  } else {
    model_rmse <- apply(model_current_wider.df,
                        MARGIN = 1,
                        FUN = myrmse,
                        target_features,
                        weight_vector)
  }
  return(model_rmse)
}
