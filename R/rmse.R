
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
#' @param data The observed data being fit to
#' @param local_epi_start_date_0 The date on which the model results start
#' @param local_obs_start_date_0 The data on which the observed data start
#' @param data Data frame The observed data
#' @param report_lag_0 integer The lag between model data and observed data
#' @param weights Vector Leght 2 vector containing weights indicating the relative weighting to give cases and deaths
#' @param under_report_factor fraction Estimated fraction of deaths that get reported as covid-19 deaths
#'
#' @export
#'
modelrmse <- function(modelOutput = modelOutput,
                      local_epi_start_date_0 = local_epi_start_date,
                      local_obs_start_date_0 = local_obs_start_date,
                      data = region_data,
                      report_lag_0 = report_lag,
                      weights = region_weights,
                      under_report_factor = 1){

  #model_end_date <- local_epi_start_date_0 + max(modelOutput$time) + report_lag_0

  country_data <- data #%>%
  #filter(date >= local_obs_start_date_0 & date <= model_end_date)

  target_features <- c(country_data$cum_cases / tail(country_data$cum_cases, 1),
                       country_data$cum_est_deaths / tail(country_data$cum_est_deaths, 1))

  target_features[is.na(target_features)] <- 0

  model.df_current <- modelOutput %>%
    dplyr::mutate(AllDeaths = D_s + D_h + D_c,
           local_epi_start_date = local_epi_start_date_0,
           date = local_epi_start_date_0 + time + report_lag_0)

  model.df_current <- model.df_current %>%
    dplyr::filter(date >= local_obs_start_date_0 & date <= max(country_data$date))

  model.df_current <- model.df_current %>%
    dplyr::mutate(ConfirmedCasesRescaled = ConfirmedCases / tail(country_data$cum_cases, 1),
           AllDeathsRescaled = AllDeaths / tail(country_data$cum_est_deaths, 1))


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
