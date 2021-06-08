
se <- function(actual, predicted){
  #min_length = min(length(actual), length(predicted))
  #return((actual[1:min_length] - predicted[1:min_length]) ^ 2)
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
                      weights = c(1.0, 0.5, 0.8),
                      under_report_factor = 1)
{
  #TODO this function needs a rewrite badly

  # keep the passed in start_date as reference to time = 0 equivalent
  date_zero <- start_date

  # Calculate the date-range where the predicted and target data intersects.
  # Hacky but, "date" used for target data in date format
  # "time" also used to indicate date but in "number of days"
  # from the start of the simulation format used in the model output
  end_date <- pmin(
    start_date + max(modelOutput$time),
    max(target_data$date)
  )
  end_time <- as.integer(end_date - date_zero)

  start_date <- pmax(
    start_date,
    min(target_data$date)
  )
  start_time <- as.integer(start_date - date_zero)

  if (start_time > end_time)
  {stop("Cannot calculate rmse, the predicted and actual values to not intersect in time.")}

  # Filter the observed data so only looking at the part that overlaps with with model prediction
  target_data <- target_data %>%
  dplyr::filter(date >= start_date & date <= end_date)

  modelOutput <- dplyr::filter(modelOutput, time >= start_time & time <= end_time)

  target_features <- c(target_data$cum_cases / tail(target_data$cum_cases, 1),
                       target_data$cum_est_deaths / tail(target_data$cum_est_deaths, 1),
                       target_data$cum_vaccinations / tail(target_data$cum_vaccinations, 1))

  target_features[is.na(target_features)] <- 0

  model.df_current <- modelOutput %>%
    dplyr::mutate(AllDeaths = D_s1 + D_h1 + D_c1 + D_s2 + D_h2 + D_c2 + D_s3 + D_h3 + D_c3,
#           local_epi_start_date = local_epi_start_date_0,
           date = date_zero + time)

  model.df_current <- model.df_current %>%
    dplyr::filter(date >= start_date & date <= max(target_data$date))

  model.df_current <- model.df_current %>%
    dplyr::mutate(ConfirmedCases = ConfirmedCases1 + ConfirmedCases2 + ConfirmedCases3,
                  Dose1Vaccinated = Vaccination_dose1_flow1 + Vaccination_dose1_flow2 + Vaccination_dose1_flow3,
                  FullyVaccinated = Vaccination_fully_flow1 + Vaccination_fully_flow2 + Vaccination_fully_flow3,
                  AllVaccinations = Dose1Vaccinated + FullyVaccinated,
                  ConfirmedCasesRescaled = ConfirmedCases / tail(target_data$cum_cases, 1),
                  AllDeathsRescaled = AllDeaths / tail(target_data$cum_est_deaths, 1),
                  AllVaccinationsRescaled = AllVaccinations / tail(target_data$cum_vaccinations, 1))


  if(under_report_factor != 1) {
    model.df_current <- model.df_current %>%
      mutate(AllDeathsRescaled = AllDeathsRescaled * under_report_factor)
  }

  if("experiment" %in% model.df_current)
  {
    model.df_current <- model.df_current %>%
      dplyr::arrange(experiment, time) %>%
      dplyr::group_by(experiment) %>%
      dplyr::select(experiment, time, ConfirmedCasesRescaled, AllDeathsRescaled, AllVaccinationsRescaled) %>%
      dplyr::ungroup()
  } else {
    model.df_current <- model.df_current %>%
      dplyr::arrange(time) %>%
      dplyr::select(time, ConfirmedCasesRescaled, AllDeathsRescaled, AllVaccinationsRescaled)
  }

  model_current_wider.df <- model.df_current %>%
    tidyr::pivot_wider(names_from = "time",
                values_from = c("ConfirmedCasesRescaled", "AllDeathsRescaled", "AllVaccinationsRescaled"))

  weight_vector <- c(rep(weights[1], (length(target_features) / 3)),
                     rep(weights[2], (length(target_features) / 3)),
                     rep(weights[3], (length(target_features) / 3)))

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
