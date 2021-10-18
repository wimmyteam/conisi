#'  Runs the conisi COVID-19 model
#'
#' Solves the system of differential equations that make up the conisi model
#' using the parameters and arguments provided.
#'
#' @param parm_table A data frame containing the time-varying parameters (in long format).
#' @param pop_size Integer The size of the population being modelled.
#' @param num_days Integer How many days to run the simulation for.
#' @param pop_prop Vector This vector contains population proportions for sub-populations
#' @param contact_matrix Vector The entries of the mixing matrix
#'
#' @return A data frame with the values of the various compartments over time
#'
#' @importFrom magrittr %>%
#' @export
#'
COVIDmodel <- function(parm_table, pop_size, num_days, pop_prop, contact_matrix, solver_method = "lsode"){
  # library(tidyverse)
  # library(deSolve)

  ngroups = length(pop_prop) # number of sub populations in the model

  N =  pop_size*pop_prop     # number in each group

  CM = matrix(contact_matrix, nrow = ngroups, byrow = TRUE)

  # Named vector containing, starting populations for compartments
  y <- c(S = N - c(1,3,6),
         E_u = c(1,3,6), #Number of seeds (people assumed to have entered the country with SARS-CoV-2 infection)
         I_1u = rep(0, ngroups),
         I_2u = rep(0, ngroups),
         I_mu = rep(0, ngroups),
         I_su = rep(0, ngroups),
         R_2u = rep(0, ngroups),
         R_mu = rep(0, ngroups),
         E_d = rep(0, ngroups),
         I_1d = rep(0, ngroups),
         I_2d = rep(0, ngroups),
         I_md = rep(0, ngroups),
         I_sd = rep(0, ngroups),
         R_2d = rep(0, ngroups),
         R_md = rep(0, ngroups),
         H = rep(0, ngroups),
         R_h = rep(0, ngroups),
         C = rep(0, ngroups),
         P = rep(0, ngroups),
         R_c = rep(0, ngroups),
         D_s = rep(0, ngroups),
         D_h = rep(0, ngroups),
         D_c = rep(0, ngroups),
         V_1 = rep(0, ngroups),
         V = rep(0, ngroups),
         E_uv = rep(0, ngroups),
         ConfirmedCases = rep(0, ngroups),
         ContribAll = rep(0, ngroups),
         ContribNonSympt = rep(0, ngroups),
         eta_d_cumul_flow = rep(0, ngroups),
         eta_u_cumul_flow = rep(0, ngroups),
         r_h_cumul_flow = rep(0, ngroups),
         delta_h_cumul_flow = rep(0, ngroups),
         theta_cumul_flow = rep(0, ngroups),
         Asymp_diagnozed_cumul_flow = rep(0, ngroups), # Cumul flow of diagnoses that were asymptomatic (including presymptomatic) at time of diagnosis
         Symp_diagnozed_cumul_flow = rep(0, ngroups), # Cumul flow of diagnoses that were symptomatic at time of diagnosis
         Asymp_inf_cumul_flow = rep(0, ngroups), # Cumul flow of infections that stay symptomatic
         Symp_inf_cumul_flow = rep(0, ngroups), # Cumul flow of infections that become symptomatic
         Vaccination_dose1_flow = rep(0, ngroups),
         Vaccination_fully_flow = rep(0, ngroups)
  )

  model <- function(t, y, parms){

    ncompartment = 40 #25 in actual, rest are to track other things
    ngroups = length(y)/ncompartment

    S = as.matrix(y[1:ngroups])
    E_u = as.matrix(y[(ngroups+1):(2*ngroups)])
    I_1u = as.matrix(y[(2*ngroups+1):(3*ngroups)])
    I_2u = as.matrix(y[(3*ngroups+1):(4*ngroups)])
    I_mu = as.matrix(y[(4*ngroups+1):(5*ngroups)])
    I_su = as.matrix(y[(5*ngroups+1):(6*ngroups)])
    R_2u = as.matrix(y[(6*ngroups+1):(7*ngroups)])
    R_mu = as.matrix(y[(7*ngroups+1):(8*ngroups)])
    E_d = as.matrix(y[(8*ngroups+1):(9*ngroups)])
    I_1d = as.matrix(y[(9*ngroups+1):(10*ngroups)])
    I_2d = as.matrix(y[(10*ngroups+1):(11*ngroups)])
    I_md = as.matrix(y[(11*ngroups+1):(12*ngroups)])
    I_sd = as.matrix(y[(12*ngroups+1):(13*ngroups)])
    R_2d = as.matrix(y[(13*ngroups+1):(14*ngroups)])
    R_md = as.matrix(y[(14*ngroups+1):(15*ngroups)])
    H = as.matrix(y[(15*ngroups+1):(16*ngroups)])
    R_h = as.matrix(y[(16*ngroups+1):(17*ngroups)])
    C = as.matrix(y[(17*ngroups+1):(18*ngroups)])
    P = as.matrix(y[(18*ngroups+1):(19*ngroups)])
    R_c = as.matrix(y[(19*ngroups+1):(20*ngroups)])
    D_s = as.matrix(y[(20*ngroups+1):(21*ngroups)])
    D_h = as.matrix(y[(21*ngroups+1):(22*ngroups)])
    D_c = as.matrix(y[(22*ngroups+1):(23*ngroups)])
    V_1 = as.matrix(y[(23*ngroups+1):(24*ngroups)])
    V = as.matrix(y[(24*ngroups+1):(25*ngroups)])
    E_uv = as.matrix(y[(25*ngroups+1):(26*ngroups)])
    ConfirmedCases = as.matrix(y[(26*ngroups+1):(27*ngroups)])
    ContribAll = as.matrix(y[(27*ngroups+1):(28*ngroups)])
    ContribNonSympt = as.matrix(y[(28*ngroups+1):(29*ngroups)])
    eta_d_cumul_flow = as.matrix(y[(29*ngroups+1):(30*ngroups)])
    eta_u_cumul_flow = as.matrix(y[(30*ngroups+1):(31*ngroups)])
    r_h_cumul_flow = as.matrix(y[(31*ngroups+1):(32*ngroups)])
    delta_h_cumul_flow = as.matrix(y[(32*ngroups+1):(33*ngroups)])
    theta_cumul_flow = as.matrix(y[(33*ngroups+1):(34*ngroups)])
    Asymp_diagnozed_cumul_flow = as.matrix(y[(34*ngroups+1):(35*ngroups)])
    Symp_diagnozed_cumul_flow = as.matrix(y[(35*ngroups+1):(36*ngroups)])
    Asymp_inf_cumul_flow = as.matrix(y[(36*ngroups+1):(37*ngroups)])
    Symp_inf_cumul_flow = as.matrix(y[(37*ngroups+1):(38*ngroups)])
    Vaccination_dose1_flow = as.matrix(y[(38*ngroups+1):(39*ngroups)])
    Vaccination_fully_flow = as.matrix(y[(39*ngroups+1):(40*ngroups)])

    parms <- parm_table %>%
      dplyr::filter((start_time <= t | is.na(start_time))& (t < end_time | is.na(end_time))) %>%
      dplyr::select(parameter_name, value) %>%
      tibble::deframe()

    with(as.list(c(y, parms)),{
      CM =  matrix(CM, nrow=ngroups, byrow = T)
      psi = c(psi_1,psi_2,psi_3)
      psi_d1 = c(psi_d1_1,psi_d1_2,psi_d1_3)
      psi_d2 = c(psi_d2_1,psi_d2_2,psi_d2_3)

      hdf <- ifelse(exists("hdf"), hdf, 1.0) # Uses default value of 1 if parameter is missing
      ddf <- ifelse(exists("ddf"), ddf, 1.0) # Uses default value of 1 if parameter is missing
      upsilon <- ifelse(exists("upsilon"), upsilon, 0.0) # Uses default value of 0.0 if parameter is missing

      II <- (a_1u * I_1u +
               a_2u * I_2u +
               1 * I_mu +
               a_su * I_su +
               a_1d * I_1d +
               a_2d * I_2d +
               a_md * I_md +
               a_sd * I_sd
      ) / N

      S_f_E_u <- FOIadjust * (S * q) * (CM %*% II)

      ContribNonSympt <- (
        a_1u * q * (CM %*% I_1u) +
          a_2u * q * (CM %*% I_2u) +
          a_1d * q * (CM %*% I_1d) +
          a_2d * q * (CM %*% I_2d)
      ) * FOIadjust * S / pop_size

      E_u_f_I_1u <- c_e1u * E_u

      I_1u_f_I_2u <- c_12u * I_1u

      I_1u_f_I_mu <- c_1mu * I_1u

      I_1u_f_I_su <- c_1su * I_1u

      I_1u_f_I_1d <- d_1 * I_1u

      I_2u_f_R_2u <- r_2u * I_2u

      I_2u_f_I_2d <- d_2 * I_2u

      I_mu_f_R_mu <- r_mu * I_mu

      I_mu_f_I_md <- d_m * I_mu

      eta_u_eff <- eta_u * (H < h_ceil) + (eta_u * exp(500 * (1 - H / h_ceil))) * (H >= h_ceil) # The effective rate, taking into account the ceiling v

      I_su_f_D_s <- (delta_su + eta_u - eta_u_eff) * I_su

      I_su_f_H <- eta_u_eff * I_su


      I_su_f_I_sd <- d_s * I_su


      E_d_f_I_1d <- c_e1u * E_d

      I_1d_f_I_2d <- c_12u * I_1d

      I_1d_f_I_md <- c_1mu * I_1d

      I_1d_f_I_sd <- c_1su * I_1d

      I_2d_f_R_2d <- r_2u * I_2d

      I_md_f_R_md <- r_mu * I_md

      eta_d_eff <- eta_u * (H < h_ceil) + (eta_u * exp(500 * (1 - H / h_ceil))) * (H >= h_ceil)

      I_sd_f_D_s <- (delta_su + (eta_u - eta_d_eff)) * I_sd

      I_sd_f_H <- eta_d_eff * I_sd

      H_f_R_h <- r_h * H

      theta_eff <- theta * (C < c_ceil) + (theta * exp(500 * (1 - C / c_ceil))) * (C >= c_ceil)

      H_f_D_h <- (delta_h * delta_h_adjust + (theta - theta_eff)) * H

      H_f_C <- theta_eff * H

      C_f_P <- rho * C

      C_f_D_c <- delta_c * delta_h_adjust * C

      P_f_R_c <- r_p * P


      # Adding the assumptions of waning immunity: flows from R_* to S, at rate upsilon (simplest possible flow assumptions: constant rate)
      R_2u_f_S <- upsilon * R_2u
      R_2d_f_S <- upsilon * R_2d
      R_mu_f_S <- upsilon * R_mu
      R_md_f_S <- upsilon * R_md
      R_h_f_S <- upsilon * R_h
      R_c_f_S <- upsilon * R_c

      # Adding vaccination: flows from S to V_1 and V and from V_1 to V
      S_f_V <- iota_1 * psi * S
      S_f_V_1 <- iota_2 * psi_d1 * S
      V_1_f_V <- iota_2 * psi_d2 * V_1
      V_1_f_E_uv <- iota_2 * FOIadjust * (V_1 * q) * (CM %*% II) * (1 - alpha_d1)

      c_s = dplyr::case_when(iota_1 == 1 & iota_2 == 0 ~ alpha,
                             iota_1 == 0 & iota_2 == 1 ~ alpha_d2,
                             iota_1 == 1 & iota_2 == 1 ~ (alpha + alpha_d2)/2,
                             TRUE ~ 1)

      V_f_E_uv <- FOIadjust * (V * q) * (CM %*% II) * (1 - c_s)

      c_12uv = dplyr::case_when(iota_1 == 1 | iota_2 == 1 ~ c_12uv,
                                TRUE ~ 0)

      E_uv_f_I_2u <- c_12uv * E_uv

      upsilon_v = dplyr::case_when(iota_1 == 1 | iota_2 == 1 ~ upsilon_v,
                                TRUE ~ 0)
      V_f_S <- upsilon_v * V

      # Defining the system of differential equations
      dS <- -S_f_E_u + R_2u_f_S + R_2d_f_S + R_mu_f_S + R_md_f_S + R_h_f_S + R_c_f_S - S_f_V - S_f_V_1 + V_f_S
      dE_u <- -E_u_f_I_1u + S_f_E_u
      dI_1u <- -I_1u_f_I_2u - I_1u_f_I_mu - I_1u_f_I_su - I_1u_f_I_1d + E_u_f_I_1u
      dI_2u <- -I_2u_f_R_2u - I_2u_f_I_2d + I_1u_f_I_2u + E_uv_f_I_2u
      dI_mu <- -I_mu_f_R_mu - I_mu_f_I_md + I_1u_f_I_mu
      dI_su <- -I_su_f_D_s - I_su_f_H - I_su_f_I_sd + I_1u_f_I_su
      dR_2u <- I_2u_f_R_2u - R_2u_f_S
      dR_mu <- I_mu_f_R_mu - R_mu_f_S
      dE_d <- -E_d_f_I_1d
      dI_1d <- -I_1d_f_I_2d - I_1d_f_I_md - I_1d_f_I_sd + I_1u_f_I_1d + E_d_f_I_1d
      dI_2d <- -I_2d_f_R_2d + I_2u_f_I_2d + I_1d_f_I_2d
      dI_md <- -I_md_f_R_md + I_mu_f_I_md + I_1d_f_I_md
      dI_sd <- -I_sd_f_D_s - I_sd_f_H + I_su_f_I_sd + I_1d_f_I_sd
      dR_2d <- I_2d_f_R_2d - R_2d_f_S
      dR_md <- I_md_f_R_md - R_md_f_S
      dH <- -H_f_R_h - H_f_D_h - H_f_C + I_su_f_H + I_sd_f_H
      dR_h <- H_f_R_h - R_h_f_S
      dC <- -C_f_P - C_f_D_c + H_f_C
      dP <- -P_f_R_c + C_f_P
      dR_c <- P_f_R_c - R_c_f_S
      dD_s <- I_su_f_D_s + I_sd_f_D_s
      dD_h <- H_f_D_h
      dD_c <- C_f_D_c
      dV_1 <- S_f_V_1 - V_1_f_V - V_1_f_E_uv
      dV <- S_f_V + V_1_f_V - V_f_E_uv - V_f_S
      dE_uv <- V_1_f_E_uv + V_f_E_uv - E_uv_f_I_2u
      dConfirmedCases <- I_1u_f_I_1d + I_2u_f_I_2d + I_mu_f_I_md + I_su_f_I_sd + I_su_f_H * hdf + I_su_f_D_s * ddf
      dContributionAll <- S_f_E_u
      dContributionNonSymptomatics <- ContribNonSympt # new infections caused by non-symptomatics
      eta_d_flow <- I_sd_f_H
      eta_u_flow <- I_su_f_H
      r_h_flow <- H_f_R_h
      delta_h_flow <- H_f_D_h
      theta_flow <- H_f_C
      Asymp_diagnozed_flow <- I_1u_f_I_1d + I_2u_f_I_2d
      Symp_diagnozed_flow <- I_mu_f_I_md + I_su_f_I_sd + I_su_f_H * hdf + I_su_f_D_s * ddf
      Asymp_inf_flow <- I_1u_f_I_2u + I_1d_f_I_2d
      Symp_inf_flow <- I_1u_f_I_mu + I_1u_f_I_su + I_1d_f_I_md + I_1d_f_I_sd
      dVaccination_dose1_flow <- S_f_V_1
      dVaccination_fully_flow <- S_f_V + V_1_f_V

      return(list(c(dS,
                    dE_u,
                    dI_1u,
                    dI_2u,
                    dI_mu,
                    dI_su,
                    dR_2u,
                    dR_mu,
                    dE_d,
                    dI_1d,
                    dI_2d,
                    dI_md,
                    dI_sd,
                    dR_2d,
                    dR_md,
                    dH,
                    dR_h,
                    dC,
                    dP,
                    dR_c,
                    dD_s,
                    dD_h,
                    dD_c,
                    dV_1,
                    dV,
                    dE_uv,
                    dConfirmedCases,
                    dContributionAll,
                    dContributionNonSymptomatics,
                    eta_d_flow,
                    eta_u_flow,
                    r_h_flow,
                    delta_h_flow,
                    theta_flow,
                    Asymp_diagnozed_flow,
                    Symp_diagnozed_flow,
                    Asymp_inf_flow,
                    Symp_inf_flow,
                    dVaccination_dose1_flow,
                    dVaccination_fully_flow)))

    })
  }


  tspan <- seq(0, num_days, 1)

  model_output <- as.data.frame(
    deSolve::ode(
      y,
      tspan,
      model,
      c(100, 1000, 10000),
      #These are not actually used, just passing it because "ode" wants a vector
      method = solver_method
    )
  )
  return(model_output)
}


model_result_extend <- function(mod_result, par_table) {
  par_df_wide <- par_table %>%
    tidyr::pivot_wider(id_cols = c(experiment, start_time),
                       names_from = parameter_name,
                       values_from = value)


  mod_result <- mod_result %>%
    dplyr::left_join(par_df_wide, by = c("time" = "start_time")) %>%
    tidyr::fill(everything(), .direction = "down")

  return(mod_result)
}

#' Runs the conisi COVID-19 model and applies mutate_model to the result
#'
#' @param parm_table A data frame containing the time-varying parameters (in long format).
#' @param pop_size Integer The size of the population being modelled.
#' @param num_days Integer How many days to run the simulation for.
#' @param pop_prop Vector This vector contains population proportions for sub-populations
#' @param contact_matrix Vector The entries of the mixing matrix
#' @param start_date Date This is the start date for the local epidemic and it is used for adding a date column to the output.
#' @param  report_lag Integer This is the number of days we assume pass between the infection and when it is first reported.
#'
#' @return A data frame with the values of the various compartments, parameters and other statistics over time
#'
#' @export
#'
COVIDmodel_run_and_mutate <- function(parm_table, pop_size, num_days, pop_prop, contact_matrix, start_date = NULL, report_lag = 0)
{
  mod_result <- COVIDmodel(parm_table, pop_size, num_days, pop_prop, contact_matrix)

  if(! "experiment" %in% colnames(parm_table)){
    parm_table <- dplyr::mutate(parm_table, experiment = 1)
  }

  par_df_wide <- parm_table %>%
    tidyr::pivot_wider(id_cols = c(experiment, start_time),
                       names_from = parameter_name,
                       values_from = value)

  mod_result <- mod_result %>%
    dplyr::left_join(par_df_wide, by = c("time" = "start_time")) %>%
    tidyr::fill(everything(), .direction = "down")

  mod_result <- conisi::mutate_model_output(mod_result, pop_size, start_date, report_lag, pop_prop)

  return(mod_result)
}

#' Runs the conisi COVID-19 model on multiple parameter sets at once
#'
#' @param parm_table A data frame containing the time-varying parameters (in long format) for multiple experiments.
#' @param pop_size Integer The size of the population being modelled.
#' @param num_days Integer How many days to run the simulation for.
#' @param pop_prop A vector with population proportions for sub-population groups
#' @param contact_matrix A vector with entries of the mixing matrix
#'
#' @return A data frame with the values of the various compartments, parameters and other statistics over time
#' @importFrom foreach %dopar%
#' @export
COVIDmodel_run_many <- function(parm_table, pop_size, num_days, pop_prop, contact_matrix){

  experiments <- unique(parm_table$experiment)

  mod_results <- foreach::foreach (j = experiments, .combine = dplyr::bind_rows)  %dopar% {

    single_experiment_params <- parm_table %>%
      dplyr::filter(experiment == j)

    exp_result <- COVIDmodel(parm_table = single_experiment_params,
                              pop_size,
                              num_days,
                              pop_prop,
                              contact_matrix)

    model_result_extend(exp_result, parm_table)

  }
  return(mod_results)
}

#' Runs the conisi COVID-19 model on multiple parameter sets at once and applies mutate_model to the result
#'
#' @param parm_table A data frame containing the time-varying parameters (in long format) for multiple experiments.
#' @param pop_size Integer The size of the population being modelled.
#' @param num_days Integer How many days to run the simulation for.
#' @param pop_prop Vector This vector contains population proportions for sub-populations
#' @param contact_matrix Vector The entries of the mixing matrix
#' @param start_date Date This is the start date for the local epidemic and it is used for adding a date column to the output.
#' @param  report_lag Integer This is the number of days we assume pass between the infection and when it is first reported.
#' @return A data frame with the values of the various compartments
#'
#' @export
#'
COVIDmodel_run_and_mutate_many <- function(parm_table, pop_size, num_days, pop_prop, contact_matrix, start_date = NULL, report_lag = 0){
  mod_result <- COVIDmodel_run_many(parm_table, pop_size, num_days, pop_prop, contact_matrix)
  # conisi::mutate_model_output(mod_result, pop_size)
  mutate_model_output(mod_result, pop_size, start_date, report_lag, pop_prop)
}
