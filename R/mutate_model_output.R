
#' This function processes the model output
#' and creates new variables that are needed for metrics of interest. The
#' output from this function is another dataframe / tibble.
#'
#' @param df Dataframe Output from COVIDmodel. The input data frame must have a row that represents
#'   each experiment and
#'   time point within the experiment. The columns must contain the parameter combinations that
#'   were used for each experiment and time point, as well as basic model output showing
#'   the numbers in each compartment at a given time.
#' @param pop Integer The size of the modeled population
#'
#' @param start_date Date This is the start date for the local epidemic and it is used for adding a
#'   date column to the output.
#'
#' @param  report_lag Integer This is the number of days we assume pass between the infection and
#'   when it is first reported.
#'
#' @export
#'
#' @importFrom magrittr %>%

mutate_model_output <- function(df, pop, start_date = NULL, report_lag = 0) {

  # used the calculate fractions without creating NaN values.
  # return NA for 0/0
  fraction <- function(a, b)
  {
    ifelse(b > 1e-10, a/b, NA)
  }

  # Do we need to create a date variable
  if(!is.null(start_date)){

  df <- df %>%
    mutate(date = start_date + time + report_lag)

  }

  # Create vars that apply to each row
  df1 <- df %>%
    dplyr::mutate(r_2d = r_2u,
           r_md = r_mu,
           c_12d = c_12u,
           c_1md = c_1mu,
           c_1sd = c_1su,
           eta_d = eta_u,
           delta_sd = delta_su,
           b_a = 0,
           d_e = 0,
           alpha_2 = a_2d * d_2 + a_2u * r_2d,
           alpha_m = a_md * d_m + 1 * r_md,
           chi_e = c_e1u + d_e,
           chi_u = c_12u + c_1mu + c_1su + d_1,
           chi_d = c_12d + c_1md + c_1sd,
           xi_2 = r_2u + d_2,
           xi_m = r_mu + d_m,
           zeta_u = eta_u + d_s + delta_su,
           zeta_d = eta_d + delta_sd,
           K = a_sd * c_1sd * c_e1u * d_1 * r_2d * r_md + a_sd * c_1sd * d_e * r_2d * r_md * chi_u,
           psi = a_1d * r_2d * r_md + a_2d * c_12d * r_md + a_md * c_1md * r_2d,
           G = a_1d * d_1 * r_2d * r_md + a_1u * chi_d * r_2d * r_md + a_2d * c_12d * d_1 * r_md + a_md * c_1md * d_1 * r_2d,
           A = ((((d_e * psi * chi_u + c_e1u * G) * zeta_u +
                    a_su * c_1su * c_e1u * chi_d * r_2d * r_md) * zeta_d +
                   K * zeta_u  + a_sd * c_1su * c_e1u * chi_d * d_s * r_2d * r_md) * xi_m +
                  c_e1u * c_1mu * r_2d * alpha_m * chi_d * zeta_u * zeta_d) * xi_2 +
             c_e1u * c_12u * r_md * alpha_2 * chi_d * xi_m * zeta_u * zeta_d,
           B = r_2d * r_md * chi_e * chi_u * chi_d * xi_2 * xi_m * zeta_u * zeta_d,
           R0 = (b_a + b_b) * A / B,
           Reff = FOIadjust * R0 * (S / pop),
           AllInfections = E_u + I_1u + I_2u + I_mu + I_su + E_d + I_1d + I_2d + I_md + I_sd + H + P + C,
           ActiveInfections = I_1u + I_2u + I_mu + I_su + I_1d + I_2d + I_md + I_sd,
           SymptKnownAsymptInfections = I_mu + I_su + I_1d + I_2d + I_md + I_sd,
           SymptKnownInfections = I_md + I_sd,
           SymptUnknownInfections = I_mu + I_su,
           AsymptKnownInfections = I_1d + I_2d,
           AsymptUnknownInfections = I_1u + I_2u,
           SymptInfections = I_mu + I_md + I_su + I_sd,
           AsymptInfections = I_1u + I_2u + I_1d + I_2d,
           KnownInfections = I_1d + I_2d + I_md + I_sd,
           UnknownInfections = I_1u + I_2u + I_mu + I_su,
           SevereKnownMildInfections = I_sd + I_su + I_md,
           SevereInfections = I_sd + I_su,
           Hospitalizations = H + P + C,
           Hosp_I_sd = I_sd + Hospitalizations,
           Hosp_SevereInfections = SevereInfections + Hospitalizations,
           Hosp_SevereKnownMildInfections = SevereKnownMildInfections + Hospitalizations,
           Hosp_SymptInfections = SymptInfections + Hospitalizations,
           Hosp_SymptKnownAsymptInfections = SymptKnownAsymptInfections + Hospitalizations,
           Hosp_ActiveInfections = ActiveInfections + Hospitalizations,
           Hosp_SymptKnownInfections = SymptKnownInfections + Hospitalizations,
           hosp_nonicu = H + P,
           deaths_hosp = D_h + D_c,
           NotWorking = KnownInfections + Hospitalizations,
           ReturnWork_cumul_flow = R_2d + R_md + R_c + R_h,
           AllDeaths = D_s + D_h + D_c,
           Prevalence = ActiveInfections / pop,
           Exposure = ContribAll / pop,
           FracSymptKnown = fraction(SymptKnownInfections, ActiveInfections),
           FracSymptUnknown = fraction(SymptUnknownInfections, ActiveInfections),
           FracAsymptKnown = fraction(AsymptKnownInfections, ActiveInfections),
           FracAsymptUnknown = fraction(AsymptUnknownInfections, ActiveInfections),
           FracHospSymptKnown = fraction(Hosp_SymptKnownInfections, Hosp_SymptInfections),
           FracAsymptKnown2 = fraction(AsymptKnownInfections, AsymptInfections),
           idf = fraction(KnownInfections + H + P + C, ActiveInfections + H + P + C),
           ifr = fraction(AllDeaths, ContribAll), # All deaths at a point in time / Cumulative Incidence (all people who have been infected)
           cfr = fraction(AllDeaths, ConfirmedCases))#, # Cum deaths at a point in time / all cumulative detected cases

  # Create vars that apply to each experiment
  df2 <- df1 %>%
    dplyr::arrange(experiment, time) %>%
    dplyr::group_by(experiment) %>%
    dplyr::mutate(AllDailyInfections = ContribAll - lag(ContribAll, default = 0),
                  NonSymptDailyInfections = ContribNonSympt - lag(ContribNonSympt, default = 0),
                  RelContribNonSympt = fraction(NonSymptDailyInfections, AllDailyInfections),
                  NewCases = ConfirmedCases - lag(ConfirmedCases, default = 0),
                  NewDeaths = AllDeaths - lag(AllDeaths, default = 0),
                  eta_d_flow = eta_d_cumul_flow - lag(eta_d_cumul_flow, default = 0),
                  eta_u_flow = eta_u_cumul_flow - lag(eta_u_cumul_flow, default = 0),
                  r_h_flow = r_h_cumul_flow - lag(r_h_cumul_flow, default = 0),
                  delta_h_flow = delta_h_cumul_flow - lag(delta_h_cumul_flow, default = 0),
                  theta_flow = theta_cumul_flow - lag(theta_cumul_flow, default = 0),
                  Symp_diagnozed_flow = Symp_diagnozed_cumul_flow - lag(Symp_diagnozed_cumul_flow, default = 0),
                  Asymp_diagnozed_flow = Asymp_diagnozed_cumul_flow - lag(Asymp_diagnozed_cumul_flow, default = 0),
                  ReturnWork_flow = ReturnWork_cumul_flow - lag(ReturnWork_cumul_flow, default = 0))

  # Create vars that apply to each time point, by summarized by experiment
  df3 <- df2 %>%
    dplyr::group_by(time) %>%
    dplyr::mutate(dplyr::across(.cols = c(Exposure, ConfirmedCases, NewCases, ActiveInfections, Prevalence,
                                          SevereInfections, I_sd, Hospitalizations, C,
                                          SymptInfections, AsymptInfections, NotWorking,
                                          KnownInfections, SymptKnownInfections, AsymptKnownInfections,
                                          AllDeaths, NewDeaths, Hosp_I_sd, Hosp_SevereInfections,
                                          Hosp_SevereKnownMildInfections, Hosp_SymptInfections,
                                          Hosp_SymptKnownAsymptInfections, Hosp_ActiveInfections,
                                          Hosp_SymptKnownInfections, hosp_nonicu, deaths_hosp,
                                          eta_d_flow, eta_u_flow, r_h_flow, delta_h_flow, theta_flow,
                                          Symp_diagnozed_flow, Asymp_diagnozed_flow, ReturnWork_flow),
                  .fns = list(mean = mean, min = min, max = max),
                  .names = "{col}_{fn}")) %>%
    dplyr::ungroup()

  #Prepend parameters with a "par_"
  df4 <- df3 %>%
    dplyr::rename_with(function(x){paste("par_", x)}, a_1d:theta | r_2d:zeta_d)

  return(df4)

}
