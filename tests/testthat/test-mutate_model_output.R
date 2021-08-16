test_that("mutate_model_output works", {
  par_file <- system.file("test-data", "test123_parameters.csv", package="conisi")
  par_table <- read.csv(par_file, sep=";")

  pop <- 1000000
  pop_prop <- c(0.03,0.39,0.58)
  contact_matrix <- c(13.6, 10.69, 15.71,0.80, 4, 7.2, 0.80, 4.90, 4)

  mod_result <- COVIDmodel(par_table, 1000000, 100, pop_prop, contact_matrix)
  #mod_result <- dplyr::mutate(mod_result, experiment = 1)

  par_df_wide <- par_table %>%
    dplyr::mutate(experiment = 1) %>%
    tidyr::pivot_wider(id_cols = c(experiment, start_time),
                       names_from = parameter_name,
                       values_from = value)


  mod_result <- mod_result %>%
    dplyr::left_join(par_df_wide, by = c("time" = "start_time")) %>%
    tidyr::fill(everything(), .direction = "down")

  mod_result <- conisi::mutate_model_output(mod_result, pop_size = pop, start_date = NULL, report_lag = 0, pop_prop)

  expect_length(mod_result, 884)

})
