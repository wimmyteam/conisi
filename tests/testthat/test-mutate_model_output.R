test_that("mutate_model_output works", {
  par_file <- system.file("test-data", "test123_parameters.csv", package="conisi")
  par_table <- read.csv(par_file, sep=";")

  pop <- 1000000
  mod_result <- COVIDmodel(par_table, pop, 100)
  #mod_result <- dplyr::mutate(mod_result, experiment = 1)

  par_df_wide <- par_table %>%
    dplyr::mutate(experiment = 1) %>%
    tidyr::pivot_wider(id_cols = c(experiment, start_time),
                       names_from = parameter_name,
                       values_from = value)


  mod_result <- mod_result %>%
    dplyr::left_join(par_df_wide, by = c("time" = "start_time")) %>%
    tidyr::fill(everything(), .direction = "down")

  mod_result <- conisi::mutate_model_output(mod_result, pop)

  expect_length(mod_result, 240)

})
