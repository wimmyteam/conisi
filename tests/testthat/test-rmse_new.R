test_that("rmse works", {
  par_file <- system.file("test-data", "manually_edited_parameters.csv", package="conisi")
  par_table <- read.csv(par_file, sep=";")

  par_table <- conisi::load_parameters(par_file, sim_name = "sa_37", start_date = "2020-02-29")

  pop_prop <- c(0.03,0.39,0.58)
  contact_matrix <- c(13.6, 10.69, 15.71,0.80, 4, 7.2, 0.80, 4.90, 4)
  mod_result <- conisi::COVIDmodel(par_table, 42255960, 100, pop_prop, contact_matrix)
  mod_result <- dplyr::mutate(mod_result, experiment = 1)

  load(system.file("test-data", "observed_data.RData", package="conisi"))
  # test_data

  start_date <- test_data$date[1]

  rmse <- conisi::modelrmse(modelOutput = mod_result,
                            start_date,
                            test_data,
                            weights = c(1.0, 0.5, 0, 0),
  )

  expect_lt(abs(rmse - 0.002173373), 0.0001)

  smaller_target <- test_data %>%
    dplyr::filter(date >= "2020-04-01")

  rmse <- conisi::modelrmse(modelOutput = mod_result,
                            start_date,
                            smaller_target,
                            weights = c(1.0, 0.5, 0.25, 0.25)
  )
  expect_lt(abs(rmse - 0.001219555), 0.0001)

  #TODO write more tests with varying start and end dates and with many experiments in single mod_result
  #TODO create test data with easily verifiable RMSE such as 1

})
