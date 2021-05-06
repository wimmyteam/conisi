test_that("rmse works", {
  par_file <- system.file("test-data", "test123_parameters.csv", package="conisi")
  par_table <- read.csv(par_file, sep=";")

  mod_result <- COVIDmodel(par_table, 1000000, 100)
  mod_result <- dplyr::mutate(mod_result, experiment = 1)

  load(system.file("test-data", "observed_data.RData", package="conisi"))
  # test_data

  start_date <- test_data$date[1]

  rmse <- modelrmse(modelOutput = mod_result,
                    start_date,
                    test_data,
                    weights = c(1, 0.5)
  )

  expect_lt(abs(rmse - 0.31747), 0.0001)

  smaller_target <- test_data %>%
    dplyr::filter(date >= "2020-04-01")

  rmse <- modelrmse(modelOutput = mod_result,
                    start_date,
                    smaller_target,
                    weights = c(1, 0.5)
  )
  expect_lt(abs(rmse - 0.358956), 0.0001)

  #TODO write more tests with varying start and end dates and with many experiments in single mod_result
  #TODO create test data with easily verifiable RMSE such as 1

})
