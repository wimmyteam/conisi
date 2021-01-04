test_that("rmse works", {
  par_file <- system.file("test-data", "test123_parameters.csv", package="conisi")
  par_table <- read.csv(par_file, sep=";")

  mod_result <- COVIDmodel(par_table, 1000000, 100)
  mod_result <- dplyr::mutate(mod_result, experiment = 1)

  load(system.file("test-data", "observed_data.RData", package="conisi"))
  start_date <- test_data$date[2]

  rmse <- modelrmse(modelOutput = mod_result,
                    local_epi_start_date_0 = start_date,
                    local_obs_start_date_0 = start_date,
                    data = test_data,
                    report_lag_0 = 2,
                    weights = c(1, 0.5)
  )

  expect_gt(rmse, 0)
})
