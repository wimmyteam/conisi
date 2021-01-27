test_that("Modelling function runs", {
  par_file <- system.file("test-data", "test123_parameters.csv", package="conisi")
  par_table <- read.csv(par_file, sep=";")
  mod_result <- COVIDmodel(par_table, 1000000, 100)
  expect_length(mod_result, 36)
})

test_that("Model helper functions work", {
  par_file <- system.file("test-data", "test123_parameters.csv", package="conisi")
  par_table <- read.csv(par_file, sep=";")

  pop <- 1000000
  mod_result <- COVIDmodel_run_and_mutate(par_table, pop, 100)
  #mod_result <- dplyr::mutate(mod_result, experiment = 1)

  expect_length(mod_result, 240)

  availablecores <- 2
  doParallel::registerDoParallel(cores = availablecores)

  xparm_table <- dplyr::mutate(par_table, experiment = 1) %>%
    tibble::add_row(
      dplyr::mutate(par_table, experiment = 2)
      )

  mod_result_2 <- COVIDmodel_run_many(xparm_table, pop, 100)
  expect_length(mod_result_2, 69)

  mod_result_3 <- COVIDmodel_run_and_mutate_many(xparm_table, pop, 100)
  expect_length(mod_result, 240)
})
