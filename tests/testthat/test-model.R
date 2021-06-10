test_that("Modelling function runs", {
  par_file <- system.file("test-data", "test123_parameters.csv", package="conisi")
  par_table <- read.csv(par_file, sep=";")
  pop_prop <- c(0.03,0.39,0.58)
  contact_matrix <- c(13.6, 10.69, 15.71,0.80, 4, 7.2, 0.80, 4.90, 4)
  mod_result <- COVIDmodel(par_table, 1000000, 100, pop_prop, contact_matrix)
  expect_length(mod_result, 118)
})

test_that("Model helper functions work", {
  par_file <- system.file("test-data", "test123_parameters.csv", package="conisi")
  par_table <- read.csv(par_file, sep=";")

  pop <- 1000000
  pop_prop <- c(0.03,0.39,0.58)
  contact_matrix <- c(13.6, 10.69, 15.71,0.80, 4, 7.2, 0.80, 4.90, 4)

  mod_result <- COVIDmodel_run_and_mutate(par_table, pop, 100, pop_prop, contact_matrix)
  #mod_result <- dplyr::mutate(mod_result, experiment = 1)

  expect_length(mod_result, 875)

  availablecores <- 2
  doParallel::registerDoParallel(cores = availablecores)

  xparm_table <- dplyr::mutate(par_table, experiment = 1) %>%
    tibble::add_row(
      dplyr::mutate(par_table, experiment = 2)
      )

  mod_result_2 <- COVIDmodel_run_many(xparm_table, pop, 100, pop_prop, contact_matrix)
  expect_length(mod_result_2, 118)

  mod_result_3 <- COVIDmodel_run_and_mutate_many(xparm_table, pop, 100, pop_prop, contact_matrix)
  expect_length(mod_result3, 875)
})
