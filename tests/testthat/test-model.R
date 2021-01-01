test_that("Modelling function runs", {
  par_file <- system.file("test-data", "test123_parameters.csv", package="conisi")
  par_table <- read.csv(par_file, sep=";")
  mod_result <- COVIDmodel(par_table, 1000000, 100)
  expect_length(mod_result, 32)
})
