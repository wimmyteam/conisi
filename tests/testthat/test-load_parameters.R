

test_that("loading parameters works", {

  expected_columns <- c("sim_name", "parameter_name", "value", "start_time", "end_time")

  file_path <- system.file("test-data", "manually_edited_parameters.csv", package="conisi")
  parameters <- conisi::load_parameters(file_path, sim_name = "sa_37", start_date = "2020-03-10")

  matching_cols <- intersect(expected_columns, names(parameters))
  expect_equal(expected_columns, matching_cols)
  expect_equal(nrow(parameters), 101)

  foi_adjust <- dplyr::filter(parameters, parameter_name == "FOIadjust")

  expect_equal(foi_adjust$end_time[1], foi_adjust$start_time[2])
})
