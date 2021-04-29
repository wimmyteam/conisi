
#' Loads simulation parameters from a file
#'
#' This function reads parameters from a file and reformats them so they are ready to be used in the model
#'
#' @param file_path String The path to the file to be loaded
#' @param sim_name String Optional The name of the simulation to be extracted from the file
#' @param start_date Date Must be provided if the file contains start_date instead of start_time.
#'   It is the date that will be considered time zero.
#'
#'@export

load_parameters <- function(file_path, sim_name = NULL, start_date = NULL) {
  data <- readr::read_delim(file=file_path,
                            delim=";",
                            col_types=readr::cols()
                            ) %>%
    dplyr::mutate(start_date = lubridate::dmy(start_date))



  if (!is.null(sim_name)){
    if(!"sim_name" %in% colnames(data))
      {stop("Trying to filter on sim_name but parameter file does not have sim_name column")}
    data <- data %>%
      dplyr::filter(sim_name == {{sim_name}})
  }

  if (!is.null(start_date)){
    data <- data %>%
    dplyr::mutate(start_time = ifelse(is.na(start_date), 0, start_date - {{start_date}}))
  } else {
    stopifnot("start_time" %in% colnames(data))
  }

  data <- data %>%
    dplyr::arrange(sim_name, parameter_name, start_time) %>%
    dplyr::group_by(parameter_name) %>%
    dplyr::mutate(end_time = dplyr::lead(start_time)) %>%
    dplyr::ungroup()

  return(data)
}
