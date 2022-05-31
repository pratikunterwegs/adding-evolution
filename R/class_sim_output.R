#' Checks the validity of the `simulation_output` class.
#'
#' @param object Instance of the `simulation_output` class.
#'
#' @return Boolean or errors
#' @export
#'
check_simulation_output <- function(object) {
  errors <- character()
  length_scenario <- length(object@scenario)
  if (length_scenario != 1) {
    msg <- sprintf("scenario is length %i, should be 1.", length_scenario)
    errors <- c(errors, msg)
  }

  if (length(errors) == 0) {
    TRUE
  } else {
    errors
  }
}

#' Defines the `simulation_output` class.
#'
#' @slot parameters list.
#' @slot trait_data list.
#' @slot edge_list data.frame.
#' @slot move_data list.
#'
#' @export
setClass(
  # name of the class
  Class = "simulation_output",

  # define the types of the class
  slots = c(
    parameters = "list",
    trait_data = "list",
    edge_list = "data.frame"
  ),

  # define the default values of the slots
  prototype = list(
    parameters = list(),
    trait_data = list(),
    edge_list = data.frame()
  ),

  # check validity of class
  validity = check_simulation_output
)
