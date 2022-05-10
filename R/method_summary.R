# methods to summarise and show simulation_output class objects

#' Summarise a `simulation_output` object.
#'
#' @name summary.simulation_output
#' @docType methods
#' @rdname summary
#'
setMethod(
  "summary", signature(object = "simulation_output"),
  function(object) {
    print(
      glue::glue(
        "Scenario: {object@parameters$scenario}
        Popsize: {object@parameters$pop_size}
        "
      )
    )
  }
)


#' Show the structure of a `simulation_output` object.
#'
#' @name str.simulation_output
#' @docType methods
#' @rdname str
#'
setMethod(
  "str", signature(object = "simulation_output"),
  function(object) {
    utils::str(
      object,
      max.level = 3
    )
  }
)
