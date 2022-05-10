
#' Get movement data.
#'
#' @param object A `simulation_output` object.
#'
#' @return A `data.table` object.
#' @export
#'
get_move_data = function(object) {
    move_data = object@move_data

    move_data = data.table::rbindlist(move_data)

    move_data
}
