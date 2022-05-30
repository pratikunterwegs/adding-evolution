
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

#' Get trait data.
#'
#' @param object A `simulation_output` object.
#'
#' @return A `data.table` object.
#' @export
#'
get_trait_data = function(object) {
    trait_data = object@trait_data
    trait_data = Map(
        trait_data, seq(length(trait_data)),
        f = function(df, g) {
            df$gen = g
            df
        }
    )
    trait_data = data.table::rbindlist(trait_data)
    trait_data
}
