
#' Get network data from `simulation_output` objects.
#'
#' @param object A `simulation_output` object.
#' @param weight_threshold Weight below which edges are excluded.
#'
#' @return A `tidygraph` objects.
#' @export
make_network = function(object, weight_threshold) {
    assertthat::assert_that(
        class(object) == "simulation_output",
        msg = "pmv get traits: object is not simulation output"
    )

    # handle edges
    edgelist = object@edge_list
    colnames(edgelist) = c("from", "to", "weight")
    edgelist = data.table::as.data.table(edgelist)
    
    edgelist = edgelist[edgelist$from != edgelist$to, ]

    edgelist = edgelist[edgelist$weight >= weight_threshold, ]

    edgelist$to = edgelist$to + 1
    edgelist$from = edgelist$from + 1

    # handle nodes
    nodes = object@trait_data
    data.table::setDT(nodes)
    nodes$id = seq(nrow(nodes))

    # make tidygraph objects
    g = tidygraph::tbl_graph(
        nodes = nodes,
        edges = edgelist,
        directed = FALSE
    )

    g
}

#' Plot networks from data from `simulation_output` objects.
#'
#' @param g A `tidygraph` object.
#' @param colour_by A variable name, _not_ as a string.
#' @return A `ggplot` object.
#' @export
plot_network = function(g, colour_by) {

    colour_by = rlang::enquo(colour_by)

    ggraph::ggraph(
        g, x = xn, y = yn
    ) +
    ggraph::geom_edge_fan(
        edge_width = 0.5,
        edge_colour = "grey70",
        ggplot2::aes(
            edge_alpha = weight
        ),
        show.legend = FALSE
    )+
    ggraph::geom_node_point(
        ggplot2::aes(
          fill = !!colour_by,
          size = assoc
        ),
        shape = 21
    )+
    ggplot2::scale_size_continuous(
        range = c(0.5, 5)
    )+
    ggraph::scale_edge_alpha(
        range = c(0.3, 1)
    )+
    ggplot2::theme_test()
}
