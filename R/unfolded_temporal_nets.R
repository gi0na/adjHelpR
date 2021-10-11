# adapt ifelse to return complex objects
ifelse_v <- function(test, yes, no){
  ifelse(test, list(yes), list(no))[[1]]
}

# create unfolded node names pasting names and timestamps together
.get_unfolded_nodes_vertical <- function(node.labels, ts, flat = FALSE){
  paste(rep(node.labels,max(ts)-min(ts)+1),rep((as.numeric(flat)+min(ts)-1):(max(ts)), each=length(node.labels)), sep = '_')
}
.get_unfolded_nodes_horizontal <- function(node.labels, ts, flat = FALSE){
  nodes <- paste(rep(node.labels,each=max(ts)+1-min(ts)+!flat),(as.numeric(flat)+min(ts)-1):(max(ts)), sep = "_")
}

# internal function to return unfolded adjacency
.get_unfolded_adjacency_internal <- function(edge.list, sparse, node.labels, ts, vertical, flat = FALSE){
  if(vertical){
    nodes <- .get_unfolded_nodes_vertical(node.labels, ts, flat)
  } else{
    nodes <- .get_unfolded_nodes_horizontal(node.labels, ts, flat)
  }
  get_adjacency(edge.list, nodes = nodes, sparse = sparse, multiedge = TRUE, select_cols = 1:2)
}

# internal function construct unfolded adjacency from original edge.list
.get_unfolded_adjacency <- function(edge.list, node.labels, vertical, sparse, flat = FALSE){
  # construct unfolded edge list
  unfolded.edge.list <- as_tibble(
    t(apply(edge.list, 1, function(row)
      paste(row[1:2], c(-1,-flat)+as.integer(row[3]) + flat, sep = '_'))), .name_repair = ~c('source', 'target'))

  return(
    # construct unfolded adjacency matrix from unfolded edgelist
    .get_unfolded_adjacency_internal(edge.list = unfolded.edge.list, sparse, node.labels, ts = edge.list$timestamp, vertical = vertical, flat)
  )
}

# rename edgelist for practicality using .select_cols
.rename_edgelist <- function(edge.list, select_cols){
  col_names <- colnames(edge.list)
  col_names <- .select_cols_temporal(col_names = col_names,
                                     select_cols = select_cols, attr_cols = FALSE)
  colnames(edge.list) <- col_names
  edge.list %>%
    select(.data$source, .data$target, .data$timestamp) %>%
    mutate(pos = .data$timestamp - min(.data$timestamp) + 1)
}

#' Create time unfolded adjacency matrix
#' @inheritParams get_adjacency
#' @param edge.list data.frame or tibble containing the edge list. It needs at
#'   least three column: the column with edge sources, the edge targets, and the
#'   timestamps of each edge. The order of columns should be 'timestamp',
#'   'source', 'target'. If the edge.list columns come in different orders, use
#'   `select_cols` to specify the right order. See the example for details.
#' @param vertical unfold vertically or horizontally? Defaults to FALSE
#'   (horizontal)
#' @param ... extra parameters passed to internal methods
#'
#' @return The time unfolded (nodes x timestamps) x (nodes x timestamps)
#' adjacency matrix. The node order follows the vertical or horizontal algnment
#'
#' @export
#'
#' @examples
#' el <- data.frame(
#' from = c('A','B',  'A','B','B',   'A','C','C',    'A','B','C',    'D'),
#' to   = c('C','C',  'C','C','D',   'C','D','E',    'C','C','E',    'E'),
#' ts   = c( 1,  1,    2,  2,  2,     3,  3,  3,      4,  4,  4,      6)
#' )
#' get_unfolded_adjacency(el, select_cols = 1:3)
get_unfolded_adjacency <- function(edge.list, select_cols = NULL, nodes = NULL, vertical = FALSE, sparse = TRUE, ...){
  # rename and reorganise edge list
  edge.list <- .rename_edgelist(edge.list, select_cols)
  # get node labels
  node.labels <- ifelse_v(is.null(nodes), nodes_from_el(edge.list, select_cols = 1:2), nodes)
  # construct unfolded adjacency matrix
  .get_unfolded_adjacency(edge.list = edge.list, node.labels = node.labels, vertical = vertical, sparse = sparse, ...)
}

# internal function to construct unfolded ggraph layout
.get_unfolded_layout <- function(edge.list, node.labels, vertical, resize_ratio, enlarge_ratio, ...){
  # construct unfolded adjacency matrix
  tadj <- .get_unfolded_adjacency(edge.list = edge.list, node.labels = node.labels, vertical = vertical, sparse = FALSE, ...)

  # set default resize ratio if not passed
  # different defaults for vertical or horizontal layouts
  if(is.null(resize_ratio))
    resize_ratio <- ifelse(vertical, .6, 2.5)

  # set height and width of grid according to number of nodes and number of timestamps
  # horizontal and vertical are swapped
  width <- ifelse(vertical, length(node.labels), max(edge.list$pos))
  height <- ifelse(vertical, max(edge.list$pos), length(node.labels))

  # generate grid layout according to the dimensions defined above
  lyt <- ggraph::create_layout(graph = tadj, layout = 'on_grid', height=height, width=width)
  lyt$y <- max(lyt$y)-lyt$y
  lyt$x <- resize_ratio*lyt$x

  # enlarge by enlarge_ratio
  lyt$y <- lyt$y*enlarge_ratio
  lyt$x <- lyt$x*enlarge_ratio

  # return layout
  return(lyt)
}

#' Generate the ggraph layout of a time unfolded graph
#'
#' @inheritParams get_unfolded_adjacency
#' @param resize_ratio ratio between horizontal and vertical dimensions of the grid layout.
#' value < 1 gives a longer vertical side, >1 longer horizontal side.
#' @param enlarge_ratio enlarge both y and x axis by this parameter. Defaults is 1.
#'
#' @return ggraph layout object to plot the time unfolded network
#' @export
#' @import ggraph
#' @examples
#' el <- data.frame(
#' from = c('A','B',  'A','B','B',   'A','C','C',    'A','B','C',    'D'),
#' to   = c('C','C',  'C','C','D',   'C','D','E',    'C','C','E',    'E'),
#' ts   = c( 1,  1,    2,  2,  2,     3,  3,  3,      4,  4,  4,      6)
#' )
#' get_unfolded_adjacency(el, select_cols = 1:3)
get_unfolded_layout <- function(edge.list, select_cols = NULL, nodes = NULL, vertical = FALSE, resize_ratio = NULL, enlarge_ratio = 1, ...){
  # rename and reorganize edge.list
  edge.list <- .rename_edgelist(edge.list, select_cols)
  # get node labels
  node.labels <- ifelse_v(is.null(nodes), nodes_from_el(edge.list, select_cols = 1:2), nodes)

  # generate the ggraph layout
  .get_unfolded_layout(edge.list, node.labels, vertical, resize_ratio, enlarge_ratio, ...)

}

#' Plot the time unfolded graph based on an edge list
#'
#' The function generates a ggraph plot with some sensible defaults. For full
#' customisation, it is recommended to generated the ggraph layout from
#' `get_unfolded_layout` and then customise the plot.
#'
#' @inheritParams get_unfolded_layout
#'
#' @return ggraph plot object with directed edges whose direction follows the
#'   arrow of time, aligned nodes in the order passed originally, node labels
#'   along the nodes and timestamps labels along the direction of time.
#' @export
#'
#' @examples
#' el <- data.frame(
#' from = c('A','B',  'A','B','B',   'A','C','C',    'A','B','C',    'D'),
#' to   = c('C','C',  'C','C','D',   'C','D','E',    'C','C','E',    'E'),
#' ts   = c( 1,  1,    2,  2,  2,     3,  3,  3,      4,  4,  4,      6)
#' )
#' get_unfolded_plot(el, select_cols = 1:3)
get_unfolded_plot <- function(edge.list, select_cols = NULL, nodes = NULL, vertical = FALSE, resize_ratio = NULL, enlarge_ratio = 1, ...){

  edge.list <- .rename_edgelist(edge.list, select_cols)
  edge.list$pos <- edge.list$pos + 1
  node.labels <- ifelse_v(is.null(nodes), nodes_from_el(edge.list, select_cols = 1:2), nodes)

  lyt <- .get_unfolded_layout(edge.list, node.labels, vertical, resize_ratio, enlarge_ratio, ...)
  lyt %>%
    ggraph::ggraph() +
    ggraph::geom_edge_fan(width = grid::unit(0.5, 'mm'),
                  arrow = ggplot2::arrow(length = grid::unit(2.5, 'mm'), type = 'closed'),
                  end_cap = ggraph::circle(2, 'mm'),
                  start_cap = ggraph::circle(2, 'mm'),
                  color='gray50') +
    ggraph::geom_node_point(size=3) +
    # times
    ggraph::geom_node_text(ggplot2::aes(
      label=ifelse(gsub('_.*', '', .data$name)==node.labels[1],
                   gsub(paste0(node.labels[1],'_'), '', .data$name),
                   NA)),
      nudge_y = (!vertical)*0.6,
      nudge_x = - vertical*0.6,
      size=5
    ) +
    # names
    ggraph::geom_node_text(ggplot2::aes(
      label=ifelse(gsub('.*_', '', .data$name)==(min(edge.list$timestamp)-1), gsub('_.*', '', .data$name), NA)),
      nudge_x = - (!vertical)*0.6,
      nudge_y = vertical*0.6,
      size=5
    ) +
    ggplot2::coord_fixed(clip='off') + ggraph::theme_graph() +
    ggplot2::theme(legend.position = "none",
          panel.background = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.spacing = grid::unit(0,"null"),
          plot.margin = rep(grid::unit(0,"null"),4),
          axis.ticks = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.ticks.length = grid::unit(0,"null"))
}

#' Plot the flat unfolded graph based on an edge list
#'
#' The function generates a ggraph plot with some sensible defaults. For full
#' customisation, it is recommended to generated the ggraph layout from
#' `get_unfolded_layout` and then customise the plot.
#'
#' @inheritParams get_unfolded_layout
#'
#' @return ggraph plot object with directed edges whose direction follows the
#'   arrow of time, aligned nodes in the order passed originally, node labels
#'   along the nodes and timestamps labels along the direction of time.
#' @export
#'
#' @examples
#' el <- data.frame(
#' from = c('A','B',  'A','B','B',   'A','C','C',    'A','B','C',    'D'),
#' to   = c('C','C',  'C','C','D',   'C','D','E',    'C','C','E',    'E'),
#' ts   = c( 1,  1,    2,  2,  2,     3,  3,  3,      4,  4,  4,      6)
#' )
#' get_flat_unfolded_plot(el, select_cols = 1:3)
get_flat_unfolded_plot <- function(edge.list, select_cols = NULL, nodes = NULL, vertical = FALSE, resize_ratio = NULL, enlarge_ratio = 1, ...){

  edge.list <- .rename_edgelist(edge.list, select_cols)
  node.labels <- ifelse_v(is.null(nodes), nodes_from_el(edge.list, select_cols = 1:2), nodes)

  lyt <- .get_unfolded_layout(edge.list, node.labels, vertical, resize_ratio, enlarge_ratio, flat = TRUE)
  lyt %>%
    ggraph::ggraph() +
    ggraph::geom_edge_arc(width = grid::unit(0.5, 'mm'),
                          end_cap = ggraph::circle(2, 'mm'),
                          start_cap = ggraph::circle(2, 'mm'),
                          color='gray50') +
    ggraph::geom_node_point(size=3) +
    # times
    ggraph::geom_node_text(ggplot2::aes(
      label=ifelse(gsub('_.*', '', .data$name)==node.labels[1],
                   gsub(paste0(node.labels[1],'_'), '', .data$name),
                   NA)),
      nudge_y = (!vertical)*0.6,
      nudge_x = - vertical*0.6,
      size=5
    ) +
    # names
    ggraph::geom_node_text(ggplot2::aes(
      label=ifelse(gsub('.*_', '', .data$name)==(min(edge.list$timestamp)), gsub('_.*', '', .data$name), NA)),
      nudge_x = - (!vertical)*0.6,
      nudge_y = vertical*0.6,
      size=5
    ) +
    ggplot2::coord_fixed(clip='off') + ggraph::theme_graph() +
    ggplot2::theme(legend.position = "none",
                   panel.background = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.spacing = grid::unit(0,"null"),
                   plot.margin = rep(grid::unit(0,"null"),4),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.ticks.length = grid::unit(0,"null"))
}
