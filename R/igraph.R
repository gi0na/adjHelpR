# auxiliary function for to extract igraph properties
check_specs.igraph <- function(object, ...) {
  if (requireNamespace("igraph", quietly = TRUE) && igraph::is.igraph(object)) {
    if (is.null(directed)) {
      if (igraph::is.directed(object)) {
        directed <- FALSE
      } else {
        directed <- TRUE
      }
    }

    if (is.null(selfloops)) {
      if (igraph::is.simple(igraph::simplify(object, remove.multiple = TRUE, remove.loops = FALSE))) {
        selfloops <- FALSE
      } else {
        selfloops <- TRUE
      }
    }
  }
  return(c("directed" = directed, "selfloops" = selfloops))
}

#' Convert a list of adjacency matrices to a list of igraph graphs.
#'
#' @param adjlist a list of adjacency matrices
#' @param directed a boolean argument specifying whether object is directed or not.
#' @param selfloops a boolean argument specifying whether the model should incorporate selfloops.
#' @param weighted boolean, generate weighted graphs?
#'
#' @return
#'
#' list of igraph graphs.
#'
#' @export
#'
#' @examples
#' data("adj_karate")
#' adj_list <- list(adj_karate)
#' glist <- CreateIgGraphs(adj_list, FALSE, FALSE)
CreateIgGraphs <- function(adjlist, directed, selfloops, weighted = NULL) {
  if (directed) {
    mode <- "directed"
  }
  if (!directed) {
    mode <- "undirected"
  }

  lapply(X = adjlist, FUN = igraph::graph_from_adjacency_matrix, mode = mode, diag = selfloops, weighted = weighted)
}
