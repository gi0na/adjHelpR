#' Compute the adjacency matrix from various types of objects.
#'
#' If `x` is a tibble or a dataframe, `get_adjacency` assumes that it has been passed an edgelist.
#' If `x` is a matrix, `get_adjacency` assumes that it has been passed an incidence matrix and constructs the corresponding adjacency matrix.
#' If `x` is a different object, `get_adjacency` uses the appropriate method if it exists. Otherwise it tries to coerce `x` to a `data.frame`.
#'
#' @param x object from which generating the adjacency matrix
#' @inheritParams el2adj
#' @param ... additional parameters to and from the main method. Currently not
#'   used.
#'
#' @return A transition matrix in the same format of `adj`.
#' @export
#'
#' @examples
#' el <- data.frame(from= c('a','b','b','c','d','d'),
#'                 to  = c('b','c','d','a','b','a'),
#'                 attr= c( 12, 6, 12 , 6 , 6 , 6 ))
#' adj <- get_adjacency(el)
#'
get_adjacency <- function(x, select_cols = NULL, multiedge = NULL, aggr_expression = NULL, nodes = NULL, sparse = TRUE,
                          drop_names = FALSE, directed = NULL, selfloops = NULL, ...) {
  UseMethod("get_adjacency")
}

#' @rdname get_adjacency
#' @export
get_adjacency.tbl <- function(x, select_cols = NULL, multiedge = NULL, aggr_expression = NULL, nodes = NULL, sparse = TRUE,
                              drop_names = FALSE, directed = NULL, selfloops = NULL, ...){
  if(is.null(multiedge)) multiedge <- FALSE
  el2adj(el = x, select_cols = select_cols, multiedge = multiedge, aggr_expression = aggr_expression, nodes = nodes, sparse = sparse,
         drop_names = drop_names, directed = directed, selfloops = selfloops)
}

#' @rdname get_adjacency
#' @export
get_adjacency.data.frame <- function(x, select_cols = NULL, multiedge = NULL, aggr_expression = NULL, nodes = NULL, sparse = TRUE,
                                     drop_names = FALSE, directed = NULL, selfloops = NULL, ...){
  get_adjacency(dplyr::as.tbl(x, "tbl"), select_cols = select_cols, multiedge = multiedge, aggr_expression = aggr_expression, nodes = nodes, sparse = sparse,
                drop_names = drop_names, directed = directed, selfloops = selfloops, ...)
}

#' @rdname get_adjacency
#' @export
get_adjacency.matrix <- function(x, select_cols = NULL, multiedge = NULL, aggr_expression = NULL, nodes = NULL, sparse = TRUE,
                                 drop_names = FALSE, directed = NULL, selfloops = NULL, ...){

  y <- cbind(matrix(0, nrow = nrow(x), ncol=nrow(x)),x)
  colnames(y) <- c(rownames(x),colnames(x))
  y1 <- cbind(t(x),matrix(0, nrow = ncol(x), ncol=ncol(x)))
  z <- rbind(y,y1)
  if(isTRUE(sparse))
    z <- methods::as(z,"dsCMatrix")
  return(z)
}

#' @rdname get_adjacency
#' @export
get_adjacency.dgTMatrix <- function(x, select_cols = NULL, multiedge = NULL, aggr_expression = NULL, nodes = NULL, sparse = TRUE,
                                    drop_names = FALSE, directed = NULL, selfloops = NULL, ...){

  y <- cbind(Matrix::Matrix(0, nrow = nrow(x), ncol=nrow(x)),x)
  colnames(y) <- c(rownames(x),colnames(x))
  y1 <- cbind(Matrix::t(x),Matrix::Matrix(0, nrow = ncol(x), ncol=ncol(x)))
  z <- methods::as(rbind(y,y1),"symmetricMatrix")
  if(isFALSE(sparse))
    z <- as.matrix(z)
  return(z)
}

#' @rdname get_adjacency
#' @export
get_adjacency.dgCMatrix <- function(x, select_cols = NULL, multiedge = NULL, aggr_expression = NULL, nodes = NULL, sparse = TRUE,
                                    drop_names = FALSE, directed = NULL, selfloops = NULL, ...){
  get_adjacency(methods::as(x, "dgTMatrix"), select_cols = select_cols, multiedge = multiedge, aggr_expression = aggr_expression, nodes = nodes, sparse = sparse,
                drop_names = drop_names, directed = directed, selfloops = selfloops, ...)
}


#' @rdname get_adjacency
#' @export
get_adjacency.sparseMatrix <- function(x, select_cols = NULL, multiedge = NULL, aggr_expression = NULL, nodes = NULL, sparse = TRUE,
                                       drop_names = FALSE, directed = NULL, selfloops = NULL, ...){
  get_adjacency(methods::as(x, "dgTMatrix"), select_cols = select_cols, multiedge = multiedge, aggr_expression = aggr_expression, nodes = nodes, sparse = sparse,
                drop_names = drop_names, directed = directed, selfloops = selfloops, ...)
}

#' @rdname get_adjacency
#' @export
get_adjacency.default <- function(x, select_cols = NULL, multiedge = NULL, aggr_expression = NULL, nodes = NULL, sparse = TRUE,
                                     drop_names = FALSE, directed = NULL, selfloops = NULL, ...){
  get_adjacency(methods::as(x, "tbl"), select_cols = select_cols, multiedge = multiedge, aggr_expression = aggr_expression, nodes = nodes, sparse = sparse,
                drop_names = drop_names, directed = directed, selfloops = selfloops, ...)
}

#' @rdname get_adjacency
#' @export
get_adjacency.igraph <- function(x, select_cols = NULL, multiedge = NULL, aggr_expression = NULL, nodes = NULL, sparse = TRUE,
                                  drop_names = FALSE, directed = NULL, selfloops = NULL, ...){
  if(requireNamespace("igraph", quietly = TRUE) && igraph::is.igraph(x)){
    return(igraph::as_adjacency_matrix(x, sparse = sparse, names = !drop_names, attr = aggr_expression))
  }
}
