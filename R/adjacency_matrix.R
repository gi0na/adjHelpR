#' Compute the adjacency matrix from various types of objects.
#'
#' If `x` is a tibble or a dataframe, `get_adjacency` assumes that it has been
#' passed an edgelist. If `x` is a matrix, `get_adjacency` assumes that it has
#' been passed an incidence matrix and constructs the corresponding adjacency
#' matrix. If `x` is a different object, `get_adjacency` uses the appropriate
#' method if it exists. Otherwise it tries to coerce `x` to a `data.frame`.
#'
#' @param x object from which generating the adjacency matrix
#' @inheritParams el2adj
#' @param edgelist optional boolean. If `x` is a matrix, should treat `x` as an
#'   edgelist? If not provided, tries to infer if `x` is an edgelist or an
#'   incidence matrix from its shape.
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
                          drop_names = FALSE, directed = NULL, selfloops = NULL, edgelist = NULL, ...) {
  UseMethod("get_adjacency")
}

#' @rdname get_adjacency
#' @export
get_adjacency.tbl <- function(x, select_cols = NULL, multiedge = NULL, aggr_expression = NULL, nodes = NULL, sparse = TRUE,
                              drop_names = FALSE, directed = NULL, selfloops = NULL, edgelist = NULL, ...){
  if(is.null(multiedge)) multiedge <- FALSE
  el2adj(edge.list = x, select_cols = select_cols, multiedge = multiedge, aggr_expression = aggr_expression, nodes = nodes, sparse = sparse,
         drop_names = drop_names, directed = directed, selfloops = selfloops)
}

#' @rdname get_adjacency
#' @export
get_adjacency.tbl_graph <- function(x, select_cols = NULL, multiedge = NULL, aggr_expression = NULL, nodes = NULL, sparse = TRUE,
                              drop_names = FALSE, directed = NULL, selfloops = NULL, edgelist = NULL, ...){
  if(is.null(multiedge)) multiedge <- FALSE
  x %>% tidygraph::activate('edges') %>% dplyr::as_tibble() %>%
    el2adj(select_cols = select_cols, multiedge = multiedge, aggr_expression = aggr_expression, nodes = nodes, sparse = sparse,
         drop_names = drop_names, directed = directed, selfloops = selfloops)
}

#' @rdname get_adjacency
#' @export
get_adjacency.data.frame <- function(x, select_cols = NULL, multiedge = NULL, aggr_expression = NULL, nodes = NULL, sparse = TRUE,
                                     drop_names = FALSE, directed = NULL, selfloops = NULL, edgelist = NULL, ...){
  get_adjacency(dplyr::as_tibble(x), select_cols = select_cols, multiedge = multiedge, aggr_expression = aggr_expression, nodes = nodes, sparse = sparse,
                drop_names = drop_names, directed = directed, selfloops = selfloops, ...)
}

#' @rdname get_adjacency
#' @export
get_adjacency.matrix <- function(x, select_cols = NULL, multiedge = NULL, aggr_expression = NULL, nodes = NULL, sparse = TRUE,
                                 drop_names = FALSE, directed = NULL, selfloops = NULL, edgelist = NULL, ...){

  if(is.null(edgelist)){
    if(is.null(c(select_cols,multiedge,aggr_expression,nodes))){
      if(ncol(x)==3 & ncol(x)<nrow(x)){
        warning('Treating x as an edgelist. If x was an incidence matrix, set edgelist=FALSE.')
        edgelist <- TRUE
      } else{
        edgelist <- FALSE
      }
    } else{
      message('Treating x as an edgelist.')
      edgelist <- TRUE
    }
  }

  if(isTRUE(edgelist)){
    z <- get_adjacency(as.data.frame(x), select_cols = select_cols, multiedge = multiedge, aggr_expression = aggr_expression, nodes = nodes, sparse = sparse,
                       drop_names = drop_names, directed = directed, selfloops = selfloops, ...)
    return(z)
  }
  if(isFALSE(edgelist)){
    warning('Treating x as an incidence matrix. If x was an edgelist matrix, set edgelist=TRUE')
    if(is.null(rownames(x))){
      rownames(x) <- 1:nrow(x)
    }
    if(is.null(colnames(x))){
      colnames(x) <- (nrow(x)+1):(nrow(x)+ncol(x))
    }
    y <- cbind(matrix(0, nrow = nrow(x), ncol=nrow(x)),x)
    colnames(y) <- c(rownames(x),colnames(x))
    y1 <- cbind(t(x),matrix(0, nrow = ncol(x), ncol=ncol(x)))
    z <- rbind(y,y1)
    if(isTRUE(sparse))
      z <- methods::as(z,"dsCMatrix")
    return(z)
  }
  stop('Wrong Format for x')
}

#' @rdname get_adjacency
#' @export
get_adjacency.dgTMatrix <- function(x, select_cols = NULL, multiedge = NULL, aggr_expression = NULL, nodes = NULL, sparse = TRUE,
                                    drop_names = FALSE, directed = NULL, selfloops = NULL, ...){

  warning('Treating x as an incidence matrix.')
  if(is.null(rownames(x))){
    rownames(x) <- 1:nrow(x)
  }
  if(is.null(colnames(x))){
    colnames(x) <- 1:ncol(x)
  }
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
  get_adjacency(dplyr::as_tibble(x), select_cols = select_cols, multiedge = multiedge, aggr_expression = aggr_expression, nodes = nodes, sparse = sparse,
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
