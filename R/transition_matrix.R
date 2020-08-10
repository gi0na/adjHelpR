#' Compute a transition matrix from an adjacency matrix.
#'
#' This method take an adjacency matrix (dense or sparse), and returns the
#' corresponding transition matrix in the same format. If the optional argument
#' `directed` is specified to be FALSE, then before computing the transition
#' matrix the adjacency matrix is forced be symmetric.
#'
#' @param adj adjacency matrix, either dense or sparse
#' @param directed boolean, optional argument specifying if the matrix refers to
#'   a directed graph.
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
#' adj <- el2adj(el)
#' t_mat <- get_transition(adj)
get_transition <- function(adj, directed=NULL, ...) {
  UseMethod("get_transition")
}

#' @rdname get_transition
#' @export
get_transition.matrix <- function(adj, directed=NULL, ...){
  if(isFALSE(directed) && !isSymmetric(adj))
    adj <- adj + t(adj)
  t_mat <- adj/rowSums(adj)
  return(t_mat)
}

#' @rdname get_transition
#' @export
get_transition.dgTMatrix <- function(adj, directed=NULL, ...){
  if(isFALSE(directed) && !Matrix::isSymmetric(adj))
    adj <- adj + Matrix::t(adj)
  t_mat <- adj/Matrix::rowSums(adj)
  return(t_mat)
}

#' @rdname get_transition
#' @export
get_transition.dgCMatrix <- function(adj, directed=NULL, ...){
  get_transition(methods::as(adj, "dgTMatrix"), directed=directed, ...)
}


#' @rdname get_transition
#' @export
get_transition.sparseMatrix <- function(adj, directed=NULL, ...){
  get_transition(methods::as(adj, "dgTMatrix"), directed=directed, ...)
}

#' @rdname get_transition
#' @export
get_transition.default <- function(adj, directed=NULL, ...){
  warning('Coercing adj to matrix')
  get_transition(methods::as(adj, "matrix"), directed=directed, ...)
}
