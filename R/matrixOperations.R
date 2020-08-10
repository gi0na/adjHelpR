#' Auxiliary function, gives mask for matrix for directed, undirected etc.
#'
#' `x` can either be the matrix, the number of nodes or a vector of nodes. If x
#' is a list, it is assumed that the first element of the list is the vector of
#' row-nodes, and the second element the vector of column-nodes. If x is a
#' numeric vector of length 2, it is assumed that x is the dimensions of the
#' matrix.
#'
#' @param x either a matrix, the number of nodes, or the vector containing the
#'   nodes
#' @param directed  a boolean argument specifying whether object is directed or
#'   not. When absent it is inferred from `x`.
#' @param selfloops  a boolean argument specifying whether the model should
#'   incorporate selfloops. When absent it is inferred from `x`.
#' @param arr.ind boolean, return vector of indices or tuple with array indices?
#'   defaults to FALSE
#'
#' @return a boolean matrix that can be used to mask adjacency matrices.
#'
#' @export
#' @importFrom methods new
#' @examples
#' el <- data.frame(from= c('a','b','b','c','d','d'),
#'                 to  = c('b','c','d','a','b','a'),
#'                 attr= c( 12, 6, 12 , 6 , 6 , 6 ))
#' adj <- el2adj(el)
#' mat2vec_ix(adj, directed = TRUE, selfloops=FALSE)
#' mat2vec_ix(adj)
#' mat2vec_ix(x=nrow(adj), directed = TRUE, selfloops=FALSE)
#' mat2vec_ix(x=rownames(adj), directed = TRUE, selfloops=FALSE)
mat2vec_ix <- function(x=NULL, directed=NULL, selfloops=NULL,
                       arr.ind = FALSE) {
  # Returns the indices to
  # vectorise adjacency matrices
  # removing unused entries in the
  # case of undirected or no
  # selfloops graphs

  n <- nodes <- NULL
  if(!(is.matrix(x) | "Matrix"%in%attributes(class(x)))){
    if(is.list(x) | length(x)>2 | (length(x)==2 & !is.numeric(x))){
      nodes <- x
      x <- NULL
    } else{
      if(!is.list(x) & length(x)<=2 & is.numeric(x)){
        n <- x
        x <- NULL
      } else{
        x <- NULL
      }
    }
  }

  if(is.null(x) & is.null(n) & is.null(nodes))
    stop('Specify x in a compatible format.')
  if(any(c(is.null(directed), is.null(selfloops)))){
    if(is.null(x))
      stop('directed or selfloops were not specified and cannot be inferred.')
    out <- check_specs(x)
    directed <- out['directed']
    selfloops <- out['selfloops']
  }

  # get nrow and ncol. x has priority over n over nodes
  if(!is.null(x)){
    dims <- dim(x)
  } else{
    if(!is.null(n)){
      if(length(n)==2){
        dims <- c(n[1], n[2])
      } else{
        dims <- rep(n[1],2)
      }
    } else{
      if(!is.null(nodes)){
        if(!is.list(nodes)){
          dims <- rep(length(nodes),2)
        } else{
          if(is.list(nodes)){
            dims <- c(length(nodes[[1]]), length(nodes[[2]]))
          } else{
            stop('Wrong nodes dimensions.')
          }
        }
      }
    }
  }
  if(dims[1]!=dims[2]){
    directed <- TRUE
    selfloops <- TRUE
  }

  dat <- .mat2vec_ix(x, dims[1], dims[2], directed, selfloops)
  if(arr.ind)
    return(dat)
  y <- methods::new("ngTMatrix", i=as.integer(dat$row-1), j=as.integer(dat$col-1), Dim=as.integer(dims))
  return(Matrix::which(y))
}

.mat2vec_ix <- function(x, n_row, n_col, directed, selfloops){
  # Returns the indices to
  # vectorise adjacency matrices
  # removing unused entries in the
  # case of undirected or no
  # selfloops graphs

  dplyr::left_join(
    dplyr::tibble(row = 1:n_row, hash = 0),
    dplyr::tibble(col = 1:n_col, hash = 0), by='hash') %>%
    dplyr::select(row,col) -> ix

  if(isFALSE(directed))
    ix %>% dplyr::filter(col>=row) -> ix
  if(isFALSE(selfloops))
    ix %>% dplyr::filter(col!=row) -> ix

  return(ix)
}


#' Auxiliary function, produces matrix from vector
#'
#' The number of elements of vec are the number of non-zero elements in the
#' adjacency matrix.
#' It performs the opposite operation of `mat2vec_ix`.
#'
#'
#' @param vec  vector to be put in matrix form
#' @param directed  a boolean argument specifying whether object is directed or not.
#' @param selfloops  a boolean argument specifying whether the model should
#' incorporate selfloops.
#' @param n vector. if length(n)==1, n is the number of vertices. If length(n)==3
#' first element is number of vertices, second and third elements are number of
#' vertices for row and column of bipartite matrix.
#' @return
#' matrix nxn generated from vector.
#' @export
#'
#' @examples
#' el <- data.frame(from= c('a','b','b','c','d','d'),
#'                 to  = c('b','c','d','a','b','a'),
#'                 attr= c( 12, 6, 12 , 6 , 6 , 6 ))
#' adj <- el2adj(el)
#' ix <- mat2vec_ix(adj, FALSE, FALSE)
#' vec <- adj[ix]
#' vec2mat(vec, FALSE, FALSE, nrow(adj))
vec2mat <- function(vec, directed, selfloops, n) {
  if (length(n) > 1) {
    mat <- matrix(0, n[2], n[3])
  } else {
    mat <- matrix(0, n, n)
  }

  idx <- mat2vec_ix(mat, directed, selfloops)
  mat[idx] <- vec
  if (!directed) mat <- mat + t(mat)
  return(mat)
}
