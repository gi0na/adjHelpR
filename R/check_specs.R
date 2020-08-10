#' Check an adjacency matrix object for directedness and presence of selfloops
#'
#' This method allows different adjacency matrix classes to be checked.
#'
#' @param object matrix to be analyzed. Different classes are allowed.
#' @param ... additional parameters to and from internal functions, currently
#'   not used.
#'
#' @return named boolean vector of length 2 whose first element is TRUE for
#'   directed `objects` and whose second element is TRUE for `object` with
#'   selfloops.
#' @export
#'
#' @examples
#' el <- data.frame(from= c('a','b','b','c','d','d'),
#'                 to  = c('b','c','d','a','b','a'),
#'                 attr= c( 12, 6, 12 , 6 , 6 , 6 ))
#' adj <- el2adj(el)
#' check_specs(adj)
check_specs <- function(object, ...) {
  UseMethod("check_specs")
}

#' @rdname check_specs
#' @export
check_specs.matrix <- function(object, ...) {
  if (isSymmetric(object)) {
    directed <- FALSE
  } else {
    directed <- TRUE
  }

  if (all(diag(object) == 0)) {
    selfloops <- FALSE
  } else {
    selfloops <- TRUE
  }
  return(c("directed" = directed, "selfloops" = selfloops))
}
#' @rdname check_specs
#' @export
check_specs.dgTMatrix <- function(object, ...) {
  if (Matrix::isSymmetric(object)) {
    directed <- FALSE
  } else {
    directed <- TRUE
  }

  if (all(Matrix::diag(object) == 0)) {
    selfloops <- FALSE
  } else {
    selfloops <- TRUE
  }
  return(c("directed" = directed, "selfloops" = selfloops))
}

#' @rdname check_specs
#' @export
check_specs.dgCMatrix <- function(object, ...) {
  check_specs(methods::as(object, "dgTMatrix"))
}


#' @rdname check_specs
#' @export
check_specs.sparseMatrix <- function(object, ...) {
  check_specs(methods::as(object, "dgTMatrix"))
}

#' @rdname check_specs
#' @export
check_specs.default <- function(object, ...) {
  warning('Coercing x to matrix.')
  check_specs(methods::as(object, "matrix"))
}
