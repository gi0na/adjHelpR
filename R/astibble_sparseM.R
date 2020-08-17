#' Method to convert 'Matrix' sparse matrices to 'tibble's
#' Method for 'as_tibble' to coerce sparse matrix objects
#' to a 'tibble'
#'
#' @param x A Matrix object
#' @param ... Extra arguments, not used
#'
#' @export
as_tibble.dgTMatrix <- function(x, ...) {
  s <- Matrix::summary(x)

  row <- s$i
  if (!is.null(rownames(x))) {
    row <- rownames(x)[row]
  }
  col <- s$j
  if (!is.null(colnames(x))) {
    col <- colnames(x)[col]
  }

  return(
    tibble(row = row, column = col, value = s$x)
  )
}

#' @rdname as_tibble.dgTMatrix
#' @export
as_tibble.dgCMatrix <- function(x, ...) {
  as_tibble(methods::as(x, "dgTMatrix"))
}


#' @rdname as_tibble.dgTMatrix
#' @export
as_tibble.sparseMatrix <- function(x, ...) {
  as_tibble(methods::as(x, "dgTMatrix"))
}
