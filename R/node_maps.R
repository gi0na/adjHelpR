#' Extract node names from an edgelist
#'
#' This function allows to extract the node names for an edgelist. When the
#' parameter select_cols is absent, existing source and target columns are used.
#' If no columns named source or target are present, and the parameter is not
#' specified, the first and second column are assumed to be source and target.
#'
#' @param edge.list edgelist, either in matrix or dataframe/tibble format
#' @param select_cols optional vector specifying source and target columns of
#'   the edgelist. When absent, existing source and target columns are used. If
#'   no columns named source or target are present, and the parameter is not
#'   specified, the first and second column are assumed to be source and target.
#'   If `select_cols` is a numeric vector, its first two elements are assumed to
#'   be the source and target columns indices. If `select_cols` is a character
#'   vector, its first two elements are assumed to be the names of the source
#'   and target columns.
#'
#' @return vector containing unique node names
#' @export
#'
#' @examples
#' el <- data.frame(from= c('a','b','b','c','d','d'),
#'                 to  = c('b','c','d','a','b','a'),
#'                 attr= c( 12, 6, 12 , 6 , 6 , 6 ))
#'
#' nodes_from_el(el, select_cols = c('from','to'))
#' nodes_from_el(el, select_cols = 1:2)
#' nodes_from_el(el)
#'
nodes_from_el <- function(edge.list, select_cols = NULL){
  col_names <- colnames(edge.list)
  colnames(edge.list) <- .select_cols(col_names, select_cols, FALSE)
  return(unique(c(edge.list$source, edge.list$target)))
}


# return the columns corresponding to source, target and attr
.select_cols <- function(col_names, select_cols, attr_cols){
  # check sizes of objects
  if ((!is.null(select_cols)) & length(select_cols) < 2) stop("select_cols have the wrong size")


  # if select_cols is not specified and source and targets do not exist set them to col 1 and 2
  if(is.null(select_cols)){
    col_names <- .map_colnames(col_names)
    if(!all(c('source','target') %in% col_names)){
      select_cols <- 1:2
      if(isTRUE(attr_cols))
        select_cols <- c(select_cols,3)
    }
  }

  # if select_cols was specified and
  if(!is.null(select_cols)){
    if(is.character(select_cols))
      select_cols <- match(select_cols, col_names)
    col_names[col_names %in% c('source','target')] <- 1:sum(col_names %in% c('source','target'))
    col_names[select_cols[1:2]] <- c('source','target')
    if(attr_cols & length(select_cols)>2){
      col_names[col_names %in% 'attr'] <- 3
      col_names[select_cols[3]] <- 'attr'
    }
  }
  return(col_names)
}

# return the columns corresponding to source, target and attr
.select_cols_temporal <- function(col_names, select_cols, attr_cols){
  # check sizes of objects
  if ((!is.null(select_cols)) & length(select_cols) < 3) stop("select_cols have the wrong size")


  # if select_cols is not specified and source and targets do not exist set them to col 1 and 2
  if(is.null(select_cols)){
    col_names <- .map_colnames(col_names)
    if(!all(c('timestamp','source','target') %in% col_names)){
      select_cols <- 1:3
      if(isTRUE(attr_cols))
        select_cols <- c(select_cols,4)
    }
  }

  # if select_cols was specified and
  if(!is.null(select_cols)){
    if(is.character(select_cols))
      select_cols <- match(select_cols, col_names)
    col_names[col_names %in% c('timestamp','source','target')] <- 1:sum(col_names %in% c('timestamp','source','target'))
    col_names[select_cols[1:3]] <- c('source','target','timestamp')
    if(attr_cols & length(select_cols)>3){
      col_names[select_cols[4]] <- 'attr'
    }
  }
  return(col_names)
}

.map_colnames <- function(col_names){
  col_names[!is.na(match(col_names, c('time','t')))] <- 'timestamp'
  col_names[!is.na(match(col_names, c('from','i')))] <- 'source'
  col_names[!is.na(match(col_names, c('to','j')))] <- 'target'
  return(col_names)
}



