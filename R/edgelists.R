#' Maps adjacency matrix to edgelist
#'
#' @param adj matrix or sparse matrix object, the adjacency matrix
#' @param directed boolean, is the graph directed?
#' @param names optional vector with node names. When `NULL` row- and column
#'   names of `adj` are used.
#'
#' @return a dataframe containing the edgelist
#'
#' @examples
#' el <- data.frame(from= c('a','b','b','c','d','d'),
#'                 to  = c('b','c','d','a','b','a'),
#'                 attr= c( 12, 6, 12 , 6 , 6 , 6 ))
#' adj <- el2adj(el)
#' el1 <- adj2el(adj, directed = TRUE)
#' @export
#' @import dplyr
#' @importFrom rlang .data
#'
adj2el <- function(adj, directed = TRUE, names = NULL) {
  if((!is.null(names)) && length(names) != nrow(adj))
    stop('length(names) is different from nrow(adj)')
  if(!is.null(names))
    rownames(adj) <- colnames(adj) <- names
  if (isFALSE(directed)) adj[lower.tri(adj, F)] <- 0
  if (is.matrix(adj)) {
    reshape2::melt(adj, value.name = "value") %>%
      as_tibble() %>%
      filter(.data$value != 0) %>%
      rename(source = "Var1", target = "Var2", attr = "value") %>%
      mutate_at(c("source", "target"), as.character) -> el
  } else {
    if (attributes(class(adj))$package == "Matrix") {
      as_tibble(adj) %>%
        rename(source = "row", target = "column", attr = "value") -> el
    } else {
      stop("Wrong adjacency format")
    }
  }
  return(el)
}


#' Maps edgelist to adjacency matrix
#'
#' @param edge.list dataframe containing a (weighted) edgelist.
#' @param nodes optional vector containing all node names in case disconnected
#'   nodes should be included.
#' @param select_cols optional vector of 3 (2 for multi-graphs) elements
#'   specifying which columns are the source,target, and attributes from which
#'   building the graph. Otherwise Column 1 is assumed to be the source, column
#'   2 the target, column 3 the attribute. In the case of multi-graphs, the
#'   third element is not needed and the number of edges between each pair of
#'   vertices is computed according to '\code{aggr_expression}'.
#' @param multiedge boolean, are there multiple edges? defaults to FALSE.
#' @param aggr_expression string, the expression used to compute the aggregated
#'   value in the adjacency matrix in the presence of multiedges. It defaults to
#'   '\code{dplyr::n()}'. If `edge.list` contains other columns such as `weight`,
#'   examples of `aggr_expression` could be `"mean(weight)"` to report the
#'   average of the `weight` column of the multiedges, `"min(weight)"` to report
#'   the minimum of the `weight` column of the multiedges.
#' @param sparse boolean, return sparse matrix? default to TRUE
#' @param drop_names boolean, drop names from matrix and return a vector of
#'   names together with the matrix. This option saves considerable memory when
#'   dealing with graphs with long node names.
#' @param directed boolean, optional parameter, if FALSE, forces to undirected
#'   upper triangular adjacency matrix
#' @param selfloops boolean, optional parameter, if ignores selfloops.
#'
#' @return the (weighted) adjacency matrix corresponding the edgelist passed. If
#'   `drop_names` is TRUE, returns a list with the adjacency matrix and a vector
#'   with node names.
#' @export
#'
#' @examples
#' el <- data.frame(from= c('a','b','b','c','d','d'),
#'                 to  = c('b','c','d','a','b','a'),
#'                 attr= c( 12, 6, 12 , 6 , 6 , 6 ))
#' adj <- el2adj(el)
#'
el2adj <- function(edge.list, select_cols = NULL, multiedge = FALSE, aggr_expression = NULL, nodes = NULL, sparse = TRUE,
                   drop_names = FALSE, directed = NULL, selfloops = NULL) {
  dat <- .multi2weight(edge.list = edge.list, select_cols = select_cols, multiedge = multiedge, aggr_expression = aggr_expression)

  if(isFALSE(selfloops))
    dat %>% filter(.data$source!=.data$target) -> dat

  if (is.null(nodes)) {
    nodes <- nodes_from_el(dat, 1:2)
  } else{
    dat %<>%
      filter(.data$source %in% nodes, .data$target %in% nodes)
  }
# TODO: add subset for nodes if nodes are not all
  if(isFALSE(drop_names))
    adj <- Matrix::sparseMatrix(
      i = match(dat$source,nodes), j = match(dat$target, nodes),
      x = dat$attr, dims = c(length(nodes), length(nodes)), dimnames = list(nodes, nodes)
    )
  if(isTRUE(drop_names))
    adj <- Matrix::sparseMatrix(
      i = match(dat$source,nodes), j = match(dat$target,nodes),
      x = dat$attr, dims = c(length(nodes), length(nodes))
    )

  if(isFALSE(directed))
    adj <- methods::as(Matrix::t(adj) + adj, 'symmetricMatrix')

  if (isFALSE(sparse)) adj <- as.matrix(adj)
  if(isFALSE(drop_names)) return(adj)
  if(isTRUE(drop_names)) return(list(adj=adj,nodes=nodes))
}

#' Generate weighted edge list from an edge list with multiedges.
#'
#' The attribute reporting the aggregated weight of the multiedges is computed according to the expression specified in `aggr_expression`.
#' When no expression is specified, the 'attr' column returns the multiedge multiplicity.
#'
#' @param select_cols optional vector of 2 elements
#'   specifying which columns are the source,target for the edges from which
#'   building the weighted edge list. Otherwise column 1 is assumed to be the source and column
#'   2 the target.
#'
#' @inheritParams el2adj
#'
#' @return a tibble with 3 columns: source, target, and attr, where attr has been computed according to `aggr_expression`.
#' @export
#' @import magrittr
#' @examples
#' el <- data.frame(from= c('a','b','b','c','d','d','d','d'),
#'              to  = c('b','c','d','a','b','a','a','a'),
#'              attr= c( 1,  1,  1 , 1 , 1 , 1,  2 , 3 ),
#'             attr2= c( TRUE,  FALSE,  TRUE , TRUE , FALSE , TRUE,  TRUE , FALSE ))
#' # edgelist with multiedge multiplicities
#' multi2weight(el)
#'
#' # edgelist with sum of attr
#' multi2weight(el, aggr_expression='sum(attr)')
#'
#' # edgelist with mean of attr
#' multi2weight(el, aggr_expression='mean(attr)')
#'
#' # edgelist with min of attr
#' multi2weight(el, aggr_expression='min(attr)')
#'
#' # edgelist with max of attr
#' multi2weight(el, aggr_expression='max(attr)')
#'
#' # edgelist with any() of attr2
#' multi2weight(el, aggr_expression='any(attr2)')
#'
#' # edgelist with all() of attr2
#' multi2weight(el, aggr_expression='all(attr2)')
#'
multi2weight <- function(edge.list, select_cols = NULL, aggr_expression = NULL) {
  .multi2weight(edge.list = edge.list, select_cols = select_cols, multiedge = TRUE, aggr_expression = aggr_expression)
}

# internal function to convert multiedge edgelists to weighted edge lists
.multi2weight <- function(edge.list, select_cols, multiedge, aggr_expression) {
  if (ncol(edge.list) < 2) stop("Not enough columns.")

  attr_cols <- isFALSE(multiedge) & ncol(edge.list) > 2 & length(select_cols) > 2
  col_names <- colnames(edge.list)
  col_names <- .select_cols(col_names = col_names,
                            select_cols = select_cols, attr_cols = attr_cols)
  colnames(edge.list) <- col_names


  if (isFALSE(multiedge)) {
    if (ncol(edge.list) == 2 | length(select_cols) == 2) {
      edge.list$attr <- 1
    }
  }
  if (isTRUE(multiedge)) {
    if (is.null(aggr_expression)) {
      aggr_expression <- parse(text = "dplyr::n()")
    } else {
      aggr_expression <- parse(text = aggr_expression)
    }
  }

  dplyr::as_tibble(edge.list) %>% mutate_at(c("source", "target"), as.character) -> dat

  if (isFALSE(multiedge) && dat %>% count(.data$source, .data$target) %$% any(n > 1)) {
    stop('Multiedges detected. Set "multiedge=TRUE"')
  }

  if (isTRUE(multiedge)) {
    dat %>%
      group_by(.data$source, .data$target) %>%
      summarise(
        attr = eval(aggr_expression)
      ) %>%
      select(.data$source, .data$target, .data$attr) -> dat
  }
  if (isFALSE(multiedge)) {
    dat %>%
      mutate_at(c("source", "target"), as.character) %>%
      select(.data$source, .data$target, .data$attr) -> dat
  }
  return(dat)
}
