% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjacency_matrix.R
\name{get_adjacency}
\alias{get_adjacency}
\alias{get_adjacency.tbl}
\alias{get_adjacency.tbl_graph}
\alias{get_adjacency.data.frame}
\alias{get_adjacency.matrix}
\alias{get_adjacency.dgTMatrix}
\alias{get_adjacency.dgCMatrix}
\alias{get_adjacency.sparseMatrix}
\alias{get_adjacency.default}
\alias{get_adjacency.igraph}
\title{Compute the adjacency matrix from various types of objects.}
\usage{
get_adjacency(
  x,
  select_cols = NULL,
  multiedge = NULL,
  aggr_expression = NULL,
  nodes = NULL,
  sparse = TRUE,
  drop_names = FALSE,
  directed = NULL,
  selfloops = NULL,
  edgelist = NULL,
  ...
)

\method{get_adjacency}{tbl}(
  x,
  select_cols = NULL,
  multiedge = NULL,
  aggr_expression = NULL,
  nodes = NULL,
  sparse = TRUE,
  drop_names = FALSE,
  directed = NULL,
  selfloops = NULL,
  edgelist = NULL,
  ...
)

\method{get_adjacency}{tbl_graph}(
  x,
  select_cols = NULL,
  multiedge = NULL,
  aggr_expression = NULL,
  nodes = NULL,
  sparse = TRUE,
  drop_names = FALSE,
  directed = NULL,
  selfloops = NULL,
  edgelist = NULL,
  ...
)

\method{get_adjacency}{data.frame}(
  x,
  select_cols = NULL,
  multiedge = NULL,
  aggr_expression = NULL,
  nodes = NULL,
  sparse = TRUE,
  drop_names = FALSE,
  directed = NULL,
  selfloops = NULL,
  edgelist = NULL,
  ...
)

\method{get_adjacency}{matrix}(
  x,
  select_cols = NULL,
  multiedge = NULL,
  aggr_expression = NULL,
  nodes = NULL,
  sparse = TRUE,
  drop_names = FALSE,
  directed = NULL,
  selfloops = NULL,
  edgelist = NULL,
  ...
)

\method{get_adjacency}{dgTMatrix}(
  x,
  select_cols = NULL,
  multiedge = NULL,
  aggr_expression = NULL,
  nodes = NULL,
  sparse = TRUE,
  drop_names = FALSE,
  directed = NULL,
  selfloops = NULL,
  ...
)

\method{get_adjacency}{dgCMatrix}(
  x,
  select_cols = NULL,
  multiedge = NULL,
  aggr_expression = NULL,
  nodes = NULL,
  sparse = TRUE,
  drop_names = FALSE,
  directed = NULL,
  selfloops = NULL,
  ...
)

\method{get_adjacency}{sparseMatrix}(
  x,
  select_cols = NULL,
  multiedge = NULL,
  aggr_expression = NULL,
  nodes = NULL,
  sparse = TRUE,
  drop_names = FALSE,
  directed = NULL,
  selfloops = NULL,
  ...
)

\method{get_adjacency}{default}(
  x,
  select_cols = NULL,
  multiedge = NULL,
  aggr_expression = NULL,
  nodes = NULL,
  sparse = TRUE,
  drop_names = FALSE,
  directed = NULL,
  selfloops = NULL,
  ...
)

\method{get_adjacency}{igraph}(
  x,
  select_cols = NULL,
  multiedge = NULL,
  aggr_expression = NULL,
  nodes = NULL,
  sparse = TRUE,
  drop_names = FALSE,
  directed = NULL,
  selfloops = NULL,
  ...
)
}
\arguments{
\item{x}{object from which generating the adjacency matrix}

\item{select_cols}{optional vector of 3 (2 for multi-graphs) elements
specifying which columns are the source,target, and attributes from which
building the graph. Otherwise Column 1 is assumed to be the source, column
2 the target, column 3 the attribute. In the case of multi-graphs, the
third element is not needed and the number of edges between each pair of
vertices is computed according to '\code{aggr_expression}'.}

\item{multiedge}{boolean, are there multiple edges? defaults to FALSE.}

\item{aggr_expression}{string, the expression used to compute the aggregated
value in the adjacency matrix in the presence of multiedges. It defaults to
'\code{dplyr::n()}'. If `edge.list` contains other columns such as `weight`,
examples of `aggr_expression` could be `"mean(weight)"` to report the
average of the `weight` column of the multiedges, `"min(weight)"` to report
the minimum of the `weight` column of the multiedges.}

\item{nodes}{optional vector containing all node names in case disconnected
nodes should be included.}

\item{sparse}{boolean, return sparse matrix? default to TRUE}

\item{drop_names}{boolean, drop names from matrix and return a vector of
names together with the matrix. This option saves considerable memory when
dealing with graphs with long node names.}

\item{directed}{boolean, optional parameter, if FALSE, forces to undirected
upper triangular adjacency matrix}

\item{selfloops}{boolean, optional parameter, if ignores selfloops.}

\item{edgelist}{optional boolean. If `x` is a matrix, should treat `x` as an
edgelist? If not provided, tries to infer if `x` is an edgelist or an
incidence matrix from its shape.}

\item{...}{additional parameters to and from the main method. Currently not
used.}
}
\value{
A transition matrix in the same format of `adj`.
}
\description{
If `x` is a tibble or a dataframe, `get_adjacency` assumes that it has been
passed an edgelist. If `x` is a matrix, `get_adjacency` assumes that it has
been passed an incidence matrix and constructs the corresponding adjacency
matrix. If `x` is a different object, `get_adjacency` uses the appropriate
method if it exists. Otherwise it tries to coerce `x` to a `data.frame`.
}
\examples{
el <- data.frame(from= c('a','b','b','c','d','d'),
                to  = c('b','c','d','a','b','a'),
                attr= c( 12, 6, 12 , 6 , 6 , 6 ))
adj <- get_adjacency(el)

}
