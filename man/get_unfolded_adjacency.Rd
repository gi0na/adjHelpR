% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unfolded_temporal_nets.R
\name{get_unfolded_adjacency}
\alias{get_unfolded_adjacency}
\title{Create time unfolded adjacency matrix}
\usage{
get_unfolded_adjacency(
  edge.list,
  select_cols = NULL,
  nodes = NULL,
  vertical = FALSE,
  sparse = TRUE,
  ...
)
}
\arguments{
\item{edge.list}{data.frame or tibble containing the edge list. It needs at
least three column: the column with edge sources, the edge targets, and the
timestamps of each edge. The order of columns should be 'timestamp',
'source', 'target'. If the edge.list columns come in different orders, use
`select_cols` to specify the right order. See the example for details.}

\item{select_cols}{optional vector of 3 (2 for multi-graphs) elements
specifying which columns are the source,target, and attributes from which
building the graph. Otherwise Column 1 is assumed to be the source, column
2 the target, column 3 the attribute. In the case of multi-graphs, the
third element is not needed and the number of edges between each pair of
vertices is computed according to '\code{aggr_expression}'.}

\item{nodes}{optional vector containing all node names in case disconnected
nodes should be included.}

\item{vertical}{unfold vertically or horizontally? Defaults to FALSE
(horizontal)}

\item{sparse}{boolean, return sparse matrix? default to TRUE}

\item{...}{extra parameters passed to internal methods}
}
\value{
The time unfolded (nodes x timestamps) x (nodes x timestamps)
adjacency matrix. The node order follows the vertical or horizontal algnment
}
\description{
Create time unfolded adjacency matrix
}
\examples{
el <- data.frame(
from = c('A','B',  'A','B','B',   'A','C','C',    'A','B','C',    'D'),
to   = c('C','C',  'C','C','D',   'C','D','E',    'C','C','E',    'E'),
ts   = c( 1,  1,    2,  2,  2,     3,  3,  3,      4,  4,  4,      6)
)
get_unfolded_adjacency(el, select_cols = 1:3)
}
