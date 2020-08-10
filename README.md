# adjHelpR
R tools to build and work with adjacency matrices.
This package provides a set of efficient functions to deal with adjacency matrices in R.
It allows to generate adjacency matrices from weighted edgelists, create weighted edgelists from edgelists with repeated edges, and handle these objects.
It exploits `dplyr` and `tibble`s to deal efficiently with large edgelists, and `Matrix` sparse matrices to store adjacency matrices in sparse format.

```
# Install the development version from GitHub:
devtools::install_github("gi0na/adjHelpR")
```
