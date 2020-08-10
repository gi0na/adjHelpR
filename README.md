
<!-- README.md is generated from README.Rmd. Please edit that file -->

# adjHelpR

<!-- badges: start -->

<!-- badges: end -->

The goal of adjHelpR is to provide tools to generate adjacency matrices
from weighted edgelists, create weighted edgelists from edgelists with
repeated edges, and handle these objects. It exploits `dplyr` and
`tibble`s to deal efficiently with large edgelists, and `Matrix` sparse
matrices to store adjacency matrices in sparse format.

## Installation

You can install the released version of adjHelpR from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("adjHelpR")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gi0na/adjHelpR")
```

## Example

This is a basic example on how to generate a weighted adjacency matrix
from an edge list:

``` r
library(adjHelpR)

el <- data.frame(from= c('a','b','b','c','d','d'),
                 to  = c('b','c','d','a','b','a'),
                 attr= c( 12, 6, 12 , 6 , 6 , 6 ))
(adj <- get_adjacency(el))
#> 4 x 4 sparse Matrix of class "dgCMatrix"
#>   a  b c  d
#> a . 12 .  .
#> b .  . 6 12
#> c 6  . .  .
#> d 6  6 .  .

(el_unweighted <- tibble::as_tibble(do.call(rbind, 
                        apply(el, 1, function(row) matrix(rep(row[1:2], each=as.integer(row[3])), ncol = 2))), 
                        .name_repair='minimal'))
#> # A tibble: 48 x 2
#>    ``    ``   
#>    <chr> <chr>
#>  1 a     b    
#>  2 a     b    
#>  3 a     b    
#>  4 a     b    
#>  5 a     b    
#>  6 a     b    
#>  7 a     b    
#>  8 a     b    
#>  9 a     b    
#> 10 a     b    
#> # … with 38 more rows
(adj <- get_adjacency(el_unweighted, multiedge = TRUE))
#> 4 x 4 sparse Matrix of class "dgCMatrix"
#>   a  b c  d
#> a . 12 .  .
#> b .  . 6 12
#> c 6  . .  .
#> d 6  6 .  .
```
