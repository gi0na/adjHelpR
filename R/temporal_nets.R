################################################################################
### Time-window Extraction
################################################################################

#' Builds rolling multiedge time-window networks from the edge_list.
#'
#' @param edge.list A data-frame, tibble, with at least 3 columns: "source",
#'   "target", "timestamp", or a graph object that can be coerced to such table,
#'   where "from" and "to" are vertex identifiers, and "time" is a
#'   Unix-timestamp when the edge occurred. Alternatively, the user can select
#'   which column is which passing `select_cols`
#' @param window_size The duration of each time-window network given as the
#'   number of seconds.
#' @param step_size The shift between consecutive time-window network given as
#'   the number of seconds.
#' @param start_time (optional) If given, the rolling time-windows are built
#'   starting at this timestamp. Else, they are started at min(edgelist$time).
#' @param end_time (optional) If given, rolling time-windows are created until
#'   they exceed this value. Else, this value defaults to max(edgelist$time).
#' @param flush (optional) If "earliest" or unset, the first time-window starts
#'   at `start_time`, possibly excluding an incomplete time-window at
#'   `end_time`. If "latest", the last time-window ends at `end_time`, and
#'   possibly an incomplete time-window at start_time is discarded.
#' @param select_cols an (optional) vector specifying which columns to use. The
#'   first column should contain either the name or index for the sources of
#'   links, the second for the targets of links. The third entry should contain
#'   either the name or index for the timestamp column. When not passed, it will
#'   be assumed that the first three columns are `c(timestamp,source,target)`.
#'   If select_cols has a fourth entry, it will use this as edgeweight.
#' @param as_date (optional) boolean identifying if timestamps are in date
#'   formats, or unix seconds.
#' @param out_format character vector specifying the format of the output:
#'   'edgelist' for edgelists, 'adjacency' for adjacency matrices, 'tbl_graph' for `tbl_graph`s, 'igraph' for
#'   `igraph` graphs. Defaults to 'edgelist'. In the case of igraph and tbl_graph graphs, adding
#'   '_weighted' (as in 'igraph_weighted') allows to generate weighted graphs.
#' @param ... additional parameters passed to internal constructors. E.g., to
#'   `get_adjacency`.
#' @param ncores integer, number of cores to use. Defaults to 1.
#'
#' @return A tibble with two columns "network" and "time", where "network" is
#'   the time-window and "time" the corresponding time at which the time-window
#'   starts.
#' @export
#' @author CZ, GC
#'
get_rolling_windows <- function(edge.list,
                                window_size,
                                step_size = NULL,
                                start_time = NULL,
                                end_time = NULL,
                                out_format = 'edgelist',
                                flush = "earliest",
                                select_cols = NULL,
                                as_date = NULL,
                                ncores = NULL,
                                ...) {
  UseMethod("get_rolling_windows")
}

#' @rdname get_rolling_windows
#' @param nodes_tbl (optional) tibble containing node attributes to be added to
#'   tbl_graph objects when tbl_graph output is chosen. The tibble should have
#'   a node_key column matching node names used in the edge.list.
#' @export
get_rolling_windows.tbl <- function(edge.list,
                                    window_size,
                                    step_size = NULL,
                                    start_time = NULL,
                                    end_time = NULL,
                                    out_format = 'edgelist',
                                    flush = "earliest",
                                    select_cols = NULL,
                                    as_date = NULL,
                                    ncores = NULL,
                                    nodes_tbl = NULL,
                                    ...) {
  # Builds rolling multiedge time-window networks from the edge_list.
  #
  # Args:
  #     edge_list: A data-frame with 3 columns: "from", "to", "time",
  #         where "from" and "to" are vertex identifiers, and "time" is a
  #         Unix-timestamp when the edge occurred.
  #     window_size: The duration of each time-window network given as the
  #         number of seconds.
  #     step_size: The shift between consecutive time-window network given as
  #         the number of seconds.
  #     directed: Whether to create directed or undirected edges.
  #     start_time: (optional) If given, the rolling time-windows are built
  #         starting at this timestamp. Else, they are started at
  #         min(edgelist$time).
  #     end_time: (optional) If given, rolling time-windows are created until they
  #         exceed this value. Else, this value defaults to max(edgelist$time).
  #     flush: (optional) If "earliest" or unset, the first time-window starts at
  #         `start_time`, possibly excluding an incomplete time-window at `end_time`.
  #         If "latest", the last time-window ends at `end_time`, and possibly an
  #         incomplete time-window at start_time is discarded.
  #
  # Returns:
  #     A tibble with two columns "network" and "time", where "network"
  #     is the time-window and "time" the corresponding time at which the
  #     time-window starts.


  # Fix metadata
  if (ncol(edge.list) < 3) stop("Not enough columns.")

  col_names <- colnames(edge.list)
  col_names <- .select_cols_temporal(col_names = col_names,
                                     select_cols = select_cols, attr_cols = length(select_cols)==4)
  colnames(edge.list) <- col_names
  if(length(select_cols)>3){
    edge.list %>% select('source','target','timestamp','attr') -> edge.list
  } else{
    edge.list %>% select('source','target','timestamp') -> edge.list
  }

  if(is.null(start_time)) start_time <- min(edge.list$timestamp)
  if(is.null(end_time)) end_time <- max(edge.list$timestamp)

  if(is.null(as_date)){
    as_date <- suppressWarnings(is.na(as.numeric(edge.list$timestamp[1])))
  }

  if(isTRUE(as_date)){
    date_format <- detect_date_format(edge.list$timestamp[1])
    if(is.null(date_format)) stop('Wrong timestamp format.')
    edge.list %>% mutate(timestamp = datestring_to_unix(datestring = .data$timestamp, in_format = date_format))
    date_format <- detect_date_format(start_time)
    if(is.null(date_format)) stop('Wrong start_time format.')
    start_time <- datestring_to_unix(datestring = start_time, in_format = date_format)
    if(!is.null(end_time)){
      date_format <- detect_date_format(end_time)
      if(is.null(date_format)) stop('Wrong end_time format.')
      end_time <- datestring_to_unix(datestring = end_time, in_format = date_format)
    }
  }

  if(is.null(step_size)) step_size <- window_size
  # Build window timestamps
  if (flush == "earliest") {
    starting_times <- seq(start_time, end_time - window_size, by = step_size)
    ending_times <- starting_times + window_size
  } else if (flush == "latest") {
    seq(end_time - window_size, start_time, by = -step_size) %>%
      rev() ->
      starting_times
    ending_times <- starting_times + window_size
  } else {
    stop("Unknown flush-option given:  ", flush)
  }

  # Create corresponding networks
  if(is.null(ncores)) ncores <- 1
  windows <-
  # for (i in seq_along(starting_times)) {
    pbmcapply::pbmclapply(X = seq_along(starting_times), FUN = function(i){
      lower <- starting_times[i]
      upper <- ending_times[i]
      if(grepl('edge', out_format)){
        edge.list %>%
          .el2slice(start_time = lower, end_time = upper, index = FALSE) ->
          curr_window
      }
      if(grepl('graph', out_format)){
        if(requireNamespace("tidygraph", quietly = TRUE)){
          if(is.null(nodes_tbl)){
          edge.list %>%
            .el2slice(start_time = lower, end_time = upper, index = FALSE) %>%
            conditional_makeweighted(select_cols = select_cols, out_format=out_format) %>%
            tidygraph::as_tbl_graph() ->
            curr_window
          } else{
            edge.list %>%
              .el2slice(start_time = lower, end_time = upper, index = FALSE) %>%
              conditional_makeweighted(select_cols = select_cols, out_format=out_format) %>%
              tidygraph::as_tbl_graph() %>%
              tidygraph::activate('nodes') %>%
              rename(node_key = 'name') %>%
              inner_join(nodes_tbl) %>%
              dplyr::select(-'node_key') ->
              curr_window
          }
        }
      }
      if(grepl('adj', out_format)){
        edge.list %>%
          .el2slice(start_time = lower, end_time = upper, index = FALSE) %>%
          get_adjacency(select_cols = select_cols, multiedge = TRUE, ...) ->
          curr_window
      }
      # if(grepl('igraph', out_format)){
      #   weighted <- NULL
      #   if(grepl('weight', out_format)) weighted <- TRUE
      #   if(requireNamespace("igraph", quietly = TRUE)){
      #     edge.list %>%
      #       .el2slice(start_time = lower, end_time = upper, index = FALSE) %>%
      #       igraph::graph_from_adjacency_matrix(weighted = weighted) ->
      #       curr_window
      #   }
      # }
      return(curr_window)
    }, mc.cores = ncores)

  # Return
  return(dplyr::tibble("time" = starting_times, "network" = windows))
}

conditional_makeweighted <- function(x, select_cols, out_format){
  if(grepl('weight', out_format)){
    x %>% multi2weight(select_cols) -> x
  }
  return(x)
}

#' @rdname get_rolling_windows
#' @export
get_rolling_windows.tbl_graph <- function(edge.list,
                                           window_size,
                                           step_size = NULL,
                                           start_time = NULL,
                                           end_time = NULL,
                                           out_format = 'edgelist',
                                           flush = "earliest",
                                           select_cols = NULL,
                                           as_date = NULL,
                                           ncores = NULL,
                                           ...) {
  edge.list %>% tidygraph::activate('edges') %>% dplyr::as_tibble() %>%
    get_rolling_windows(window_size, step_size,
                        start_time, end_time, out_format,
                        flush, select_cols, as_date, ncores,
                        nodes_tbl = edge.list %>%
                          tidygraph::activate('nodes') %>%
                          dplyr::as_tibble() %>%
                          dplyr::mutate(node_key = as.character(1:length(.data$name))), ...)

}

#' @rdname get_rolling_windows
#' @export
get_rolling_windows.igraph <- function(edge.list,
                                          window_size,
                                          step_size = NULL,
                                          start_time = NULL,
                                          end_time = NULL,
                                          out_format = 'edgelist',
                                          flush = "earliest",
                                          select_cols = NULL,
                                          as_date = NULL,
                                          ncores = NULL,
                                          ...) {
  edge.list %>% tidygraph::as_tbl_graph() %>%
    get_rolling_windows(edge.list, window_size, step_size,
                        start_time, end_time, out_format,
                        flush, select_cols, as_date, ncores,
                        tbl_graph = edge.list %>% tidygraph::as_tbl_graph(), ...)

}

#' @rdname get_rolling_windows
#' @export
get_rolling_windows.default <- function(edge.list,
                                          window_size,
                                          step_size = NULL,
                                          start_time = NULL,
                                          end_time = NULL,
                                          out_format = 'edgelist',
                                          flush = "earliest",
                                          select_cols = NULL,
                                          as_date = NULL,
                                          ncores = NULL,
                                          ...) {
  edge.list %>% dplyr::as_tibble() %>%
    get_rolling_windows(edge.list, window_size, step_size,
                        start_time, end_time, out_format,
                        flush, select_cols, as_date, ncores,...)

}

#' Filter edgelist for a time slice
#'
#' @inheritParams get_adjacency
#' @param edge.list the edgelist from which to compute the time slices
#' @param start_time the timestamp a which start the timewindow slice. It can be an integer defining either the index of the start point or the unit of time at which to cut, or a datestring string.
#' @param end_time (optional) the timestamp a which to end the timewindow slice. It can be an integer defining either the index of the end point or the unit of time at which to cut, or a datestring string.
#' @param duration (optional) the duration of the timewindow. If `end_time` is not provided it defines how long after `start_time` the window should be cut.
#' @param index boolean, set to TRUE if start_time and end_time are supplied as indices. Defaults to FALSE.
#' @param as_date (optional) boolean, are the timestamps supplied as datestrings?
#'
#' @return a tibble with the edgelist cut such that timestamps are between `start_time` and `end_time`
#' @export
#'
#' @examples
#' el <- data.frame(from= c('a','b','b','c','d','d'),
#'                 to  = c('b','c','d','a','b','a'),
#'                 t = 1:6 )
#' slice <- get_slice_edgelist(el, select_cols = 1:3, start_time = 2, duration = 3)
get_slice_edgelist <- function(edge.list, select_cols = NULL, start_time, end_time=NULL, duration=NULL, index=FALSE, as_date=NULL){
  if (ncol(edge.list) < 3) stop("Not enough columns.")

  col_names <- colnames(edge.list)
  col_names <- .select_cols_temporal(col_names = col_names,
                                     select_cols = select_cols, attr_cols = FALSE)
  colnames(edge.list) <- col_names

  if(is.null(as_date)){
    as_date <- suppressWarnings(is.na(as.numeric(edge.list$timestamp[1])))
  }

  if(isTRUE(as_date)){
    date_format <- detect_date_format(edge.list$timestamp[1])
    if(is.null(date_format)) stop('Wrong timestamp format.')
    edge.list %>% mutate(timestamp = datestring_to_unix(datestring = .data$timestamp, in_format = date_format))
    date_format <- detect_date_format(start_time)
    if(is.null(date_format)) stop('Wrong start_time format.')
    start_time <- datestring_to_unix(datestring = start_time, in_format = date_format)
    if(!is.null(end_time)){
      date_format <- detect_date_format(end_time)
      if(is.null(date_format)) stop('Wrong end_time format.')
      end_time <- datestring_to_unix(datestring = end_time, in_format = date_format)
    }
  }

  if(is.null(c(end_time,duration)))
    stop('Specify either end_time or duration')

  if(is.null(end_time))
    end_time <- start_time+duration

  .el2slice(edge.list, start_time, end_time, index)
}

.el2slice <- function(edge.list, start_time, end_time, index){

  if(!index){
    edge.list %>%
      dplyr::filter(start_time <= .data$timestamp & .data$timestamp < end_time) ->
      curr_slice
  } else{
    start_index <- start_time
    end_index <- end_time-1
    dplyr::as_tibble(edge.list) %>%
      slice(start_index:end_index) ->
      curr_slice
  }
  return(curr_slice)
}
