################################################################################
### Time-window Extraction
################################################################################

#' Builds rolling multiedge time-window networks from the edge_list.
#'
#' @param el A data-frame with at least 3 columns: "source", "target",
#'   "timestamp", where "from" and "to" are vertex identifiers, and "time" is a
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
#'   at `earliest`, possibly excluding an incomplete time-window at `latest`. If
#'   "latest", the last time-window ends at `latest`, and possibly an incomplete
#'   time-window at earliest is discarded.
#' @param select_cols an (optional) vector specifying which columns to use. The
#'   first entry should contain either the name or index for the timestamp
#'   column, the second for the sources of links, the third for the targets of
#'   links. When not passed, it will be assumed that the first three colummns
#'   are `c(timestamp,source,target)`.
#' @param as_date (optional) boolean identifying if timestamps are in date
#'   formats, or unix seconds.
#' @param out_format character vector specifying the format of the output:
#'   'edgelist' for edgelists, 'adjacency' for adjacency matrices, `igraph` for
#'   igraph graphs. Defaults to 'edgelist'. In the case of igraph graphs, adding
#'   `igraph_weighted` allows to generated weighted graphs.
#' @param ... additional parameters passed to internal constructors. E.g., to
#'   `get_adjacency`.
#'
#' @return A tibble with two columns "network" and "time", where "network" is
#'   the time-window and "time" the corresponding time at which the time-window
#'   starts.
#' @export
#' @author CZ, GC
#' @examples
get_rolling_windows <- function(el,
                                window_size,
                                step_size = NULL,
                                start_time = NULL,
                                end_time = NULL,
                                out_format = 'edgelist',
                                flush = "earliest",
                                select_cols = NULL,
                                as_date=NULL,
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
  #     earliest: (optional) If given, the rolling time-windows are built
  #         starting at this timestamp. Else, they are started at
  #         min(edgelist$time).
  #     latest: (optional) If given, rolling time-windows are created until they
  #         exceed this value. Else, this value defaults to max(edgelist$time).
  #     flush: (optional) If "earliest" or unset, the first time-window starts at
  #         `earliest`, possibly excluding an incomplete time-window at `latest`.
  #         If "latest", the last time-window ends at `latest`, and possibly an
  #         incomplete time-window at earliest is discarded.
  #
  # Returns:
  #     A tibble with two columns "network" and "time", where "network"
  #     is the time-window and "time" the corresponding time at which the
  #     time-window starts.


  # Fix metadata
  if (ncol(el) < 3) stop("Not enough columns.")

  col_names <- colnames(el)
  col_names <- .select_cols_temporal(col_names = col_names,
                                     select_cols = select_cols, attr_cols = FALSE)
  colnames(el) <- col_names

  if(is.null(start_time)) start_time <- min(el$timestamp)
  if(is.null(end_time)) end_time <- max(el$timestamp)

  if(is.null(as_date)){
    as_date <- suppressWarnings(is.na(as.numeric(el$timestamp[1])))
  }

  if(isTRUE(as_date)){
    date_format <- detect_date_format(el$timestamp[1])
    if(is.null(date_format)) stop('Wrong timestamp format.')
    el %>% mutate(timestamp = datestring_to_Unix(datestring = timestamp, in_format = date_format))
    date_format <- detect_date_format(start_time)
    if(is.null(date_format)) stop('Wrong start_time format.')
    start_time <- datestring_to_Unix(datestring = start_time, in_format = date_format)
    if(!is.null(end_time)){
      date_format <- detect_date_format(end_time)
      if(is.null(date_format)) stop('Wrong end_time format.')
      end_time <- datestring_to_Unix(datestring = end_time, in_format = date_format)
    }
  }

  if(is.null(step_size)) step_size <- window_size
  # Build window timestamps
  if (flush == "earliest") {
    starting_times <- seq(start_time, end_time - window_size, by = step_size)
    ending_times <- starting_times + window_size
  } else if (flush == "latest") {
    seq(latest - window_size, earliest, by = -step_size) %>%
      rev() ->
      starting_times
    ending_times <- starting_times + window_size
  } else {
    stop("Unknown flush-option given:  ", flush)
  }

  # Create corresponding networks
  windows <- list()
  for (i in seq_along(starting_times)) {
    lower <- starting_times[i]
    upper <- ending_times[i]
    if(grepl('edge', out_format)){
      el %>%
        .el2slice(start_time = lower, end_time = upper, index = FALSE) ->
        curr_window
    }
    if(grepl('adj', out_format)){
      el %>%
        .el2slice(start_time = lower, end_time = upper, index = FALSE) %>%
        get_adjacency(select_cols = select_cols, ...) ->
        curr_window
    }
    if(grepl('igraph', out_format)){
      weighted <- NULL
      if(grepl('weight', out_format)) weighted <- TRUE
      if(requireNamespace("igraph", quietly = TRUE)){
        el %>%
          .el2slice(start_time = lower, end_time = upper, index = FALSE) %>%
          get_adjacency(select_cols = select_cols, ...) %>%
          igraph::graph_from_adjacency_matrix(weighted = weighted) ->
          curr_window
      }
    }
    windows <- append(windows, list(curr_window))
  }

  # Return
  return(dplyr::tibble("time" = starting_times, "network" = windows))
}

#' Filter edgelist for a time slice
#'
#' @param el edgelist
#' @param select_cols
#' @param start_time
#' @param end_time
#' @param duration
#' @param index
#' @param as_date
#' @param attr_cols
#'
#' @return
#' @export
#'
#' @examples
get_slice_edgelist <- function(el, select_cols = NULL, start_time, end_time=NULL, duration=NULL, index=FALSE, as_date=NULL){
  if (ncol(el) < 3) stop("Not enough columns.")

  col_names <- colnames(el)
  col_names <- .select_cols_temporal(col_names = col_names,
                                     select_cols = select_cols, attr_cols = FALSE)
  colnames(el) <- col_names

  if(is.null(as_date)){
    as_date <- suppressWarnings(is.na(as.numeric(el$timestamp[1])))
  }

  if(isTRUE(as_date)){
    date_format <- detect_date_format(el$timestamp[1])
    if(is.null(date_format)) stop('Wrong timestamp format.')
    el %>% mutate(timestamp = datestring_to_Unix(datestring = timestamp, in_format = date_format))
    date_format <- detect_date_format(start_time)
    if(is.null(date_format)) stop('Wrong start_time format.')
    start_time <- datestring_to_Unix(datestring = start_time, in_format = date_format)
    if(!is.null(end_time)){
      date_format <- detect_date_format(end_time)
      if(is.null(date_format)) stop('Wrong end_time format.')
      end_time <- datestring_to_Unix(datestring = end_time, in_format = date_format)
    }
  }

  if(is.null(c(end_time,duration)))
    stop('Specify either end_time or duration')

  if(is.null(end_time))
    end_time <- start_time+duration

  .el2slice(el, start_time, end_time, index)
}

.el2slice <- function(el, start_time, end_time, index){

  if(!index){
    el %>%
      dplyr::filter(start_time <= timestamp & timestamp < end_time) ->
      curr_slice
  } else{
    start_index <- start_time
    end_index <- end_time-1
    dplyr::as_tibble(el) %>%
      slice(start_index:end_index) ->
      curr_slice
  }
  return(curr_slice)
}
