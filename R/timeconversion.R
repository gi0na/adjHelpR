################################################################################
### General
################################################################################

#' Converts the given unix_timestamp into a date-string.
#'
#' @param unix_timestamp The unix-timestamp to convert. A unix_timestamp is the
#'       number of seconds since the 1st January 1970 (as an integer).
#' @param out_format The format code of the return value. For the specification
#'       see section "Details" in help(strptime)
#'
#' @return String representing the date corresponding to the unix_timestamp
#'       formatted according to out_format.
#' @export
#'
#' @examples
#'       >>> UnixToDatestring(0)
#'       1970-01-01
#'       >>> UnixToDatestring(86400)
#'       1970-01-02
unix_to_datestring <- function(unix_timestamp, out_format = "%Y-%m-%d"){
  # Converts the given unix_timestamp into a date-string.
  #
  # Examples:
      # >>> UnixToDatestring(0)
      # 1970-01-01
      # >>> UnixToDatestring(86400)
      # 1970-01-02
  #
  # Args:
  #     unix_timestamp: The unix-timestamp to convert. A unix_timestamp is the
  #         number of seconds since the 1st January 1970 (as an integer).
#     out_format: The format code of the return value. For the specification
#         see section "Details" in help(strptime)
#
# Returns:
#     String representing the date corresponding to the unix_timestamp
#         formatted according to out_format.
#
# Source:
#     This function is built around a code-snippet in gentoo.r (i.e. in
#         lib_external.R)
unix_timestamp %>%
  as.POSIXlt(origin = "1970-01-01 00:00.00 UTC") %>%
  format(out_format)
}


#' Converts the given datestring into a unix-timestamp.
#'
#' @param datestring The date to convert given as a string of the specified
#'        in_format.
#' @param in_format The format code of the given datestring. For the
#'         specification see section "Details" in help(strptime)
#'
#' @return Integer representing the date as a unix_timestamp.
#' @export
#'
#' @examples
#'       >>> DatestringToUnix("1970-01-01")
#'       0
#'       >>> DatestringToUnix("1970-01-02")
#'       86400
datestring_to_unix <- function(datestring, in_format = "%Y-%m-%d"){
  # Converts the given datestring into a unix-timestamp.
  #
  # Examples:
      # >>> DatestringToUnix("1970-01-01")
      # 0
      # >>> DatestringToUnix("1970-01-02")
      # 86400
  #
  # Args:
  #     datestring: The date to convert given as a string of the specified
  #         in_format.
#     in_format: The format code of the given datestring. For the
#         specification see section "Details" in help(strptime)
#
# Returns:
#     Number representing the date as a unix_timestamp.
datestring %>%
  strptime(format = in_format, tz = "GMT") %>%
  as.POSIXct(origin = "1970-01-01 00:00.00 UTC") %>%
  as.numeric()
}

#' Detect date format of a string
#'
#' @param datestring string representing a date
#'
#' @return a string with the format code of the given datestring
#' @export
#'
#' @examples
#' detect_date_format("1970-01-01")
#'
detect_date_format <- function(datestring){
  format <- NULL
  if(all(grepl(pattern = '[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}', datestring)))
    format <- '%Y-%m-%d %H:%M:%S'
  if(all(grepl(pattern = '[0-9]{2}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}', datestring)))
    format <- '%y-%m-%d %H:%M:%S'
  if(all(grepl(pattern = '[0-9]{2}-[a-zA-Z]{3}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}', datestring)))
    format <- '%y-%h-%d %H:%M:%S'
  if(all(grepl(pattern = '[0-9]{4}-[a-zA-Z]{3}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}', datestring)))
    format <- '%Y-%h-%d %H:%M:%S'
  if(all(grepl(pattern = '[0-9]{2}/[0-9]{2}/[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}', datestring)))
    format <- '%D %H:%M:%S'
  return(format)
}
