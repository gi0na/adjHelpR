% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/timeconversion.R
\name{unix_to_datestring}
\alias{unix_to_datestring}
\title{Converts the given unix_timestamp into a date-string.}
\usage{
unix_to_datestring(unix_timestamp, out_format = "\%Y-\%m-\%d")
}
\arguments{
\item{unix_timestamp}{The unix-timestamp to convert. A unix_timestamp is the
number of seconds since the 1st January 1970 (as an integer).}

\item{out_format}{The format code of the return value. For the specification
see section "Details" in help(strptime)}
}
\value{
String representing the date corresponding to the unix_timestamp
      formatted according to out_format.
}
\description{
Converts the given unix_timestamp into a date-string.
}
\examples{
      unix_to_datestring(0)
      # 1970-01-01
      unix_to_datestring(86400)
      # 1970-01-02
}
