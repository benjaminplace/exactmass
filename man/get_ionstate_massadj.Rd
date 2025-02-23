% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exactmass_calc.R
\name{get_ionstate_massadj}
\alias{get_ionstate_massadj}
\title{Get mass adjustment for compound's ion state}
\usage{
get_ionstate_massadj(ionstate = "none")
}
\arguments{
\item{ionstate}{character string containing the ion state of the formula in `[M+H]+` format}
}
\value{
numeric monoisotopic mass
}
\description{
Get mass adjustment for compound's ion state
}
\examples{
get_ionstate_massadj("[M-H]+")
}
