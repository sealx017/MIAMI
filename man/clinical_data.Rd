% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{clinical_data}
\alias{clinical_data}
\title{Clinical outcomes data.}
\format{
A data frame with 190276 rows and 6 variables:
\describe{
  \item{ID}{Subject IDs}
  \item{Recurrence}{The censoring indicator for recurrence, = 0 implies event}
  \item{Recurrence_time}{Time to recurrence}
  \item{Survival}{The censoring indicator for survival, = 0 implies event}
  \item{Survival_time}{Time to death}
  \item{Age}{Age of the subjects}
}
}
\source{
\url{https://www.angelolab.com/mibi-data}
}
\usage{
clinical_data
}
\description{
A dataset two clinical outcomes, time to death and time to recurrence along with the censoring indicators
and one covariate, Age for the same set of subjects as marker_data.
}
\keyword{datasets}
