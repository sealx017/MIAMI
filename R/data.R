#' Marker intensity data.
#'
#' A dataset with expression data of five markers, HLA-DR, CD45RO, H3K27me3, H3K9ac
#' and HLA_Class_1 in different cells of 39 subjects.
#'
#' @format A data frame with 190276 rows and 6 variables:
#' \describe{
#'   \item{ID}{Subject IDs}
#'   \item{HLA.DR}{Scaled intensity of HLA-DR marker}
#'   \item{CD45RO}{Scaled intensity of CD45RO marker}
#'   \item{HLA_Class_1}{Scaled intensity of HLA_Class_1 marker}
#'   \item{H3K27me3}{Scaled intensity of H3K27me3 marker}
#'   \item{H3K9ac}{Scaled intensity of H3K9ac marker}
#' }
#' @source \url{https://www.angelolab.com/mibi-data}
"marker_data"



#' Clinical outcomes data.
#'
#' A dataset two clinical outcomes, time to death and time to recurrence along with the censoring indicators
#' and one covariate, Age for the same set of subjects as marker_data.
#'
#' @format A data frame with 190276 rows and 6 variables:
#' \describe{
#'   \item{ID}{Subject IDs}
#'   \item{Recurrence}{The censoring indicator for recurrence, = 0 implies event}
#'   \item{Recurrence_time}{Time to recurrence}
#'   \item{Survival}{The censoring indicator for survival, = 0 implies event}
#'   \item{Survival_time}{Time to death}
#'   \item{Age}{Age of the subjects}
#' }
#' @source \url{https://www.angelolab.com/mibi-data}
"clinical_data"
