#'@title Standardizing a marker expression vector between 0 and 1.
#' @param x is the vector of expression of a marker.
#' @return It returns a standardized expression vector.

#' @export
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

#'@title Compute the PDF of a univariate normal distribution.
#' @param x is the point where the density will be evaluated at.
#' @param mean is the mean of the normal distribution.
#' @param sd is the standard deviation of the normal distribution.
#' @return It returns the density  of a univariate normal distribution with specified mean and sd at point x.

#' @export
dnorm_mine<-function(x, mean, sd){
  return(1/sqrt(2*pi)/sd*exp(-(x-mean)^2/2/sd^2))
}
