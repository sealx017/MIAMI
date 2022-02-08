#' @title Estimate the marginal density of a marker (random variable)
#' @param dat is a vector with intensity data of a single marker or a matrix of intensity of multiple markers.
#' @param ngrids number of grid-points to be used in estimating the kernel density. The bandwidth is chosen
#' using Silverman's rule.
#' @return It returns a list whose first element is the estimated density and the second element is the grid-
#' points where the density is estimated.

#' @export

univ_dens <- function(dat, ngrids = 1024){
  x = as.matrix(dat)
  p = dim(x)[2]
  min_coef = 0
  max_coef = 1
  den = matrix(0, nrow = p, ncol = ngrids)
  den_grid = matrix(0, nrow = p, ncol = ngrids)
  for(i in 1:p){
    vec = x[,i]
    temp_vec = c(vec[!is.na(vec)])
    if(length(temp_vec)==0){
      den[i,] = rep(0,ngrids)
    }
    else{
      if(range(temp_vec)[1]==range(temp_vec)[2]){
        den[i,1] = 1
      }
      else if(length(temp_vec)<5){
        s = density(temp_vec, from = min_coef, to = max_coef,
                    n = ngrids, bw = mean(range(temp_vec)))
        den[i,] = s$y
        den_grid[i,] = s$x
      }
      else{
        s = density(temp_vec, from = min_coef, to = max_coef,
        n = ngrids, bw='nrd0')
        den[i,] = s$y
        den_grid[i,] = s$x
      }}}
  return(list(den, den_grid))
}

#' @title Plot the estimated marginal density of a random variable (marker)
#' @param f is the output of the function univ_dens.
#' @return Plots the estimated marginal density.

#' @export

univ_dens_plot<- function(f)
{
  x = f[[2]]
  y = f[[1]]
  S = plot(x, y, xlab = "x", ylab = "Estimated Density",
           pch = 19, cex = 0.2)
  return(S)
}

#' @title Estimate the joint density of two random variables (markers)
#' @param dat is a matrix of intensity of two random variables.
#' @param ngrids is the number of grid-points to be used in estimating the kernel density. The diagonal
#' multivariate plug-in bandwidth estimator is used.
#' @param maxs is either a character named "default" or a vector of two values which will respectively be
#' the maximum values considered for the two variables. If maxs = "default", the maximum of each of the data
#' vectors are considered.
#' @return It returns a kde object with estimated density and the grid-points from the ks package.

#' @export

biv_dens <- function(dat, ngrids = 512, maxs = "default"){
  x = dat
  if(maxs == "default"){
    maxs = c(max(x[,1]), max(x[,2]))
  }
  else{maxs = maxs}
  min_coef = c(0, 0); max_coef = maxs;
  H = ks::Hpi.diag(x=x)
  fit = ks::kde(x = x, H = H, gridsize = rep(ngrids, dim(x)[2]),
                xmin = min_coef, xmax = max_coef)
  return(fit)
}

#' @title Plot the estimated joint density of two random variables (markers)
#' @param fit is the output of the function biv_dens.
#' @return Plots the estimated bivariate joint density.

#' @export

biv_dens_plot<- function(fit, maxs = c(0.5, 0.5))
{
  q1 = maxs[1]; q2 = maxs[2]
  S = image(fit$eval.points[[1]][which(fit$eval.points[[1]]<=q1)],
            fit$eval.points[[2]][which(fit$eval.points[[2]]<=q2)],
            fit$estimate[which(fit$eval.points[[1]]<=q1),
                         which(fit$eval.points[[2]]<=q2)],
            col = viridis::inferno(50), xlab = "X_1", ylab = "X_2")
  return(S)
}



Checking_Pairwise_EQMI<- function(data, rand_size = 10000, measure = "EQMI_star"){
  n = dim(data)[1];
  p = dim(data)[2]-1;
  prelim_check = matrix(0, p, p)
  sampled = sample(1:n, rand_size)
  for(k1 in 1:(p-1)){
    for(k2 in (k1+1):p){
      list_all = NULL
      list_all[[1]]  = data[sampled,(1+k1)]
      list_all[[2]]  = data[sampled,(1+k2)]
      QMI_comp = QMI(list_all, bandwidth = "HPI",measure = "EQMI_star")
      prelim_check[k1, k2] = QMI_comp[[1]][2]
    }
  }
  prelim_check = prelim_check + t(prelim_check)
  colnames(prelim_check) = rownames(prelim_check) = colnames(data)[-1]
  return(prelim_check)
}
