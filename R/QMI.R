#' @title Compute different mutual information (MI) theoretic measures such as EQMI*, EQMI and CSQMI
#' in a general dataset
#' @param dat is the data matrix whose p columns correspond to p different random variables (markers in a
#' multiplex imaging dataset) and whose rows are samples (usually, different cells from different subjects in a
#' multiplex imaging dataset).
#' @param bandwidth is the bandwidth selection procedure user wants to use. The default option "HPI"
#' corresponds to multivariate plug-in bandwidth estimates. For faster computation, use "Ind" which would
#' independently select bandwidth for every r.v. by using Silverman's rule.
#' @param measure corresponds to the MI measure one wants to compute. For example, measure = "EQMI_star"
#' will output the measure EQMI_star. Other options are, "EQMI, "CSQMI" and "ALL", where the last one
#' outputs all the three measures.
#' @param var_names can be \code{TRUE} (T) or \code{FALSE} (F). If it is True the output list will have name of the
#' measures.
#' @return It returns a list whose first element is the vector of estimated MI between the following combination
#' of the markers, (1, 2), (1, 2, 3), (1, 2, 3, 4), ..., and (1, 2, ..., p), and the second element is the
#' vector of bandiwidth parameters used. If measure = "All", the first element is a list of three elements
#' corresponding to \code{EQMI_star}, \code{EQMI} and \code{CSQMI}.

#' @export

QMI <- function(dat, bandwidth = "HPI", measure = "EQMI_star", var_names = F){
  n_markers = p = dim(dat)[2]
  min_coef = apply(dat, 2, min)
  max_coef = apply(dat, 2, max)
  mat = dat+.Machine$double.xmin
  n = dim(mat)[1]
  if(bandwidth=="HPI"){
    H =  ks::Hpi.diag(x=mat)
    for(k in 1:n_markers){
      assign(paste0("h_",k), H[k,k])
    }
  }
  else if(bandwidth=="Ind"){
    for(k in 1:n_markers){
      assign(paste0("h_",k),(1.06*n^(-0.2)*sd(mat[,k]))^2)
    }
  }
  V_J = 1/n^2; V_C = 1/n; V_M = 1
  I_array_ED = I_array_ED_raw = I_array_CS = matrix(0,n_markers,1)
  for(k in 1:n_markers){
    S = matrix(0,n,n)
    for (i in 1:(n-1)){
      S[i,(i+1):n] = dnorm_mine(mat[(i+1):n,k]-mat[i,k], mean = 0, sd = sqrt(2*get(paste0("h_",k))))
    }
    #fastComp = mclapply(1:(n-1), function(i){dnorm_mine(mat[(i+1):n,k]-mat[i,k], mean = 0,
    #sd = sqrt(2*get(paste0("h_",k))))}, mc.cores=detectCores()-1, mc.silent = T)
    #for (i in 1:(n-1)){
    #  S[i,(i+1):n] = fastComp[[i]]
    #}
    S = S + t(S)
    diag(S) = dnorm_mine(0, mean = 0, sd = sqrt(2*get(paste0("h_",k))))
    S_vec = 1/n*rowSums(S)
    S_cons = 1/n*sum(S_vec)
    V_J = V_J*S
    V_C = V_C*S_vec
    V_M = V_M*S_cons
    I_array_ED[k,1] = (sum(V_J) - 2*sum(V_C) + V_M)/(sum(V_J) +  V_M)
    I_array_ED_raw[k,1] = (sum(V_J) - 2*sum(V_C) + V_M)
    I_array_CS[k,1] = log(sum(V_J)) - 2*log(sum(V_C)) + log(V_M)
  }
  h_s = NULL

  for(k in 1:n_markers){
    h_s = c(h_s, get(paste0("h_",k)))
  }
  if(var_names == T){h_s = as.data.frame(t(h_s)); colnames(h_s) = paste0("h", 1:n_markers)}

  #---Setting names of EQMI_star-----

  names_all = NULL
  name = "EQMI*_1"
  for(i in 2:p){
    names_all = c(names_all,  paste0(name, i))
    name = paste0(name, i)
  }
  res1 = as.matrix(I_array_ED[-1,])
  if(var_names == T){rownames(res1) = names_all; colnames(res1) = "Estimate"}

  #---Setting names of EQMI-----

  names_all = NULL
  name = "EQMI_1"
  for(i in 2:p){
    names_all = c(names_all,  paste0(name, i))
    name = paste0(name, i)
  }
  res2 = as.matrix(I_array_ED_raw[-1,])
  if(var_names == T){rownames(res2) = names_all; colnames(res2) = "Estimate"}

  #---Setting names of CSQMI-----
  names_all = NULL
  name = "CSQMI_1"
  for(i in 2:p){
    names_all = c(names_all,  paste0(name, i))
    name = paste0(name, i)
  }
  res3 = as.matrix(I_array_CS[-1,])
  if(var_names == T){rownames(res3) = names_all; colnames(res3) = "Estimate"}

  #--------------------------------

  if(measure == "EQMI_star"){
    final_list = list(EQMI_star = res1, Bandwidth_parameters = h_s)
  }

  else if (measure == "EQMI"){
    final_list = list(EQMI = res2, Bandwidth_parameters = h_s)
  }

  else if (measure == "CSQMI"){
    final_list = list(CSQMI = res3, Bandwidth_parameters = h_s)
  }

  else{
    final_list = list(EQMI_star = res1, EQMI = res2,
                      CSQMI = res3, Bandwidth_parameters = h_s)
  }
  return(final_list)
}
