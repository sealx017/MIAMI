#'@title Fit CoxPH model to study association between survival/recurrence outcome and estimated EQMi*
#' @param clinical_dat is a data.frame whose first column is "ID" denoting the subject IDs. The next two columns
#' respectively correspond to the time (time to death or time to recurrence) and censoring indicator
#' ( = 0 imples death or recurrence).
#' @param cov is a data.frame whose first column is "ID" denoting the subject IDs. The other columns correspond
#' to different covariates.
#' @param EQMI_dat corresponds to the output of the function QMI_of_all_images. It is a matrix of estimated
#' EQMI* of the sets of markers, (1, 2), (1, 2, 3), (1, 2, 3, 4), ..., and (1, 2, ..., p), where p is the
#' number of available markers.
#' @param degree corresponds to the degree of polynomial to use in the CoxPH model. If degree > 1,
#' a non-linear relationship between the clinical outcomes and EQMI* is tested.
#' @return It returns the p-values of the association tests between the clinical outcome and EQMI* of
#' the different sets of markers.
#' @export

Cox_PH<-function(clinical_dat, cov, EQMI_dat, degree = 1)
{
  n_cov = dim(cov)[2]-1
  surv_full = merge(merge(clinical_dat, cov), EQMI_dat)
  colnames(surv_full)[2:3] = c("Censored", "Time")
  p = dim(EQMI_dat)[2]-1
  pval_surv = NULL
  for(i in 1:p){
    surv = surv_full[,c(1:(n_cov+3),(i+4))]
    colnames(surv)[(n_cov+4)] = "QMI"
    res.cox_noQMI<-survival::coxph(survival::Surv(`Time`,Censored) ~ surv[,c(4:(n_cov+3))],
                         data = surv)
    res.cox_QMI<-survival::coxph(survival::Surv(`Time`,Censored) ~ poly(QMI,degree) + surv[,c(4:(n_cov+3))],
                       data = surv)
    myLRT = -2*(res.cox_noQMI$loglik[2] - res.cox_QMI$loglik[2])
    pval_surv = c(pval_surv, pchisq(myLRT, degree, lower.tail = F))
  }
  pval_surv = data.frame(pval_surv)
  colnames(pval_surv) = "p-value"
  pval_surv$Variable = colnames(EQMI_dat)[-1]
  pval_surv = pval_surv[,c(2,1)]
  df = pval_surv
  #sign_formatter <- formatter("span",
  #style = x ~ style(color = ifelse(x <0.01, "green",
  #ifelse(x < 0.05, "green", "black"))))
  #sign_formatter(c(-1, 0, 1))
  #print(formattable(df, list("Variable" = allblack_formatter,
  #                  "p-value" = sign_formatter), align = c('l', 'r')))
  return(pval_surv)
}

