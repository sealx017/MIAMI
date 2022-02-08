#' @title Compute different mutual information (MI) theoretic measures such as EQMI*, EQMI and CSQMI
#' for every image (subject) in a multiplex imaging dataset
#' @param data is the data.frame of marker intensity observed in different cells of different subjects.
#' The first column is "ID" which corresponds to the subject ID i.e., the subject to which the row
#' (cell) belongs to. The rest of the p columns correspond to p available markers.
#' @param bandwidth is the bandwidth selection procedure user wants to use. The default option "HPI"
#' corresponds to multivariate plug-in bandwidth estimates. For faster computation, use "Ind" which would
#' independently select bandwidth for every r.v. by using Silverman's rule.
#' @param measure corresponds to the MI measure one wants to compute. For example, measure = \code{EQMI_star}
#' will output the measure EQMI_star. Other options are, \code{EQMI} and \code{CSQMI}. We recommend using
#' \code{EQMI_star} for association analysis with clinical outcomes down the line.
#' @param progress_bar if \code{TRUE} will show a progress bar.
#' @return It returns a data.frame whose first column is the vector of subject ID's. The other (p-1) columns
#' respectively correspond to estimated MI between the sets of markers, (1, 2), (1, 2, 3), (1, 2, 3, 4), ..., and
#' (1, 2, ..., p).

#' @export


QMI_all<-function(data, bandwidth = "HPI", measure = "EQMI_star", progress_bar = "True")
  {
  QMI_vec = NULL
  p = dim(data)[2]-1
  k = 1
  unique_IDS = unique(data$ID)
  n = length(unique_IDS)
  if(progress_bar == "True"){
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
       max = n,                      # Maximum value of the progress bar
       style = 3,                    # Progress bar style (also available style = 1 and style = 2)
       width = 50,                   # Progress bar width. Defaults to getOption("width")
       char = "=", initial = 0)}
  for(subject in unique_IDS){
    subject_level_data = data[data$ID==subject,-1]
    EQMI_star = QMI(subject_level_data, bandwidth, measure);
    QMI_vec = rbind(QMI_vec, c(subject, as.vector(EQMI_star[[1]])));
    if(progress_bar == "True"){
    setTxtProgressBar(pb, k);
    k = k+1;}
  }
  QMI_vec = data.frame(QMI_vec)
  names_all = NULL
  if(measure == "EQMI_star"){
  name = "EQMI*_1"}
  else if(measure == "EQMI"){
  name = "EQMI_1"}
  else if(measure == "CSQMI"){
  name = "CSQMI_1"}
  for(i in 2:p){
  names_all = c(names_all,  paste0(name, i))
  name = paste0(name, i)
  }
  colnames(QMI_vec) = c("ID", names_all)
  if(progress_bar == "True"){close(pb)}
  return(QMI_vec)
}
