
if(FALSE){
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(LaplacesDemon))
suppressMessages(library(ggplot2))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(igraph))
suppressMessages(library(PerformanceAnalytics))
suppressMessages(require(survival))
suppressMessages(require(survminer))
suppressMessages(source("/Users/seals/Desktop/CSPH/Jun9/alex_funcs.R"))
suppressMessages(source("/Users/seals/Desktop/CSPH/May4/density_finder.R"))
suppressMessages(source("/Users/seals/Documents/GitHub/MIAMI/R/density_and_EQMI.R"))

my_working_dir = '/Users/seals/Desktop/Debashish_data'
python_path = "/usr/bin/python"
cell_data = read.csv(paste0(my_working_dir,"/TNBC_shareCellData/cellData.csv"))
p_status = read.csv(paste0(my_working_dir,"/TNBC_shareCellData/patient_class.csv"),header = F)
clinical_data = readxl::read_excel(paste0(my_working_dir,"/TNBC_shareCellData/tnbc_clinical.xlsx"),col_names = T,skip = 1)
recur_data = read.csv("/Users/seals/Downloads/rasp-mibi-main/rawdata/clinical_data.csv")
clinical_data = clinical_data[1:39,c("InternalId", "AGE_AT_DX")]
clinical_data$InternalId = as.numeric(clinical_data$InternalId)
colnames(clinical_data) = c("ID", "Age")
marker_names = names(cell_data)[4:52]
cell_data = na.omit(cell_data)
recur_data = merge(recur_data, clinical_data)

IRP_names = c( "HLA.DR", "CD45RO",  "HLA_Class_1", "H3K27me3" , "H3K9ac")
group_covariate_names = c("SampleID")
#selected_data = cell_data[cell_data$Group %in% groups_pruned,c(group_covariate_names, IRP_names)]
selected_data = cell_data[,c(group_covariate_names, IRP_names)]
check = selected_data  %>% count(SampleID)
range_std_data =  selected_data
range_std_data[,IRP_names] = apply(selected_data[,IRP_names],2,range01)
colnames(range_std_data)[1] = "ID"
all_patients = intersect(unique(range_std_data$ID),recur_data$ID)
sel_images = as.numeric(all_patients)
range_std_data = range_std_data[range_std_data$ID %in% sel_images,]
write.csv(range_std_data[,1:6], "/Users/seals/Documents/GitHub/MIAMI/Data/Marker_Data.csv", row.names = F)


recur_data$Recurrence = 1-recur_data$Recurrence
recur_data$Survival = 1-recur_data$Survival
write.csv(recur_data, "/Users/seals/Documents/GitHub/MIAMI/Data/Clinical_Data.csv", row.names = F)


ngrids = 1024
n = dim(range_std_data)[1];
p = dim(range_std_data)[2]-1;

RParray = array(NA, dim = c(n, 2))
QMI_array_ind = QMI_array_hpi_diag  = NULL
x = y = NULL
which_images = 1
rand_size = 20000

Pairwise_EQMI_check = Checking_Pairwise_EQMI(range_std_data, 30000)

png(paste0("/Users/seals/Desktop/CSPH/Mutual_information/MIBI/Plots/",
           "heatmap_of_correlation", ".png"),
    width = 800, height = 800, units = "px", res = 150)
prelim_check = cor(range_std_data[,-1])
print(pheatmap::pheatmap(prelim_check))
#IRP_names = c( "HLA.DR", "CD45RO", "HLA_Class_1", "H3K27me3" , "H3K9ac")
dev.off()

IRP_names = c( "HLA.DR", "CD45RO", "HLA_Class_1", "H3K27me3" , "H3K9ac")

for(images in sel_images){
  list_all = NULL ; m = 1
  for(k in IRP_names){list_all[[m]]  = range_std_data[range_std_data$SampleID==images,k]; m = m+1}
  QMI = QMI_principe_general(list_all, bandwidth = "HPI")
  QMI_array_hpi_diag = rbind(QMI_array_hpi_diag, unlist(QMI[[3]]))
  which_images = which_images + 1
  print(images)
}
saveRDS(QMI_array_hpi_diag,
        file = paste0("/Users/seals/Desktop/CSPH/Mutual_information/MIBI/QMI_array_with_order.RData"))
load("/Users/seals/Desktop/CSPH/Mutual_information/MIBI/QMI_array_with_order.RData")
deg = 1
all_pval_surv = NULL
all_pval_recur = NULL
band_method = c("Hpi_diag")

for(band in band_method){
  if(band == "Hpi_diag"){mymat = QMI_array_hpi_diag}
  else {mymat = QMI_array_ind}
  pval_surv = pval_recur = NULL
  for(k in c(2:5)){
    #for(k in c(6:8)){
    each_qmi = as.data.frame(mymat[,k]); rownames(each_qmi) = sel_images
    transformed = na.omit(each_qmi)
    hier_groups <- as.matrix(ifelse(transformed<mean(transformed[,1]),1,2))
    hier_groups <- cbind(hier_groups, transformed)
    surv_mat <- na.omit(clinical_data[,c("InternalId", "Survival_days_capped*", "Censored","AGE_AT_DX")])

    surv_mat$Censored <- 1+surv_mat$Censored
    new_vec <- cbind(as.matrix(surv_mat$`Survival_days_capped*`),as.matrix(surv_mat$`Censored`),
                     as.matrix(surv_mat$AGE_AT_DX))
    rownames(new_vec) <- surv_mat$InternalId
    new_vec <- as.data.frame(new_vec)
    surv <- as.data.frame(na.omit(cbind(new_vec, hier_groups[match(rownames(new_vec), rownames(hier_groups)),])))
    colnames(surv) <- c("Survival", "Censored","Age","Clus_group", "QMI")
    surv$Clus_group <- as.factor(surv$Clus_group)


    recur_mat <- merge(recur_data, clinical_data, by.x = "ID", by.y = "InternalId")
    recur_mat <- recur_mat[,c(1,2,3,9)]
    recur_mat$Recurrence <- 1-recur_mat$Recurrence
    rownames(recur_mat) <- recur_mat$ID
    recur <- as.data.frame(na.omit(cbind(recur_mat, hier_groups[match(rownames(recur_mat), rownames(hier_groups)),])))
    colnames(recur)[4] <- "Age"
    colnames(recur)[5] <- "Clus_group"
    colnames(recur)[6] <- "QMI"


    res.cox_noQMI<-coxph(Surv(`Survival`,Censored) ~ Age,
                         data = surv)

    res.cox_QMI<-coxph(Surv(`Survival`,Censored) ~ poly(QMI,deg)+Age,
                       data = surv)
    myLRT = -2*(res.cox_noQMI$loglik[2] - res.cox_QMI$loglik[2])
    pval_surv = c(pval_surv,pchisq(myLRT, deg, lower.tail = F))
    res.cox_noQMI<-coxph(Surv(`Recurrence_time`,Recurrence) ~ Age,
                         data = recur)

    res.cox_QMI<-coxph(Surv(`Recurrence_time`,Recurrence) ~ poly(QMI,deg)+Age,
                       data = recur)
    myLRT = -2*(res.cox_noQMI$loglik[2] - res.cox_QMI$loglik[2])
    pval_recur = c(pval_recur,pchisq(myLRT, deg, lower.tail = F))
  }
  pval_surv = cbind(as.data.frame(t(pval_surv)), band)
  pval_recur = cbind(as.data.frame(t(pval_recur)), band)
  all_pval_surv = rbind(all_pval_surv, pval_surv)
  all_pval_recur = rbind(all_pval_recur, pval_recur)
}

print(all_pval_surv)
print(all_pval_recur)
colnames(all_pval) = c("CSQMI", "EQMI", "Bandwidth_selection_method")
}
