## define path ----
path_mri <- "E:/PostD/Project/UKB/Database/MRI/"
path_PUFA <- "E:/PostD/Project//PUFAs/MRI2/PUFA/"
path_cov <- "E:/PostD/Project//PUFAs/MRI2/PUFA/"
path_function <- "E:/PostD/Project/PUFAs/MRI2/Function/"
path_output <- "E:/PostD/Project/PUFAs/MRI2/"
## imput PUFA ----
library(data.table)
PUFA <- as.data.frame(fread(paste0(path_PUFA,"PUFA_composition.csv")))
PUFA <- PUFA[1:6]
## imput MRI data ----
library(data.table)
aparc_volume <- as.data.frame(fread(paste0(path_mri,"UKB_2_0_aparc_GrayVol.csv")))
aparc_area <- as.data.frame(fread(paste0(path_mri,"UKB_2_0_aparc_SurfArea.csv")))
aparc_thickness <- as.data.frame(fread(paste0(path_mri,"UKB_2_0_aparc_ThickAvg.csv")))
aseg_volume <- as.data.frame(fread(paste0(path_mri,"UKB_2_0_aseg_Volume.csv")))
source(paste0(path_function,"qc_col2022b.R"))
aseg_volume <- qc_col2022b(aseg_volume)
names(aparc_volume)[1] <- "eid"
names(aparc_area)[1] <- "eid"
names(aparc_thickness)[1] <- "eid"
names(aseg_volume)[1] <- "eid"
## imput covariates ----
library(data.table)
cov_int2 <- as.data.frame(fread(paste0(path_cov,"PUFA_composition.csv")))
cov_int2 <- cov_int2[c(1,7:ncol(cov_int2))]
site_cov <- as.data.frame(fread(paste0(path_mri,"UKB_site_cov.txt")))
cov_site_int2 <- merge(cov_int2,site_cov,by="eid") 
## lm: PUFA + aparc_volume + cov_site_int2 ----
PUFA_aparc_volume_Cov <- Reduce(function(x,y) merge(x,y,by="eid"),list(PUFA,aparc_volume,cov_site_int2),accumulate = F)
## compute ----
source(paste0(path_function,"lm2022d.R"))
PUFA_aparc_volume_list <- NULL
for (i in names(PUFA)[-1]){
  for (j in names(aparc_volume)[-1]){
    PUFA_aparc_volume_list[[i]][[j]] <- lm2022d(PUFA_aparc_volume_Cov[i],PUFA_aparc_volume_Cov[j],PUFA_aparc_volume_Cov)
  }
}
PUFA_aparc_volume_Pt <- NULL
for (i in names(PUFA)[-1]){PUFA_aparc_volume_Pt[[i]] <- Reduce(rbind.data.frame,PUFA_aparc_volume_list[[i]],accumulate = F)}
PUFA_aparc_volume_results <- Reduce(rbind.data.frame,PUFA_aparc_volume_Pt,accumulate = F)
## lm: PUFA + aparc_area + cov_site_int2 ----
PUFA_aparc_area_Cov <- Reduce(function(x,y) merge(x,y,by="eid"),list(PUFA,aparc_area,cov_site_int2),accumulate = F)
## compute ----
source(paste0(path_function,"lm2022d.R"))
PUFA_aparc_area_list <- NULL
for (i in names(PUFA)[-1]){
  for (j in names(aparc_area)[-1]){
    PUFA_aparc_area_list[[i]][[j]] <- lm2022d(PUFA_aparc_area_Cov[i],PUFA_aparc_area_Cov[j],PUFA_aparc_area_Cov)
  }
}
PUFA_aparc_area_Pt <- NULL
for (i in names(PUFA)[-1]){PUFA_aparc_area_Pt[[i]] <- Reduce(rbind.data.frame,PUFA_aparc_area_list[[i]],accumulate = F)}
PUFA_aparc_area_results <- Reduce(rbind.data.frame,PUFA_aparc_area_Pt,accumulate = F)
## lm: PUFA + aparc_thickness + cov_site_int2 ----
PUFA_aparc_thickness_Cov <- Reduce(function(x,y) merge(x,y,by="eid"),list(PUFA,aparc_thickness,cov_site_int2),accumulate = F)
## compute ----
source(paste0(path_function,"lm2022d.R"))
PUFA_aparc_thickness_list <- NULL
for (i in names(PUFA)[-1]){
  for (j in names(aparc_thickness)[-1]){
    PUFA_aparc_thickness_list[[i]][[j]] <- lm2022d(PUFA_aparc_thickness_Cov[i],PUFA_aparc_thickness_Cov[j],PUFA_aparc_thickness_Cov)
  }
}
PUFA_aparc_thickness_Pt <- NULL
for (i in names(PUFA)[-1]){PUFA_aparc_thickness_Pt[[i]] <- Reduce(rbind.data.frame,PUFA_aparc_thickness_list[[i]],accumulate = F)}
PUFA_aparc_thickness_results <- Reduce(rbind.data.frame,PUFA_aparc_thickness_Pt,accumulate = F)
## lm: PUFA + aseg_volume + cov_site_int2 ----
PUFA_aseg_volume_Cov <- Reduce(function(x,y) merge(x,y,by="eid"),list(PUFA,aseg_volume,cov_site_int2),accumulate = F)
## compute ----
source(paste0(path_function,"lm2022d.R"))
PUFA_aseg_volume_list <- NULL
for (i in names(PUFA)[-1]){
  for (j in names(aseg_volume)[-1]){
    PUFA_aseg_volume_list[[i]][[j]] <- lm2022d(PUFA_aseg_volume_Cov[i],PUFA_aseg_volume_Cov[j],PUFA_aseg_volume_Cov)
  }
}
PUFA_aseg_volume_Pt <- NULL
for (i in names(PUFA)[-1]){PUFA_aseg_volume_Pt[[i]] <- Reduce(rbind.data.frame,PUFA_aseg_volume_list[[i]],accumulate = F)}
PUFA_aseg_volume_results <- Reduce(rbind.data.frame,PUFA_aseg_volume_Pt,accumulate = F)
aseg_ROI <- c("Hippocampus","Amygdala","Thalamus","Accumbens","Pallidum","Putamen")
ROI_list <- NULL 
for (i in aseg_ROI){
  ROI_list[[i]] <- grep(i,PUFA_aseg_volume_results$outcome)
}
ROI_index <- Reduce(c,ROI_list,accumulate = F)
ROI_index <- unique(ROI_index)
PUFA_aseg_volume_ROI <- PUFA_aseg_volume_results[ROI_index,]
## results: rbind.data.frame ----
PUFA_volume_Results <- rbind.data.frame(PUFA_aparc_volume_results,PUFA_aseg_volume_ROI)
write.csv(PUFA_volume_Results,file = paste0(path_output,"PUFA_Volume_Results.csv"),row.names = F,quote = F)
write.csv(PUFA_aparc_volume_results,file = paste0(path_output,"PUFA_aparc_volume.csv"),row.names = F,quote = F)
write.csv(PUFA_aparc_area_results,file = paste0(path_output,"PUFA_aparc_area.csv"),row.names = F,quote = F)
write.csv(PUFA_aparc_thickness_results,file = paste0(path_output,"PUFA_aparc_thickness.csv"),row.names = F,quote = F)
write.csv(PUFA_aseg_volume_results,file = paste0(path_output,"PUFA_aseg_volume.csv"),row.names = F,quote = F)
