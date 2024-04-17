## define path ----
path_mri <- "E:/PostD/Project/UKB/Database/MRI/"
path_PUFA <- "E:/PostD/Project//PUFAs/MRI2/PUFA/"
path_cov <- "E:/PostD/Project//PUFAs/MRI2/PUFA/"
path_function <- "E:/PostD/Project/PUFAs/MRI2/Function/"
path_output <- "E:/PostD/Project/PUFAs/MRI2/"
## imput Fishoil ----
library(data.table)
Fishoil <- as.data.frame(fread(paste0(path_PUFA,"fishoil_use.csv")))
Fishoil <- Fishoil[1:2]
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
cov_int2 <- as.data.frame(fread(paste0(path_cov,"fishoil_use.csv")))
cov_int2 <- cov_int2[c(1,3:ncol(cov_int2))]
site_cov <- as.data.frame(fread(paste0(path_mri,"UKB_site_cov.txt")))
cov_site_int2 <- merge(cov_int2,site_cov,by="eid") 
## lm: Fishoil + aparc_volume + cov_site_int2 ----
Fishoil_aparc_volume_Cov <- Reduce(function(x,y) merge(x,y,by="eid"),list(Fishoil,aparc_volume,cov_site_int2),accumulate = F)
## compute ----
source(paste0(path_function,"lm2022c.R"))
Fishoil_aparc_volume_list <- NULL
for (i in names(Fishoil)[-1]){
  for (j in names(aparc_volume)[-1]){
    Fishoil_aparc_volume_list[[i]][[j]] <- lm2022c(Fishoil_aparc_volume_Cov[i],Fishoil_aparc_volume_Cov[j],Fishoil_aparc_volume_Cov)
  }
}
Fishoil_aparc_volume_Pt <- NULL
for (i in names(Fishoil)[-1]){Fishoil_aparc_volume_Pt[[i]] <- Reduce(rbind.data.frame,Fishoil_aparc_volume_list[[i]],accumulate = F)}
Fishoil_aparc_volume_results <- Reduce(rbind.data.frame,Fishoil_aparc_volume_Pt,accumulate = F)
## lm: Fishoil + aparc_area + cov_site_int2 ----
Fishoil_aparc_area_Cov <- Reduce(function(x,y) merge(x,y,by="eid"),list(Fishoil,aparc_area,cov_site_int2),accumulate = F)
## compute ----
source(paste0(path_function,"lm2022c.R"))
Fishoil_aparc_area_list <- NULL
for (i in names(Fishoil)[-1]){
  for (j in names(aparc_area)[-1]){
    Fishoil_aparc_area_list[[i]][[j]] <- lm2022c(Fishoil_aparc_area_Cov[i],Fishoil_aparc_area_Cov[j],Fishoil_aparc_area_Cov)
  }
}
Fishoil_aparc_area_Pt <- NULL
for (i in names(Fishoil)[-1]){Fishoil_aparc_area_Pt[[i]] <- Reduce(rbind.data.frame,Fishoil_aparc_area_list[[i]],accumulate = F)}
Fishoil_aparc_area_results <- Reduce(rbind.data.frame,Fishoil_aparc_area_Pt,accumulate = F)
## lm: Fishoil + aparc_thickness + cov_site_int2 ----
Fishoil_aparc_thickness_Cov <- Reduce(function(x,y) merge(x,y,by="eid"),list(Fishoil,aparc_thickness,cov_site_int2),accumulate = F)
## compute ----
source(paste0(path_function,"lm2022c.R"))
Fishoil_aparc_thickness_list <- NULL
for (i in names(Fishoil)[-1]){
  for (j in names(aparc_thickness)[-1]){
    Fishoil_aparc_thickness_list[[i]][[j]] <- lm2022c(Fishoil_aparc_thickness_Cov[i],Fishoil_aparc_thickness_Cov[j],Fishoil_aparc_thickness_Cov)
  }
}
Fishoil_aparc_thickness_Pt <- NULL
for (i in names(Fishoil)[-1]){Fishoil_aparc_thickness_Pt[[i]] <- Reduce(rbind.data.frame,Fishoil_aparc_thickness_list[[i]],accumulate = F)}
Fishoil_aparc_thickness_results <- Reduce(rbind.data.frame,Fishoil_aparc_thickness_Pt,accumulate = F)
## lm: Fishoil + aseg_volume + cov_site_int2 ----
Fishoil_aseg_volume_Cov <- Reduce(function(x,y) merge(x,y,by="eid"),list(Fishoil,aseg_volume,cov_site_int2),accumulate = F)
## compute ----
source(paste0(path_function,"lm2022c.R"))
Fishoil_aseg_volume_list <- NULL
for (i in names(Fishoil)[-1]){
  for (j in names(aseg_volume)[-1]){
    Fishoil_aseg_volume_list[[i]][[j]] <- lm2022c(Fishoil_aseg_volume_Cov[i],Fishoil_aseg_volume_Cov[j],Fishoil_aseg_volume_Cov)
  }
}
Fishoil_aseg_volume_Pt <- NULL
for (i in names(Fishoil)[-1]){Fishoil_aseg_volume_Pt[[i]] <- Reduce(rbind.data.frame,Fishoil_aseg_volume_list[[i]],accumulate = F)}
Fishoil_aseg_volume_results <- Reduce(rbind.data.frame,Fishoil_aseg_volume_Pt,accumulate = F)
aseg_ROI <- c("Hippocampus","Amygdala","Thalamus","Accumbens","Pallidum","Putamen")
ROI_list <- NULL 
for (i in aseg_ROI){
  ROI_list[[i]] <- grep(i,Fishoil_aseg_volume_results$outcome)
}
ROI_index <- Reduce(c,ROI_list,accumulate = F)
ROI_index <- unique(ROI_index)
Fishoil_aseg_volume_ROI <- Fishoil_aseg_volume_results[ROI_index,]
## results: rbind.data.frame ----
Fishoil_volume_Results <- rbind.data.frame(Fishoil_aparc_volume_results,Fishoil_aseg_volume_ROI)
write.csv(Fishoil_volume_Results,file = paste0(path_output,"Fishoil_Volume_Results.csv"),row.names = F,quote = F)
write.csv(Fishoil_aparc_volume_results,file = paste0(path_output,"Fishoil_aparc_volume.csv"),row.names = F,quote = F)
write.csv(Fishoil_aparc_area_results,file = paste0(path_output,"Fishoil_aparc_area.csv"),row.names = F,quote = F)
write.csv(Fishoil_aparc_thickness_results,file = paste0(path_output,"Fishoil_aparc_thickness.csv"),row.names = F,quote = F)
write.csv(Fishoil_aseg_volume_results,file = paste0(path_output,"Fishoil_aseg_volume.csv"),row.names = F,quote = F)
