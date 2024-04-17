## define path ----
path_dti <- "E:/PostD/Project/UKB/Database/DTI/"
path_PUFA <- "E:/PostD/Project//PUFAs/DTI2/PUFA/"
path_cov <- "E:/PostD/Project//PUFAs/DTI2/PUFA/"
path_function <- "E:/PostD/Project/PUFAs/DTI2/Function/"
path_output <- "E:/PostD/Project//PUFAs/DTI2/"
## imput PUFA ----
library(data.table)
PUFA <- as.data.frame(fread(paste0(path_PUFA,"PUFA_concentration.csv")))
PUFA <- PUFA[1:6]
## imput DTI data ----
library(data.table)
dti_data <- as.data.frame(fread(paste0(path_dti,"UKB_6_DTI.txt")))
dti_name <- as.data.frame(fread(paste0(path_dti,"UKB_DTI_name.txt")))
dti_name[1] <- paste0("X",dti_name$FieldID)
names(dti_data) <- c("eid",dti_name$FieldID)
FA_index <- 1+grep("mean FA",dti_name$Description,ignore.case = T)
MD_index <- 1+grep("mean MD",dti_name$Description,ignore.case = T)
dti_data <- dti_data[c(1,FA_index,MD_index)]
## imput covariates ----
library(data.table)
cov_int2 <- as.data.frame(fread(paste0(path_cov,"PUFA_concentration.csv")))
cov_int2 <- cov_int2[c(1,7:ncol(cov_int2))]
site_cov <- as.data.frame(fread(paste0(path_dti,"UKB_site_cov.txt")))
cov_site_int2 <- merge(cov_int2,site_cov,by="eid") 
## lm: PUFA + dti_data + cov_site_int2 ----
PUFA_dti_data_Cov <- Reduce(function(x,y) merge(x,y,by="eid"),list(PUFA,dti_data,cov_site_int2),accumulate = F)
## compute ----
source(paste0(path_function,"lm2022c.R"))
PUFA_dti_data_list <- NULL
for (i in names(PUFA)[-1]){
  for (j in names(dti_data)[-1]){
    PUFA_dti_data_list[[i]][[j]] <- lm2022c(PUFA_dti_data_Cov[i],PUFA_dti_data_Cov[j],PUFA_dti_data_Cov)
  }
}
PUFA_dti_data_Pt <- NULL
for (i in names(PUFA)[-1]){PUFA_dti_data_Pt[[i]] <- Reduce(rbind.data.frame,PUFA_dti_data_list[[i]],accumulate = F)}
PUFA_dti_data_results <- Reduce(rbind.data.frame,PUFA_dti_data_Pt,accumulate = F)
PUFA_dti_results <- merge(PUFA_dti_data_results,dti_name,by.x = "outcome",by.y = "FieldID")
## output ----
PUFA_DTI_Results <-PUFA_dti_results[c("predictor","Description","ES","p","N")]
names(PUFA_DTI_Results)[2] <- "outcome"
FA_index <- grep("Mean FA",PUFA_DTI_Results$outcome,ignore.case = F)
MD_index <- grep("Mean MD",PUFA_DTI_Results$outcome,ignore.case = F)
PUFA_FA_Results <- PUFA_DTI_Results[FA_index,]
PUFA_MD_Results <- PUFA_DTI_Results[MD_index,]
write.csv(PUFA_DTI_Results,file = paste0(path_output,"PUFA_DTI_Results.csv"),row.names = F,quote = F)
write.csv(PUFA_FA_Results,file = paste0(path_output,"PUFA_FA_Results.csv"),row.names = F,quote = F)
write.csv(PUFA_MD_Results,file = paste0(path_output,"PUFA_MD_Results.csv"),row.names = F,quote = F)


