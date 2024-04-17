## define path ----
path_dti <- "E:/PostD/Project/UKB/Database/DTI/"
path_Fishoil <- "E:/PostD/Project/PUFAs/DTI2/PUFA/"
path_cov <- "E:/PostD/Project/PUFAs/DTI2/PUFA/"
path_function <- "E:/PostD/Project/PUFAs/DTI2/Function/"
path_output <- "E:/PostD/Project/PUFAs/DTI2/"
## imput Fishoil ----
library(data.table)
Fishoil <- as.data.frame(fread(paste0(path_Fishoil,"fishoil_use.csv")))
Fishoil <- Fishoil[1:2]
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
cov_int2 <- as.data.frame(fread(paste0(path_cov,"fishoil_use.csv")))
cov_int2 <- cov_int2[c(1,3:ncol(cov_int2))]
site_cov <- as.data.frame(fread(paste0(path_dti,"UKB_site_cov.txt")))
cov_site_int2 <- merge(cov_int2,site_cov,by="eid") 
## lm: Fishoil + dti_data + cov_site_int2 ----
Fishoil_dti_data_Cov <- Reduce(function(x,y) merge(x,y,by="eid"),list(Fishoil,dti_data,cov_site_int2),accumulate = F)
## compute ----
source(paste0(path_function,"lm2022c.R"))
Fishoil_dti_data_list <- NULL
for (i in names(Fishoil)[-1]){
  for (j in names(dti_data)[-1]){
    Fishoil_dti_data_list[[i]][[j]] <- lm2022c(Fishoil_dti_data_Cov[i],Fishoil_dti_data_Cov[j],Fishoil_dti_data_Cov)
  }
}
Fishoil_dti_data_Pt <- NULL
for (i in names(Fishoil)[-1]){Fishoil_dti_data_Pt[[i]] <- Reduce(rbind.data.frame,Fishoil_dti_data_list[[i]],accumulate = F)}
Fishoil_dti_data_results <- Reduce(rbind.data.frame,Fishoil_dti_data_Pt,accumulate = F)
Fishoil_dti_results <- merge(Fishoil_dti_data_results,dti_name,by.x = "outcome",by.y = "FieldID")
## output ----
Fishoil_DTI_Results <-Fishoil_dti_results[c("predictor","Description","ES","p","N")]
names(Fishoil_DTI_Results)[2] <- "outcome"
FA_index <- grep("Mean FA",Fishoil_DTI_Results$outcome,ignore.case = F)
MD_index <- grep("Mean MD",Fishoil_DTI_Results$outcome,ignore.case = F)
Fishoil_FA_Results <- Fishoil_DTI_Results[FA_index,]
Fishoil_MD_Results <- Fishoil_DTI_Results[MD_index,]
write.csv(Fishoil_DTI_Results,file = paste0(path_output,"Fishoil_DTI_Results.csv"),row.names = F,quote = F)
write.csv(Fishoil_FA_Results,file = paste0(path_output,"Fishoil_FA_Results.csv"),row.names = F,quote = F)
write.csv(Fishoil_MD_Results,file = paste0(path_output,"Fishoil_MD_Results.csv"),row.names = F,quote = F)


