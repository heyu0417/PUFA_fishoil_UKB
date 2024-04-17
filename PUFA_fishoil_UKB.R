#Cox analysis#
Exposures_P1 <- read.csv("Circulating PUFAs, fish oil supplementaiton, and risk of incident dementia/Circulating PUFAs, fish oil supplementation, and risk of incident dementia/Data_new/Exposures_new.csv")
fishoil_new <- Exposures_P1[,c(1,15)]
library(dplyr)
fishoil_new <- mutate(fishoil_new,fishoil=NA)
fishoil_new$fishoil <- 0
fishoil_new$fishoil[which(fishoil_new$X6179.0.0==""|fishoil_new$X6179.0.0=="[-3]")] <- NA
fishoil_new$fishoil[which(fishoil_new$X6179.0.0=="[1, 2, 3, 4, 5, 6]")] <- 1
fishoil_new <- mutate(fishoil_new,mineral=fishoil)
fishoil_new$mineral[which(fishoil_new$fishoil==0)] <- 1
fishoil_new$mineral[which(fishoil_new$X6179.0.0=="[1]")] <- 0
P1_incom <- Dementia_new[,c(1,9,12)]
P1_incom <- merge(P1_incom,P1_covar,by="eid")
P1_incom$oilyfish <- factor(P1_incom$oilyfish)
library(mice)
tail(md.pattern(P1_incom))
x <- P1_incom[,c(1:3)]
imp <- mice(x,seed=1234)
P1_com <- complete(imp,action=3)
P1_com <- cbind(x,P1_com)
write.csv(P1_com,file = "P1_com.csv")
X23456 <- Exposures_P1[,c(1,2,13)]
X23456 <- na.omit(X23456)
X23456 <- mutate(X23456,X23456_z=(X23456.0.0-mean(X23456.0.0))/sd(X23456.0.0))
X23456 <- merge(X23456,P1_com1,by="eid")
library(survival)
cox0 <- coxph(Surv(dementia_years,dementia_status)~X23456_z+X23442.0.0+age+sex+edu_level+APOE4,data = PUFA_Dementia)
cox0 <- coxph(Surv(dementia_years,dementia_status)~X23456_z+X23442.0.0+age+sex+edu_level+APOE4+eth_white+TDI+BMI+smoking+alcohol+TPA+CVD,data = PUFA_Dementia)
cox0 <- coxph(Surv(dementia_years,dementia_status)~X23456_z+X23442.0.0+age+sex+edu_level+APOE4+eth_white+TDI+BMI+smoking+alcohol+TPA+CVD+fruit+vegetable+oilyfish+meat_processed+red_meat+vitamin+mineral,data = PUFA_Dementia)
summary(cox0)
cox.zph(cox0)
library(rms)
library(ggplot2)
dd <- datadist(X23446)
options(datadist="dd")
fit <- cph(Surv(dementia_years,dementia_status)~rcs(X23446.0.0,5)+age+sex+edu_level+APOE4,data = X23446)
fit <- cph(Surv(dementia_years,dementia_status)~rcs(X23446.0.0,5)+age+sex+edu_level+APOE4+eth_white+TDI+BMI+smoking+alcohol+TPA+CVD,data = X23446)
anova(fit)
HR <- Predict(fit,X23446.0.0,fun=exp,ref.zero = TRUE)
P1 <- ggplot(HR)
P1
P2 <- ggplot()+geom_line(data = HR,aes(X23449.0.0,yhat),linetype="solid",size=1,alpha=0.7,colour="#069AF3")+geom_ribbon(data = HR,aes(X23449.0.0,ymin=lower,ymax=upper),alpha=0.1,fill="#069AF3")+geom_hline(yintercept = 1,linetype=2,size=0.5)+theme_classic()+labs(title = "Dementia",x="LA (mmol/L)",y="Hazard Ratio (95% CI)")
P2
fishoil <- fishoil_new[,c(1,3)]
fishoil <- na.omit(fishoil)
fishoil$fishoil <- factor(fishoil$fishoil)
fishoil <- merge(fishoil,P1_com1,by="eid")
cox0 <- coxph(Surv(dementia_years,dementia_status)~fishoil+age+sex+edu_level+APOE4,data = fishoil_Dementia)
cox0 <- coxph(Surv(VD_years,VD_status)~fishoil+age+sex+edu_level+APOE4+eth_white+TDI+BMI+smoking+alcohol+TPA+CVD,data = fishoil_VD)
cox0 <- coxph(Surv(dementia_years,dementia_status)~fishoil+age+sex+edu_level+APOE4+eth_white+TDI+BMI+smoking+alcohol+TPA+CVD+fruit+vegetable+oilyfish+meat_processed+red_meat+vitamin+mineral,data = fishoil_Dementia)

#Table 1#
P1_table1 <- Exposures_P1[,-c(4,9,15)]
x <- fishoil_new[,c(1,3)]
x$fishoil <- factor(x$fishoil)
P1_table1 <- merge(P1_table1,x,by="eid")
P1_table1 <- merge(P1_table1,P1_incom,by="eid")
library(compareGroups)
restab <- descrTable(fishoil~.-eid,data = P1_tableS4)
restab
export2word(restab,file = "Table1.docx")
mean(x1$X23456.0.0,na.rm = TRUE)
sd(x0$BMI,na.rm = TRUE)

#Mediation analysis#
fishoil_mediation <- merge(fishoil_Dementia,PUFA_z,by="eid")
library(mediation)
library(survival)
med.fit <- lm(X23446_z~fishoil+age+sex+edu_level+APOE4,data = fishoil_mediation)
out.fit <- survreg(Surv(dementia_years,dementia_status)~fishoil+X23446_z+age+sex+edu_level+APOE4,data = fishoil_mediation,dist = "exponential")
med.out <- mediate(med.fit,out.fit,treat = "fishoil",mediator = "X23446_z",robustSE = TRUE)
summary(med.out)

#Inflammatory analysis#
inflammatory_markers <- read.csv("Circulating PUFAs, fish oil supplementaiton, and risk of incident dementia/Circulating PUFAs, fish oil supplementation, and risk of incident dementia/Data_new/Inflammatory markers.csv")
inflammatory_markers <- mutate(inflammatory_markers,LMR=Lymphocyte/Monocyte)
inflammatory_markers$LMR[which(inflammatory_markers$Monocyte==0)] <- NA
LMR <- inflammatory_markers[,c(1,9)]
LMR <- LMR[which(LMR$LMR>0),]
LMR <- mutate(LMR,LMR_log=log(LMR))
LMR <- mutate(LMR,LMR_z=(LMR_log-mean(LMR_log))/sd(LMR_log))
hist(LMR_fishoil$LMR_z)
LMR_fishoil <- fishoil[,c(1:3,6:8)]
LMR_fishoil <- merge(LMR,LMR_fishoil,by="eid")
fit <- lm(LMR_z~fishoil+age+sex+edu_level+APOE4,data = LMR_fishoil)
summary(fit)

#Neuroimaging_heatmap#
MD_l <- read.csv("Circulating PUFAs, fish oil supplementaiton, and risk of incident dementia/Circulating PUFAs, fish oil supplementation, and risk of incident dementia/Data_new/Heatmap/MD_left.csv",header = T,row.names = 1)
MD_r <- read.csv("Circulating PUFAs, fish oil supplementaiton, and risk of incident dementia/Circulating PUFAs, fish oil supplementation, and risk of incident dementia/Data_new/Heatmap/MD_right.csv",header = T,row.names = 1)
library(ComplexHeatmap)
a <- as.matrix(MD_l)
b <- as.matrix(MD_r)
library(circlize)
ht_r <- Heatmap(b,col = colorRamp2(c(0.001,0.05,0.10),c("#E50000", "#FFFFFF", "#0343DF")),rect_gp = gpar(col="white",lty=1,lwd=1.5),cluster_columns = F,cluster_rows = F,column_labels = c("Total PUFAs","Omega-6 PUFAs","Linoleic Acid"),column_names_rot = 45,row_names_gp = gpar(fontsize=10),column_names_gp = gpar(fontsize=10),width = unit(2,"cm"),height = unit(8,"cm"),heatmap_legend_param = list(title = "Bonferroni-corrected P value",at=c(0.001,0.05,0.10),labels = c("0.001","0.050","0.100"),legend_height = unit(4,"cm")))
ht_l <- Heatmap(a,col = colorRamp2(c(0.001,0.05,0.10),c("#E50000", "#FFFFFF", "#0343DF")),rect_gp = gpar(col="white",lty=1,lwd=1.5),cluster_columns = F,cluster_rows = F,column_labels = c("Linoleic Acid","Omega-6 PUFAs","Total PUFAs"),column_names_rot = 45,row_names_side = "left",row_names_gp = gpar(fontsize=10),column_names_gp = gpar(fontsize=10),width = unit(2,"cm"),height = unit(8,"cm"),show_heatmap_legend = F)
ht = ht_l + ht_r
draw(ht,ht_gap=unit(1,"mm"))

#Subgroup analysis#
forestplot_fishoil <- read.csv("Circulating PUFAs, fish oil supplementaiton, and risk of incident dementia/Circulating PUFAs, fish oil supplementation, and risk of incident dementia/Data_new/fishoil_stratified.csv",header = FALSE)
library(forestplot)
attach(forestplot_fishoil)
forestplot(labeltext=as.matrix(forestplot_fishoil[,1:4]),mean=forestplot_fishoil$V5,lower=forestplot_fishoil$V6,upper=forestplot_fishoil$V7,is.summary=c(T,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F),graph.pos=3,hrzl_lines=list("2"=gpar(lwd=1.5),"5"=gpar(lty=5,lwd=1),"8"=gpar(lty=5,lwd=1),"11"=gpar(lty=5,lwd=1),"14"=gpar(lty=5,lwd=1),"17"=gpar(lty=5,lwd=1),"20"=gpar(lty=5,lwd=1),"23"=gpar(lty=5,lwd=1)),clip=c(0.7,1.2),xlab="Hazard Ratio (HR)",zero=1,graphwidth=unit(100,"mm"),colgap=unit(5,"mm"),lineheight="auto",line.margin=unit(10,"mm"),col=fpColors(box = "#069AF3",lines = "#069AF3",zero = "gray"),xlog=FALSE,xticks=c(0.7,0.8,0.9,1,1.1,1.2),lwd.xaxis=2,lwd.zero=2,lwd.ci=2,ci.vertices= TRUE,ci.vertices.height=0.2,boxsize=0.3)

#Forestplot#
forestplot_main <- read.csv("Forestplot.csv",header = FALSE)
library(forestplot)
attach(forestplot_main)
forestplot(labeltext=as.matrix(forestplot_main[,1:3]),mean=forestplot_main$V4,lower=forestplot_main$V5,upper=forestplot_main$V6,is.summary=c(T,T,F,F,F,F,F,T,F,F,F,F,F,T,F),graph.pos=2,hrzl_lines=list("2"=gpar(lwd=1.5),"8"=gpar(lty=5,lwd=1),"14"=gpar(lty=5,lwd=1)),clip=c(0.8,1.1),xlab="Hazard Ratio (HR)",zero=1,graphwidth=unit(100,"mm"),colgap=unit(10,"mm"),lineheight="auto",line.margin=unit(20,"mm"),col=fpColors(box = "#069AF3",lines = "#069AF3",zero = "gray"),xlog=FALSE,xticks=c(0.8,0.9,1,1.1),lwd.xaxis=2,lwd.zero=2,lwd.ci=2,ci.vertices= TRUE,ci.vertices.height=0.2,boxsize=0.3)
