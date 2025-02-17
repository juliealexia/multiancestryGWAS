library(data.table)
library(dplyr)

meta_analysis <- function(trait){
  AFR_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/plink2_step2_firth_out_AFR_",trait,".regenie")))
  AFR_res$P <- 10^(-AFR_res$LOG10P)
  AFR_temp <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_AFR.out.",trait,".glm.linear")))
  AFR_temp$BETA <- NULL
  AFR_temp$SE <- NULL
  AFR_temp$T_STAT <- NULL
  AFR_temp$P <- NULL
  AFR_res <- left_join(AFR_temp,AFR_res[,c("ID","BETA","SE","P")],by="ID")
 
  AMR_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/plink2_step2_firth_out_AMR_",trait,".regenie")))
  AMR_res$P <- 10^(-AMR_res$LOG10P)
  AMR_temp <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_AMR.out.",trait,".glm.linear")))
  AMR_temp$BETA <- NULL
  AMR_temp$SE <- NULL
  AMR_temp$T_STAT <- NULL
  AMR_temp$P <- NULL
  AMR_res <- left_join(AMR_temp,AMR_res[,c("ID","BETA","SE","P")],by="ID")

  EAS_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/plink2_step2_firth_out_EAS_",trait,".regenie")))
  EAS_res$P <- 10^(-EAS_res$LOG10P)
  EAS_temp <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_EAS.out.",trait,".glm.linear")))
  EAS_temp$BETA <- NULL
  EAS_temp$SE <- NULL
  EAS_temp$T_STAT <- NULL
  EAS_temp$P <- NULL
  EAS_res <- left_join(EAS_temp,EAS_res[,c("ID","BETA","SE","P")],by="ID")

  EUR_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/plink2_step2_firth_out_EUR_",trait,".regenie")))
  EUR_res$P <- 10^(-EUR_res$LOG10P)
  EUR_temp <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_EUR.out.",trait,".glm.linear")))
  EUR_temp$BETA <- NULL
  EUR_temp$SE <- NULL
  EUR_temp$T_STAT <- NULL
  EUR_temp$P <- NULL
  EUR_res <- left_join(EUR_temp,EUR_res[,c("ID","BETA","SE","P")],by="ID")

  SAS_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/plink2_step2_firth_out_SAS_",trait,".regenie")))
  SAS_res$P <- 10^(-SAS_res$LOG10P)
  SAS_temp <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_SAS.out.",trait,".glm.linear")))
  SAS_temp$BETA <- NULL
  SAS_temp$SE <- NULL
  SAS_temp$T_STAT <- NULL
  SAS_temp$P <- NULL
  SAS_res <- left_join(SAS_temp,SAS_res[,c("ID","BETA","SE","P")],by="ID")
 
  ALL_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/plink2_step2_firth_out_ALL_",trait,".regenie")))
  colnames(ALL_res)[7] <- "OBS_CT"
  
  met_AFR <- data.frame("ID"=AFR_res$ID, "N_AFR"=AFR_res$OBS_CT, "BETA_AFR"=AFR_res$BETA, "VAR_AFR"=(AFR_res$SE)^2)   
  met_AMR <- data.frame("ID"=AMR_res$ID, "N_AMR"=AMR_res$OBS_CT, "BETA_AMR"=AMR_res$BETA, "VAR_AMR"=(AMR_res$SE)^2)
  met_EUR <- data.frame("ID"=EUR_res$ID, "N_EUR"=EUR_res$OBS_CT, "BETA_EUR"=EUR_res$BETA, "VAR_EUR"=(EUR_res$SE)^2)   
  met_EAS <- data.frame("ID"=EAS_res$ID, "N_EAS"=EAS_res$OBS_CT, "BETA_EAS"=EAS_res$BETA, "VAR_EAS"=(EAS_res$SE)^2)   
  met_SAS <- data.frame("ID"=SAS_res$ID, "N_SAS"=SAS_res$OBS_CT, "BETA_SAS"=SAS_res$BETA, "VAR_SAS"=(SAS_res$SE)^2)   
  
  all_SNPs1 <- union(met_AFR$ID, met_EUR$ID)
  all_SNPs2 <- union(all_SNPs1, met_AMR$ID)  
  all_SNPs3 <- union(all_SNPs2, met_EAS$ID)  
  all_SNPs <- union(all_SNPs3, met_SAS$ID)  
  length(all_SNPs)    
  
  meta_tab <- data.frame("ID"=all_SNPs)  
  meta_tab <- left_join(meta_tab, met_AFR, by = "ID")  
  meta_tab <- left_join(meta_tab, met_EUR, by = "ID")  
  meta_tab <- left_join(meta_tab, met_AMR, by = "ID")
  meta_tab <- left_join(meta_tab, met_EAS, by = "ID")  
  meta_tab <- left_join(meta_tab, met_SAS, by = "ID")  
  
  temp_afr <- meta_tab$BETA_AFR/meta_tab$VAR_AFR  
  temp_eur <- meta_tab$BETA_EUR/meta_tab$VAR_EUR  
  temp_amr <- meta_tab$BETA_AMR/meta_tab$VAR_AMR
  temp_eas <- meta_tab$BETA_EAS/meta_tab$VAR_EAS 
  temp_sas <- meta_tab$BETA_SAS/meta_tab$VAR_SAS  
  
  temp_afr2 <- 1/meta_tab$VAR_AFR
  temp_eur2 <- 1/meta_tab$VAR_EUR  
  temp_amr2 <- 1/meta_tab$VAR_AMR
  temp_eas2 <- 1/meta_tab$VAR_EAS
  temp_sas2 <- 1/meta_tab$VAR_SAS
  
  temp_beta <- rowSums( cbind (temp_afr,temp_eur,temp_amr,temp_eas,temp_sas), na.rm=TRUE)    
  temp_var <- rowSums( cbind (temp_afr2,temp_eur2,temp_amr2,temp_eas2,temp_sas2), na.rm=TRUE)    
  
  meta_tab$BETA_ALL <- temp_beta/temp_var  
  meta_tab$SE_ALL <- sqrt(1/temp_var)  
  meta_tab$T_STAT_ALL <- meta_tab$BETA_ALL/meta_tab$SE_ALL  
  meta_tab$P <- 2*pnorm(abs(meta_tab$T_STAT_ALL), 0, 1, lower = F)
  meta_tab$OBS_CT <- rowSums( meta_tab[,c("N_AFR","N_EUR","N_AMR","N_EAS","N_SAS")], na.rm=TRUE)   
  
  meta_tab <- left_join(meta_tab[,c("ID","T_STAT_ALL","P","OBS_CT")], ALL_res[,c("ID","A1FREQ","CHROM","GENPOS")], by = "ID")
  colnames(meta_tab)[5] <- "A1_FREQ"
  colnames(meta_tab)[6] <- "CHR"
  colnames(meta_tab)[7] <- "POS"
  
  ALL_res <- ALL_res[ALL_res$ID %in% meta_tab$ID,]
  ALL_res$P <- 10^(-ALL_res$LOG10P)
  meta_tab$P <- as.numeric(meta_tab$P)
  ALL_res$EXTRA<-NULL
  colnames(ALL_res)[1] <- "CHR"
  colnames(ALL_res)[2] <- "POS"
  colnames(ALL_res)[6] <- "A1_FREQ"
  ALL_res$TEST <- NULL
  ALL_res$ALLELE0 <- NULL
  ALL_res$ALLELE1 <- NULL
  ALL_res$BETA <- NULL
  ALL_res$SE <- NULL
  ALL_res$CHISQ <- NULL
  
  mrmegares <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/mrmega_",trait,"_UKB.rg.result")))
  mrmegares <- mrmegares[is.na(mrmegares$Comments),]
  ndf <- as.integer(unique(mrmegares$ndf_association))
  mrmegares <- mrmegares[,c("MarkerName","Chromosome","Position","EAF","Nsample","chisq_association")]
  mrmegares$P <- pchisq(mrmegares$chisq_association,ndf, lower.tail = F)
  colnames(mrmegares)[1] <- "ID"
  colnames(mrmegares)[2] <- "CHR"
  colnames(mrmegares)[3] <- "POS"
  colnames(mrmegares)[4] <- "A1_FREQ"
  colnames(mrmegares)[5] <- "OBS_CT"
  
  fwrite(mrmegares,file=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/totransfer/regenie_",trait,"_MRMEGA.txt"), sep="\t")
  fwrite(ALL_res,file=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/totransfer/regenie_",trait,"_MEGA.txt"), sep="\t")
  fwrite(meta_tab,file=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/totransfer/regenie_",trait,"_META.txt"), sep="\t")
}

meta_analysis_bin <- function(trait){
  AFR_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/plink2_step2_firth_out_AFR_",trait,".regenie")))

  AMR_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/plink2_step2_firth_out_AMR_",trait,".regenie")))

  EAS_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/plink2_step2_firth_out_EAS_",trait,".regenie")))

  EUR_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/plink2_step2_firth_out_EUR_",trait,".regenie")))

  SAS_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/plink2_step2_firth_out_SAS_",trait,".regenie")))

  ALL_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/plink2_step2_firth_out_ALL_",trait,".regenie")))
  colnames(ALL_res)[7] <- "OBS_CT"
  
  met_AFR <- data.frame("ID"=AFR_res$ID, "N_AFR"=AFR_res$N, "BETA_AFR"=AFR_res$BETA, "VAR_AFR"=(AFR_res$SE)^2)   
  met_AMR <- data.frame("ID"=AMR_res$ID, "N_AMR"=AMR_res$N, "BETA_AMR"=AMR_res$BETA, "VAR_AMR"=(AMR_res$SE)^2)
  met_EUR <- data.frame("ID"=EUR_res$ID, "N_EUR"=EUR_res$N, "BETA_EUR"=EUR_res$BETA, "VAR_EUR"=(EUR_res$SE)^2)   
  met_EAS <- data.frame("ID"=EAS_res$ID, "N_EAS"=EAS_res$N, "BETA_EAS"=EAS_res$BETA, "VAR_EAS"=(EAS_res$SE)^2)   
  met_SAS <- data.frame("ID"=SAS_res$ID, "N_SAS"=SAS_res$N, "BETA_SAS"=SAS_res$BETA, "VAR_SAS"=(SAS_res$SE)^2)   
  
  all_SNPs1 <- union(met_AFR$ID, met_EUR$ID)
  all_SNPs2 <- union(all_SNPs1, met_AMR$ID)  
  all_SNPs3 <- union(all_SNPs2, met_EAS$ID)  
  all_SNPs <- union(all_SNPs3, met_SAS$ID)  
  length(all_SNPs)    
  
  meta_tab <- data.frame("ID"=all_SNPs)  
  meta_tab <- left_join(meta_tab, met_AFR, by = "ID")  
  meta_tab <- left_join(meta_tab, met_EUR, by = "ID")  
  meta_tab <- left_join(meta_tab, met_AMR, by = "ID")
  meta_tab <- left_join(meta_tab, met_EAS, by = "ID")  
  meta_tab <- left_join(meta_tab, met_SAS, by = "ID")  
  
  temp_afr <- meta_tab$BETA_AFR/meta_tab$VAR_AFR  
  temp_eur <- meta_tab$BETA_EUR/meta_tab$VAR_EUR  
  temp_amr <- meta_tab$BETA_AMR/meta_tab$VAR_AMR
  temp_eas <- meta_tab$BETA_EAS/meta_tab$VAR_EAS 
  temp_sas <- meta_tab$BETA_SAS/meta_tab$VAR_SAS  
  
  temp_afr2 <- 1/meta_tab$VAR_AFR
  temp_eur2 <- 1/meta_tab$VAR_EUR  
  temp_amr2 <- 1/meta_tab$VAR_AMR
  temp_eas2 <- 1/meta_tab$VAR_EAS
  temp_sas2 <- 1/meta_tab$VAR_SAS
  
  temp_beta <- rowSums( cbind (temp_afr,temp_eur,temp_amr,temp_eas,temp_sas), na.rm=TRUE)    
  temp_var <- rowSums( cbind (temp_afr2,temp_eur2,temp_amr2,temp_eas2,temp_sas2), na.rm=TRUE)    
  
  meta_tab$BETA_ALL <- temp_beta/temp_var  
  meta_tab$SE_ALL <- sqrt(1/temp_var)  
  meta_tab$T_STAT_ALL <- meta_tab$BETA_ALL/meta_tab$SE_ALL  
  meta_tab$P <- 2*pnorm(abs(meta_tab$T_STAT_ALL), 0, 1, lower = F)
  meta_tab$OBS_CT <- rowSums( meta_tab[,c("N_AFR","N_EUR","N_AMR","N_EAS","N_SAS")], na.rm=TRUE)   
  
  meta_tab <- left_join(meta_tab[,c("ID","T_STAT_ALL","P","OBS_CT")], ALL_res[,c("ID","A1FREQ","CHROM","GENPOS")], by = "ID")
  colnames(meta_tab)[5] <- "A1_FREQ"
  colnames(meta_tab)[6] <- "CHR"
  colnames(meta_tab)[7] <- "POS"
  
  ALL_res <- ALL_res[ALL_res$ID %in% meta_tab$ID,]
  ALL_res$P <- 10^(-ALL_res$LOG10P)
  meta_tab$P <- as.numeric(meta_tab$P)
  ALL_res$EXTRA<-NULL
  colnames(ALL_res)[1] <- "CHR"
  colnames(ALL_res)[2] <- "POS"
  colnames(ALL_res)[6] <- "A1_FREQ"
  ALL_res$TEST <- NULL
  ALL_res$ALLELE0 <- NULL
  ALL_res$ALLELE1 <- NULL
  ALL_res$BETA <- NULL
  ALL_res$SE <- NULL
  ALL_res$CHISQ <- NULL
  
  mrmegares <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/mrmega_",trait,"_UKB.rg.result")))
  mrmegares <- mrmegares[is.na(mrmegares$Comments),]
  ndf <- as.integer(unique(mrmegares$ndf_association))
  mrmegares <- mrmegares[,c("MarkerName","Chromosome","Position","EAF","Nsample","chisq_association")]
  mrmegares$P <- pchisq(mrmegares$chisq_association,ndf, lower.tail = F)
  colnames(mrmegares)[1] <- "ID"
  colnames(mrmegares)[2] <- "CHR"
  colnames(mrmegares)[3] <- "POS"
  colnames(mrmegares)[4] <- "A1_FREQ"
  colnames(mrmegares)[5] <- "OBS_CT"
  
  fwrite(mrmegares,file=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/totransfer/regenie_",trait,"_MRMEGA.txt"), sep="\t")
  fwrite(ALL_res,file=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/totransfer/regenie_",trait,"_MEGA.txt"), sep="\t")
  fwrite(meta_tab,file=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/totransfer/regenie_",trait,"_META.txt"), sep="\t")
}

meta_analysis("Height")
meta_analysis("LDL")
meta_analysis("HDL")
meta_analysis("TC")
meta_analysis("Waist")
meta_analysis("Calcium")
meta_analysis("Creatinine")
meta_analysis("EGFR")
meta_analysis_bin("Asthma")
meta_analysis_bin("Breast")
meta_analysis_bin("Prostate")
meta_analysis_bin("T2D")
meta_analysis_bin("CAD")