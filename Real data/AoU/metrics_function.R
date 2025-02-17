metrics_function <- function(trait){
  data_MEGA <- as.data.frame(fread(paste0("summary_ALL_step2_",trait,".regenie")))
  data_MEGA$EXTRA <- NULL
  data_MEGA <- na.omit(data_MEGA)
  colnames(data_MEGA)[2] <- "POS"
  colnames(data_MEGA)[4] <- "REF"
  colnames(data_MEGA)[5] <- "ALT"
  colnames(data_MEGA)[6] <- "A1_FREQ"
  colnames(data_MEGA)[7] <- "OBS_CT"
  data_MEGA$P <- 10^(-data_MEGA$LOG10P)
  
  data_META <- as.data.frame(fread(paste0(trait,"_META_regenie.txt")))
  data_META <- na.omit(data_META)
  data_META$P <- 2*pnorm(abs(data_META$T_STAT), 0, 1, lower = F)
  data_MRMEGA <- as.data.frame(fread(paste0(trait,"_MEGA_regenie.txt")))
  colnames(data_MEGA)[1]<-"CHR"
  meta_tab <- left_join(data_META[,c("ID","P")], data_MEGA[,c("ID","A1_FREQ","CHR","POS","OBS_CT")], by = "ID")
  ALL_res <- data_MEGA[data_MEGA$ID %in% meta_tab$ID,]
  ALL_res$P <- as.numeric(ALL_res$P)
  meta_tab$P <- as.numeric(meta_tab$P)
  lambda_val_pooled <- lambda_function(ALL_res)
  lambda_val_meta <- lambda_function(meta_tab)
  mrmegares <- data_MRMEGA[is.na(data_MRMEGA$Comments),]
  mrmegares <- mrmegares[,c("MarkerName","Nsample","chisq_association")]
  mrmegares$P <- pchisq(mrmegares$chisq_association,4, lower.tail = F)
  colnames(mrmegares)[2] <- "OBS_CT"
  lambda_val_mrmega <- lambda_function(mrmegares)
  fwrite(data.frame("ID"=ALL_res$ID[ALL_res$P<5e-8]),paste0(trait,"_sig_snps_MEGA.txt"),sep = "\t", col.names = F)
  fwrite(data.frame("ID"=meta_tab$ID[meta_tab$P<5e-8]),paste0(trait,"_sig_snps_META.txt"),sep = "\t", col.names = F)
  fwrite(data.frame("ID"=mrmegares$MarkerName[mrmegares$P<5e-8]),paste0(trait,"_sig_snps_MRMEGA.txt"),sep = "\t", col.names = F)
  system(paste0("plink --bfile ALL_filt --extract ",trait,"_sig_snps_MEGA.txt --indep-pairwise 500 kb 1 0.1 --out MEGA_prunedsnps"))
  system(paste0("plink --bfile ALL_filt --extract ",trait,"_sig_snps_META.txt --indep-pairwise 500 kb 1 0.1 --out META_prunedsnps"))
  system(paste0("plink --bfile ALL_filt --extract ",trait,"_sig_snps_MRMEGA.txt --indep-pairwise 500 kb 1 0.1 --out MRMEGA_prunedsnps"))
  remainingSNPSMEGA <- fread("MEGA_prunedsnps.prune.in")
  remainingSNPSMETA <- fread("META_prunedsnps.prune.in")
  remainingSNPSMRMEGA <- fread("MRMEGA_prunedsnps.prune.in")
  print(paste0("Number of subjects: ", max(ALL_res$OBS_CT)))
  
  print(paste0("Pooled analysis lambda_1000: ", lambda_val_pooled[2]," lambda_GC: ", lambda_val_pooled[1], " Number of sig. SNPs: ", sum(ALL_res$P<5e-8)," Number of ind.sig. SNPs: ", nrow(remainingSNPSMEGA)))
  print(paste0("Meta analysis lambda_1000: ", lambda_val_meta[2]," lambda_GC: ", lambda_val_meta[1], " Number of sig. SNPs: ", sum(meta_tab$P<5e-8)," Number of ind.sig. SNPs: ", nrow(remainingSNPSMETA)))
  print(paste0("MRMEGA lambda_1000: ", lambda_val_mrmega[2]," lambda_GC: ", lambda_val_mrmega[1], " Number of sig. SNPs: ", sum(mrmegares$P<5e-8)," Number of ind.sig. SNPs: ", nrow(remainingSNPSMRMEGA)))
  
}
