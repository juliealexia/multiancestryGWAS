args <- commandArgs(trailingOnly = T)
#i is the number of EUR subjects in simulation
i <- as.numeric(args[[1]])
#i is the number of AFR subjects in simulation
j <- as.numeric(args[[2]])
#i is the number of AMR subjects in simulation
k <- as.numeric(args[[3]])
#i is the number of EAS subjects in simulation
l <- as.numeric(args[[4]])
#i is the number of SAS subjects in simulation
m <- as.numeric(args[[5]])

library(data.table)
library(dplyr)

name<- paste0(i/1000,"_",j/1000,"_",k/1000,"_",l/1000,"_",m/1000)

idskeep<-c()

if(i>0){
  idskeep<-c(idskeep,paste0("AFR_",c(1:j)))
}
if(j>0){
  idskeep<-c(idskeep,paste0("AMR_",c(1:k)))
}
if(k>0){
  idskeep<-c(idskeep,paste0("EAS_",c(1:l)))
}
if(l>0){
  idskeep<-c(idskeep,paste0("EUR_",c(1:i)))
}
if(m>0){
  idskeep<-c(idskeep,paste0("SAS_",c(1:m)))
}

tot<-c(i,j,k,l,m)

setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/IDlists")
iids<-data.frame(idskeep,idskeep)
colnames(iids)<-c("#FID","IID")
fwrite(iids, file = paste0("subID_keep_",name,".txt"), sep = "\t", na="NA", quote=FALSE)

pheno<-as.data.frame(fread("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/ALL_pcpheno.phen"))
pheno <- pheno[pheno$IID%in%idskeep,]
setwd(paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/subsets"))
fwrite(pheno, file = paste0("pheno_",name,".phen"), sep = "\t", na="NA", quote=FALSE)

eth <- c("AFR","AMR","EAS","EUR","SAS")
for(i in 1:5){
  pheno<-as.data.frame(fread(paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/",eth[i],"_pcpheno.phen")))
  idskeep<- c(1:tot[i])
  pheno <- pheno[pheno$IID%in%idskeep,]
  setwd(paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/subsets"))
  fwrite(pheno, file = paste0(eth[i],"_pheno_",tot[i],".phen"), sep = "\t", na="NA", quote=FALSE)
  system(paste0("/n/home12/jdias/plink2 --threads 4 --pfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/",eth[i],"/all_chr --out /n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/subsets/summary_",tot[i],"_",eth[i],".out --linear --vif 3000 --covar /n/holystore01/LABS/xlin/Everyone/JD_HZ/",eth[i],"/pca_005.eigenvec --pheno /n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/subsets/",eth[i],"_pheno_",tot[i],".phen"))
}

system(paste0("/n/home12/jdias/plink2 --covar-variance-standardize --threads 4 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/ALL/all_chr --out /n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/subsets/summary_",name,".out --linear --vif 3000 --covar /n/holystore01/LABS/xlin/Everyone/JD_HZ/ALL/pca_005.eigenvec --pheno /n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/subsets/pheno_",name,".phen"))

setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/")
load("cau.snp.infor.rdata")

ncausal<-nrow(cau.snp.infor)

load("snp.infor.match.hm.rdata")
chr_rsid<-snp.infor.match.hm[,c(12,13)]

setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/meta_analysis")

meta_res<-data.frame("Ethnicity"=rep(name,10),"Phenotype"=c(1:10),"Val1"=rep(0,10),"Val2"=rep(0,10), "Val3"=rep(0,10),"NoCauSNPs"=rep(0,10))

for(j in 1:10){
  meta_tab <- as.data.frame(fread(paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/subsets/summary_",name,".out.PHENO",j,".glm.linear")))
  meta_tab <- meta_tab[meta_tab$TEST=="ADD",]
  meta_tab <- na.omit(meta_tab)
  meta_tab$P <- as.numeric(meta_tab$P)
  
  sighits <- meta_tab[meta_tab$P<5e-08,c(3,13,14,15)]
  
  meta_res[j,6] <- nrow(sighits)
  
  meta_res[j,3] <- sum(sighits$ID %in% cau.snp.infor$id)/ncausal
  
  x <- strsplit(sighits$ID,":")
  sighitsid <- data.frame("rs_id"=sapply(x,"[",1),"position"=sapply(x,"[",2),"P"=sighits$P)
  sighitsid <- left_join(sighitsid,chr_rsid,by="rs_id")
  sighitsid$CHR[1]<-1
  for(k in 1:nrow(sighitsid)){
    if(is.na(sighitsid$CHR[k])){
      sighitsid$CHR[k] <- sighitsid$CHR[k-1]
    }
  }
  
  cau_snp_res<-data.frame("ID"=cau.snp.infor$id, "rs_id"=cau.snp.infor$rs_id, "CHR"= cau.snp.infor$CHR, "position"=cau.snp.infor$position)
  cau_snp_res <- left_join(cau_snp_res,meta_tab, by="ID")
  cau_snp_res <- cau_snp_res[,c(1,2,3,4,18)]
  cau_snp_res$MOST_SIG_500KB <- NA
  cau_snp_res$P_MOST_SIG_500KB <- NA
  
  sighitsincausalreg<-c()
  
  for (l in 1:nrow(cau_snp_res)) {
    if(l%%10000 == 0){
      print(paste0(l," causal SNPs checked"))
    }
    chr<-cau.snp.infor$CHR[l]
    pos<-cau.snp.infor$position[l]
    #Within 500kb
    in500kb<-subset(sighitsid, CHR==chr & position>=(pos-500000) & position<=(pos+500000))
    if(nrow(in500kb)>0){
      cau_snp_res$P_MOST_SIG_500KB[l]<-min(in500kb$P)
      cau_snp_res$MOST_SIG_500KB[l]<-in500kb[in500kb$P==min(in500kb$P),1]
      sighitsincausalreg <- union(sighitsincausalreg,in500kb$rs_id)
    }
  }
  
  meta_res[j,4]<-(ncausal-sum(is.na(cau_snp_res$P_MOST_SIG_500KB)))/ncausal
  
  meta_res[j,5]<-(nrow(sighits)-length(sighitsincausalreg))/nrow(sighits)
}

pooledresults <- colMeans(meta_res[,c(3,4,5,6)])

meta_res<-data.frame("Ethnicity"=rep("MRMEGA_ALL",10),"Phenotype"=c(1:10),"Val1"=rep(0,10),"Val2"=rep(0,10), "Val3"=rep(0,10),"NoCauSNPs"=rep(0,10))

for(j in 1:10){
  
  #Convert Plink output files into MR-MEGA files
  for(i in 1:5){
    dat_assoc <- as.data.frame(fread(paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/subsets/summary_",tot[i],"_",eth[i],".out.PHENO",j,".glm.linear")))
    dat_assoc <- dat_assoc[dat_assoc$TEST=="ADD",]
    dat_assoc <- na.omit(dat_assoc)
    head(dat_assoc)
    print(sum(dat_assoc$P<5e-08))
    
    out.frame <- data.frame(MARKERNAME=dat_assoc$ID,EA=dat_assoc$ALT,NEA=dat_assoc$REF,BETA=dat_assoc$BETA,SE=dat_assoc$SE,EAF=dat_assoc$A1_FREQ,N=dat_assoc$OBS_CT,CHROMOSOME=dat_assoc[,1],POSITION=dat_assoc$POS,P=dat_assoc$P)
    
    setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/MRMEGA")
    write.table(out.frame, paste0(eth[i],"_",tot[i],"_MRMEGA.txt"),  sep="\t", quote=F, row.names=F, col.names=T, append=F)
  }
  
  setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/MRMEGA")
  mrmega <- data.frame("name"=paste0(eth,"_",tot,"_MRMEGA.txt"))
  
  fwrite(mrmega, file = paste0(name,"_mrmega.in"), sep = "\t", na="NA", quote=FALSE, col.names = F)
  
  #Run MR-MEGA (5 cohorts so 2 PCs max)
  system(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/MRMEGA/MR-MEGA -i ",name,"_mrmega.in --qt --pc 2 --out /n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/subsets/mrmega_",name,"_pheno_",j))
  
  mrmegares <- as.data.frame(fread(paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/subsets/mrmega_",name,"_pheno_",j,".result")))
  head(mrmegares)
  
  keep <- mrmegares[is.na(mrmegares$Comments),]
  sighits <- keep[,c(1,18)]
  colnames(sighits) <- c("ID","P")
  sighits <- sighits[sighits$P<5e-08,]
  meta_res[j,6] <- nrow(sighits)
  
  meta_res[j,3] <- sum(sighits$ID %in% cau.snp.infor$id)/ncausal
  
  x <- strsplit(sighits$ID,":")
  sighitsid <- data.frame("rs_id"=sapply(x,"[",1),"position"=sapply(x,"[",2),"P"=sighits$P)
  sighitsid <- left_join(sighitsid,chr_rsid,by="rs_id")
  sighitsid$CHR[1]<-1
  for(k in 1:nrow(sighitsid)){
    if(is.na(sighitsid$CHR[k])){
      sighitsid$CHR[k] <- sighitsid$CHR[k-1]
    }
  }
  cau_snp_res<-data.frame("ID"=cau.snp.infor$id, "rs_id"=cau.snp.infor$rs_id, "CHR"= cau.snp.infor$CHR, "position"=cau.snp.infor$position_GRCh37)
  colnames(keep)[1]<-"ID"
  cau_snp_res <- left_join(cau_snp_res,keep, by="ID")
  cau_snp_res <- cau_snp_res[,c(1,2,3,4,21)]
  cau_snp_res$MOST_SIG_500KB <- NA
  cau_snp_res$P_MOST_SIG_500KB <- NA
  
  sighitsincausalreg<-c()
  
  for (l in 1:nrow(cau_snp_res)) {
    if(l%%10000 == 0){
      print(paste0(l," causal SNPs checked"))
    }
    chr<-cau.snp.infor$CHR[l]
    pos<-cau.snp.infor$position[l]
    #Within 500kb
    in500kb<-subset(sighitsid, CHR==chr & position>=(pos-500000) & position<=(pos+500000))
    if(nrow(in500kb)>0){
      cau_snp_res$P_MOST_SIG_500KB[l]<-min(in500kb$P)
      cau_snp_res$MOST_SIG_500KB[l]<-in500kb[in500kb$P==min(in500kb$P),1]
      sighitsincausalreg <- union(sighitsincausalreg,in500kb$rs_id)
    }
  }
  
  meta_res[j,4]<-(ncausal-sum(is.na(cau_snp_res$P_MOST_SIG_500KB)))/ncausal
  
  meta_res[j,5]<-(nrow(sighits)-length(sighitsincausalreg))/nrow(sighits)
  
  print(meta_res)
}

mrmegaresults <- colMeans(meta_res[,c(3,4,5,6)])

for(j in 1:10){
  dat_AFR <- as.data.frame(fread(paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/subsets/summary_",tot[1],"_AFR.out.PHENO",j,".glm.linear")))
  dat_AFR <- dat_AFR[dat_AFR$TEST=="ADD",]
  met_AFR <- data.frame("ID"=dat_AFR$ID, "N_AFR"=dat_AFR$OBS_CT, "BETA_AFR"=dat_AFR$BETA, "VAR_AFR"=(dat_AFR$SE)^2)
  
  dat_AMR <- as.data.frame(fread(paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/subsets/summary_",tot[2],"_AMR.out.PHENO",j,".glm.linear")))
  dat_AMR <- dat_AMR[dat_AMR$TEST=="ADD",]
  met_AMR <- data.frame("ID"=dat_AMR$ID, "N_AMR"=dat_AMR$OBS_CT, "BETA_AMR"=dat_AMR$BETA, "VAR_AMR"=(dat_AMR$SE)^2)
  
  dat_EAS <- as.data.frame(fread(paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/subsets/summary_",tot[3],"_EAS.out.PHENO",j,".glm.linear")))
  dat_EAS <- dat_EAS[dat_EAS$TEST=="ADD",]
  met_EAS <- data.frame("ID"=dat_EAS$ID, "N_EAS"=dat_EAS$OBS_CT, "BETA_EAS"=dat_EAS$BETA, "VAR_EAS"=(dat_EAS$SE)^2)
  
  dat_EUR <- as.data.frame(fread(paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/subsets/summary_",tot[4],"_EUR.out.PHENO",j,".glm.linear")))
  dat_EUR <- dat_EUR[dat_EUR$TEST=="ADD",]
  met_EUR <- data.frame("ID"=dat_EUR$ID, "N_EUR"=dat_EUR$OBS_CT, "BETA_EUR"=dat_EUR$BETA, "VAR_EUR"=(dat_EUR$SE)^2)
  
  dat_SAS <- as.data.frame(fread(paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/subsets/summary_",tot[5],"_SAS.out.PHENO",j,".glm.linear")))
  dat_SAS <- dat_SAS[dat_SAS$TEST=="ADD",]
  met_SAS <- data.frame("ID"=dat_SAS$ID, "N_SAS"=dat_SAS$OBS_CT, "BETA_SAS"=dat_SAS$BETA, "VAR_SAS"=(dat_SAS$SE)^2)
  
  sum(dat_AFR$P<5e-08)
  sum(dat_AMR$P<5e-08)
  sum(dat_EAS$P<5e-08)
  sum(dat_EUR$P<5e-08)
  sum(dat_SAS$P<5e-08)
  
  all_SNPs1 <- union(dat_EUR$ID, dat_AFR$ID)
  all_SNPs2 <- union(dat_AMR$ID, dat_EAS$ID)
  all_SNPs3 <- union(all_SNPs1, all_SNPs2)
  all_SNPs <- union(all_SNPs3, dat_SAS$ID)
  
  meta_tab <- data.frame("ID"=all_SNPs)
  meta_tab <- left_join(meta_tab, met_AFR, by = "ID")
  meta_tab <- left_join(meta_tab, met_AMR, by = "ID")
  meta_tab <- left_join(meta_tab, met_EAS, by = "ID")
  meta_tab <- left_join(meta_tab, met_EUR, by = "ID")
  meta_tab <- left_join(meta_tab, met_SAS, by = "ID")
  
  temp_afr <- meta_tab$BETA_AFR/meta_tab$VAR_AFR
  temp_amr <- meta_tab$BETA_AMR/meta_tab$VAR_AMR
  temp_eas <- meta_tab$BETA_EAS/meta_tab$VAR_EAS
  temp_eur <- meta_tab$BETA_EUR/meta_tab$VAR_EUR
  temp_sas <- meta_tab$BETA_SAS/meta_tab$VAR_SAS
  
  temp_afr2 <- 1/meta_tab$VAR_AFR
  temp_amr2 <- 1/meta_tab$VAR_AMR
  temp_eas2 <- 1/meta_tab$VAR_EAS
  temp_eur2 <- 1/meta_tab$VAR_EUR
  temp_sas2 <- 1/meta_tab$VAR_SAS
  
  temp_beta<-rowSums( cbind (temp_afr,temp_amr,temp_eas,temp_eur,temp_sas), na.rm=TRUE)
  
  temp_var <- rowSums( cbind (temp_afr2,temp_amr2,temp_eas2,temp_eur2,temp_sas2), na.rm=TRUE)
  
  meta_tab$BETA_ALL <- temp_beta/temp_var
  meta_tab$SE_ALL <- sqrt(1/temp_var)
  meta_tab$T_STAT_ALL <- meta_tab$BETA_ALL/meta_tab$SE_ALL
  meta_tab$P_ALL <- 2*pnorm(abs(meta_tab$T_STAT_ALL), 0, 1, lower = F)
  
  sighits<-meta_tab[meta_tab$P_ALL<5e-08,c(1,18,19,20)]
  
  meta_res[j,6]<-nrow(sighits)
  
  meta_res[j,3]<-sum(sighits$ID %in% cau.snp.infor$id)/ncausal
  
  x<-strsplit(sighits$ID,":")
  sighitsid<-data.frame("rs_id"=sapply(x,"[",1),"position"=sapply(x,"[",2),"P"=sighits$P_ALL)
  sighitsid<-left_join(sighitsid,chr_rsid,by="rs_id")
  sighitsid$CHR[1]<-1
  for(k in 1:nrow(sighitsid)){
    if(is.na(sighitsid$CHR[k])){
      sighitsid$CHR[k]<-sighitsid$CHR[k-1]
    }
  }
  
  cau_snp_res<-data.frame("ID"=cau.snp.infor$id, "rs_id"=cau.snp.infor$rs_id, "CHR"= cau.snp.infor$CHR, "position"=cau.snp.infor$position)
  cau_snp_res <- inner_join(cau_snp_res,meta_tab, by="ID")
  cau_snp_res <- cau_snp_res[,c(1,2,3,4,23)]
  cau_snp_res$MOST_SIG_500KB <- NA
  cau_snp_res$P_MOST_SIG_500KB <- NA
  
  sighitsincausalreg<-c()
  
  for (l in 1:nrow(cau_snp_res)) {
    if(l%%10000 == 0){
      print(paste0(l," causal SNPs checked"))
    }
    chr<-cau.snp.infor$CHR[l]
    pos<-cau.snp.infor$position[l]
    #Within 500kb
    in500kb<-subset(sighitsid, CHR==chr & position>=(pos-500000) & position<=(pos+500000))
    if(nrow(in500kb)>0){
      cau_snp_res$P_MOST_SIG_500KB[l]<-min(in500kb$P)
      cau_snp_res$MOST_SIG_500KB[l]<-in500kb[in500kb$P==min(in500kb$P),1]
      sighitsincausalreg <- union(sighitsincausalreg,in500kb$rs_id)
    }
  }
  
  meta_res[j,4]<-(ncausal-sum(is.na(cau_snp_res$P_MOST_SIG_500KB)))/ncausal
  
  meta_res[j,5]<-(nrow(sighits)-length(sighitsincausalreg))/nrow(sighits)
}

metaresults <- colMeans(meta_res[,c(3,4,5,6)])


print(pooledresults)
print(mrmegaresults)
print(metaresults)