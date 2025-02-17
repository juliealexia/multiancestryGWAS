library(data.table)
#library(ggplot2)
library(dplyr)

#Genotypes : /n/holystore01/LABS/xlin/Lab/tonychen/UKB/Genotypes/ukb_hm3_mega_pca.bed
#Cov and Pheno : /n/holystore01/LABS/xlin/Lab/tonychen/2.Multi-Ancestry/ukb_multi.cov and /n/holystore01/LABS/xlin/Lab/tonychen/2.Multi-Ancestry/ukb_multi.pheno

cov <- as.data.frame(fread("/n/holystore01/LABS/xlin/Lab/tonychen/2.Multi-Ancestry/ukb_multi.cov"))
head(cov)

# AFR 6409
# EAS 599
# EUR 318043
# SAS 7520
# MIX 2335
# UNK 5248

ancestry <- as.data.frame(fread("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/ukb_multi.anc"))
head(ancestry)

pheno <- as.data.frame(fread("/n/holystore01/LABS/xlin/Lab/tonychen/2.Multi-Ancestry/ukb_multi.pheno"))
head(pheno)

setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/")

id_EUR <- ancestry$IID[ancestry$predicted=="EUR" & ancestry$P.EUR>=0.9]
length(id_EUR)
#311061
id_EAS <- ancestry$IID[ancestry$predicted=="EAS" & ancestry$P.EAS>=0.9]
length(id_EAS)
#585
id_AFR <- ancestry$IID[ancestry$predicted=="AFR" & ancestry$P.AFR>=0.9]
length(id_AFR)
#6863
id_SAS <- ancestry$IID[ancestry$predicted=="SAS" & ancestry$P.SAS>=0.9]
length(id_SAS)
#5734
id_AMR <- ancestry$IID[ancestry$predicted=="AMR" & ancestry$P.AMR>=0.9]
length(id_AMR)
#589

bloodtraits <- as.data.frame(fread("/n/holystore01/LABS/xlin/Lab/tonychen/2.Multi-Ancestry/UKB/ukb_blood_pheno.txt"))
for(eth in c("AFR","AMR","EAS","EUR","SAS","ALL")){
  for(i in c(5,6,7,8,12,13,14,15,16)){
    if(eth == "AFR"){
      ids <- id_AFR
    }
    if(eth == "AMR"){
      ids <- id_AMR
    }
    if(eth == "EAS"){
      ids <- id_EAS
    }
    if(eth == "EUR"){
      ids <- id_EUR
    }
    if(eth == "SAS"){
      ids <- id_SAS
    }
    if(eth == "ALL"){
      ids <- c(id_AFR,id_AMR,id_EAS,id_EUR,id_SAS)
    }
    phenos <- pheno[pheno$IID %in% ids, c(1,2,i)]
    phenos <- na.omit(phenos)
    fwrite(phenos, paste0(colnames(pheno)[i],"_",eth,".pheno"))
  }
}

#binary traits
for(eth in c("AFR","AMR","EAS","EUR","SAS","ALL")){
  for(i in c(35,32,31,25,27)){
    if(eth == "AFR"){
      ids <- id_AFR
    }
    if(eth == "AMR"){
      ids <- id_AMR
    }
    if(eth == "EAS"){
      ids <- id_EAS
    }
    if(eth == "EUR"){
      ids <- id_EUR
    }
    if(eth == "SAS"){
      ids <- id_SAS
    }
    if(eth == "ALL"){
      ids <- c(id_AFR,id_AMR,id_EAS,id_EUR,id_SAS)
    }
    phenos <- pheno[pheno$IID %in% ids, c(1,2,i)]
    phenos <- na.omit(phenos)
    fwrite(phenos, paste0(colnames(pheno)[i],"_",eth,".pheno"))
  }
}

height_all <- pheno[pheno$IID%in%c(id_AFR,id_AMR,id_EAS,id_EUR,id_SAS),c(1,2,3)]
height_all <-na.omit(height_all)
height_all <- height_all[height_all$Height>130,]
fwrite(height_all,file="height_ALL.pheno",sep="\t")
BMI_all <- pheno[pheno$IID%in%c(id_AFR,id_AMR,id_EAS,id_EUR,id_SAS),c(1,2,4)]
BMI_all <-na.omit(BMI_all)
fwrite(BMI_all,file="BMI_ALL.pheno",sep="\t")
id_ALL <- height_all[,c(1:2)]
fwrite(id_ALL,file="id_ALL.txt",sep="\t",quote=F,col.names = T, row.names = F)

height_AFR <- pheno[pheno$IID%in%id_AFR,c(1:3)]
height_AFR <- na.omit(height_AFR)
height_AFR <- height_AFR[height_AFR$Height>130,]
fwrite(height_AFR,file="height_AFR.pheno",sep="\t")
BMI_AFR <- pheno[pheno$IID%in%c(id_AFR),c(1,2,4)]
BMI_AFR <-na.omit(BMI_AFR)
fwrite(BMI_AFR,file="BMI_AFR.pheno",sep="\t")
id_AFR <- pheno[pheno$IID%in%id_AFR,c(1:2)]
fwrite(id_AFR,file="id_AFR.txt",sep="\t",quote=F,col.names = T, row.names = F)

height_AMR <- pheno[pheno$IID%in%id_AMR,c(1:3)]
height_AMR <- na.omit(height_AMR)
height_AMR <- height_AMR[height_AMR$Height>130,]
fwrite(height_AMR,file="height_AMR.pheno",sep="\t")
BMI_AMR <- pheno[pheno$IID%in%c(id_AMR),c(1,2,4)]
BMI_AMR <-na.omit(BMI_AMR)
fwrite(BMI_AMR,file="BMI_AMR.pheno",sep="\t")
id_AMR <- pheno[pheno$IID%in%id_AMR,c(1:2)]
fwrite(id_AMR,file="id_AMR.txt",sep="\t",quote=F,col.names = T, row.names = F)

height_EAS <- pheno[pheno$IID%in%id_EAS,c(1:3)]
height_EAS <- na.omit(height_EAS)
height_EAS <- height_EAS[height_EAS$Height>130,]
fwrite(height_EAS,file="height_EAS.pheno",sep="\t")
BMI_EAS <- pheno[pheno$IID%in%c(id_EAS),c(1,2,4)]
BMI_EAS <-na.omit(BMI_EAS)
fwrite(BMI_EAS,file="BMI_EAS.pheno",sep="\t")
id_EAS <- pheno[pheno$IID%in%id_EAS,c(1:2)]
fwrite(id_EAS,file="id_EAS.txt",sep="\t",quote=F,col.names = T, row.names = F)

height_EUR <- pheno[pheno$IID%in%id_EUR,c(1:3)]
height_EUR <- na.omit(height_EUR)
height_EUR <- height_EUR[height_EUR$Height>130,]
fwrite(height_EUR,file="height_EUR.pheno",sep="\t")
BMI_EUR <- pheno[pheno$IID%in%c(id_EUR),c(1,2,4)]
BMI_EUR <-na.omit(BMI_EUR)
fwrite(BMI_EUR,file="BMI_EUR.pheno",sep="\t")
id_EUR <- pheno[pheno$IID%in%id_EUR,c(1:2)]
fwrite(id_EUR,file="id_EUR.txt",sep="\t",quote=F,col.names = T, row.names = F)

height_SAS <- pheno[pheno$IID%in%id_SAS,c(1:3)]
height_SAS <- na.omit(height_SAS)
height_SAS <- height_SAS[height_SAS$Height>130,]
fwrite(height_SAS,file="height_SAS.pheno",sep="\t")
SAS_BMI <- pheno[pheno$IID%in%c(id_SAS),c(1,2,4)]
SAS_BMI <-na.omit(SAS_BMI)
fwrite(SAS_BMI,file="BMI_SAS.pheno",sep="\t")
id_SAS <- pheno[pheno$IID%in%id_SAS,c(1:2)]
fwrite(id_SAS,file="id_SAS.txt",sep="\t",quote=F,col.names = T, row.names = F)

# AFR
system("/n/home12/jdias/plink2 --bfile /n/holystore01/LABS/xlin/Lab/tonychen/UKB/Genotypes/ukb_hm3_mega_pca --keep /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/id_AFR.txt --maf 0.01 --make-bed --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/AFR")
# AMR
system("/n/home12/jdias/plink2 --bfile /n/holystore01/LABS/xlin/Lab/tonychen/UKB/Genotypes/ukb_hm3_mega_pca --keep /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/id_AMR.txt --maf 0.01 --make-bed --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/AMR")
# EAS
system("/n/home12/jdias/plink2 --bfile /n/holystore01/LABS/xlin/Lab/tonychen/UKB/Genotypes/ukb_hm3_mega_pca --keep /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/id_EAS.txt --maf 0.01 --make-bed --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/EAS")
# EUR
system("/n/home12/jdias/plink2 --bfile /n/holystore01/LABS/xlin/Lab/tonychen/UKB/Genotypes/ukb_hm3_mega_pca --keep /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/id_EUR.txt --maf 0.01 --make-bed --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/EUR")
# SAS
system("/n/home12/jdias/plink2 --bfile /n/holystore01/LABS/xlin/Lab/tonychen/UKB/Genotypes/ukb_hm3_mega_pca --keep /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/id_SAS.txt --maf 0.01 --make-bed --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/SAS")
# ALL
system("/n/home12/jdias/plink2 --bfile /n/holystore01/LABS/xlin/Lab/tonychen/UKB/Genotypes/ukb_hm3_mega_pca --keep /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/id_ALL.txt --make-bed --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/ALL")

# Generate SNP list for PC
system("/n/home12/jdias/plink --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/ALL --indep-pairwise 500 kb 1 0.05 --maf 0.01 --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/prunedsnps")

# Extract SNPs for PC
system("/n/home12/jdias/plink2 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/AFR --extract /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/prunedsnps.prune.in --export vcf --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/subset_AFR")
system("/n/home12/jdias/plink2 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/EAS --extract /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/prunedsnps.prune.in --export vcf --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/subset_EAS")
system("/n/home12/jdias/plink2 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/EUR --extract /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/prunedsnps.prune.in --export vcf --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/subset_EUR")
system("/n/home12/jdias/plink2 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/SAS --extract /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/prunedsnps.prune.in --export vcf --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/subset_SAS")
system("/n/home12/jdias/plink2 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/AMR --extract /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/prunedsnps.prune.in --export vcf --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/subset_AMR")

system("/n/home12/jdias/plink2 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/ALL --extract /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/prunedsnps.prune.in --export vcf --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/subset_ALL")

# Generate PCs
system("/n/home12/jdias/plink2 --vcf /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/subset_AFR.vcf --make-bed --pca approx --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_AFR")
system("/n/home12/jdias/plink2 --vcf /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/subset_EAS.vcf --make-bed --pca approx --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_EAS")
system("/n/home12/jdias/plink2 --vcf /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/subset_EUR.vcf --make-bed --pca approx --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_EUR")
system("/n/home12/jdias/plink2 --vcf /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/subset_SAS.vcf --make-bed --pca approx --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_SAS")
system("/n/home12/jdias/plink2 --vcf /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/subset_AMR.vcf --make-bed --pca approx --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_AMR")

system("/n/home12/jdias/plink2 --vcf /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/subset_ALL.vcf --make-bed --pca approx --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_ALL")

# Reformulate PC IDs and add sex+age info for covariates
ethnicities <- c("AFR","EAS","EUR","SAS","AMR","ALL")
agesex <- cov[,c(2,4,5)]
agesex$IID <- as.character(agesex$IID)
for(eth in ethnicities){
  pcs <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_",eth,".eigenvec")))
  pcstemp <- pcs[,c(2:11)]
  nametemp <- substr(pcs[,1],1,nchar(pcs[,1])/2)
  iids <- data.frame(nametemp,nametemp)
  colnames(iids) <- c("FID","IID")
  pcs <- cbind(iids,pcstemp)
  pcs <-  left_join(pcs, agesex, by="IID")
  pcs$female <- ifelse(pcs$female==1,"Female","Male")
  fwrite(pcs, file=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_",eth,".eigenvec"), sep="\t", quote = F, row.names = F)
}

system("/n/home12/jdias/plink2 --threads 64 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/AFR --covar-variance-standardize --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_AFR.out --linear hide-covar --vif 3000 --covar /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_AFR.eigenvec --covar-name PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 --pheno /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/height_AFR.pheno")
system("/n/home12/jdias/plink2 --threads 64 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/AMR --covar-variance-standardize --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_AMR.out --linear hide-covar --vif 3000 --covar /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_AMR.eigenvec --pheno /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/height_AMR.pheno")
system("/n/home12/jdias/plink2 --threads 64 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/EAS --covar-variance-standardize --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_EAS.out --linear hide-covar --vif 3000 --covar /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_EAS.eigenvec --pheno /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/height_EAS.pheno")
system("/n/home12/jdias/plink2 --threads 64 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/EUR --covar-variance-standardize --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_EUR.out --linear hide-covar --vif 3000 --covar /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_EUR.eigenvec --pheno /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/height_EUR.pheno")
system("/n/home12/jdias/plink2 --threads 64 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/SAS --covar-variance-standardize --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_SAS.out --linear hide-covar --vif 3000 --covar /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_SAS.eigenvec --pheno /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/height_SAS.pheno")
system("/n/home12/jdias/plink2 --threads 64 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/ALL --covar-variance-standardize --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_ALL.out --linear hide-covar --vif 3000 --covar /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_ALL.eigenvec --pheno /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/height_ALL.pheno")


lambda_function <- function(tab){
  x = tab$P
  z = qnorm(x / 2)
  lambda = round(median(z^2) / qchisq(0.5,1), 3)
  N.effect = median(tab$OBS_CT)
  lambda_1000 = round(1+1000*(lambda-1)/N.effect  ,3)
  return(c(lambda,lambda_1000))
}

meta_analysis <- function(trait){
  AFR_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_AFR.out.",trait,".glm.linear")))
  AFR_res$P <- as.numeric(AFR_res$P)
  lambda_val <- lambda_function(AFR_res)
  print(paste0("AFR lambda_GC: ", lambda_val[1]," lambda_1000: ", lambda_val[2]))
  out.frame <- data.frame(MARKERNAME=AFR_res$ID,EA=AFR_res$ALT,NEA=AFR_res$REF,BETA=AFR_res$BETA,SE=AFR_res$SE,EAF=AFR_res$A1_FREQ,N=AFR_res$OBS_CT,CHROMOSOME=AFR_res[,1],POSITION=AFR_res$POS,P=AFR_res$P)
  write.table(out.frame, file="/n/holystore01/LABS/xlin/Everyone/JD_HZ/MRMEGA/AFR_MRMEGA.txt",  sep="\t", quote=F, row.names=F, col.names=T, append=F)
  
  AMR_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_AMR.out.",trait,".glm.linear")))
  AMR_res$P <- as.numeric(AMR_res$P)
  lambda_val <- lambda_function(AMR_res)
  print(paste0("AMR lambda_GC: ", lambda_val[1]," lambda_1000: ", lambda_val[2]))
  out.frame <- data.frame(MARKERNAME=AMR_res$ID,EA=AMR_res$ALT,NEA=AMR_res$REF,BETA=AMR_res$BETA,SE=AMR_res$SE,EAF=AMR_res$A1_FREQ,N=AMR_res$OBS_CT,CHROMOSOME=AMR_res[,1],POSITION=AMR_res$POS,P=AMR_res$P)
  write.table(out.frame, file="/n/holystore01/LABS/xlin/Everyone/JD_HZ/MRMEGA/AMR_MRMEGA.txt",  sep="\t", quote=F, row.names=F, col.names=T, append=F)
  
  EAS_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_EAS.out.",trait,".glm.linear")))
  EAS_res$P <- as.numeric(EAS_res$P)
  lambda_val <- lambda_function(EAS_res)
  print(paste0("EAS lambda_GC: ", lambda_val[1]," lambda_1000: ", lambda_val[2]))
  out.frame <- data.frame(MARKERNAME=EAS_res$ID,EA=EAS_res$ALT,NEA=EAS_res$REF,BETA=EAS_res$BETA,SE=EAS_res$SE,EAF=EAS_res$A1_FREQ,N=EAS_res$OBS_CT,CHROMOSOME=EAS_res[,1],POSITION=EAS_res$POS,P=EAS_res$P)
  write.table(out.frame, file="/n/holystore01/LABS/xlin/Everyone/JD_HZ/MRMEGA/EAS_MRMEGA.txt",  sep="\t", quote=F, row.names=F, col.names=T, append=F)
  
  EUR_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_EUR.out.",trait,".glm.linear")))
  EUR_res$P <- as.numeric(EUR_res$P)
  lambda_val <- lambda_function(EUR_res)
  print(paste0("EUR lambda_GC: ", lambda_val[1]," lambda_1000: ", lambda_val[2]))
  out.frame <- data.frame(MARKERNAME=EUR_res$ID,EA=EUR_res$ALT,NEA=EUR_res$REF,BETA=EUR_res$BETA,SE=EUR_res$SE,EAF=EUR_res$A1_FREQ,N=EUR_res$OBS_CT,CHROMOSOME=EUR_res[,1],POSITION=EUR_res$POS,P=EUR_res$P)
  write.table(out.frame, file="/n/holystore01/LABS/xlin/Everyone/JD_HZ/MRMEGA/EUR_MRMEGA.txt",  sep="\t", quote=F, row.names=F, col.names=T, append=F)
  
  SAS_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_SAS.out.",trait,".glm.linear")))
  lambda_val <- lambda_function(SAS_res)
  print(paste0("SAS lambda_GC: ", lambda_val[1]," lambda_1000: ", lambda_val[2]))
  out.frame <- data.frame(MARKERNAME=SAS_res$ID,EA=SAS_res$ALT,NEA=SAS_res$REF,BETA=SAS_res$BETA,SE=SAS_res$SE,EAF=SAS_res$A1_FREQ,N=SAS_res$OBS_CT,CHROMOSOME=SAS_res[,1],POSITION=SAS_res$POS,P=SAS_res$P)
  write.table(out.frame, file="/n/holystore01/LABS/xlin/Everyone/JD_HZ/MRMEGA/SAS_MRMEGA.txt",  sep="\t", quote=F, row.names=F, col.names=T, append=F)
  
  system(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/MRMEGA/MR-MEGA -i mr-mega.in --qt --pc 2 --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/mrmega_",trait,"_UKB.result"))
  
  ALL_res <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_ALL.out.",trait,".glm.linear")))
  ALL_res <- na.omit(ALL_res)
  
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
  
  colnames(ALL_res)[1]<-"CHR"
  meta_tab <- left_join(meta_tab[,c("ID","P","OBS_CT")], ALL_res[,c("ID","A1_FREQ","CHR","POS")], by = "ID")
  
  ALL_res <- ALL_res[ALL_res$ID %in% meta_tab$ID,]
  ALL_res$P <- as.numeric(ALL_res$P)
  lambda_val_pooled <- lambda_function(ALL_res)
  meta_tab$P <- as.numeric(meta_tab$P)
  lambda_val_meta <- lambda_function(meta_tab)
  
  mrmegares <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/mrmega_",trait,"_UKB.result.result")))
  mrmegares <- mrmegares[is.na(mrmegares$Comments),]
  mrmegares <- mrmegares[,c("MarkerName","Nsample","chisq_association")]
  mrmegares$P <- pchisq(mrmegares$chisq_association,3, lower.tail = F)
  colnames(mrmegares)[2] <- "OBS_CT"
  lambda_val_mrmega <- lambda_function(mrmegares)
  
  fwrite(data.frame("ID"=ALL_res$ID[ALL_res$P<5e-8]),paste0(trait,"_sig_snps_MEGA.txt"),sep = "\t", col.names = F)
  fwrite(data.frame("ID"=meta_tab$ID[meta_tab$P<5e-8]),paste0(trait,"_sig_snps_META.txt"),sep = "\t", col.names = F)
  fwrite(data.frame("ID"=mrmegares$MarkerName[mrmegares$P<5e-8]),paste0(trait,"_sig_snps_MRMEGA.txt"),sep = "\t", col.names = F)
  system(paste0("/n/home12/jdias/plink --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/ALL --extract ",trait,"_sig_snps_MEGA.txt --indep-pairwise 500 kb 1 0.1 --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/MEGA_prunedsnps"))
  system(paste0("/n/home12/jdias/plink --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/ALL --extract ",trait,"_sig_snps_META.txt --indep-pairwise 500 kb 1 0.1 --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/META_prunedsnps"))
  system(paste0("/n/home12/jdias/plink --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/ALL --extract ",trait,"_sig_snps_MRMEGA.txt --indep-pairwise 500 kb 1 0.1 --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/MRMEGA_prunedsnps"))
  remainingSNPSMEGA <- fread("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/MEGA_prunedsnps.prune.in")
  remainingSNPSMETA <- fread("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/META_prunedsnps.prune.in")
  remainingSNPSMRMEGA <- fread("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/MRMEGA_prunedsnps.prune.in")
  
  print(paste0("Number of subjects: ", max(ALL_res$OBS_CT)))
  
  print(paste0("Pooled analysis lambda_1000: ", lambda_val_pooled[2]," lambda_GC: ", lambda_val_pooled[1], " Number of sig. SNPs: ", sum(ALL_res$P<5e-8)," Number of ind.sig. SNPs: ", nrow(remainingSNPSMEGA)))
  print(paste0("Meta analysis lambda_1000: ", lambda_val_meta[2]," lambda_GC: ", lambda_val_meta[1], " Number of sig. SNPs: ", sum(meta_tab$P<5e-8)," Number of ind.sig. SNPs: ", nrow(remainingSNPSMETA)))
  print(paste0("MRMEGA lambda_1000: ", lambda_val_mrmega[2]," lambda_GC: ", lambda_val_mrmega[1], " Number of sig. SNPs: ", sum(mrmegares$P<5e-8)," Number of ind.sig. SNPs: ", nrow(remainingSNPSMRMEGA)))
}

meta_analysis("Height")
meta_analysis("LDL")
meta_analysis("HDL")
meta_analysis("TC")
meta_analysis("BMI")
meta_analysis("Waist")
meta_analysis("Calcium")
meta_analysis("Creatinine")
meta_analysis("Cystatin")
meta_analysis("EGFR")
meta_analysis("HbA1C")



ancestry <- as.data.frame(fread("/Users/juliealexia/Downloads/ukb_multi.anc"))

test <- as.data.frame(fread("/Users/juliealexia/Documents/A_WORKTODO/pca_ALL.eigenvec"))
test <- left_join(test, ancestry[,c("IID","predicted")], by="IID")
ggplot(test, aes(x=PC1,y=PC2,color=predicted)) +geom_point(size=1)+theme_Publication()+scale_colour_Publication() + guides(color=guide_legend("Predicted ancestry"))

test <- as.data.frame(fread("pca_005.eigenvec"))
test$ancestry <- gsub("\\_.*", "", test$IID)
test$ancestry <- ifelse(test$ancestry=="CEU-YRI","Admixed 1",test$ancestry)
test$ancestry <- ifelse(test$ancestry=="YRI-CEU","Admixed 2",test$ancestry)
ggplot(test, aes(x=PC1,y=PC2,color=ancestry)) +geom_point(size=1)+theme_Publication()+scale_colour_Publication() + guides(color=guide_legend("Ancestry"))

test2 <- as.data.frame(fread("/Users/juliealexia/Documents/A_WORKTODO/ukb_het_cov.txt"))
test3 <- as.data.frame(fread("/Users/juliealexia/Documents/A_WORKTODO/ukb_blood_pheno.txt"))
test3 <- left_join(test2, test3[,c("IID","predicted")])
ggplot(test3, aes(x=pc1,y=pc2,color=predicted)) + geom_point()

ggplot(sub, aes(x=PC1,y=PC2,color=predicted)) + geom_point()
ancestry_sub <- ancestry[(ancestry$P.AFR>=0.9 | ancestry$P.EUR>=0.9 | ancestry$P.EAS>=0.9 | ancestry$P.AMR>=0.9 | ancestry$P.SAS>=0.9),]
ancestry_sub <- ancestry[(ancestry$P.EUR>=0.9),]
ggplot(ancestry_sub, aes(x=PC1,y=PC2,color=predicted)) + geom_point() 
ggplot(test, aes(x=PC1,y=PC2)) + geom_point() + geom_point(data=test[test$IID%in%c(1901042, 3987438, 4129760),], aes(x=PC1,y=PC2),color='red')

testgwas <- as.data.frame(fread("/Users/juliealexia/Documents/A_WORKTODO/summary_EUR.out.Height.glm.linear"))
colnames(testgwas)[1] <- "CHR"
testgwas <- testgwas[testgwas$P<5e-4,]
manhattan(testgwas, chr="CHR", bp="POS", p="P", snp="ID")

