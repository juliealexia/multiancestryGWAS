library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(magrittr)

setwd("/lus/eagle/projects/prs-atlas/JulieOutput")

eth <- c("AFR","AMR","EAS","EUR","SAS")

## Load SNP information file
load("/lus/eagle/projects/prs-atlas/JulieOutput/snp.infor.match.hm.rdata")
head(snp.infor.match.hm)

ALL_SNP <- as.data.frame(read.table("/lus/eagle/projects/prs-atlas/JulieOutput/ALL_CC.bim"))
head(ALL_SNP)
ALL_SNP_id <- sub("\\:.*", "", ALL_SNP$V2)
sum(ALL_SNP_id %in% snp.infor.match.hm$rs_id)
#1205353 HapMap SNPs
#1205289 HapMap SNPs present in our data

nrow(ALL_SNP)
#2026482 total SNPs in our data

#If we want 1% of the HapMap SNPs to be causal SNPs that means we want 12053 causal SNPs
potentialSNPs <- ALL_SNP_id[ALL_SNP_id %in% snp.infor.match.hm$rs_id]
set.seed(3007)
keep <- sort(sample(1:1205289,12053, replace=F))
causalSNPs_id <- potentialSNPs[keep]

cau.snp.infor <- snp.infor.match.hm[snp.infor.match.hm$rs_id%in%causalSNPs_id,]
nrow(cau.snp.infor)
head(cau.snp.infor)

# Save causal SNPs (1% of HapMap3 SNPs)
save(cau.snp.infor,file = "/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/cau.snp.infor.rdata")
n.total.snp <- nrow(cau.snp.infor)

## Create variant file to extract the causal SNPs from the pgen files
causnps<-data.frame(cau.snp.infor$rs_id,cau.snp.infor$CHR,cau.snp.infor$position)
colnames(causnps)<-c("ID","CHR","pos")
write.table(causnps, file = "select.cau.snp.txt", sep = "\t", na="NA", quote=FALSE, row.names = F, col.names = F)
## Some rs_id have the form 1/2/3/.../22 so remove those because they are not part of our causal SNPs
exclude <- list(c(1:22))
write.table(exclude, file = "exclude.txt", sep = "\t", na="NA", quote=FALSE, row.names = F, col.names = F)

## Extract causal SNPs from ancestry datasets
for(i in 1:5){
  system(paste0("/lus/eagle/projects/covid-xray/arodriguez/tools/plink2/plink2 --pfile /lus/eagle/projects/prs-atlas/JulieOutput/",eth[i],"_ALL --extract /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/select.cau.snp.txt --exclude /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/exclude.txt --make-pgen --out /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/",eth[i],"_CauSNPs"))
}

## Generate effect sizes i.e. beta ~ N(0,sigma^2)
set.seed(666)
beta <- rnorm(n.total.snp,mean=0,sd=0.0123)

## Do beta_hat=beta*sqrt(2*af*(1-af))
## h^2=n_c*Var(beta_hat)=sum(beta_hat^2) where n_c = number of causal snps i.e. n.total.snp
## We want var(G*beta)=0.5683 where var(G*beta)=h^2
## Given the breast cancer GWAS paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5798588/pdf/emss-74191.pdf), roughly 41% of the FRR can be explained by GWAS chips. Then the h^2 for binary outcome is 0.41*2*log(2) = 0.5683807.

setwd("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal")

#Save original effect size for each ancestry and make recode effect size according to minor/major allele
for (i in 1:5){
  select.cau <- data.frame("snpid"=cau.snp.infor$id,"effect_size"=beta)
  #head(select.cau)
  #plink format used minor allele as coding allele
  #the fifth column is minor allele
  #need to match the minor allele with the coding allele
  snp.infor <- read.table(paste0("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/",eth[i],"_CauSNPs.pvar"))
  #head(snp.infor)
  colnames(snp.infor) <- c("chr","position","rsid","major","minor")
  
  tempsplit <- noquote(strsplit(select.cau$snpid, split = ":"))
  
  select.cau_split <- as.data.frame(matrix(unlist(tempsplit),ncol=4,byrow=T))
  
  select.cau$rsid <- select.cau_split$V1
  select.cau$noncoding <- select.cau_split$V3
  select.cau$coding <- select.cau_split$V4
  
  select.cau.infor <- select.cau[which(select.cau$rsid %in% snp.infor$rsid),]
  select.cau.infor$major <- snp.infor$major
  select.cau.infor$minor <- snp.infor$minor
  
  whicheth <- which(colnames(cau.snp.infor)==eth[i])
  af <- cau.snp.infor[which(cau.snp.infor$rs_id %in% select.cau.infor$rsid),whicheth]
  beta_hat <- select.cau.infor$effect_size*sqrt(2*af*(1-af))
  print(paste0(eth[i]," heritability: ", sum(beta_hat^2)))
  prs_prep <- data.frame(select.cau.infor$rsid,select.cau.infor$coding,select.cau.infor$effect_size)
  
  write.table(prs_prep, file = paste0("prs_prep",eth[i],".txt"), sep = "\t", na="NA", quote=FALSE, row.names = F, col.names = F)
}

## With sigma^2 = 0.0123 we have var(G*beta) to be 0.558, 0.582, 0.524, 0.574, 0.577 for AFR, AMR, EAS, EUR and SAS respectively.


## Do the PRS within each ancestry
for(i in 1:5){
  system(paste0("/lus/eagle/projects/covid-xray/arodriguez/tools/plink2/plink2 --threads 2 --score /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/prs_prep",eth[i],".txt cols=+scoresums,-scoreavgs header no-mean-imputation --pfile /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/",eth[i],"_CauSNPs --out /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/",eth[i],"_PRS"))
}

## 11672 variants for AFR
## 12012 variants for AMR
## 10812 variants for EAS
## 11774 variants for EUR
## 11848 variants for SAS

alphvec <- c(-4.981,-4.448,-4.781,-3.06,-4.406)
for(i in 1:5){
  score <- as.data.frame(read.table(paste0("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/",eth[i],"_PRS.sscore")))
  score <- score[,c(1,2,6)]
  
  alpha <- alphvec[i]
  probs <- exp(alpha + score$V6)/(1+exp(alpha + score$V6))
  set.seed(666)
  pheno <- rbinom(length(probs),1,probs)
  print(paste0(eth[i]," prevalence: ",mean(pheno)))
  
  phenodf <- data.frame(score$V1,score$V1,pheno)
  colnames(phenodf)<-c("#FID","IID","PHENO1")
  write.table(phenodf,file=paste0("pheno_bin_",eth[i],".txt"),sep = "\t", na="NA", quote=FALSE, row.names = F)
}

## For AFR alpha=-4.981 gives prevalence of 0.0139*3 = 0.0417
## For AMR alpha=-4.448 gives prevalence of 0.0091*3 = 0.0273
## For EAS alpha=-4.781 gives prevalence of 0.0154*3 = 0.0462
## For EUR alpha=-3.06 gives prevalence of 0.0252*3 = 0.0756
## For SAS alpha=-4.406 gives prevalence of 0.0154*3 = 0.0462

## Pick-out our case-controls and create new pgen files for each population of the matched case-controls

setwd("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/")
for(i in 1:5){
  pheno <- as.data.frame(read.table(paste0("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/pheno_bin_",eth[i],".txt")))
  colnames(pheno)<-c("#FID","IID","PHENO1")
  caseID <- which(pheno$PHENO1==1)
  set.seed(666)
  controlID <- sample(setdiff(1:length(pheno$IID),caseID),length(caseID))
  allID <- c(caseID,controlID)
  allID <- sort(allID)
  idskeep <- pheno$IID[allID]
  iids<-data.frame(idskeep,idskeep)
  colnames(iids)<-c("#FID","IID")
  write.table(iids, file = paste0("CC_",eth[i],".txt"), sep = "\t", na="NA", quote=FALSE, col.names = T, row.names = F)
  
  system(paste0("/lus/eagle/projects/covid-xray/arodriguez/tools/plink2/plink2 --pfile /lus/eagle/projects/prs-atlas/JulieOutput/",eth[i],"_ALL --keep /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/CC_",eth[i],".txt --make-pgen --out /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/",eth[i],"_CC"))
}


## Convert to bed/bim/fam in order to merge using plink1.9
for(i in 1:5){
  system(paste0("/lus/eagle/projects/covid-xray/arodriguez/tools/plink2/plink2 --pfile /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/",eth[i],"_CC --make-bed --out /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/",eth[i],"_CC"))
}

## Reformat bim files
for(i in 1:5){
  oldbim <- as.data.frame(read.table(paste0("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/",eth[i],"_CC.bim")))
  oldbim$V2 <- paste0(oldbim$V2,":",oldbim$V4,":",oldbim$V5,":",oldbim$V6)
  write.table(oldbim,file=paste0(eth[i],"_CC.bim"),quote = F, col.names = F, row.names = F, sep = "\t")
}

setwd("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/")
## Merge all case-controls together
pgen <- c(paste0(eth,"_CC.bed"))
pvar <- c(paste0(eth,"_CC.bim"))
psam <- c(paste0(eth,"_CC.fam"))
mergedf <- data.frame(pgen, pvar, psam)
write.table(mergedf, file = "mergeCC.txt", sep = "\t", na="NA", quote=FALSE ,col.names = FALSE, row.names = FALSE)

system(paste0("/lus/eagle/projects/prs-atlas/JulieOutput/plink --merge-list /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/mergeCC.txt --out /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/CC_ALL"))

## Generate pruned SNPs for PCs
system(paste0("/lus/eagle/projects/prs-atlas/JulieOutput/plink --bfile /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/CC_ALL --indep-pairwise 500 kb 1 0.05 --out /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/prunedsnps" ))

## Generate PCs
#For ALL
system(paste0("/lus/eagle/projects/covid-xray/arodriguez/tools/plink2/plink2 --bfile /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/CC_ALL --extract /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/prunedsnps.prune.in --export vcf --out /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/subset_005_ALL" ))

system(paste0("/lus/eagle/projects/covid-xray/arodriguez/tools/plink2/plink2 --vcf /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/subset_005_ALL.vcf --make-bed --pca approx --out /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/pca_ALL" ))

#For each ancestry individually
for(i in 1:5){
  system(paste0("/lus/eagle/projects/covid-xray/arodriguez/tools/plink2/plink2 --bfile /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/",eth[i],"_CC --extract /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/prunedsnps.prune.in --export vcf --out /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/subset_005_",eth[i] ))
  
  system(paste0("/lus/eagle/projects/covid-xray/arodriguez/tools/plink2/plink2 --vcf /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/subset_005_",eth[i],".vcf --make-bed --pca approx --out /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/pca_",eth[i] ))
  
}

#Reformulate ID for PCs - need first two columns to match the sample ID format of the bed/bim/fam files
eth2 <- c("ALL","AFR","AMR","EAS","EUR","SAS")
setwd(paste0("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal"))
for(i in 1:6){
  pcs<-as.data.frame(read.table(paste0("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/pca_",eth2[i],".eigenvec")))
  pcstemp<-pcs[,c(2:11)]
  nametemp<-substr(pcs[,1], 1, nchar(pcs[,1])/2)
  iids<-data.frame(nametemp,nametemp)
  colnames(iids)<-c("#FID","IID")
  pcs<-cbind(iids,pcstemp)
  write.table(pcs, file = paste0("pca_",eth2[i],".eigenvec"), sep = "\t", na="NA", quote=FALSE, row.names = F)
}


#Write phenotype CC files for each ancestry
eth <- c("AFR","AMR","EAS","EUR","SAS")
for(i in 1:5){
  iids <- as.data.frame(read.table(paste0("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/CC_",eth[i],".txt")))
  pheno <- as.data.frame(read.table(paste0("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/pheno_bin_",eth[i],".txt")))
  phenotemp <- pheno[pheno$V1 %in% iids$V1,]
  colnames(phenotemp)<-c("#FID","IID","binpheno")
  write.table(phenotemp, file = paste0(eth[i],"_CC.phen"), sep = "\t", na="NA", quote=FALSE, row.names = F)
}

#Write phenotype CC file for pooled
cc_afr <- as.data.frame(read.table("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/AFR_CC.phen"))
cc_amr <- as.data.frame(read.table("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/AMR_CC.phen"))
cc_eas <- as.data.frame(read.table("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/EAS_CC.phen"))
cc_eur <- as.data.frame(read.table("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/EUR_CC.phen"))
cc_sas <- as.data.frame(read.table("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/SAS_CC.phen"))
cc_all <- rbind(cc_afr,cc_amr,cc_eas,cc_eur,cc_sas)
colnames(cc_all)<-c("#FID","IID","binpheno")
write.table(cc_all, file = "ALL_CC.phen", sep = "\t", na="NA", quote=FALSE, row.names = F)

#AFR
system(paste0("/lus/eagle/projects/covid-xray/arodriguez/tools/plink2/plink2 --threads 64 --pfile /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/AFR_CC --out /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/results/summary_pc_AFR.out --logistic hide-covar --covar /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/pca_AFR.eigenvec --pheno /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/AFR_CC.phen --no-psam-pheno --1"))
#AMR
system(paste0("/lus/eagle/projects/covid-xray/arodriguez/tools/plink2/plink2 --threads 64 --pfile /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/AMR_CC --out /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/results/summary_pc_AMR.out --logistic hide-covar --covar /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/pca_AMR.eigenvec --pheno /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/AMR_CC.phen --no-psam-pheno --1 "))
#EAS
system(paste0("/lus/eagle/projects/covid-xray/arodriguez/tools/plink2/plink2 --threads 64 --pfile /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/EAS_CC --out /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/results/summary_pc_EAS.out --logistic hide-covar --covar /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/pca_EAS.eigenvec --pheno /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/EAS_CC.phen --no-psam-pheno --1"))
#EUR
system(paste0("/lus/eagle/projects/covid-xray/arodriguez/tools/plink2/plink2 --threads 64 --pfile /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/EUR_CC --out /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/results/summary_pc_EUR.out --logistic hide-covar --covar /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/pca_EUR.eigenvec --pheno /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/EUR_CC.phen --no-psam-pheno --1"))
#SAS
system(paste0("/lus/eagle/projects/covid-xray/arodriguez/tools/plink2/plink2 --threads 64 --pfile /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/SAS_CC --out /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/results/summary_pc_SAS.out --logistic hide-covar --covar /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/pca_SAS.eigenvec --pheno /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/SAS_CC.phen --no-psam-pheno --1"))
#ALL
system(paste0("/lus/eagle/projects/covid-xray/arodriguez/tools/plink2/plink2 --threads 64 --bfile /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/CC_ALL --out /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/results/summary_pc_ALL.out --logistic hide-covar --covar /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/pca_ALL.eigenvec --pheno /lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/ALL_CC.phen --no-psam-pheno --1 --vif 3000"))


# #CHROM  POS     ID      REF     ALT     A1      FIRTH?  TEST    OBS_CT  OR      LOG(OR)_SE      Z_STAT  P       ERRCODE
afr_res <- as.data.frame(read.table("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/results/summary_pc_AFR.out.binpheno.glm.logistic.hybrid"))
min(afr_res$V13)
nrow(afr_res)
sum(afr_res$V13<5e-08)

amr_res <- as.data.frame(read.table("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/results/summary_pc_AMR.out.binpheno.glm.logistic.hybrid"))
min(amr_res$V13)
nrow(amr_res)
sum(amr_res$V13<5e-08)

eas_res <- as.data.frame(read.table("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/results/summary_pc_EAS.out.binpheno.glm.logistic.hybrid"))
min(eas_res$V13)
nrow(eas_res)
sum(eas_res$V13<5e-08)

eur_res <- as.data.frame(read.table("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/results/summary_pc_EUR.out.binpheno.glm.logistic.hybrid"))
min(eur_res$V13)
nrow(eur_res)
sum(eur_res$V13<5e-08)

sas_res <- as.data.frame(read.table("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/results/summary_pc_SAS.out.binpheno.glm.logistic.hybrid"))
min(sas_res$V13)
nrow(sas_res)
sum(sas_res$V13<5e-08)

setwd("/lus/eagle/projects/prs-atlas/JulieOutput")
load("snp.infor.match.hm.rdata")

setwd("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/results")
load("/lus/eagle/projects/prs-atlas/JulieOutput/onepercentcausal/cau.snp.infor.rdata")


met_AFR <- data.frame("ID"=afr_res$V3, "N_AFR"=afr_res$V9, "BETA_AFR"=log(afr_res$V10), "VAR_AFR"=(afr_res$V11)^2)

met_AMR <- data.frame("ID"=amr_res$V3, "N_AFR"=amr_res$V9, "BETA_AFR"=log(amr_res$V10), "VAR_AFR"=(amr_res$V11)^2)

met_EAS <- data.frame("ID"=eas_res$V3, "N_AFR"=eas_res$V9, "BETA_AFR"=log(eas_res$V10), "VAR_AFR"=(eas_res$V11)^2)

met_EUR <- data.frame("ID"=eur_res$V3, "N_AFR"=eur_res$V9, "BETA_AFR"=log(eur_res$V10), "VAR_AFR"=(eur_res$V11)^2)

met_SAS <- data.frame("ID"=sas_res$V3, "N_AFR"=sas_res$V9, "BETA_AFR"=log(sas_res$V10), "VAR_AFR"=(sas_res$V11)^2)

all_SNPs1 <- union(met_AFR$ID, met_AFR$ID)
all_SNPs2 <- union(met_AMR$ID, met_EAS$ID)
all_SNPs3 <- union(all_SNPs1, all_SNPs2)
all_SNPs <- union(all_SNPs3, met_SAS$ID)
#length(all_SNPs)

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

setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/meta_analysis")
fwrite(meta_tab, file = paste0("meta_pheno_005_afcor",j,".txt"), sep = "\t")

sighits<-meta_tab[meta_tab$P_ALL<5e-08,c(1,18,19,20)]

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
setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/meta_analysis")

meta_res[j,4]<-(ncausal-sum(is.na(cau_snp_res$P_MOST_SIG_500KB)))/ncausal

meta_res[j,5]<-(nrow(sighits)-length(sighitsincausalreg))/nrow(sighits)

