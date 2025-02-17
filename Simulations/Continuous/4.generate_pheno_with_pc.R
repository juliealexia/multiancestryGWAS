eth <- c("AFR","AMR","EAS","EUR","SAS")

#Do standardization
for(i in 1:5){
  select.cau<-as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/",eth[i],"/select.cau.snp.bim")))
  select.cau2<-as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/",eth[i],"/select.cau.GCTA_2")))
  colnames(select.cau)<-c("chr","snpid","nonthing","position","minor","major")
  colnames(select.cau2)<-c("snpid","effectsize")
  select.cau<-left_join(select.cau,select.cau2,by="snpid")
  load("/n/holystore01/LABS/xlin/Everyone/JD_HZ/cau.snp.infor.rdata")
  colnames(cau.snp.infor)[1]<-"snpid"
  select.cau.infor <- left_join(select.cau,cau.snp.infor,by="snpid")
  
  j<-which(colnames(select.cau.infor)==eth[i])
  
  EF_AF<-select.cau.infor[,c(2,7,j)]
  print(paste0(eth[i]," h2 before AF correction:",sum(EF_AF$effectsize^2)))
  #EF_AF has snp_id, original effect size and allele frequency
  
  betas <- EF_AF$effectsize*sqrt(2*EF_AF[,3]*(1-EF_AF[,3]))
  
  EF_AF$effectsize2<-betas
  
  herit <- sum(betas^2)
  
  print(paste0("Actual heritability: ",herit))
  
  prs_prep <- data.frame(select.cau$snpid,select.cau$minor,select.cau$effectsize)
  
  fwrite(prs_prep,file = paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/prs_prep",eth[i],".txt"),row.names = F,col.names = F,quote=F, sep="\t")

}

pheno_var <- rep(0,5)
for(i in 1:5){
  system(paste0("/n/home12/jdias/plink2 --threads 2 --score /n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/prs_prep",eth[i],".txt cols=+scoresums,-scoreavgs header no-mean-imputation --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/",eth[i],"/select.cau.snp --out /n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/",eth[i],"_PRS"))
  score <- as.data.frame(fread(paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/",eth[i],"_PRS.sscore")))
  pheno_var[i] <- var(score$SCORE1_SUM)
  print(paste0(eth[i], " phenotype variance: ", pheno_var[i]))
}

set.seed(666)
fixdefx<-rnorm(5, mean=0, sd=1)

for(i in 1:5){
  pcs <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/",eth[i],"/pca_005.eigenvec")))
  #target variance combined from first 2 PCs is 1.180602% (~1.2%) of original phenotypic variance from --score
  tarvar <- 0.01180602
  varpc1 <- var(pcs$PC1)
  beta1 <- sqrt((tarvar/2)/varpc1)
  varpc2 <- var(pcs$PC2)
  beta2 <- sqrt((tarvar/2)/varpc2)
  pc_part <- beta1*pcs$PC1 + beta2*pcs$PC2
  tot_var <- pheno_var[i] + var(beta1*pcs$PC1 + beta2*pcs$PC2)
  score <- as.data.frame(fread(paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/",eth[i],"_PRS.sscore")))
  pheno_withpc <- score[,c(1,2)]
  score_val <- score$SCORE1_SUM
  
  for(j in 1:10){
    noise <- rnorm(nrow(pcs),0,sqrt(1-tot_var))
    print(paste0("Variance from PCs: ",round(var(beta1*pcs$PC1 + beta2*pcs$PC2),4), " , variance from G*beta: ", round(pheno_var[i],4), " , variance from noise: ", round(var(noise),4)))
    pheno_temp <- score_val + pc_part + noise + fixdefx[i]
    pheno_withpc <- cbind(pheno_withpc, pheno_temp)
    print(paste0("Variance of overall phenotype: ", var(pheno_temp)))
  }
  
  colnames(pheno_withpc)[3:12]<-paste0("PHENO",c(1:10))
  
  fwrite(pheno_withpc, file=paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/",eth[i],"_pcpheno.phen"),row.names = F,quote=F, sep="\t")
}

#Make pooled file

temppheno<-data.frame()

for(i in 1:5){
  temppheno2<-as.data.frame(fread(paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/",eth[i],"_pcpheno.phen")))
  temppheno<-rbind(temppheno,temppheno2)
}

names<-c(paste0("AFR_",c(1:120000)),paste0("AMR_",c(1:120000)),paste0("EAS_",c(1:120000)),paste0("EUR_",c(1:120000)),paste0("SAS_",c(1:120000)))
temppheno[,1]<-names
temppheno[,2]<-names
head(temppheno)

fwrite(temppheno, file="/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/ALL_pcpheno.phen",sep = "\t")