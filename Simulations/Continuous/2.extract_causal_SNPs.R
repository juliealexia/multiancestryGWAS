#goal extract causal SNPs
#since GCTA has limited read-in speed compared to plink
#first extract the causal SNPs
#then simulate phenotypes
args <- commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
load("/n/holystore01/LABS/xlin/Everyone/JD_HZ/cau.snp.infor.rdata")
eth <- c("EUR","AFR","AMR","EAS","SAS")
cau.snp.infor <- cau.snp.infor.list[[1]]
library(data.table)
library(dplyr)
MAF <- cau.snp.infor %>% select(EUR,AFR,AMR,EAS,SAS)
MAF <- as.data.frame(MAF)

##############
# start here #
##############
for(i in 1:5){
  write.table(cau.snp.infor,file = paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/",eth[i],"/select.cau.snp.txt"),row.names = F,col.names = F,quote=F)
}

cur.dir <- "/n/holystore01/LABS/xlin/Everyone/JD_HZ/"
eth <- c("EUR","AFR","AMR","EAS","SAS")


#### with all chr
for(i in 1:5){
  system(paste0("/n/home12/jdias/plink2 --pfile ",cur.dir,eth[i],"/all_chr --extract ",cur.dir,"/",eth[i],"/select.cau.snp.txt --make-bed --out ",cur.dir,"onepercentcausal/",eth[i],"/select.cau.snp"))
}

#### only if doing chr by chr. Need to change ancestry before running!
for(i in 1:22){
  system(paste0("/n/home12/jdias/plink2 --bfile /n/holystore01/LABS/xlin/Everyone/multi_ancestry_simu/EUR_mega/chr",i," --extract /n/holystore01/LABS/xlin/Everyone/JD_HZ/EUR/select.cau.snp.txt --make-bed --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/EUR/select.cau.snp.chr",i))
}



