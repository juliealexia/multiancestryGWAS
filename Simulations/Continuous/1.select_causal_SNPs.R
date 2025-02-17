#Goal: select the causal SNPs
setwd("/n/holystore01/LABS/xlin/Everyone/multi_ancestry_simu")
#load the SNPs information
load("snp.infor.match37_38.rdata")
head(snp.infor.match)
nrow(snp.infor.match)
# 19239052 snps
hm3rsid <- read.table("hm3rsid.txt", header = TRUE, sep = "", dec = ".")
head(hm3rsid)
nrow(hm3rsid)
# 1217311 snps

library(dplyr)
library(magrittr)

snp.infor.match.hm <- snp.infor.match %>% inner_join(hm3rsid, by=c('rs_id'='SNP'))
head(snp.infor.match.hm) 
nrow(snp.infor.match.hm)
# 1205353 snps

save(snp.infor.match.hm,file = "/n/holystore01/LABS/xlin/Everyone/JD_HZ/snp.infor.match.hm.rdata")

setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/")
load("snp.infor.match.hm.rdata")
head(snp.infor.match.hm)

set.seed(1)
#Uncomment for 5% causal
#cau.idx <- sort(sample(1:nrow(snp.infor.match.hm),ceiling(nrow(snp.infor.match.hm)*0.05), replace=F))
#Uncomment for 1% causal
#cau.idx <- sort(sample(1:nrow(snp.infor.match.hm),ceiling(nrow(snp.infor.match.hm)*0.01), replace=F))
cau.snp.infor <- snp.infor.match.hm[cau.idx,]

save(cau.snp.infor,file = "/n/holystore01/LABS/xlin/Everyone/JD_HZ/onepercentcausal/cau.snp.infor.rdata")
