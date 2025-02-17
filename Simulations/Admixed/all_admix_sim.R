library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)
library(mvtnorm)

for(i in 1:22){
  system(paste0("/n/home12/jdias/plink2  --pfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/admix-kit/data/CEU-YRI.chrom",i," --pca approx --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/CEU-YRI.chr",i))
}
for(i in 1:22){
  system(paste0("scp jdias@login.rc.fas.harvard.edu:/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/CEU-YRI.chr",i,".eigenvec ."))
}

for(i in 1:22){
  PCS <- as.data.frame(fread(paste0("/Users/juliealexia/Documents/Summer Project Pete+Haoyu/CEU-YRI.chr",i,".eigenvec")))
  
  ggplot(PCS,aes(x=PC1,y=PC2))+geom_point()+ggtitle(paste0("Chromosome ",i))
}

#Merge all chromosomes for admix 80-20 i.e. 80 YRI 20 CEU
system(paste0("/n/home12/jdias/plink2  --pmerge-list /n/home12/jdias/mergelistadmix.txt --pmerge-list-dir /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/admix-kit/data --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/CEU-YRI_all_chr"))

#Merge all chromosomes for admix 50-50
system(paste0("/n/home12/jdias/plink2  --pmerge-list /n/home12/jdias/mergelistadmix2.txt --pmerge-list-dir /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/admix-kit/data --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/YRI-CEU_all_chr"))

#Merge all chromosomes for CEU
system(paste0("/n/home12/jdias/plink2  --pmerge-list /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/admix-kit/data/ancestry/mergelistCEU.txt --pmerge-list-dir /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/admix-kit/data/ancestry --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/CEU_all_chr"))

#Merge all chromosomes for YRI
system(paste0("/n/home12/jdias/plink2  --pmerge-list /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/admix-kit/data/ancestry/mergelistYRI.txt --pmerge-list-dir /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/admix-kit/data/ancestry --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/YRI_all_chr"))

#Rename subjects
eth <- c("CEU","YRI","CEU-YRI","YRI-CEU")

for(i in 1:4){
  psamfile <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/",eth[i],"_all_chr.psam")))
  psamfile[,1] <- paste0(eth[i],"_",c(1:60000))
  psamfile[,2] <- paste0(eth[i],"_",c(1:60000))
  colnames(psamfile)<-c("#FID","IID")
  psamfile$SEX<-NA
  setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/")
  fwrite(psamfile, file = paste0(eth[i],"_all_chr.psam"), sep = "\t", na="NA", quote=FALSE)
}

#Merge all data i.e. CEU YRI and CEU-YRI/YRI-CEU(admixed)
for(i in 1:4){
  system(paste0("/n/home12/jdias/plink2 --pfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/",eth[i],"_all_chr --make-bed --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/",eth[i],"_all_chr"))
}

system(paste0("/n/home12/jdias/plink --merge-list /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/mergelistadmixall.txt --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/ALL_chr"))

#Need to update causal snp list as we are transferring to hg38 reference instead of hg37
setwd("/n/holystore01/LABS/xlin/Everyone/multi_ancestry_simu")
#load the SNPs information
load("snp.infor.match37_38.rdata")
head(snp.infor.match)
select.snp<-paste(snp.infor.match$CHR,snp.infor.match$position_GRCh38,snp.infor.match$a0,snp.infor.match$a1,sep=":")
snp.infor.match$x<-select.snp
sim.snp<-fread("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/ALL_all_chr.bim")
match<-intersect(select.snp,sim.snp$V2)
select.snp.filt<-snp.infor.match[snp.infor.match$x%in%match,]
nrow(select.snp.filt)
hm3rsid <- read.table("hm3rsid.txt", header = TRUE, sep = "", dec = ".")
head(hm3rsid)
snp.infor.match.hm <- select.snp.filt %>% inner_join(hm3rsid, by=c('rs_id'='SNP'))
nrow(snp.infor.match.hm) 

save(snp.infor.match.hm,file = "/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/snp.infor.match.hm.rdata")

set.seed(1)
cau.idx <- sort(sample(1:nrow(snp.infor.match.hm),ceiling(nrow(snp.infor.match.hm)*0.05), replace=F))
cau.snp.infor <- snp.infor.match.hm[cau.idx,]

save(cau.snp.infor,file = "/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/cau.snp.infor.rdata")
select.cau.snp<-cau.snp.infor$x
write.table(cau.snp.infor,file = paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.snp.txt"),row.names = F,col.names = F,quote=F)

setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/")

#Extract causal SNPs
system(paste0("/n/home12/jdias/plink2 --pfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/CEU-YRI_all_chr --extract /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.snp.txt --make-bed --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.snp"))
system(paste0("/n/home12/jdias/plink2 --pfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/YRI-CEU_all_chr --extract /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.snp.txt --make-bed --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/YRI-CEU.select.cau.snp"))
system(paste0("/n/home12/jdias/plink2 --pfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/YRI_all_chr --extract /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.snp.txt --make-bed --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/YRI.select.cau.snp"))
system(paste0("/n/home12/jdias/plink2 --pfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/CEU_all_chr --extract /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.snp.txt --make-bed --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/CEU.select.cau.snp"))


load("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/cau.snp.infor.rdata")
#Generate effect size
GenSigma <- function(sigma,n1,n2,n3,n4,n5,
                     gr12,gr13,gr14,gr15,
                     gr23,gr24,gr25,
                     gr34,gr35,
                     gr45){
  vecsigma <- sigma*c(1/n1,
                      gr12/sqrt(n1*n2),
                      gr13/sqrt(n1*n3),
                      gr14/sqrt(n1*n4),
                      gr15/sqrt(n1*n5),
                      gr12/sqrt(n1*n2),
                      1/n2,
                      gr23/sqrt(n2*n3),
                      gr24/sqrt(n2*n4),
                      gr25/sqrt(n2*n5),
                      gr13/sqrt(n1*n3),
                      gr23/sqrt(n2*n3),
                      1/n3,
                      gr34/sqrt(n3*n4),
                      gr35/sqrt(n3*n5),
                      gr14/sqrt(n1*n4),
                      gr24/sqrt(n2*n4),
                      gr34/sqrt(n3*n4),
                      1/n4,
                      gr45/sqrt(n4*n5),
                      gr15/sqrt(n1*n5),
                      gr25/sqrt(n2*n5),
                      gr35/sqrt(n3*n5),
                      gr45/sqrt(n4*n5),
                      1/n5)
  Sigma <- matrix(vecsigma,5,5)
  return(Sigma)
  
}
n.total.snp <- nrow(cau.snp.infor)

sigma = 1.25
gr12 = 1
gr13 = 1
gr14 = 1
gr15 = 1
gr23 = 1
gr24 = 1
gr25 = 1
gr34 = 1
gr35 = 1
gr45 = 1
n = as.numeric(n.total.snp)
n1 = n
n2 = n
n3 = n
n4 = n
n5 = n
Sigma <- GenSigma(sigma,n1,n2,n3,n4,n5,
                  gr12,gr13,gr14,gr15,
                  gr23,gr24,gr25,
                  gr34,gr35,
                  gr45)

set.seed(666)
beta =  rmvnorm(n.total.snp,c(0,0,0,0,0),sigma=Sigma)

select.cau <- cbind(cau.snp.infor[,16],beta[,1])
write.table(select.cau,file = paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau_rho_5"),row.names = F,col.names = F,quote=F)
  
select.cau <- read.table(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau_rho_5"),header=F)
#head(select.cau)
colnames(select.cau) <- c("snpid","effect_size")
#plink format used minor allele as coding allele
#the fifth column is minor allele
#need to match the minor allele with the coding allele
snp.infor <- read.table("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.snp.bim")
#head(snp.infor)
colnames(snp.infor) <- c("chr","snpid","nonthing","position","minor","major")

select.cau.infor <- left_join(select.cau,snp.infor,by="snpid")

select.cau.infor.split <- select.cau.infor %>% separate(snpid,into=c("rsid","position2","noncoding","coding"),sep=":")
select.cau.infor.split$rsid<-select.cau.infor$snpid
idx <- which(select.cau.infor.split$coding!=select.cau.infor.split$minor)
coding_effect_size <- select.cau$effect_size
minor_effect_size <- coding_effect_size
minor_effect_size[idx] <- -coding_effect_size[idx]
herit <- nrow(select.cau)*var(minor_effect_size)
select.cau.GCTA <- select.cau
select.cau.GCTA$effect_size = minor_effect_size
write.table(select.cau.GCTA,file = paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.GCTA"),row.names = F,col.names = F,quote=F)

select.cau<-as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.snp.bim")))
select.cau2<-as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.GCTA")))
colnames(select.cau)<-c("chr","snpid","nonthing","position","minor","major")
colnames(select.cau2)<-c("snpid","effectsize")
select.cau<-left_join(select.cau,select.cau2,by="snpid")
load("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/cau.snp.infor.rdata")
colnames(cau.snp.infor)[16]<-"snpid"
select.cau.infor <- left_join(select.cau,cau.snp.infor,by="snpid")

#CEU
j<-which(colnames(select.cau.infor)=="EUR")

EF_AF<-select.cau.infor[,c(2,7,j)]
print(paste0("CEU h2 before AF correction: ",sum(EF_AF$effectsize^2)))
#EF_AF has snp_id, original effect size and allele frequency

betas<-EF_AF$effectsize*sqrt(2*EF_AF[,3]*(1-EF_AF[,3]))

EF_AF$effectsize2<-betas

herit<-nrow(EF_AF)*var(betas)

print(paste0("Actual heritability: ",herit))

temp<-left_join(select.cau2,EF_AF,by="snpid")
temp<-temp[,c(1,5)]
colnames(temp)<-c("snpid","effect_size")
temp[is.na(temp$effect_size),2]<-0

write.table(temp,file = paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.GCTA_afcor_CEU"),row.names = F,col.names = F,quote=F)

res <- system(paste0("/n/home12/jdias/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/CEU.select.cau.snp --simu-qt --simu-causal-loci /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.GCTA_afcor_CEU --simu-hsq ",herit," --simu-rep 10 --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_CEU"))


#YRI
j<-which(colnames(select.cau.infor)=="AFR")

EF_AF<-select.cau.infor[,c(2,7,j)]
print(paste0("YRI h2 before AF correction: ",sum(EF_AF$effectsize^2)))
#EF_AF has snp_id, original effect size and allele frequency

betas<-EF_AF$effectsize*sqrt(2*EF_AF[,3]*(1-EF_AF[,3]))

EF_AF$effectsize2<-betas

herit<-nrow(EF_AF)*var(betas)

print(paste0("Actual heritability: ",herit))

temp<-left_join(select.cau2,EF_AF,by="snpid")
temp<-temp[,c(1,5)]
colnames(temp)<-c("snpid","effect_size")
temp[is.na(temp$effect_size),2]<-0

write.table(temp,file = paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.GCTA_afcor_YRI"),row.names = F,col.names = F,quote=F)

res <- system(paste0("/n/home12/jdias/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/YRI.select.cau.snp --simu-qt --simu-causal-loci /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.GCTA_afcor_YRI --simu-hsq ",herit," --simu-rep 10 --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_YRI"))

select.cau<-as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.snp.bim")))
select.cau2<-as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.GCTA")))
colnames(select.cau)<-c("chr","snpid","nonthing","position","minor","major")
colnames(select.cau2)<-c("snpid","effectsize")
select.cau<-left_join(select.cau,select.cau2,by="snpid")
load("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/cau.snp.infor.rdata")
colnames(cau.snp.infor)[16]<-"snpid"
select.cau.infor <- left_join(select.cau,cau.snp.infor,by="snpid")
meanaf<-(0.2*select.cau.infor$EUR+0.8*select.cau.infor$AFR)
select.cau.infor<-cbind(select.cau.infor,meanaf)
colnames(select.cau.infor)[25]<-"EUR-AFR"

#CEU-YRI
j<-which(colnames(select.cau.infor)=="EUR-AFR")

EF_AF<-select.cau.infor[,c(2,7,j)]
print(paste0("YRI h2 before AF correction: ",sum(EF_AF$effectsize^2)))
#EF_AF has snp_id, original effect size and allele frequency

betas<-EF_AF$effectsize*sqrt(2*EF_AF[,3]*(1-EF_AF[,3]))

EF_AF$effectsize2<-betas

herit<-nrow(EF_AF)*var(betas)

print(paste0("Actual heritability: ",herit))

temp<-left_join(select.cau2,EF_AF,by="snpid")
temp<-temp[,c(1,5)]
colnames(temp)<-c("snpid","effect_size")
temp[is.na(temp$effect_size),2]<-0

write.table(temp,file = paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.GCTA_afcor_CEU-YRI"),row.names = F,col.names = F,quote=F)

res <- system(paste0("/n/home12/jdias/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.snp --simu-qt --simu-causal-loci /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.GCTA_afcor_CEU-YRI --simu-hsq ",herit," --simu-rep 10 --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_CEU-YRI"))

#YRI-CEU
meanaf2<-(0.5*select.cau.infor$EUR+0.5*select.cau.infor$AFR)
select.cau.infor<-cbind(select.cau.infor,meanaf2)
colnames(select.cau.infor)[26]<-"AFR-EUR"

#CEU-YRI
j<-which(colnames(select.cau.infor)=="AFR-EUR")

EF_AF<-select.cau.infor[,c(2,7,j)]
print(paste0("YRI-CEU h2 before AF correction: ",sum(EF_AF$effectsize^2)))
#EF_AF has snp_id, original effect size and allele frequency

betas<-EF_AF$effectsize*sqrt(2*EF_AF[,3]*(1-EF_AF[,3]))

EF_AF$effectsize2<-betas

herit<-nrow(EF_AF)*var(betas)

print(paste0("Actual heritability: ",herit))
#0.4228632

temp<-left_join(select.cau2,EF_AF,by="snpid")
temp<-temp[,c(1,5)]
colnames(temp)<-c("snpid","effect_size")
temp[is.na(temp$effect_size),2]<-0

write.table(temp,file = paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.GCTA_afcor_YRI-CEU"),row.names = F,col.names = F,quote=F)

res <- system(paste0("/n/home12/jdias/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/YRI-CEU.select.cau.snp --simu-qt --simu-causal-loci /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/select.cau.GCTA_afcor_YRI-CEU --simu-hsq ",herit," --simu-rep 10 --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_YRI-CEU"))


set.seed(666)
fixdefx<-rnorm(5, mean=0, sd=1)
#eth <- c("EUR","AFR","AMR","EAS","SAS")
pheno <- as.data.frame(fread("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_CEU.phen"))
fixdefx_CEU<-fixdefx[1]
for(i in 3:12){
  pheno[,i]<-pheno[,i]+fixdefx_CEU
}
head(pheno)
write.table(pheno,file = paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_CEU.phen"),row.names = F,col.names = F,quote=F)
pheno <- as.data.frame(fread("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_YRI.phen"))
fixdefx_YRI<-fixdefx[2]
for(i in 3:12){
  pheno[,i]<-pheno[,i]+fixdefx_YRI
}
head(pheno)
write.table(pheno,file = paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_YRI.phen"),row.names = F,col.names = F,quote=F)
pheno <- as.data.frame(fread("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_CEU-YRI.phen"))
fixdefx_admix<-0.2*fixdefx[1]+0.8*fixdefx[2]
for(i in 3:12){
  pheno[,i]<-pheno[,i]+fixdefx_admix
}
head(pheno)
write.table(pheno,file = paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_CEU-YRI.phen"),row.names = F,col.names = F,quote=F)
pheno <- as.data.frame(fread("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_YRI-CEU.phen"))
fixdefx_admix<-0.5*fixdefx[1]+0.5*fixdefx[2]
for(i in 3:12){
  pheno[,i]<-pheno[,i]+fixdefx_admix
}
head(pheno)
write.table(pheno,file = paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_YRI-CEU.phen"),row.names = F,col.names = F,quote=F)

#Merge phenotypes together in one pheno.phen file
phenoCEUYRI<-as.data.frame(fread("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_CEU-YRI.phen"))
phenoYRICEU<-as.data.frame(fread("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_YRI-CEU.phen"))
phenoCEU<-as.data.frame(fread("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_CEU.phen"))
phenoYRI<-as.data.frame(fread("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_YRI.phen"))
phenoALL<-rbind(phenoCEUYRI,phenoCEU)
phenoALL<-rbind(phenoALL,phenoYRICEU)
phenoALL<-rbind(phenoALL,phenoYRI)
write.table(phenoALL,file = paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_ALL.phen"),row.names = F,col.names = F,quote=F)

#Prune SNPs for PCs
system(paste0("/n/home12/jdias/plink --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/ALL_chr --indep-pairwise 500 kb 1 0.05 --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/prunedsnps" ))

#FOR ALL TOGETHER
#Generate PCs
system(paste0("/n/home12/jdias/plink2 -bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/ALL_chr --extract /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/prunedsnps.prune.in --export vcf --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/subset_005" ))

system(paste0("/n/home12/jdias/plink2 --vcf /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/subset_005.vcf --make-bed --pca approx --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/pca_005" ))

#Reformulate ID for PCs - need first two columns to match the sample ID format of the bed/bim/fam files
pcs<-as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/pca_005.eigenvec")))
pcstemp<-pcs[,c(2:11)]
nametemp<-substr(pcs[,1], 1, nchar(pcs[,1])/2)
iids<-data.frame(nametemp,nametemp)
colnames(iids)<-c("#FID","IID")
pcs<-cbind(iids,pcstemp)
setwd(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop"))
fwrite(pcs, file = paste0("pca_005.eigenvec"), sep = "\t", na="NA", quote=FALSE)

#run GWAS
system(paste0("/n/home12/jdias/plink2 --threads 2 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/ALL_chr --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/summary_pc_ALL.out --linear --covar /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/pca_005.eigenvec --pheno /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_ALL.phen"))

setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop")
load("cau.snp.infor.rdata")

ncausal<-nrow(cau.snp.infor)

meta_res<-data.frame("Ethnicity"=rep("Admixed-ALL",10),"Phenotype"=c(1:10),"Val1"=rep(0,10),"Val2"=rep(0,10), "Val3"=rep(0,10))

for(j in 1:10){
  meta_tab <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/summary_pc_ALL.out.PHENO",j,".glm.linear")))
  meta_tab <- meta_tab[meta_tab$TEST=="ADD",]
  meta_tab <- na.omit(meta_tab)
  meta_tab$P<-as.numeric(meta_tab$P)
  sighits <- meta_tab[meta_tab$P<5e-08,c(3,13,14,15)]
  
  meta_res[j,3] <- sum(sighits$ID %in% cau.snp.infor$x)/ncausal
  
  x <- strsplit(sighits$ID,":")
  sighitsid <- data.frame("ID"=sighits$ID,"CHR"=sapply(x,"[",1),"position"=sapply(x,"[",2),"P"=sighits$P)
  cau_snp_res<-data.frame("ID"=cau.snp.infor$x, "rs_id"=cau.snp.infor$rs_id, "CHR"= cau.snp.infor$CHR, "position"=cau.snp.infor$position_GRCh37)
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
      sighitsincausalreg <- union(sighitsincausalreg,in500kb$ID)
    }
  }
  
  setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop")
  fwrite(cau_snp_res, file = paste0("cau_snp_res_pooled_pheno",j,".txt"), sep = "\t")
  
  meta_res[j,4]<-(ncausal-sum(is.na(cau_snp_res$P_MOST_SIG_500KB)))/ncausal
  
  meta_res[j,5]<-(nrow(sighits)-length(sighitsincausalreg))/nrow(sighits)
}

setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop")
fwrite(meta_res, file = paste0("pooled_res.txt"), sep = "\t")

colMeans(meta_res[,c(3,4,5)])

#SINGLE ANCESTRY

eth <- c("CEU","YRI","CEU-YRI","YRI-CEU")
for(i in 1:4){
  system(paste0("/n/home12/jdias/plink2 -bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/",eth[i],"_all_chr --extract /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/prunedsnps.prune.in --export vcf --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/",eth[i],"_subset_005" ))

  system(paste0("/n/home12/jdias/plink2 --vcf /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/",eth[i],"_subset_005.vcf --make-bed --pca approx --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/",eth[i],"_pca_005" ))
}

for(i in 1:4){
  #Reformulate ID for PCs - need first two columns to match the sample ID format of the bed/bim/fam files
  pcs<-as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/",eth[i],"_pca_005.eigenvec")))
  pcstemp<-pcs[,c(2:11)]
  nametemp<-substr(pcs[,1], 1, nchar(pcs[,1])/2)
  iids<-data.frame(nametemp,nametemp)
  colnames(iids)<-c("#FID","IID")
  pcs<-cbind(iids,pcstemp)
  setwd(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop"))
  fwrite(pcs, file = paste0(eth[i],"_pca_005.eigenvec"), sep = "\t", na="NA", quote=FALSE)
}

for(i in 1:4){
  system(paste0("/n/home12/jdias/plink2 --threads 2 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/",eth[i],"_all_chr --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/summary_pc_",eth[i],".out --linear --covar /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/",eth[i],"_pca_005.eigenvec --pheno /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/phenotypes_rho_5_afcor_",eth[i],".phen"))
}

meta_res<-data.frame("Ethnicity"=rep("ALL_meta",10),"Phenotype"=c(1:10),"Val1"=rep(0,10),"Val2"=rep(0,10), "Val3"=rep(0,10))

for(j in 1:10){
  dat_CEU <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/summary_pc_CEU.out.PHENO",j,".glm.linear")))
  dat_CEU <- dat_CEU[dat_CEU$TEST=="ADD",]
  dat_CEU <- na.omit(dat_CEU)
  met_CEU <- data.frame("ID"=dat_CEU$ID, "N_CEU"=dat_CEU$OBS_CT, "BETA_CEU"=dat_CEU$BETA, "VAR_CEU"=(dat_CEU$SE)^2)
  
  dat_YRI <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/summary_pc_YRI.out.PHENO",j,".glm.linear")))
  dat_YRI <- dat_YRI[dat_YRI$TEST=="ADD",]
  met_YRI <- data.frame("ID"=dat_YRI$ID, "N_YRI"=dat_YRI$OBS_CT, "BETA_YRI"=dat_YRI$BETA, "VAR_YRI"=(dat_YRI$SE)^2)
  
  dat_CEU_YRI <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/summary_pc_CEU-YRI.out.PHENO",j,".glm.linear")))
  dat_CEU_YRI <- dat_CEU_YRI[dat_CEU_YRI$TEST=="ADD",]
  met_CEU_YRI <- data.frame("ID"=dat_CEU_YRI$ID, "N_CEUYRI"=dat_CEU_YRI$OBS_CT, "BETA_CEUYRI"=dat_CEU_YRI$BETA, "VAR_CEUYRI"=(dat_CEU_YRI$SE)^2)
  
  dat_YRI_CEU <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/summary_pc_YRI-CEU.out.PHENO",j,".glm.linear")))
  dat_YRI_CEU <- dat_YRI_CEU[dat_YRI_CEU$TEST=="ADD",]
  met_YRI_CEU <- data.frame("ID"=dat_YRI_CEU$ID, "N_YRICEU"=dat_YRI_CEU$OBS_CT, "BETA_YRICEU"=dat_YRI_CEU$BETA, "VAR_YRICEU"=(dat_YRI_CEU$SE)^2)
  
  
  all_SNPs1 <- union(dat_CEU$ID, dat_YRI$ID)
  all_SNPs2 <- union(all_SNPs1, dat_CEU_YRI$ID)
  all_SNPs <- union(all_SNPs2, dat_YRI_CEU$ID)
  #length(all_SNPs)
  
  meta_tab <- data.frame("ID"=all_SNPs)
  meta_tab <- left_join(meta_tab, met_CEU, by = "ID")
  meta_tab <- left_join(meta_tab, met_YRI, by = "ID")
  meta_tab <- left_join(meta_tab, met_CEU_YRI, by = "ID")
  meta_tab <- left_join(meta_tab, met_YRI_CEU, by = "ID")
  
  temp_ceu <- meta_tab$BETA_CEU/meta_tab$VAR_CEU
  temp_yri <- meta_tab$BETA_YRI/meta_tab$VAR_YRI
  temp_ceuyri <- meta_tab$BETA_CEUYRI/meta_tab$VAR_CEUYRI
  temp_yriceu <- meta_tab$BETA_YRICEU/meta_tab$VAR_YRICEU
  
  temp_ceu2 <- 1/meta_tab$VAR_CEU
  temp_yri2 <- 1/meta_tab$VAR_YRI
  temp_ceuyri2 <- 1/meta_tab$VAR_CEUYRI
  temp_yriceu2 <- 1/meta_tab$VAR_YRICEU
  
  temp_beta<-rowSums( cbind (temp_ceu,temp_yri,temp_ceuyri,temp_yriceu), na.rm=TRUE)
  
  temp_var <- rowSums( cbind (temp_ceu2,temp_yri2,temp_ceuyri2,temp_yriceu2), na.rm=TRUE)
  
  meta_tab$BETA_ALL <- temp_beta/temp_var
  meta_tab$SE_ALL <- sqrt(1/temp_var)
  meta_tab$T_STAT_ALL <- meta_tab$BETA_ALL/meta_tab$SE_ALL
  meta_tab$P_ALL <- 2*pnorm(abs(meta_tab$T_STAT_ALL), 0, 1, lower = F)
  
  setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop")
  fwrite(meta_tab, file = paste0("meta_pheno_005",j,".txt"), sep = "\t")
  
  meta_tab<-meta_tab[!is.na(meta_tab$P_ALL),]
  
  sighits<-meta_tab[meta_tab$P_ALL<5e-08,c(1,15,16,17)]
  
  meta_res[j,3]<-sum(sighits$ID %in% cau.snp.infor$x)/ncausal
  
  x <- strsplit(sighits$ID,":")
  sighitsid <- data.frame("ID"=sighits$ID,"CHR"=sapply(x,"[",1),"position"=sapply(x,"[",2),"P"=sighits$P)
  cau_snp_res<-data.frame("ID"=cau.snp.infor$x, "rs_id"=cau.snp.infor$rs_id, "CHR"= cau.snp.infor$CHR, "position"=cau.snp.infor$position_GRCh37)
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
      sighitsincausalreg <- union(sighitsincausalreg,in500kb$ID)
    }
  }
  
  setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop")
  fwrite(cau_snp_res, file = paste0("cau_snp_res_meta_pheno_",j,".txt"), sep = "\t")
  
  meta_res[j,4]<-(ncausal-sum(is.na(cau_snp_res$P_MOST_SIG_500KB)))/ncausal
  
  meta_res[j,5]<-(nrow(sighits)-length(sighitsincausalreg))/nrow(sighits)
}

setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop")
fwrite(meta_res, file = "meta_res.txt", sep = "\t")

colMeans(meta_res[,c(3,4,5)])

##MRMEGA

setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop")
load("cau.snp.infor.rdata")

ncausal<-nrow(cau.snp.infor)

eth <- c("CEU","YRI","CEU-YRI","YRI-CEU")

meta_res<-data.frame("Ethnicity"=rep("ALL_MRMEGA",10),"Phenotype"=c(1:10),"Val1"=rep(0,10),"Val2"=rep(0,10), "Val3"=rep(0,10))

for(j in 1:10){
  
  #Convert Plink output files into MR-MEGA files
  for(i in 1:4){
    dat_assoc <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/summary_pc_",eth[i],".out.PHENO",j,".glm.linear")))
    dat_assoc <- dat_assoc[dat_assoc$TEST=="ADD",]
    dat_assoc <- na.omit(dat_assoc)
    head(dat_assoc)
    
    out.frame <- data.frame(MARKERNAME=dat_assoc$ID,EA=dat_assoc$ALT,NEA=dat_assoc$REF,BETA=dat_assoc$BETA,SE=dat_assoc$SE,EAF=dat_assoc$A1_FREQ,N=dat_assoc$OBS_CT,CHROMOSOME=dat_assoc[,1],POSITION=dat_assoc$POS,P=dat_assoc$P)
    
    setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/MRMEGA")
    write.table(out.frame, paste0(eth[i],"_MRMEGA.txt"),  sep="\t", quote=F, row.names=F, col.names=T, append=F)
  }
  
  #Run MR-MEGA (4 cohorts so 1 PCs max)
  system(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/MRMEGA/MR-MEGA -i mr-mega2.in --qt --pc 1 --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/mrmega_results_005_pheno_",j))
  
  mrmegares <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop/mrmega_results_005_pheno_",j,".result")))
  head(mrmegares)
  
  keep <- mrmegares[is.na(mrmegares$Comments),]
  sighits <- keep[,c(1,16)]
  colnames(sighits) <- c("ID","P")
  sighits <- sighits[sighits$P<5e-08,]
  
  meta_res[j,3] <- sum(sighits$ID %in% cau.snp.infor$x)/ncausal
  
  # x <- strsplit(sighits$ID,":")
  # sighitsid <- data.frame("ID"=sighits$ID,"CHR"=sapply(x,"[",1),"position"=sapply(x,"[",2),"P"=sighits$P)
  # cau_snp_res<-data.frame("ID"=cau.snp.infor$x, "rs_id"=cau.snp.infor$rs_id, "CHR"= cau.snp.infor$CHR, "position"=cau.snp.infor$position_GRCh37)
  # colnames(keep)[1]<-"ID"
  # cau_snp_res <- left_join(cau_snp_res,keep, by="ID")
  # cau_snp_res <- cau_snp_res[,c(1,2,3,4,17)]
  # cau_snp_res$MOST_SIG_500KB <- NA
  # cau_snp_res$P_MOST_SIG_500KB <- NA
  # 
  # sighitsincausalreg<-c()
  # 
  # for (l in 1:nrow(cau_snp_res)) {
  #   if(l%%10000 == 0){
  #     print(paste0(l," causal SNPs checked"))
  #   }
  #   chr<-cau.snp.infor$CHR[l]
  #   pos<-cau.snp.infor$position[l]
  #   #Within 500kb
  #   in500kb<-subset(sighitsid, CHR==chr & position>=(pos-500000) & position<=(pos+500000))
  #   if(nrow(in500kb)>0){
  #     cau_snp_res$P_MOST_SIG_500KB[l]<-min(in500kb$P)
  #     cau_snp_res$MOST_SIG_500KB[l]<-in500kb[in500kb$P==min(in500kb$P),1]
  #     sighitsincausalreg <- union(sighitsincausalreg,in500kb$ID)
  #   }
  # }
  # setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop")
  # fwrite(cau_snp_res, file = paste0("cau_snp_res_005_MRMEGA_pheno",j,".txt"), sep = "\t")
  # 
  # meta_res[j,4]<-(ncausal-sum(is.na(cau_snp_res$P_MOST_SIG_500KB)))/ncausal
  # 
  # meta_res[j,5]<-(nrow(sighits)-length(sighitsincausalreg))/nrow(sighits)
  # 
  # print(meta_res)
}

setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/admixpop")
fwrite(meta_res, file = paste0("pan_res_MRMEGA.txt"), sep = "\t")

colMeans(meta_res[,c(3,4,5)])
