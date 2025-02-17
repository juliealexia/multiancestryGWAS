
library(data.table)

for(eth in c("AFR","AMR","EAS","EUR","SAS","ALL")){
  system(paste0("/n/home12/jdias/plink2 --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/",eth," q --geno 0.05 --hwe 0.0000001 --make-bed --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/",eth,"_step1"))
}

#Binary
pheno_of_int <- c("Asthma","CAD","T2D","Breast","Prostate")

for(pheno1 in pheno_of_int){
  for(eth in c("ALL","AMR","AFR","EAS","EUR","SAS")){
    pheno <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/",pheno1,"_",eth,".pheno")))
    fwrite(pheno, paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/",pheno1,"_",eth,".pheno"),sep="\t")
  }
}

for(pheno1 in pheno_of_int){
  for(eth in c("AFR","AMR","EAS","EUR","SAS","ALL")){
    
    phenofile=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/",pheno_of_int,"_",eth,".pheno")
    
    infolder="/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/"
    runstep1=function(eth="ALL")
    {
      bedfile=paste0(infolder,eth,"_step1")
      covarfile=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_",eth,".cov")
      outprefix=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_",eth,".out")
      cmd=paste0("regenie --step 1 --force-step1 --bed ",bedfile," --covarFile ",covarfile,
                 " --phenoFile ",phenofile," --bsize 1000 --exclude /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/snp_rmv.txt --bt --lowmem --lowmem-prefix /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/",eth,"_rg --loocv --out ",outprefix)
      system(cmd)
    }
    runstep1(eth)
    
    runstep2=function(eth="ALL")
    {
      covarfile=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_",eth,".cov")
      
      predprefix=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_",eth,".out")
      
      bedfile=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/",eth)
      outprefix=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/plink2_step2_firth_out_",eth) 
      cmd=paste0("regenie --step 2 --bed ",bedfile," --covarFile ",covarfile,
                 " --phenoFile ",phenofile," --bsize 1000 --bt --firth --approx --firth-se --pThresh 0.01 --pred ",predprefix,"_pred.list --out ",outprefix)
      system(cmd)
    }
    runstep2(eth)
  }
}

#Continuous
pheno_of_int <- c("height","Waist","HDL","LDL","TC","Calcium","Creatinine","EGFR")

for(pheno1 in pheno_of_int){
  for(eth in c("AFR","AMR","EAS","EUR","SAS","ALL")){
    
    phenofile=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/",pheno_of_int,"_",eth,".pheno")
    
    infolder="/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/"
    runstep1=function(eth="ALL")
    {
      bedfile=paste0(infolder,eth,"_step1")
      covarfile=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_",eth,".cov")
      outprefix=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_",eth,".out")
      cmd=paste0("regenie --step 1 --force-step1 --bed ",bedfile," --covarFile ",covarfile,
                 " --phenoFile ",phenofile," --bsize 1000 --exclude /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/snp_rmv.txt --qt --lowmem --lowmem-prefix /n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/",eth,"_rg --loocv --out ",outprefix)
      system(cmd)
    }
    runstep1(eth)
    
    runstep2=function(eth="ALL")
    {
      covarfile=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/pca_",eth,".cov")
      
      predprefix=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/summary_",eth,".out")
      
      bedfile=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/",eth)
      outprefix=paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/UKB/Results/plink2_step2_firth_out_",eth) 
      cmd=paste0("regenie --step 2 --bed ",bedfile," --covarFile ",covarfile,
                 " --phenoFile ",phenofile," --bsize 1000 --bt --firth --approx --firth-se --pThresh 0.01 --pred ",predprefix,"_pred.list --out ",outprefix)
      system(cmd)
    }
    runstep2(eth)
  }
}
