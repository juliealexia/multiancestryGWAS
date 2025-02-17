library(data.table)

eth <- c("EUR","AFR","AMR","EAS","SAS")

psamfile <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/ALL/all_chr_AFR.psam")))

for(i in 1:5){
  psamfile <- as.data.frame(fread(paste0("/n/holystore01/LABS/xlin/Everyone/JD_HZ/ALL/all_chr_",eth[i],".psam")))
  psamfile[,1] <- paste0(eth[i],"_",psamfile[,1])
  psamfile[,2] <- paste0(eth[i],"_",psamfile[,2])
  setwd("/n/holystore01/LABS/xlin/Everyone/JD_HZ/ALL/")
  fwrite(psamfile, file = paste0("all_chr_",eth[i],".psam"), sep = "\t", na="NA", quote=FALSE)
}

for(i in 1:5){
  system(paste0("/n/home12/jdias/plink2 --pfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/ALL/all_chr_",eth[i]," --make-bed --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/ALL/all_chr_",eth[i]))
}

system(paste0("/n/home12/jdias/plink --merge-list /n/home12/jdias/mergelisteth.txt --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/ALL/all_chr"))

system(paste0("/n/home12/jdias/plink2 --pca approx --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/ALL/all_chr --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/ALL/pca" ))


library(data.table)

eth <- c("EUR","AFR","AMR","EAS","SAS")

for (i in 1:5) {
  for(j in 1:22){
    system(paste0("cp /n/holystore01/LABS/xlin/Everyone/multi_ancestry_simu/",eth[i],"_mega/chr",j,".bed /n/holystore01/LABS/xlin/Everyone/JD_HZ/ALL"))
    system(paste0("cp /n/holystore01/LABS/xlin/Everyone/multi_ancestry_simu/",eth[i],"_mega/chr",j,".bim /n/holystore01/LABS/xlin/Everyone/JD_HZ/ALL"))
    system(paste0("mv -v chr",j,".bed chr",j,"_",eth[i],".bed"))
    system(paste0("mv -v chr",j,".bim chr",j,"_",eth[i],".bim"))
  }
}

for(j in 1:22){
  system(paste0("/n/home12/jdias/plink  --merge-list /n/home12/jdias/mergelist",j,".txt --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/ALL/chr",j,"_ALL"))
  for(i in 1:5){
    system(paste0("mv -v chr",j,"_",eth[i],".fam chr",j+1,"_",eth[i],".fam"))
  }
}

for(j in 1:22){
  system(paste0("/n/home12/jdias/plink2 --pca approx --bfile /n/holystore01/LABS/xlin/Everyone/JD_HZ/ALL/chr1_ALL --out /n/holystore01/LABS/xlin/Everyone/JD_HZ/ALL/pca_chr",j ))
}

