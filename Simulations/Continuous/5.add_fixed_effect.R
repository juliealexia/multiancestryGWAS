library(data.table)

eth <- c("EUR","AFR","AMR","EAS","SAS")

set.seed(666)
fixdefx<-rnorm(5, mean=0, sd=1)

#Adding fixed effect by ethnicity
for(i in 1:5){
  phendata <- as.data.frame(fread(paste0("/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/",eth[i],"_pcpheno.phen")))
  newphen <- cbind(phendata[1:2],phendata[,3:12]+fixdefx[i])
  fwrite(newphen, file = "/n/holylfs05/LABS/kraft_lab/Lab/julied/pcpheno/",eth[i],"_pcpheno.phen", sep = "\t")
}

