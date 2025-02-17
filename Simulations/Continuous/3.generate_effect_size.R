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

load("/n/holystore01/LABS/xlin/Everyone/JD_HZ/cau.snp.infor.rdata")

n.total.snp <- nrow(cau.snp.infor)

cur.dir <- "/n/holystore01/LABS/xlin/Everyone/JD_HZ/"

sigma = 0.4
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

save(Sigma,file = paste0(cur.dir,"causal_Sigma.rdata"))

load("/n/holystore01/LABS/xlin/Everyone/JD_HZ/causal_Sigma.rdata")

library(mvtnorm)
#beta represent standarize scale effect-size
set.seed(666)
beta =  rmvnorm(n.total.snp,c(0,0,0,0,0),sigma=Sigma)

#diagonal and off diag are the same i..e symmetric with 0.4/number of causal snps
colnames(beta) <- paste0("beta_",c("EUR","AFR","AMR","EAS","SAS"))

eth <- c("EUR","AFR","AMR","EAS","SAS")
for (i in 1:5){
  select.cau <- cbind(cau.snp.infor[,1],beta[,i])
  write.table(select.cau,file = paste0(cur.dir,eth[i],"/select.cau_rho"),row.names = F,col.names = F,quote=F)
}



