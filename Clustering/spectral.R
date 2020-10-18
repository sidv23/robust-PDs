# Preliminaries
root <- "/path/to/root/directory/"
setwd(root)

source('./RKDE.R')
source('./ph-functions.R')
source('./Clustering/spectral-functions.R')

####################################################################################

dims <- c(0,1,2)

load('./Clustering/data/X.RData')
load('./Clustering/data/cls.RData')
load('./Clustering/data/dist-dgms.RData')
load('./Clustering/data/rkde-dgms.RData')

ad <- lapply(dist.dgms,function(x)x$diagram)
rd <- lapply(rkde.dgms,function(x)x$diagram)

####################################################################################

# Persistence Image
Delta.imgs <- pimg.dist.matrix(ad)
AD <- Delta.imgs[[1]]
AD1 <- Delta.imgs[[2]]
AD2 <- Delta.imgs[[3]]

A.Dist  <-  pmax(AD[[1]],AD[[2]],AD[[3]])
A1.Dist <- pmax(AD1[[1]],AD1[[2]],AD1[[3]])
A2.Dist <- pmax(AD2[[1]],AD2[[2]],AD2[[3]])

save(AD,   file='./distances/AD.RData')
save(AD1,  file='./distances/AD1.RData')
save(AD2,  file='./distances/AD2.RData')


####################################################################################

# Persistence Diagrams
Delta.dgms <- dgm.dist.matrix(rd)
RD  <- Delta.dgms[[1]]
RD1 <- Delta.dgms[[2]]
RD2 <- Delta.dgms[[3]]

R.Dist  <-  pmax(RD[[1]],RD[[2]],RD[[3]])
R1.Dist <- pmax(RD1[[1]],RD1[[2]],RD1[[3]])
R2.Dist <- pmax(RD2[[1]],RD2[[2]],RD2[[3]])

save(RD,   file='./distances/RD.RData')
save(RD1,  file='./distances/RD1.RData')
save(RD2,  file='./distances/RD2.RData')


#####  H0   #######################################################################

numclasses <- c(1:6)
names(numclasses) <- c('noise','sphere','torus','circle','3clust','9clust')

set.seed(2020)
n <- 100

scc <- clust_rand_index(n,cls,numclasses,Delta.dgms, Delta.imgs)
write.csv(cbind(RD.results,AD.results),file = './pimg-alt.csv')