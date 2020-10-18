
homedir <- '/Users/vishwanathgl/Box/TDA/NeurIPS/'
setwd(paste(homedir,'/circles',sep=""))

source(paste(homedir,'ph-functions.R',sep=""))
source('count-functions.R')


packages <- c("plotrix","dplyr","SyNet","spatstat",
              "TDA","hitandrun","functional","Rfast",
              'snow','doParallel','foreach','e1071','wesanderson')
sapply(packages,require,character.only=T)

# B <- 10; k <- 0.2
# sq <- seq(0.001,0.21,0.02)

B <- 15; k <- 0.2
sq <- seq(0.001,0.32,0.02)

wes_pals <- c('BottleRocket2','Darjeeling1','Moonrise3')
pal <- wes_palette(wes_pals[2],n=3,type='discrete')


#Final Plots

for(K in c(5,7,10)){
  print(K)
  load(paste('N.RData'  ,sep=""))
  load(paste('data-k',K,'/kde_dgms.RData' ,sep=""))
  load(paste('data-k',K,'/rkde_dgms.RData',sep=""))
  load(paste('data-k',K,'/dtm_dgms.RData' ,sep=""))

  set.seed(2020)
  l.kde  <- sapply(sq,function(x)topsvm(kde_dgms ,N=N,h=x,bins=B,k=k,type='eps',runs=25))
  l.rkde <- sapply(sq,function(x)topsvm(rkde_dgms,N=N,h=x,bins=B,k=k,type='eps',runs=25))
  lines.plot2(sq,l.kde,l.rkde,pal=pal)
  
  filename=paste('./plots/lines-k',K,'.pdf',sep="")
  Lines.plot2(filename,sq,l.kde,l.rkde,pal=pal)
  
  l.dtm  <- sapply(sq,function(x)topsvm(dtm_dgms ,N=N,h=x,bins=B,k=k,type='eps',runs=25))
  lines.plot(sq,l.kde,l.rkde,l.dtm,pal=pal)
  
  filename=paste('./plots/3lines-k',K,'.pdf',sep="")
  Lines.plot(filename,sq,l.kde,l.rkde,l.dtm,pal=pal)
  
  filename=paste('./plots/lines-k',K,sep="")
  save(l.kde,l.rkde,l.dtm,file = paste(filename,'.RData',sep=''))
}

