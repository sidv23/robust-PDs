homedir <- '/storage/work/s/suv87/tda/NeurIPS/NeurIPS/'
setwd(homedir)
source(paste(homedir,'ph-functions.R',sep=""))
source(paste(homedir,'Circles/circles-svm-functions.R',sep=""))

packages <- c('e1071')
sapply(packages,require,character.only=T)


#Standard
B <- 20; k <- 0.2
sq <- seq(0.001,0.23,0.01)

wes_pals <- c('BottleRocket2','Darjeeling1','Moonrise3')
pal <- c('firebrick1','dodgerblue','grey40')


#Final Plots
for(K in c(5,7)){
  print(K)
  load(paste(homedir,'Circles/','N.RData'  ,sep=""))
  load(paste(homedir,'Circles/','data-k',K,'/kde_dgms.RData' ,sep=""))
  load(paste(homedir,'Circles/','data-k',K,'/rkde_dgms.RData',sep=""))
  load(paste(homedir,'Circles/','data-k',K,'/dtm_dgms.RData' ,sep=""))
  
  set.seed(2020)
  l.kde  <- sapply(sq,function(x)topsvm(kde_dgms ,N=N,h=x,bins=B,k=k,type='eps',runs=200))
  l.rkde <- sapply(sq,function(x)topsvm(rkde_dgms,N=N,h=x,bins=B,k=k,type='eps',runs=200))
  lines.plot(sq,l.kde,l.rkde,pal=pal)
  
  filename=paste(homedir,'Circles/plots/lines-k',K,'.pdf',sep="")
  lines.save(filename,sq,l.kde,l.rkde,pal=pal)
  
  filename=paste('Circles/plots/lines-k',K,sep="")
  df <- data.frame(kde=c(l.kde[1,]),sd.kde=c(l.kde[2,]),
                   rkde=c(l.rkde[1,]),sd.rkde=c(l.rkde[2,]))
  write.csv(df,file = paste(filename,'.csv',sep=''))
}

