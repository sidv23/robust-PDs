# influence-plots
setwd('/storage/work/s/suv87/tda/NeurIPS/NeurIPS')
source('./ph-functions.R')
source('./Influence/influence-functions.R')

packages <- c('snow','doParallel','foreach')
sapply(packages,require,character.only=T)
plot <- graphics::plot


pal <- c('dodgerblue','firebrick1','grey40')
load('z.RData')

R <- seq(10,100,10)

temp1 <- temp2 <- list()
for(i in 1:10){
  temp1[[i]] <- z$b1[[i]] + z$B1[[i]]
  temp2[[i]] <- z$b2[[i]] + z$B2[[i]]
}

filename = 'plots-new/bottleneck.pdf'
influence.plot(filename,R,temp1,temp2,ylab='Influence',main='Bottleneck Influence')
# dev.off()

filename = 'plots-new/supnorm.pdf'
influence.plot(filename,R,z$l1,z$l2,ylab='Influence',main='Sup norm influence')
# dev.off()

filename = 'plots-new/H0-bottleneck.pdf'
influence.plot(filename,R,z$b1,z$b2,ylab='Influence',main='H0: Bottleneck Influence')
# dev.off()

filename = 'plots-new/H1-bottleneck.pdf'
influence.plot(filename,R,z$B1,z$B2,ylab='Influence',main='H1: Bottleneck Influence')
# dev.off()

filename = 'plots-new/H0-W1.pdf'
influence.plot(filename,R,z$w1,z$w2,ylab='Influence',main='H0: W1 Influence')
# dev.off()

filename = 'plots-new/H1-W1.pdf'
influence.plot(filename,R,z$W1,z$W2,ylab='Influence',main='H1: W1 Influence')





######### DTM #############

influence.plot <- function(filename,domain,list1,list2,
                           ylab,main,type='sd',q=1){
  t1 <- sapply(list1,function(x)mean(x))
  t2 <- sapply(list2,function(x)mean(x))
  if(type=='sd'){
    s1 <- sapply(list1,function(x)sd(x))
    s2 <- sapply(list2,function(x)sd(x))
  } else if(type=='quantile'){
    s1 <- sapply(list1,function(x)quantile(x,q))
    s2 <- sapply(list2,function(x)quantile(x,q))
  }
  
  kludge <- max(range(c(s1,s2)))
  limy <- range(c(t1,t2))+c(-kludge,kludge)
  
  pdf(filename,height=4,width = 4)
  par(oma=c(0,0,0,0),mar=c(4,4,1,1))
  plot(domain,t1,col='grey40',type="l",ylim=limy,lwd=2,
       xlab="Distance from Support",ylab = ylab, main=main)
  polygon(x=c(domain,rev(domain)),y=c(t1+s1,rev(t1-s1)),
          col=alpha('grey40',0.05),border = NA)
  lines(domain,t2,col='firebrick1',type="l",ylim=limy,lwd=2)
  polygon(x=c(domain,rev(domain)),y=c(t2+s2,rev(t2-s2)),
          col=alpha('firebrick1',0.05),border = NA)
  
  legend('topleft',legend=c('KDE','DTM'),col=c('firebrick1','grey40'),lty=c(1,1))
  dev.off()
}

R <- seq(10,100,10)

temp3 <- list()
for(i in 1:10){
  temp3[[i]] <- z$b3[[i]] + z$B3[[i]]
}

filename = 'plots-new/dtm-bottleneck.pdf'
influence.plot(filename,R,temp3,temp2,ylab='Influence',main='Bottleneck Influence')
# dev.off()

filename = 'plots-new/dtm-supnorm.pdf'
influence.plot(filename,R,z$l3,z$l2,ylab='Influence',main='Sup norm influence')
# dev.off()

filename = 'plots-new/dtm-H0-bottleneck.pdf'
influence.plot(filename,R,z$b3,z$b2,ylab='Influence',main='H0: Bottleneck Influence')
# dev.off()

filename = 'plots-new/dtm-H1-bottleneck.pdf'
influence.plot(filename,R,z$B3,z$B2,ylab='Influence',main='H1: Bottleneck Influence')
# dev.off()

filename = 'plots-new/dtm-H0-W1.pdf'
influence.plot(filename,R,z$w3,z$w2,ylab='Influence',main='H0: W1 Influence')
# dev.off()

filename = 'plots-new/dtm-H1-W1.pdf'
influence.plot(filename,R,z$W3,z$W2,ylab='Influence',main='H1: W1 Influence')


