# setwd('/Users/vishwanathgl/Box/TDA/NeurIPS/Examples')
source('./ph-functions.R')

cls <- c('dodgerblue','firebrick1')


ph.distfun <- function(X, by, l=NA){
  Xlim <- c(range(X[,1])[1]-1/20*(range(X[,1])[2]-range(X[,1])[1]), range(X[,1])[2]+1/20*(range(X[,1])[2]-range(X[,1])[1]))
  Ylim <- c(range(X[,2])[1]-1/20*(range(X[,2])[2]-range(X[,2])[1]), range(X[,2])[2]+1/20*(range(X[,2])[2]-range(X[,2])[1]))
  if(is.na(l)){
    lims <- cbind(Xlim, Ylim)
  }
  else{
    lims <- cbind(c(-l,l),c(-l,l))
  }
  Diag <- gridDiag(X = X, FUN = distFct,
                   lim = lims,
                   by = by, sublevel = TRUE, library = "Dionysus",
                   printProgress = TRUE, location=TRUE)
  return(Diag)
}

# Signal
set.seed(2020)
lims <- 6
signal <- circleUnif(40,3)*rnorm(40,1,0.1)
d1 <- ph.distfun(signal,by=0.1,l=lims)

# par(mfrow=c(1,2))

plot(signal,asp=1,pch=20,col=alpha('black',0.7),cex=1,
     ylim=c(-lims,lims),xlim=c(-lims,lims))

plot.diagram(d1$diagram,diagLim = c(0,4),
             dimension = 0,pch=20, col=alpha(cls[1],0.8))
plot.diagram(d1$diagram,add = T,
             dimension = 1,pch=17, col=alpha(cls[2],0.8))
legend('bottomright',
       c('0th Order Homology','1st Order Homology'),
       pch=c(20,17),
       col=cls,bty = 'n')


# Observed
set.seed(2020)
noise <- rbind(c(0,0),c(4,-4),matrix(runif(2*12,-lims,lims),ncol=2))
signal <- circleUnif(40,3)
X <- rbind(signal,noise)
d2 <- ph.distfun(X,by=0.1,l=lims)


plot(signal,asp=1,pch=20,col=alpha('black',0.7),cex=1,
     ylim=c(-lims,lims),xlim=c(-lims,lims))
points(noise,asp=1,pch=20,col=alpha('firebrick1',0.7),cex=1.5,
     ylim=c(-lims,lims),xlim=c(-lims,lims))

plot.diagram(d2$diagram,diagLim = c(0,4),
             dimension = 0,pch=20, col=alpha(cls[1],0.8))
plot.diagram(d2$diagram,add = T,
             dimension = 1,pch=17, col=alpha(cls[2],0.8))
legend('bottomright',
       c('0th Order Homology','1st Order Homology'),
       pch=c(20,17),
       col=cls,bty = 'n')

# RKDE PD
d3 <- ph.kde(X,by=0.1,H=bw(X,5))

plot.diagram(d3$diagram,
             dimension = 0,pch=20, col=alpha(cls[1],0.8))
plot.diagram(d3$diagram,add = T,
             dimension = 1,pch=17, col=alpha(cls[2],0.8))
legend('bottomright',
       c('0th Order Homology','1st Order Homology'),
       pch=c(20,17),
       col=cls,bty = 'n')

# DTM PD
d4 <- ph.dtm(X,by=0.1,m0=M0(X,5))

plot.diagram(d4$diagram,
             dimension = 0,pch=20, col=alpha(cls[1],0.8))
plot.diagram(d4$diagram,add = T,
             dimension = 1,pch=17, col=alpha(cls[2],0.8))
legend('bottomright',
       c('0th Order Homology','1st Order Homology'),
       pch=c(20,17),
       col=cls,bty = 'n')

