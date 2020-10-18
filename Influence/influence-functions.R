# Influence Functions

gen.sim <- function(signal,r,F1,F2,F3,D1,D2,D3){
  l1 <- l2 <- l3 <- w1 <- w2 <- w3 <- c()
  b1 <- b2 <- b3 <- W1 <- W2 <- W3 <- B1 <- B2 <- B3 <- c()
  for(i in 1:50){
    X <- rbind(signal,sample.annulus(100,r,r+1))
    d1 <- ph.rkde2(X,by=0.5,H)
    d2 <- ph.kde(X,by=0.5,H)
    d3 <- ph.dtm(X,by=0.5,M)
    
    f1 <- kde(X,G,H,weight=(1e-10+rkde.w(X,H)[[1]]))
    f2 <- kde(X,G,H)
    f3 <- dtm(X,G,M)
    
    l1[i] <- max(abs(f1-F1))
    l2[i] <- max(abs(f2-F2))
    l3[i] <- max(abs(f3-F3))
    
    w1[i] <- TDA::wasserstein(d1$diagram,D1$diagram,dimension = 0,p=1)
    w2[i] <- TDA::wasserstein(d2$diagram,D2$diagram,dimension = 0,p=1)
    w3[i] <- TDA::wasserstein(d3$diagram,D3$diagram,dimension = 0,p=1)
    
    b1[i] <- TDA::bottleneck(d1$diagram,D1$diagram,dimension = 0)
    b2[i] <- TDA::bottleneck(d2$diagram,D2$diagram,dimension = 0)
    b3[i] <- TDA::bottleneck(d3$diagram,D3$diagram,dimension = 0)
    
    W1[i] <- TDA::wasserstein(d1$diagram,D1$diagram,dimension = 1,p=1)
    W2[i] <- TDA::wasserstein(d2$diagram,D2$diagram,dimension = 1,p=1)
    W3[i] <- TDA::wasserstein(d3$diagram,D3$diagram,dimension = 1,p=1)
    
    B1[i] <- TDA::bottleneck(d1$diagram,D1$diagram,dimension = 1)
    B2[i] <- TDA::bottleneck(d2$diagram,D2$diagram,dimension = 1)
    B3[i] <- TDA::bottleneck(d3$diagram,D3$diagram,dimension = 1)
    
  }
  
  return(list(
    l1=l1,
    l2=l2,
    l3=l3,
    w1=w1,
    w2=w2,
    w3=w3,
    b1=b1,
    b2=b2,
    b3=b3,
    W1=W1,
    W2=W2,
    W3=W3,
    B1=B1,
    B2=B2,
    B3=B3
  )
  )
  
}

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
  plot(domain,t1,col='dodgerblue',type="l",ylim=limy,lwd=2,
       xlab="Distance from Support",ylab = ylab, main=main)
  polygon(x=c(domain,rev(domain)),y=c(t1+s1,rev(t1-s1)),
          col=alpha('dodgerblue',0.05),border = NA)
  lines(domain,t2,col='firebrick1',type="l",ylim=limy,lwd=2)
  polygon(x=c(domain,rev(domain)),y=c(t2+s2,rev(t2-s2)),
          col=alpha('firebrick1',0.05),border = NA)
  
  legend('topleft',legend=c('KDE','RKDE'),col=c('firebrick1','dodgerblue'),lty=c(1,1))
  dev.off()
}

