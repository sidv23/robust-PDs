bottleneck_packages <- c('snow','doParallel','foreach')
sapply(bottleneck_packages,require,character.only=T)

# Function to simulate Data
simulate_data <- function(n=500,p=0.3,plt=F){
  rad <- 4
  n1 <- n;         n2 <- n/2;         n3 <- 1
  n4 <- p*(n1+n2)
  m1 <- c(-1,-1);    m2 <- c(2,2); m3 <- c(0.15,0.15)
  
  x1 <- m1+circleUnif(n1,rad-2)*rnorm(n1,1,0.1)
  x2 <- m2+circleUnif(n2,rad-3)*rnorm(n2,1,0.07)
  
  signal <- rbind(x1,x2)
  noise <- cbind(runif(n4,-rad,rad),runif(n4,-rad,rad))
  
  X <- rbind(signal,noise)
  
  if(plt==T){
    plot(X,asp=1,pch=".")
    points(signal,pch='*',col="blue")
    points(noise,pch='*',col="red")
  }
  return(list(X,signal))
}





