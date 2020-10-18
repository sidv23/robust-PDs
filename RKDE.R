rkde_packages <- c("dplyr","plotrix","spatstat","TDA","hitandrun","functional","Rfast")
sapply(rkde_packages, require, character.only=T)


bw <- function(x,k=3){
  return(x %>% nndist(k=k) %>% median())
}

med_phi <- function(x){
  return((((x>0)*2)-1)*(1/x))
}

med_rho <- function(x){
  return(abs(x))
}

hampel_rho <- function(x,a,b,c){
  return(
    case_when(
      x>=0&x<a ~ 0.5*(as.numeric(x)^2),
      x>=a&x<b ~ (a*as.numeric(x)) - (0.5*(a^2)),
      x>=b&x<c ~ ((0.5*(a/(b-c)))*as.numeric(x-c)^2) + (0.5*a*(b+c-a)),
      x>=c     ~ (0.5*a*(b+c-a)),
      # x<0      ~ stop("Argument must be non-negative"),
    )
  )
}

hampel_phi <- function(x,a,b,c){
  return(
    case_when(
      x>=0&x<a ~ 1,
      x>=a&x<b ~ a/as.numeric(x),
      # x>=b&x<c ~ (a/(c-b))*as.numeric(c-x),
      x>=b&x<c ~ (a/((c-b)*as.numeric(x)))*as.numeric(c-x),
      x>=c     ~ 0,
      # x<0      ~ stop("Argument must be non-negative"),
    )
  )
}

huber_rho <- function(x,a){
  return(
    case_when(
      x>=0&x<=a ~ 0.5*(as.numeric(x)^2),
      x>a       ~ (a*as.numeric(x)) - (0.5*(a^2)),
      # x<0      ~ stop("Argument must be non-negative"),
    )
  )
}

huber_phi <- function(x,a){
  return(
    case_when(
      x>=0&x<=a ~ 1,
      x>a      ~ a/as.numeric(x),
      # x<0      ~ stop("Argument must be non-negative"),
    )
  )
}


norm <- function(x){
  return(sqrt(sum(x^2)))
}

gram <- function(X,h=1){
  return(Rfast::Dist(X,method="euclidean", square=F) %>% dnorm(mean=0,sd=h))
}

rkhs_norm <- function(X,h=1,w=rep(1/nrow(X),nrow(X))){
  n <- nrow(X)
  K <- gram(X,h)
  a <- K%*%w
  b <- t(w)%*%a
  norm <- c(1:n) %>% sapply(function(i)(K[i,i]+b-(2*a[i])))
  return(sqrt(norm))
}

rkde_loss <- function(X,w,h,fun){
  return( sum(fun(rkhs_norm(X,h,w)))/nrow(X) )
}

rkde_w <- function(w,X,h,loss,phi,tolerance=1e-6){
  l.old <- rkde_loss(X,w,h,loss)
  ratio <- 1
  iter <- 0
  while(ratio>tolerance  & iter < 100){
    w.phi <- phi(rkhs_norm(X,h,w))
    w <- w.phi/sum(w.phi)
    l.new <- rkde_loss(X,w,h,loss)
    ratio <- abs((l.new-l.old)/l.old)
    # print(paste("Loss = ",l.new))
    # print(paste("Ratio = ",ratio))
    l.old <- l.new
    iter <- iter+1
  }
  return(w)
}

rkde.w <- function(X,h){
  nx <- nrow(X)
  w <- c(simplex.sample(n=nx,N=1)$samples)
  w.med <- rkde_w(w,X,h,loss=med_rho,phi=med_phi)
  d <- rkhs_norm(X,h,w.med)
  # q <- quantile(d,c(0.5,0.85,0.9))
  q <- quantile(d,c(0.5,0.9,0.95))
  a <- q[1];b <- q[2]; c <- q[3]
  h.rho <- Curry(hampel_rho,a=a,b=b,c=c)
  h.phi <- Curry(hampel_phi,a=a,b=b,c=c)
  w.new <- rkde_w(w.med,X,h,loss=h.rho,phi=h.phi)
  return(list(w.new,w.med))
}

rkde <- function(X,Grid,h){
  w <- rkde.w(X,h)[[1]]
  n <- nrow(Grid)
  nx <- nrow(X)
  normalize <- (sqrt(2*pi)*h)^(ncol(X)-1)
  B <- Rfast::dista(Grid,X,type="euclidean",square=F) %>% dnorm(mean=0,sd=h)
  # return(B%*%w)
  return((B%*%w)/normalize)
}
