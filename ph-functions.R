ph_packages <- c("TDA","plotly","viridis","plot3D")
sapply(ph_packages,require,character.only=T)

source('./RKDE.R')

sample.annulus <- function(num.points, inner.radius, outer.radius){
  theta <- runif(num.points) * 2 * pi
  radius <- sqrt(runif(num.points, inner.radius^2, outer.radius^2))
  x <- radius * cos(theta)
  y <- radius * sin(theta)
  cbind(x,y)
}

M0 <- function(X,k=7){
  return(k/nrow(X))
}

rescale_dgm <- function(dgm){
  d <- dgm
  m0 <- max(dgm$diagram[,2:3])
  d$diagram[,2:3] <- d$diagram[,2:3]/m0
  return(d)
}

total.pers <- function(Diag){
  d <- rescale_dgm(Diag)$diagram
  x <- abs(d[,3]-d[,1])
  # y <- ifelse(x>mean(x),x,0)
  pers <- sum(x)
  return(pers)
}

# Plot filter functions

make.grid <- function(X,by){
  Xlim <- c(range(X[,1])[1]-1/20*(range(X[,1])[2]-range(X[,1])[1]), range(X[,1])[2]+1/20*(range(X[,1])[2]-range(X[,1])[1]))
  Ylim <- c(range(X[,2])[1]-1/20*(range(X[,2])[2]-range(X[,2])[1]), range(X[,2])[2]+1/20*(range(X[,2])[2]-range(X[,2])[1]))
  gx <- seq.int(Xlim[1L], Xlim[2L], by=by)
  gy <- seq.int(Ylim[1L], Ylim[2L], by=by)
  Grid <- expand.grid(gx,gy)
  return(Grid)
}

plot.dtm <- function(X,by,m0){
  Grid <- make.grid(X,by)
  xseq <- unique(Grid[,1]) %>% sort()
  yseq <- unique(Grid[,2]) %>% sort()
  z <- matrix(dtm(X,Grid,m0=m0),nrow = length(xseq),ncol=length(yseq))
  kludge <- median(z)
  p <- plot_ly(x=xseq,y=xseq,z=z,type="surface",scene='scene1',
               colors = viridis(n=20,option="C")) %>%
    layout(title=paste("DTM, m0=",m0))
  p %>% add_markers(x=X[,1],y=X[,2],z=rep(-kludge,nrow(X)),size=0.1)
}

plot.kdist <- function(X,by,m0){
  Grid <- make.grid(X,by)
  xseq <- unique(Grid[,1]) %>% sort()
  yseq <- unique(Grid[,2]) %>% sort()
  z <- matrix(dtm(X,Grid,m0=m0),nrow = length(xseq),ncol=length(yseq))
  kludge <- median(z)
  p <- plot_ly(x=xseq,y=xseq,z=z,type="surface",scene='scene1',
               colors = viridis(n=20,option="C")) %>%
    layout(title=paste("DTM, m0=",m0))
  p %>% add_markers(x=X[,1],y=X[,2],z=rep(-kludge,nrow(X)),size=0.1)
}

plot.rkdist <- function(X,by,m0){
  Grid <- make.grid(X,by)
  xseq <- unique(Grid[,1]) %>% sort()
  yseq <- unique(Grid[,2]) %>% sort()
  z <- matrix(dtm(X,Grid,m0=m0),nrow = length(xseq),ncol=length(yseq))
  kludge <- median(z)
  p <- plot_ly(x=xseq,y=xseq,z=z,type="surface",scene='scene1',
               colors = viridis(n=20,option="C")) %>%
    layout(title=paste("DTM, m0=",m0))
  p %>% add_markers(x=X[,1],y=X[,2],z=rep(-kludge,nrow(X)),size=0.1)
}

plot.kde <- function(X,by,H){
  Grid <- make.grid(X,by)
  xseq <- unique(Grid[,1]) %>% sort()
  yseq <- unique(Grid[,2]) %>% sort()
  z <- matrix(kde(X,Grid,h=H),nrow = length(xseq),ncol=length(yseq))
  kludge <- max(z)
  p <- plot_ly(x=xseq,y=xseq,z=z,type="surface",scene='scene1',
               colors = viridis(n=20,option="C")) %>%
    layout(title=paste("KDE, h=",H))
  p %>% add_markers(x=X[,1],y=X[,2],z=rep(-kludge,nrow(X)),size=0.1)
}

plot.rkde <- function(X,by,H){
  Grid <- make.grid(X,by)
  xseq <- unique(Grid[,1]) %>% sort()
  yseq <- unique(Grid[,2]) %>% sort()
  z <- matrix(kde(X,Grid,h=H,weight = (1e-10+rkde.w(X,H)[[1]])),nrow = length(xseq),ncol=length(yseq))
  kludge <- max(z)
  p <- plot_ly(x=xseq,y=xseq,z=z,type="surface",scene='scene1',
               colors = viridis(n=20,option="C")) %>%
    layout(title=paste("RKDE, h=",H))
  p %>% add_markers(x=X[,1],y=X[,2],z=rep(-kludge,nrow(X)),size=0.1)
}



# 2-dimensional persistent homolgy

ph.dtm <- function(X, by, m0){
  Xlim <- c(range(X[,1])[1]-1/20*(range(X[,1])[2]-range(X[,1])[1]), range(X[,1])[2]+1/20*(range(X[,1])[2]-range(X[,1])[1]))
  Ylim <- c(range(X[,2])[1]-1/20*(range(X[,2])[2]-range(X[,2])[1]), range(X[,2])[2]+1/20*(range(X[,2])[2]-range(X[,2])[1]))
  Diag <- gridDiag(X = X, FUN = dtm, m0=m0, lim = cbind(Xlim, Ylim),
                   by = by, sublevel = TRUE, library = "Dionysus",
                   printProgress = TRUE, location=TRUE)
  return(Diag)
}

ph.kde <- function(X, by, H){
  Xlim <- c(range(X[,1])[1]-1/20*(range(X[,1])[2]-range(X[,1])[1]), range(X[,1])[2]+1/20*(range(X[,1])[2]-range(X[,1])[1]))
  Ylim <- c(range(X[,2])[1]-1/20*(range(X[,2])[2]-range(X[,2])[1]), range(X[,2])[2]+1/20*(range(X[,2])[2]-range(X[,2])[1]))
  Diag <- gridDiag(X = X, FUN = kde, h=H, lim = cbind(Xlim, Ylim),
                   by = by, sublevel = FALSE, library = "Dionysus",
                   printProgress = TRUE, location=TRUE)
  return(Diag)
}

ph.rkde <- function(X, by, H){
  Xlim <- c(range(X[,1])[1]-1/20*(range(X[,1])[2]-range(X[,1])[1]), range(X[,1])[2]+1/20*(range(X[,1])[2]-range(X[,1])[1]))
  Ylim <- c(range(X[,2])[1]-1/20*(range(X[,2])[2]-range(X[,2])[1]), range(X[,2])[2]+1/20*(range(X[,2])[2]-range(X[,2])[1]))
  Diag <- gridDiag(X = X, FUN = rkde, h=H, lim = cbind(Xlim, Ylim),
                   by = by, sublevel = FALSE, library = "Dionysus",
                   printProgress = TRUE, location=TRUE)
  return(Diag)
}

ph.rkde2 <- function(X, by, H){
  Xlim <- c(range(X[,1])[1]-1/20*(range(X[,1])[2]-range(X[,1])[1]), range(X[,1])[2]+1/20*(range(X[,1])[2]-range(X[,1])[1]))
  Ylim <- c(range(X[,2])[1]-1/20*(range(X[,2])[2]-range(X[,2])[1]), range(X[,2])[2]+1/20*(range(X[,2])[2]-range(X[,2])[1]))
  Diag <- gridDiag(X = X, FUN = kde, h=H, lim = cbind(Xlim, Ylim),weight=(1e-10+rkde.w(X,H)[[1]]),
                   by = by, sublevel = FALSE, library = "Dionysus",
                   printProgress = TRUE, location=TRUE)
  return(Diag)
}

ph.rkdist <- function(X, by, H){
  Xlim <- c(range(X[,1])[1]-1/20*(range(X[,1])[2]-range(X[,1])[1]), range(X[,1])[2]+1/20*(range(X[,1])[2]-range(X[,1])[1]))
  Ylim <- c(range(X[,2])[1]-1/20*(range(X[,2])[2]-range(X[,2])[1]), range(X[,2])[2]+1/20*(range(X[,2])[2]-range(X[,2])[1]))
  Diag <- gridDiag(X = X, FUN = kernelDist, h=H, weight=(1e-10+rkde.w(X,H)[[1]]),
                   lim = cbind(Xlim, Ylim),
                   by = by, sublevel = TRUE, library = "Dionysus",
                   printProgress = TRUE, location=TRUE)
  return(Diag)
}

ph.kdist <- function(X, by, H){
  Xlim <- c(range(X[,1])[1]-1/20*(range(X[,1])[2]-range(X[,1])[1]), range(X[,1])[2]+1/20*(range(X[,1])[2]-range(X[,1])[1]))
  Ylim <- c(range(X[,2])[1]-1/20*(range(X[,2])[2]-range(X[,2])[1]), range(X[,2])[2]+1/20*(range(X[,2])[2]-range(X[,2])[1]))
  Diag <- gridDiag(X = X, FUN = kernelDist, h=H,
                   lim = cbind(Xlim, Ylim),
                   by = by, sublevel = TRUE, library = "Dionysus",
                   printProgress = TRUE, location=TRUE)
  return(Diag)
}

ph.distfun <- function(X, by, H){
  Xlim <- c(range(X[,1])[1]-1/20*(range(X[,1])[2]-range(X[,1])[1]), range(X[,1])[2]+1/20*(range(X[,1])[2]-range(X[,1])[1]))
  Ylim <- c(range(X[,2])[1]-1/20*(range(X[,2])[2]-range(X[,2])[1]), range(X[,2])[2]+1/20*(range(X[,2])[2]-range(X[,2])[1]))
  Diag <- gridDiag(X = X, FUN = distFct,
                   lim = cbind(Xlim, Ylim),
                   by = by, sublevel = TRUE, library = "Dionysus",
                   printProgress = TRUE, location=TRUE)
  return(Diag)
}



# 3-dimensional Persistent Homology

ph3d.dtm <- function(X, by, m0){
  Xlim <- c(range(X[,1])[1]-1/20*(range(X[,1])[2]-range(X[,1])[1]), range(X[,1])[2]+1/20*(range(X[,1])[2]-range(X[,1])[1]))
  Ylim <- c(range(X[,2])[1]-1/20*(range(X[,2])[2]-range(X[,2])[1]), range(X[,2])[2]+1/20*(range(X[,2])[2]-range(X[,2])[1]))
  Zlim <- c(range(X[,3])[1]-1/20*(range(X[,3])[2]-range(X[,3])[1]), range(X[,3])[2]+1/20*(range(X[,3])[2]-range(X[,3])[1]))
  Diag <- gridDiag(X = X, FUN = dtm, m0=m0, lim = cbind(Xlim, Ylim, Zlim),
                   by = by, sublevel = TRUE, library = "Dionysus",
                   printProgress = TRUE, location=TRUE)
  return(Diag)
}

ph3d.kde <- function(X, by, H){
  Xlim <- c(range(X[,1])[1]-1/20*(range(X[,1])[2]-range(X[,1])[1]), range(X[,1])[2]+1/20*(range(X[,1])[2]-range(X[,1])[1]))
  Ylim <- c(range(X[,2])[1]-1/20*(range(X[,2])[2]-range(X[,2])[1]), range(X[,2])[2]+1/20*(range(X[,2])[2]-range(X[,2])[1]))
  Zlim <- c(range(X[,3])[1]-1/20*(range(X[,3])[2]-range(X[,3])[1]), range(X[,3])[2]+1/20*(range(X[,3])[2]-range(X[,3])[1]))
  Diag <- gridDiag(X = X, FUN = kde, h=H, lim = cbind(Xlim, Ylim, Zlim),
                   by = by, sublevel = FALSE, library = "Dionysus",
                   printProgress = TRUE, location=TRUE)
  return(Diag)
}

ph3d.rkde <- function(X, by, H){
  Xlim <- c(range(X[,1])[1]-1/20*(range(X[,1])[2]-range(X[,1])[1]), range(X[,1])[2]+1/20*(range(X[,1])[2]-range(X[,1])[1]))
  Ylim <- c(range(X[,2])[1]-1/20*(range(X[,2])[2]-range(X[,2])[1]), range(X[,2])[2]+1/20*(range(X[,2])[2]-range(X[,2])[1]))
  Zlim <- c(range(X[,3])[1]-1/20*(range(X[,3])[2]-range(X[,3])[1]), range(X[,3])[2]+1/20*(range(X[,3])[2]-range(X[,3])[1]))
  Diag <- gridDiag(X = X, FUN = rkde, h=H, lim = cbind(Xlim, Ylim, Zlim),
                   by = by, sublevel = FALSE, library = "Dionysus",
                   printProgress = TRUE, location=TRUE)
  return(Diag)
}


