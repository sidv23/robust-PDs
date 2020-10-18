spectral_packages <- c('fossil','kknn')
sapply(spectral_packages,require,character.only=T)


pers.image <- function (d1, nbins, dimension, h) 
{
  d1 = d1[d1[, 1] == dimension, 2:3, drop = F]
  
  if(!is.empty(d1)){
    d1[1,2] <- ifelse(is.infinite(d1[1,2]),2,d1[1,2])
    d1[, 2] = d1[, 2] - d1[, 1]
    maxD = max(d1)
    # maxP = min(d1[, 2])
    maxP = max(d1[, 2])
    minD = min(0, min(d1))
    dx = maxD/nbins
    x_lower = seq(minD, maxD, length.out = nbins)
    x_upper = x_lower + dx
    y_lower = seq(0, maxD, length.out = nbins)
    y_upper = y_lower + dx
    PSurface = function(point, maxP) {
      x = point[1]
      y = point[2]
      out1 = pnorm(x_upper, mean = x, sd = h) - pnorm(x_lower, 
                                                      mean = x, sd = h)
      out2 = pnorm(y_upper, mean = y, sd = h) - pnorm(y_lower, 
                                                      mean = y, sd = h)
      wgt = y/maxP * (y < maxP) + 1 * (y >= maxP)
      return(out1 %o% out2 * wgt)
    }
    Psurf_mat = apply(d1, 1, PSurface, maxP = maxP)
    out = apply(Psurf_mat, 1, sum)
    return(matrix(out, nrow = nbins))
  } else {
    return(matrix(0,nrow=nbins,ncol=nbins))
  }
}

bottle <- function(D1,D2,d=0){
  id1 <- which(D1[,1]==d)
  id2 <- which(D2[,1]==d)
  
  if(is.empty(c(id1,id2))){dist <- 0}
  else if(is.empty(id1)){dist <- max(abs(D2[id2,2]-D2[id2,3]))}
  else if(is.empty(id2)){dist <- max(abs(D1[id1,2]-D1[id1,3]))}
  else if(!is.empty(c(id1,id2))){dist <- bottleneck(D1,D2,dimension=d)}
  
  return(dist)
}


wass <- function(D1,D2,p=1,d=0){
  id1 <- which(D1[,1]==d)
  id2 <- which(D2[,1]==d)
  
  if(is.empty(c(id1,id2))){dist <- 0}
  else if(is.empty(id1)){dist <- (abs(D2[id2,2]-D2[id2,3]))^p %>% sum }
  else if(is.empty(id2)){dist <- (abs(D1[id1,2]-D1[id1,3]))^p %>% sum }
  else if(!is.empty(c(id1,id2))){dist <- wasserstein(D1,D2,p = p,dimension=d) }
  
  return(dist^(1/p))
}