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
  else if(!is.empty(c(id1,id2))){dist <- TDA::bottleneck(D1,D2,dimension=d)}
  
  return(dist)
}


wass <- function(D1,D2,p=1,d=0){
  id1 <- which(D1[,1]==d)
  id2 <- which(D2[,1]==d)
  
  if(is.empty(c(id1,id2))){dist <- 0}
  else if(is.empty(id1)){dist <- (abs(D2[id2,2]-D2[id2,3]))^p %>% sum }
  else if(is.empty(id2)){dist <- (abs(D1[id1,2]-D1[id1,3]))^p %>% sum }
  else if(!is.empty(c(id1,id2))){dist <- TDA::wasserstein(D1,D2,p = p,dimension=d) }
  
  return(dist^(1/p))
}



# Persistence Image Distance Matrix

pimg.dist.matrix <- function(dgms){
  ad <- dgms
  a.img <- list()
  AD <-  AD1 <- AD2 <- list()
  
  for(k in 1:length(dims)){
    print(paste('Computing: Homology',k-1))
    a.img[[k]] <- lapply(ad,function(x)pers.image(x,nbins = 20,dimension = dims[k],h=0.1))
    
    AD[[k]] <- outer(1:150,1:150,Vectorize(function(i,j){
      max(abs(a.img[[k]][[i]]-a.img[[k]][[j]]))
    }))
    
    AD1[[k]] <- outer(1:150,1:150,Vectorize(function(i,j){
      (abs(a.img[[k]][[i]]-a.img[[k]][[j]])) %>% abs %>% sum
    }))
    
    AD2[[k]] <- outer(1:150,1:150,Vectorize(function(i,j){
      (abs(a.img[[k]][[i]]-a.img[[k]][[j]]))^2 %>% sum %>% sqrt
    }))
    
    rownames(AD[[k]]) <- rownames(AD1[[k]]) <- rownames(AD2[[k]]) <- cls
  }
  return(list(AD,AD1,AD2))
}



# Persistence Diagram Distance Matrix
dgm.dist.matrix <- function(dgms){
  rd <- dgms
  
  RD <-  RD1 <- RD2 <- list()
  
  for(k in 1:length(dims)){
    
    print(paste('Computing: Homology',k-1))
    
    RD[[k]] <- outer(1:150,1:150,Vectorize(function(i,j){
      bottle(rd[[i]],rd[[j]],d = dims[k])
    }))
    
    RD1[[k]] <- outer(1:150,1:150,Vectorize(function(i,j){
      wass(rd[[i]],rd[[j]],p=1,d = dims[k])
    }))
    
    RD2[[k]] <- outer(1:150,1:150,Vectorize(function(i,j){
      wass(rd[[i]],rd[[j]],p=2,d = dims[k])
    }))
    
    rownames(RD[[k]]) <- rownames(RD1[[k]]) <- rownames(RD2[[k]]) <- cls
  }
  return(list(RD,RD1,RD2))
}


clust_rand_index <- function(n,cls,numclasses,Delta.dgms,Delta.imgs){
  
  # Transform Diagrams and Images
  AD <- Delta.imgs[[1]]
  AD1 <- Delta.imgs[[2]]
  AD2 <- Delta.imgs[[3]]
  
  A.Dist  <-  pmax(AD[[1]],AD[[2]],AD[[3]])
  A1.Dist <- pmax(AD1[[1]],AD1[[2]],AD1[[3]])
  A2.Dist <- pmax(AD2[[1]],AD2[[2]],AD2[[3]])
  
  RD  <- Delta.dgms[[1]]
  RD1 <- Delta.dgms[[2]]
  RD2 <- Delta.dgms[[3]]
  R.Dist  <-  pmax(RD[[1]],RD[[2]],RD[[3]])
  R1.Dist <- pmax(RD1[[1]],RD1[[2]],RD1[[3]])
  R2.Dist <- pmax(RD2[[1]],RD2[[2]],RD2[[3]])
  
  
  # Convert Class to numeric
  ncls <- numclasses[cls]
  
  noclass <- rep(0,length(cls))
  tem <- list(ad=c(),rd=c())
  temp <- list(tem,tem,tem)
  rand <- list(temp,temp,temp,temp)
  
  for(i in 1:length(rand)){
    
    for(j in 1:length(temp)){
      
      if(i!=4){
        if(j==1){
          
          for(k in 1:n){
            scc1 <- specClust(AD1[[i]],6)
            rand[[i]][[j]]$ad[k] <- rand.index(scc1$cluster,ncls)
            
            scc2 <- specClust(RD1[[i]],6)
            rand[[i]][[j]]$rd[k] <- rand.index(scc2$cluster,ncls)
          }
          
        } else if(j==2){
          
          for(k in 1:n){
            scc1 <- specClust(AD2[[i]],6)
            rand[[i]][[j]]$ad[k] <- rand.index(scc1$cluster,ncls)
            
            scc2 <- specClust(RD2[[i]],6)
            rand[[i]][[j]]$rd[k] <- rand.index(scc2$cluster,ncls)
          }
          
        } else if(j==3){
          
          for(k in 1:n){
            scc1 <- specClust(AD[[i]],6)
            rand[[i]][[j]]$ad[k] <- rand.index(scc1$cluster,ncls)
            
            scc2 <- specClust(RD[[i]],6)
            rand[[i]][[j]]$rd[k] <- rand.index(scc2$cluster,ncls)
          }
          
        }
      } else if(i==4){
        
        if(j==1){
          
          for(k in 1:n){
            scc1 <- specClust(A1.Dist,6)
            rand[[i]][[j]]$ad[k] <- rand.index(scc1$cluster,ncls)
            
            scc2 <- specClust(R1.Dist,6)
            rand[[i]][[j]]$rd[k] <- rand.index(scc2$cluster,ncls)
          }
          
        } else if(j==2){
          
          for(k in 1:n){
            scc1 <- specClust(A2.Dist,6)
            rand[[i]][[j]]$ad[k] <- rand.index(scc1$cluster,ncls)
            
            scc2 <- specClust(R2.Dist,6)
            rand[[i]][[j]]$rd[k] <- rand.index(scc2$cluster,ncls)
          }
          
        } else if(j==3){
          
          for(k in 1:n){
            scc1 <- specClust(A.Dist,6)
            rand[[i]][[j]]$ad[k] <- rand.index(scc1$cluster,ncls)
            
            scc2 <- specClust(R.Dist,6)
            rand[[i]][[j]]$rd[k] <- rand.index(scc2$cluster,ncls)
          }
          
        }
      }
    }
  }
  
  RD.results <- matrix(NA,nrow=length(rand),ncol=length(temp))
  rownames(RD.results) <- c('H0','H1','H2','C')
  colnames(RD.results) <- c('W1','W2','Winf')
  
  for(i in 1:nrow(RD.results)){
    for(j in 1:ncol(RD.results)){
      RD.results[i,j] <- rand[[i]][[j]]$rd %>% max
    }
  }
  
  AD.results <- matrix(NA,nrow=length(rand),ncol=length(temp))
  rownames(AD.results) <- c('H0','H1','H2','C')
  colnames(AD.results) <- c('W1','W2','Winf')
  
  for(i in 1:nrow(RD.results)){
    for(j in 1:ncol(RD.results)){
      AD.results[i,j] <- rand[[i]][[j]]$ad %>% max
    }
  }
  
  results <- cbind(RD.results,AD.results)
  
  return(list(rand=rand,summary=results))
}

