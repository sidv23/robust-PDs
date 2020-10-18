# Circle Functions

pers.image <- function (d1, nbins, dimension, h)
{
  d1 = d1[d1[, 1] == dimension, 2:3, drop = F]
  d1[, 2] = d1[, 2] - d1[, 1]
  maxD = max(d1)
  maxP = min(d1[, 2])
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
}


dgms_to_imgs <- function(X,h=0.05,bins=20,rescale=F){
  img0 <- img1 <- list()
  for(i in 1:length(X)){
    img0[[i]] <- X$diagram %>% pers.image(nbins=bins,dimension=0,h=h)
    if(max(X$diagram[,1]) > 0){    
      img1[[i]] <- X$diagram %>% pers.image(nbins=bins,dimension=1,h=h)
    } else{
      img1[[i]] <- matrix(0,nrow=bins,ncol=bins)
    }
  }
  return(list(img0=img0,img1=img1))
}

split_data <- function(data,N,k=0.25){
  n <- length(N)
  ID <- 1:n
  ids <- sample(1:n,round(k*n),replace=F)
  idx <- ID[-ids]
  train_in <- t(array(unlist(data[idx]), dim=c(length(unlist(data[1])),length(idx))))
  train_out <- unlist(N[idx])
  test_in <- t(array(unlist(data[ids]), dim=c(length(unlist(data[1])),length(ids))))
  test_out <- unlist(N[ids])
  
  return(list(train_in=train_in,
              train_out=train_out,
              test_in=test_in,
              test_out=test_out))
}

img_svm <- function(train,type='eps'){
  svm_model <- svm(train$train_in, train$train_out, type=type, kernel='radial')
  predict.train <- predict(svm_model, train$train_in)
  predict.test <- predict(svm_model, train$test_in)
  return(list(predict.train=predict.train,
              predict.test=predict.test))
}

topsvm <- function(dgms,cls=cls,N=N,
                   k=0.25,bins=20,h=0.1,type='C',
                   homology=c(0,1), runs=2){
  p_imgs <- lapply(dgms,function(x)dgms_to_imgs(x,h=h,bins=bins))
  
  if((0 %in% homology) & (1 %in% homology)){
    p_img <- lapply(p_imgs,function(x) cbind(x$img0[[1]],x$img1[[1]]) )
  } else if((0 %in% homology) & !(1 %in% homology)){
    p_img <- lapply(p_imgs,function(x) x$img0[[1]] )
  } else if(!(0 %in% homology) & (1 %in% homology)){
    p_img <- lapply(p_imgs,function(x) x$img1[[1]] )
  }
  
  error <- c()
  for(j in 1:runs){
    print(paste('Run = ',j))
    train <- split_data(p_img,N,k)
    ratio <- length(train$test_out) / length(train$train_out)
    results <- img_svm(train,type=type)
    if(type=='C'){
      error[j] <- sum(train$test_out != results$predict.test)/ length(train$test_out)
    } else{
      error[j] <- mean((train$test_out-results$predict.test)^2)
    }
  }
  
  return(c(mean(error),sd(error)))
  
}

lines.plot <- function(domain,l1,l2,
                       ylab="",main="",pal){
  limsy <- max(c(l1[1,],l2[1,]))
  plot(sq,l1[1,],type="l",ylim=c(0,limsy),col=pal[1],lwd=3,
       ylab="Predicted Mean Squared Error",xlab="Persistence Image Bandwidth")
  lines(sq,l2[1,],type="l",ylim=c(0,1),col=pal[2],lwd=3)
  legend('topright',
         legend = c('KDE','RKDE'),
         col=pal[1:2],
         lty= c(1,1),lwd=c(3,3)
  )
  # dev.off()
}


lines.save <- function(filename,domain,l1,l2,
                       ylab="",main="",pal){
  pdf(file=filename,height=4,width=4)
  par(oma=c(0,0,0,0),mar=c(4,4,1,1))
  limsy <- max(c(l1[1,],l2[1,]))
  plot(sq,l1[1,],type="l",ylim=c(0,limsy),col=pal[1],lwd=3,
       ylab="Predicted Mean Squared Error",xlab="Persistence Image Bandwidth")
  lines(sq,l2[1,],type="l",ylim=c(0,1),col=pal[2],lwd=3)

  legend('topright',
         legend = c('KDE','RKDE'),
         col=pal[1:2],
         lty= c(1,1),lwd=c(3,3)
  )
  dev.off()
}
