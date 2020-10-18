setwd('/storage/work/s/suv87/tda/NeurIPS/NeurIPS/')
source('RKDE.R')
source('ph-functions.R')

packages <- c('e1071','functional','wesanderson')
sapply(packages,require,character.only=T)

# Persistence Image

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


rescale_dgm <- function(dgm){
  d <- dgm
  m0 <- max(dgm$diagram[,2:3])
  d$diagram[,2:3] <- d$diagram[,2:3]/m0
  return(d)
}


dgms_to_imgs <- function(X,h=0.05,bins=20){
  X <- rescale_dgm(X)
  img0 <- img1 <- list()
  for(i in 1:length(X)){
    img0[[i]] <- X$diagram %>% pers.image(nbins=bins,dimension=0,h=h)
    img1[[i]] <- X$diagram %>% pers.image(nbins=bins,dimension=1,h=h)
  }
  return(list(img0=img0,img1=img1))
}


split_data <- function(data,img_ids,border=10){
  train_test_border <- border
  train_in <- t(array(unlist(data[img_ids < train_test_border]), dim=c(length(unlist(data[1])),sum(img_ids < train_test_border))))
  train_out <- unlist(cls[img_ids < train_test_border])
  test_in <- t(array(unlist(data[img_ids >= train_test_border]), dim=c(length(unlist(data[1])),sum(img_ids >= train_test_border))))
  test_out <- unlist(cls[img_ids >= train_test_border])
  return(list(train_in=train_in,
              train_out=train_out,
              test_in=test_in,
              test_out=test_out))
}

split_data_alt <- function(data,img_ids,k=0.25){
  n <- 20
  ID <- 1:n
  ids <- sample(1:n,round(k*n),replace=F)
  idx <- ID[-ids]
  train_in <- t(array(unlist(data[img_ids %in% idx]), dim=c(length(unlist(data[1])),sum(img_ids %in% idx))))
  train_out <- unlist(cls[img_ids %in% idx])
  test_in <- t(array(unlist(data[img_ids %in% ids]), dim=c(length(unlist(data[1])),sum(img_ids %in% ids))))
  test_out <- unlist(cls[img_ids %in% ids])
  return(list(train_in=train_in,
              train_out=train_out,
              test_in=test_in,
              test_out=test_out))
}


img_svm <- function(train,type='C'){
  svm_model <- svm(train$train_in, train$train_out, type=type, kernel='linear', scale = T)
  predict.train <- predict(svm_model, train$train_in)
  predict.test <- predict(svm_model, train$test_in)
  return(list(predict.train=predict.train,
              predict.test=predict.test))
}

topsvm <- function(dgms,cls=cls,ids=img_ids,
                   k=0.25,bins=20,h=0.1,runs=2,
                   homology=c(0,1), print.table=F){
  # print("Converting Diagrams to Persistence Images")
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
    print(paste('Runs = ',j))
    train <- split_data_alt(p_img,ids,k)
    results <- img_svm(train)
    error[j] <- sum(train$test_out != results$predict.test)/ length(train$test_out)
  }
  
  return(c(mean(error),sd(error)))
  
}

lines.plot <- function(domain,l1,l2,l3,
                       ylab="",main="",type='sd',q=1,pal=NA){
  if(is.na(pal)){
    pal <- c('firebrick1','dodgerblue','grey40')
  }
  par(oma=c(0,0,0,0),mar=c(4,4,1,1))
  limsy <- max(c(l1[1,],l2[1,],l3[1,]))
  plot(sq,l1[1,],type="l",ylim=c(0,0.65),col=pal[1],main = main,lwd=3,
       ylab='Misclassification Error',xlab="Persistence Image Bandwidth")
  lines(sq,l2[1,],type="l",ylim=c(0,1),col=pal[2],lwd=3)
  lines(sq,l3[1,],type="l",ylim=c(0,1),col=pal[3],lwd=3)

  legend('topright',
         legend = c('KDE','RKDE','DTM'),
         col=pal,
         lty= c(1,1,1),lwd=c(3,3,3) )
}


lines.save <- function(filename,domain,l1,l2,l3,
                       ylab="",main="",type='sd',q=1,pal=NA){
  if(is.na(pal)){
    pal <- c('firebrick1','dodgerblue','grey40')
  }
  limsy <- max(c(l1[1,],l2[1,],l3[1,]))
  pdf(file=filename,height=4,width=4)
  par(oma=c(0,0,0,0),mar=c(4,4,1,1))
  plot(sq,l1[1,],type="l",ylim=c(0,0.65),col=pal[1],main = main,lwd=3,
       ylab='Misclassification Error',xlab="Persistence Image Bandwidth")
  lines(sq,l2[1,],type="l",ylim=c(0,1),col=pal[2],lwd=3)
  lines(sq,l3[1,],type="l",ylim=c(0,1),col=pal[3],lwd=3)
  
  legend('topright',
         legend = c('KDE','RKDE','DTM'),
         col=pal,
         lty= c(1,1,1),lwd=c(3,3,3) )
  dev.off()
}
