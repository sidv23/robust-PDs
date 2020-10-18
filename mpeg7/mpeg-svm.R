# Analysis
setwd('/storage/work/s/suv87/tda/NeurIPS/NeurIPS/')
source('RKDE.R')
source('ph-functions.R')
source('mpeg7/mpeg7-svm-functions.R')


reload.data <- function(){
  load('mpeg7/data/cls.RData',     envir = globalenv())
  load('mpeg7/data/img_ids.RData', envir=globalenv())
  load('mpeg7/data/kde_dgms.RData', envir=globalenv())
  load('mpeg7/data/rkde_dgms.RData',envir=globalenv())
  load('mpeg7/data/dtm_dgms.RData', envir=globalenv())
}

#Load the Data
classes <- c("beetle","spring","Bone","horse","deer")


reload.data()
H <- 0.015; B <- 20; K <- 0.25

for(numcls in c(3,5)){
  
  test_classes <- classes[1:numcls]
  reload.data()
  
  img_ids <- img_ids[cls %in% test_classes]
  kde_dgms <- kde_dgms[cls %in% test_classes]
  rkde_dgms <- rkde_dgms[cls %in% test_classes]
  dtm_dgms <- dtm_dgms[cls %in% test_classes]
  cls <- cls[cls %in% test_classes]
  
  mpeg7.svm <- Curry(topsvm,bins=B,k=K,
                     ids=img_ids, cls=cls)
  
  set.seed(2020)
  sq <- seq(0.005,0.2,0.02)
  l.kde  <- sapply(sq,function(x)mpeg7.svm(kde_dgms, h=x,runs=200))
  l.rkde <- sapply(sq,function(x)mpeg7.svm(rkde_dgms,h=x,runs=200))
  l.dtm  <- sapply(sq,function(x)mpeg7.svm(dtm_dgms, h=x,runs=200))
  
  lines.plot(sq,l.kde,l.rkde,l.dtm)
  filename <- paste('mpeg7/plots/',numcls,'-class-lines.pdf',sep="")
  lines.save(filename,sq,l.kde,l.rkde,l.dtm)
  filename2 <- paste('mpeg7/plots/',numcls,'-class-lines.RData',sep="")
  save(l.kde,l.rkde,l.dtm,file=filename2)
}



