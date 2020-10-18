# Influence Analysis
setwd('/storage/work/s/suv87/tda/NeurIPS/NeurIPS')
source('./ph-functions.R')
source('./Influence/influence-functions.R')

packages <- c('snow','doParallel','foreach')
sapply(packages,require,character.only=T)
plot <- graphics::plot

# Parallel Settings
# nprocs <-10
# mp_type = "PSOCK"
# cl <- parallel::makeCluster(nprocs, type=mp_type)
# doParallel::registerDoParallel(cl)


# Main Stuff 
set.seed(2020)
signal <- sample.annulus(400,4,5)
signal <- rbind(signal,matrix(runif(300,-6,6),ncol=2))

png(filename = 'plots/points.png')
plot(signal,asp=1,pch='.')
dev.off()

k <- 5
M <- M0(signal,k)
H <- bw(signal,k)
by <- 0.5

D1 <- ph.rkde2(signal,by=0.5,H)
D2 <- ph.kde(signal,by=0.5,H)
D3 <- ph.dtm(signal,by=0.5,M)

G <- make.grid(rbind(c(-110,-110),c(110,110)),by=0.5)

F1 <- kde(signal,G,H,
          weight=(1e-10+rkde.w(signal,H)[[1]]))
F2 <- kde(signal,G,H)
F3 <- dtm(signal,G,M)

z <- list()
z$d1 <- z$d2 <- z$d3 <- list()
z$b1 <- z$b2 <- z$b3 <- list()
z$w1 <- z$w2 <- z$w3 <- list()
z$B1 <- z$B2 <- z$B3 <- list()
z$W1 <- z$W2 <- z$W3 <- list()
z$l1 <- z$l2 <- z$l3 <- list()

R <- seq(10,100,10)

for(i in 1:length(R)){
  # r <- 100
  temp <- gen.sim(signal,R[i],F1,F2,F3,D1,D2,D3)
  z$l1[[i]] <- temp$l1
  z$l2[[i]] <- temp$l2
  z$l3[[i]] <- temp$l3
  z$w1[[i]] <- temp$w1
  z$w2[[i]] <- temp$w2
  z$w3[[i]] <- temp$w3
  z$b1[[i]] <- temp$b1
  z$b2[[i]] <- temp$b2
  z$b3[[i]] <- temp$b3
  z$W1[[i]] <- temp$W1
  z$W2[[i]] <- temp$W2
  z$W3[[i]] <- temp$W3
  z$B1[[i]] <- temp$B1
  z$B2[[i]] <- temp$B2
  z$B3[[i]] <- temp$B3
}
save(z,file='z.RData')

temp1 <- temp2 <- list()
for(i in 1:length(R)){
  temp1[[i]] <- z$b1[[i]] + z$B1[[i]]
  temp2[[i]] <- z$b2[[i]] + z$B2[[i]]
}

filename = 'plots/bottleneck.pdf'
influence.plot(filename,R,temp1,temp2,ylab='Influence',main='Bottleneck Influence')
# dev.off()

filename = 'plots/supnorm.pdf'
influence.plot(filename,R,z$l1,z$l2,ylab='Influence',main='Sup norm influence')
# dev.off()

filename = 'plots/H0-bottleneck.pdf'
influence.plot(filename,R,z$b1,z$b2,ylab='Influence',main='H0: Bottleneck Influence')
# dev.off()

filename = 'plots/H1-bottleneck.pdf'
influence.plot(filename,R,z$B1,z$B2,ylab='Influence',main='H1: Bottleneck Influence')
# dev.off()

filename = 'plots/H0-W1.pdf'
influence.plot(filename,R,z$w1,z$w2,ylab='Influence',main='H0: W1 Influence')
# dev.off()

filename = 'plots/H1-W1.pdf'
influence.plot(filename,R,z$W1,z$W2,ylab='Influence',main='H1: W1 Influence')








