This page describes the method to construct robust persistence diagrams,
as implemented in our paper [Robust Persistence Diagrams using
Reproducing Kernels](https://arxiv.org/abs/2006.10012).

You will need the following dependencies for implementing the analyses
using R. Please run the following code:

``` r
pkgs <- c("dplyr","plotrix","spatstat","TDA","hitandrun","functional","Rfast","plotly","viridis","plot3D")
sapply(pkgs, install.packages, character.only=TRUE)
```

Here’s an example for computing the persistence diagrams
--------------------------------------------------------

We start by sampling points 𝕏<sub>*n*</sub> from a circle in 2D with
some uniform noise in the enclosing region.

``` r
set.seed(2020)
signal <- circleUnif(400,4)*rnorm(400,1,0.1)
noise <- matrix(runif(2*1000,-8,8),ncol=2)
X <- rbind(signal,noise)
plot(X,asp=1,pch=20,cex=0.4,col=alpha("black",0.7))
draw.circle(0,0,4,border = alpha('red',0.2),lwd=3)
```

![](README_files/figure-markdown_github/unnamed-chunk-1-1.png)

Here are the persistence diagrams.

### 1. Distance Function Persistence Diagram

``` r
dgm.distFct <- ph.distfun(X,by=0.1,H=0)$diagram
plot.diagram(dgm.distFct)
```

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)

### 2. DTM Persistence Diagram

``` r
dgm.dtm <- ph.dtm(X,by=0.1,m0 = M0(X,3))$diagram
plot.diagram(dgm.dtm)
```

![](README_files/figure-markdown_github/unnamed-chunk-3-1.png)

### 3. KDE Persistence Diagram

``` r
dgm.kde <- ph.kde(X,by=0.1,H=bw(X,3))$diagram
plot.diagram(dgm.kde)
```

![](README_files/figure-markdown_github/unnamed-chunk-5-1.png)

### 4. Robust Persistence Diagram

``` r
dgm.rkde <- ph.rkde2(X,by=0.1,H=bw(X,3))$diagram
plot.diagram(dgm.rkde)
```

![](README_files/figure-markdown_github/unnamed-chunk-7-1.png)
