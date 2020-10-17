Introduction
------------

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see
<a href="http://rmarkdown.rstudio.com" class="uri">http://rmarkdown.rstudio.com</a>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

Including Plots
---------------

You can also embed plots, for example:

![](README_files/figure-markdown_github/pressure-1.png)

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.

Here’s some TDA Stuff to test
=============================

``` r
set.seed(2020)
signal <- circleUnif(100,r=1)*rnorm(100,1,0.1)
noise  <- matrix(runif(50*2,-3,3),ncol=2)
X <- rbind(signal,noise)
plot(X,asp=1,pch=20)
```

![](README_files/figure-markdown_github/unnamed-chunk-8-1.png)

Here are the persistence diagrams.

1.  Distance Function Persistence Diagram

``` r
dgm.distFct <- ph.distfun(X,by=0.1,H=0)$diagram
plot(dgm.distFct)
```

![](README_files/figure-markdown_github/unnamed-chunk-9-1.png)

1.  DTM Persistence Diagram

``` r
dgm.dtm <- ph.dtm(X,by=0.1,m0 = M0(X,3))$diagram
plot(dgm.dtm)
```

![](README_files/figure-markdown_github/unnamed-chunk-10-1.png)

1.  KDE Persistence Diagram

``` r
dgm.kde <- ph.kde(X,by=0.1,H=bw(X,3))$diagram
plot(dgm.kde)
```

![](README_files/figure-markdown_github/unnamed-chunk-11-1.png)

1.  Robust Persistence Diagram

``` r
dgm.rkde <- ph.rkde2(X,by=0.1,H=bw(X,3))$diagram
plot(dgm.rkde)
```

![](README_files/figure-markdown_github/unnamed-chunk-12-1.png)
