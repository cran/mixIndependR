## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load package-------------------------------------------------------------
library(mixIndependR)

## ----preparation, include=FALSE-----------------------------------------------
x <- data.frame(STR1=c(12,13,13,14,15,13,14,12,14,15),
                STR1_1=c(12,14,13,15,13,14,13,12,14,15),
                SNP1=c("A","T","A","A","T","A","A","T","T","A"),
                SNP1_1=c("A","T","T","T","A","T","A","A","T","T"))
p <- AlleleFreq(x)
h <-Heterozygous(x)
H <- RxpHetero(h,p,HWE=FALSE)
AS<-AlleleShare_Table(x,replicate=TRUE)
e <-RealProAlleleShare(AS)
ObsDist_K<-FreqHetero(h)
ExpDist_K<- DistHetero(H)
ObsDist_X<-FreqAlleleShare(AS)
ExpDist_X<-DistAlleleShare(e)

## ----Simulation---------------------------------------------------------------
Simu_K <- Simulate_DistK(H,50,500)
Simu_X <- Simulate_DistX(e,25,500)

## ----Chi-square---------------------------------------------------------------
x2_K<-Dist_SimuChisq(Simu_K,ExpDist_K$Density,100)
x2_X<-Dist_SimuChisq(Simu_X,ExpDist_X$Density,100)
P1<-ecdf(x2_K)
P2<-ecdf(x2_X)

## ----Last plot, echo=FALSE,fig.show='hold'------------------------------------
plot(P1)
abline(h=0.95)
text(1,0.8,"0.95")
plot(P2)
abline(h=0.95)
text(1,0.8,"0.95")
