## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load package-------------------------------------------------------------
library(mixIndependR)

## ----preparation, include=FALSE-----------------------------------------------
x <- data.frame(STR1=c("12|12","13|14","13|13","14|15","15|13","13|14","14|13","12|12","14|14","15|15"),SNP1=c("A|A","T|T","A|T","A|T","T|A","A|T","A|A","T|A","T|T","A|T"))

p <- AlleleFreq(x)
h <-Heterozygous(x)
H <- RxpHetero(h,p,HWE=FALSE)
AS<-AlleleShare(x,replacement=FALSE)
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

