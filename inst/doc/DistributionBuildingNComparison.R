## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(mixIndependR)
library(ggplot2)

## ----include=FALSE------------------------------------------------------------
x <- data.frame(STR1=c("12|12","13|14","13|13","14|15","15|13","13|14","14|13","12|12","14|14","15|15"),SNP1=c("A|A","T|T","A|T","A|T","T|A","A|T","A|A","T|A","T|T","A|T"))

p <- AlleleFreq(x)
h <-Heterozygous(x)
H <- RxpHetero(h,p,HWE=FALSE)
AS<-AlleleShare(x,replacement =FALSE)
e <-RealProAlleleShare(AS)

## -----------------------------------------------------------------------------
ObsDist_K<-FreqHetero(h)
ExpDist_K<- DistHetero(H)
head(ObsDist_K)
head(ExpDist_K)

## -----------------------------------------------------------------------------
ObsDist_X<-FreqAlleleShare(AS)
ExpDist_X<-DistAlleleShare(e)
head(ObsDist_X)
head(ExpDist_X)

## -----------------------------------------------------------------------------
df_K <- ComposPare_K(h,ExpDist_K,trans = F)
df_X <- ComposPare_X(AS,ExpDist_X,trans = F)

## ----fig.show='hold',echo=FALSE-----------------------------------------------

ggplot(df_K, aes(x=freq, color=OvE, fill=OvE)) +
            geom_histogram(aes(y=..density..), alpha=0.5,binwidth = 1,
                 position="identity")+
            geom_density(alpha=.2,bw=T) +
            ggtitle("Number of Heterozygous Loci")
ggplot(df_X, aes(x=freq, color=OvE, fill=OvE)) +
            geom_histogram(aes(y=..density..), alpha=0.5,binwidth = 1,
                 position="identity")+
            geom_density(alpha=.2,bw=T) +
            ggtitle("Number of Shared Alleles")

