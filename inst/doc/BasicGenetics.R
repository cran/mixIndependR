## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo = FALSE------------------------------------------------------------
x <- data.frame(STR1=c("12|12","13|14","13|13","14|15","15|13","13|14","14|13","12|12","14|14","15|15"),SNP1=c("A|A","T|T","A|T","A|T","T|A","A|T","A|A","T|A","T|T","A|T"))

print(x)

## ----setup--------------------------------------------------------------------
library(mixIndependR)

## -----------------------------------------------------------------------------
AlleleFreq(x,sep = "\\|")

## -----------------------------------------------------------------------------
GenotypeFreq(x,sep = "\\|",expect = FALSE)  ####or GenotypeFreq(x)
GenotypeFreq(x,sep = "\\|",expect = TRUE) ####or GenotypeFreq(x,expect =T)

## -----------------------------------------------------------------------------
h <-Heterozygous(x,sep = "\\|") ####or Just use Heterozygous(x)
print(h)

## -----------------------------------------------------------------------------
p<-AlleleFreq(x,sep = "\\|")
H <- RxpHetero(h,p,HWE=TRUE)
head(H)

## -----------------------------------------------------------------------------
AS<-AlleleShare(x,sep = "\\|",replacement = FALSE) ###or without "sep="
head(AS)

## -----------------------------------------------------------------------------
e <-RealProAlleleShare(AS)
e0<-ExpProAlleleShare(p)
head(e)
head(e0)

## -----------------------------------------------------------------------------
g <- GenotypeFreq(x,expect=FALSE)
g0 <- GenotypeFreq(x,expect=TRUE)
HWE.Chisq(g,g0,rescale.p = T,simulate.p.value = T,B=2000)

