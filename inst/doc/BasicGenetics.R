## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)



## ----echo=FALSE---------------------------------------------------------------
x <- data.frame(STR1=c(12,13,13,14,15,13,14,12,14,15),
                STR1_1=c(12,14,13,15,13,14,13,12,14,15),
                SNP1=c("A","T","A","A","T","A","A","T","T","A"),
                SNP1_1=c("A","T","T","T","A","T","A","A","T","T"))
print(x)

## ----setup--------------------------------------------------------------------
library(mixIndependR)

## -----------------------------------------------------------------------------
AlleleFreq(x)

## -----------------------------------------------------------------------------
p <- AlleleFreq(x)
GenotypeFreq(x,p,expect = FALSE)
GenotypeFreq(x,p,expect = TRUE)

## -----------------------------------------------------------------------------
h <-Heterozygous(x)
print(h)

## -----------------------------------------------------------------------------
H <- RxpHetero(h,p,HWE=TRUE)
head(H)

## -----------------------------------------------------------------------------
AS<-AlleleShare_Table(x,replicate=TRUE)
head(AS)

## -----------------------------------------------------------------------------
e <-RealProAlleleShare(AS)
e0<-ExpProAlleleShare(p)
head(e)
head(e0)

## -----------------------------------------------------------------------------
g <- GenotypeFreq(x,p,expect=FALSE)
g0 <- GenotypeFreq(x,p,expect=TRUE)
HWE.Chisq(g,g0,rescale.p=FALSE,simulate.p.value=TRUE,2000)
HWE.Fisher(p,H,g/colSums(g))

