#'Calculate Allele Frequency
#'@details This function calculates the allele frequencies of one dataset.
#'@usage AlleleFreq(x)
#'@param x a dataset of alleles. Type needs to be Homogeneous. Each row denotes each sample. One allele in one cell.In the (2r-1)th column, there is the other allele on the same locus from that in the 2r-th column; noted: no column for ID, make row.names=1 when importing.
#'@return a matrix of allele frequencies. Each row denotes each allele; each column denotes each marker. The order of makers follows x.
#'@export
#'@examples
#'x <- data.frame(STR1=c(12,13,13,14,15,13,14,12,14,15),
#'                STR1_1=c(12,14,13,15,13,14,13,12,14,15),
#'                SNP1=c("A","T","A","A","T","A","A","T","T","A"),
#'                SNP1_1=c("A","T","T","T","A","T","A","A","T","T"))
#'AlleleFreq(x)
#'

AlleleFreq <- function(x){
    x <- as.matrix(x)
    m <- nrow(x)    #number of individuals
    n <- ncol(x)/2   #number of loci
    counta <- function(z,y){
      f <- 0
      k <- length(z)
      for (i in 1:k){
        if (is.na(z[i])){
          f <- f
        }else{
          if (z[i]==y){
            f <- f+1
          }else{
            f <- f
          }
        }
      }
      return(f)
    }
    #####funciton to count, missing value counts as 0#####
  Allele <- as.data.frame(table(x))   #####All allleles included####
  l <- nrow(Allele)   #####Number of Alleles
  p <- mat.or.vec(l,n)
  for (j in 1:n){
    a <- 2*j -1
    b <- 2*j
    for (k in 1:l){
      p[k,j] <- counta(x[,a:b],Allele[k,1])/(2*m)
    }
  }
  row.names(p) <- Allele[,1]
  M<-colnames(x)
  Ms <- rep.int(0,n)
  for (i in 1:n){
    Ms[i] <- M[2*i-1]
  }
  colnames(p)<-Ms
  return(p)
}
