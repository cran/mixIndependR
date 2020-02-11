#'Test heterozygosity at each locus
#'@details This function test the heterozygosity of each individuals at each locus.Output a table and Usually followed by write.csv(as.data.frame(y),file = "~/*.csv") to export the results.
#'@usage Heterozygous(x)
#'@param x a dataset of alleles. Each row denotes each individual.One allele in one cell.In the (2r-1)th column, there is the same locus with the 2r-th column; noted: no column for ID, make row.names=1 when importing.
#'@return a dataframe of heterozygosity.0 is homozygous;1 is heterozygous. Each row denotes each individual; Each column denotes each locus.
#'@export
#'@examples
#'x <- data.frame(STR1=c(12,13,13,14,15,13,14,12,14,15),
#'                STR1_1=c(12,14,13,15,13,14,13,12,14,15),
#'                SNP1=c("A","T","A","A","T","A","A","T","T","A"),
#'                SNP1_1=c("A","T","T","T","A","T","A","A","T","T"))
#'Heterozygous(x)
#'
Heterozygous <- function(x){
  x <- as.matrix(x)
  n <- nrow(x)
  m <- ncol(x)
  y <- mat.or.vec(n,m/2)
  y
  for (i in 1:n){
    for (r in 1:(m/2)){
      if (x[i,(2*r)]==x[i,(2*r-1)])
        y[i,r] <-0
      else
        y[i,r] <-1
    }
  }
  M<-colnames(x)
  Ms <- rep.int(0,m/2)
  l <- m/2
  for (i in 1:l){
    Ms[i] <- M[2*i-1]
  }
  colnames(y)<-Ms
  rownames(y)<-rownames(x)
  return(y)
}
