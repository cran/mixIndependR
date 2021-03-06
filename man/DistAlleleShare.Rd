% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DistAlleleShare.R
\name{DistAlleleShare}
\alias{DistAlleleShare}
\title{Build Expected Distribution of Numbers of Shared Alleles}
\usage{
DistAlleleShare(e)
}
\arguments{
\item{e}{a matrix/dataframe of probability of shared alleles; outcome of "ExpProAlleleShare" or "RealProAlleleShare". Each row denotes each locus. The first column is the case of 0 shared alleles, the second column is the case of 1 shared alleles, the third column is the case of 2 shared alleles.}
}
\value{
a dataframe of probabilities of each number of shared alleles(from 0 to 2*loci); the first column is No. of Shared Alleles; the Second Column is Expected Density
}
\description{
Build Expected Distribution of Numbers of Shared Alleles
}
\details{
This function build the expected distribution of numbers of shared alleles for known shared alleles of each pair of individuals.
}
\examples{
e0<-data.frame("P0"=runif(5,min = 0,max = 0.5),"P1"=runif(5,0,0.5))
e<-data.frame(e0,"P2"=1-rowSums(e0))
DistAlleleShare(e)
}
\references{
Chakraborty, R., Stivers, D. N., Su, B., Zhong, Y., & Budowle, B. (1999) <doi:10.1002/(SICI)1522-2683(19990101)20:8<1682::AID-ELPS1682>3.0.CO;2-Z>
}
