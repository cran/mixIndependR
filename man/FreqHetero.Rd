% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FreqHetero.R
\name{FreqHetero}
\alias{FreqHetero}
\title{Build Observed Distribution of No. of Heterozygous loci}
\usage{
FreqHetero(h)
}
\arguments{
\item{h}{a dataframe of heterozygosity, made up with 0 and 1, outcome of function "Heterozygous" Rows for individuals, and columns for markers.}
}
\value{
a dataframe of frequencies of each number of heterozygous loci(from 0 to No. of loci)
}
\description{
Build Observed Distribution of No. of Heterozygous loci
}
\details{
This function build the observed distributions from observed heterozygosity table, made up of 0,1.
}
\examples{
h<-matrix(rbinom(20,1,0.5),nrow=5)
FreqHetero(h)
}
