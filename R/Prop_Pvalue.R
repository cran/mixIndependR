#'Loop and Find distribution of P values
#'@details This function can generate a bundle of p-values for one sample or for fix-sized sub-samples. The bundle of cumulative probabilities and the proportion of p-values(1- cumulative probability, please refer to the description of Dist_SimuChisq.R) smaller than Alpha are exported.   
#'@usage Prop_Pvalue(x,N,B,t,m=NA,Part=FALSE,alpha=0.05)
#'@importFrom stats ecdf
#'@param x a dataset of alleles. Each row denotes each individual.One allele in one cell.In the (2r-1)th column, there is the same locus with the 2r-th column; noted: no column for ID, make row.names=1 when importing.
#'@param N times of running this loop, also the number of p-values in the bundle.
#'@param t times of simulation in "Simulate_DistK" and "Simulate_DistX".
#'@param B times of bootstrapping in Chi Squares Test.
#'@param m when Part is TRUE, sub-samples are chosen from x. m is the number of loci in the subsample
#'@param Part a logical variable. If TRUE, this function will calculate p-values for sub-samples with a given sample size m.
#'@param alpha 1- confidence level; if the confidence level is 95\%, alpha =0.05 
#'@return a list of bundle of cumulative probabilities for number of heterozygous loci and bundle for number of shared alleles; and the proportions of p values smaller than alpha.
#'@export
#'@examples 
#'x0 <- data.frame(STR1=c(12,13,13,14,15,13,14,12,14,15),
#'                 STR1_1=c(12,14,13,15,13,14,13,12,14,15),
#'                 SNP1=c("A","T","A","A","T","A","A","T","T","A"),
#'                 SNP1_1=c("A","T","T","T","A","T","A","A","T","T"),
#'                 STR2=c(10,12,11,9,10,12,11,12,12,10),
#'                 STR2_1=c(10,9,11,11,10,12,10,10,12,9),
#'                 SNP2=c("C","C","G","G","G","G","C","G","G","C"),
#'                 SNP2_1=c("C","C","G","G","C","G","C","C","G","G"))
#'Prop_Pvalue(x0,3,10,10,m=2,Part = TRUE,alpha = 0.05)

Prop_Pvalue<-function(x,N,B,t,m=NA,Part=FALSE,alpha=0.05){
  n <- nrow(x)
  l <- ncol(x)/2
  if (Part==T){
    r <- sample(1:l,m,replace = F)
    index <- rep.int(0,2*m)
    for (j in 1:m){
      a <- 2*j-1
      b <- 2*j
      index[a] <- (2*r-1)[j]
      index[b] <- (2*r)[j]
    }
    x <- x[,index]
  }
  else{
    x <- x
  }
  p <- AlleleFreq(x)
  p_value <- rep.int(0,N)
  p_value_AS <- rep.int(0,N)
for (i in 1:N){  
    h <- Heterozygous(x)
    H <- RxpHetero(h,p,HWE = F)
    Obs_DistHetero<-FreqHetero(h)
    Exp_DistHetero<-DistHetero(H)
    prob_K<-Exp_DistHetero$Density
    obs_K<-Obs_DistHetero$Frequency
    idx_K <-which(prob_K==0)
    if (length(idx_K)==0){
      prob_K <- prob_K
      obs_K <- obs_K
    }
    else{
      prob_K <- prob_K[-idx_K]
      obs_K <- obs_K[-idx_K]
    }
    x20_K <-chisq.test(obs_K,p=prob_K,simulate.p.value = T,B=B)
    s_K<-Simulate_DistK(H,n,t)
    x2_K<-Dist_SimuChisq(s_K,Exp_DistHetero$Density,B)
    P_K <- ecdf(x2_K)
    p_value[i]<-P_K(x20_K$statistic)
    print(paste(i,"K",p_value[i]))
    AS <- AlleleShare_Table(x,replicate = F)
    Obs_DistAlleleShare<-FreqAlleleShare(AS)
    e <- RealProAlleleShare(AS)
    Exp_DistAlleleShare <- DistAlleleShare(e)
    prob_X<-Exp_DistAlleleShare$Density
    obs_X<-Obs_DistAlleleShare$Frequency
    idx_X <-which(prob_X==0)
    if (length(idx_X)==0){
      prob_X <- prob_X
      obs_X <- obs_X
    }
    else{
      prob_X <- prob_X[-idx_X]
      obs_X <- obs_X[-idx_X]
    }
    x20_X <-chisq.test(obs_X,p=prob_X,simulate.p.value = T,B=B)
    s_X<-Simulate_DistX(e,n/2,t)
    x2_X<-Dist_SimuChisq(s_X,Exp_DistAlleleShare$Density,B)
    P_X <- ecdf(x2_X)
    p_value_AS[i]<-P_X(x20_X$statistic)
    print(paste(i,"X",p_value_AS[i]))
  }
  pdata_K <- p_value
  pdata_X <- p_value_AS
  Output_KX <- list(Pr_K=length(which(1-p_value<alpha))/N,
                    Pr_X=length(which(1-p_value_AS<alpha))/N,
                    pdata_K=p_value,
                    pdata_X=p_value_AS) 
  print(paste("K",Output_KX$Pr_K,"X",Output_KX$Pr_X))
  return(Output_KX)
}
