#'Simple count including zero###
#'@details This function counts how many the assigned elements there are in one vector.
#'@param z a vector you would like to check
#'@param y  an element you would like to count.(Even it is not included in z)
#'@return the times that y appears in z
#'@export
#'@examples
#'z <-rbinom(20,1,0.5)
#'counta(z,0)
#'


counta <- function(z,y){
  f <- 0
  n <- length(z)
  for (i in 1:n){
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
