# This function provides conditionally Poisson distributed numbers
# lambda...Poisson parameter
# n...number of conditional Poisson numbers to be generated

#' Generates conditionally Poisson distributed numbers
#'
#' @description This function provides conditionally Poisson distributed numbers
#'
#' @param lambda numeric; the MOI parameter
#' @param N integer; sample size
#'
#' @return vector; MOI for N samples
#' 
#' @export
#'
#' @examples
#' \donotrun{
#' cpoiss(2, 150)
#' }
#' 
cpoiss <- function (lambda, N){
  m <- 100 # to accelerate computation it is assumed that m<100 is generically drawn
  out <- rep(0,N)
  x <- runif(N,min=0,max=1)
  p0 <- ppois(0,lambda)
  nc <- 1/(1-exp(-lambda))
  pvec <- (ppois(1:m,lambda)-p0)*nc
  pvec <- c(pvec,1) 
  for (i in 1:N){
    k <- 1
    while(x[i] > pvec[k]){
      k <- k+1
    }
    if(k==m){ # if a m>=100 is drawn this is executed
      k <- k+1
      a <- dpois(k,lambda)*nc
      b <- pvec[m]+a
      while(x[i]>b){
        k <- k+1
        a <- a*lambda/k
        b <- b+a
      }
    }
    out[i] <- k
  }
  out
}
