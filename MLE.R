#' Derives the MLE of MOI parameter and frequency spectra
#'
#' @description derives the maximum-likelihood estimate (MLE) of the MOI
#'   parameter (Poisson parameter) and the lineage (allele) frequencies.
#'
#' @param N integer; sample size
#' @param Nk integer vector; number of lineage prevalence counts in a dataset.
#'   for a simulated data this is simply derived as \code{colSums(dataset)}. To derive
#'   the MLE and lineage prevalence counts for a real dataset please refer to
#'   the package \link[MLMOI]{moimle}.
#'
#' @return list; 
#' 1. ml...maximum log-likelihood, 
#' 2. mle_lam...the MLE of MOI parameter lambda 
#' 3. mle_psi...the MLE of mean MOI psi 
#' 4. mle_p...the MLE of lineage frequencies
#'
#' @export
#'
#' @examples
#' \donotrun{
#' m <- cpoiss(2, 150) #lambda = 2, N = 150
#' p <- c(0.6,0.4) #lineage frequencies
#' dataset <- mnom(m, p)
#' Nk <- colSums(dataset)
#' MLE(150, Nk)
#' }
#' 
MLE <- function(N, Nk){
  sel <- Nk
  Nk <- sel[sel>0]
  nk <- Nk/N
  l1 <- 2.5         # initial value
  la <- 2.5
  l0 <- 0
  eps <- 10^(-8)       # precision 
  out <- list(NA, NA,NA,NA,NA)
  k <- 1
  while(abs(l0-l1)>eps && k<50 && l1>0){
    k <- k+1
    l0 <- l1
    l1 <- l0-(l0+sum(log(1-nk*(1-exp(-l0)))))/(1-sum(nk/(exp(l0)*(1-nk)+nk)))
  }
  if(k==50 || l1<0){
    print(c(l0,l1,Nk))
    for(st in 1:10){
      print(st)
      l1 <- st
      l0 <- l1+1
      k <- 1
      while(abs(l0-l1)>eps && k<100 && l1>0){
        k <- k+1
        l0 <- l1
        l1 <- l0-(l0+sum(log(1-nk*(1-exp(-l0)))))/(1-sum(nk/(exp(l0)*(1-nk)+nk)))
      }
      if(abs(l0-l1)<eps){
        break
      }
      
    }
    if(abs(l0-l1)>eps){               # if numerical problems occur, calculations are performed with higher precision
      l1 <- mpfr(10*la,precBits=100)
      l0 <- l1+1
      while(abs(l0-l1)>eps){
        l0 <- l1
        l1=l0-(l0+sum(log(1-nk*(1-exp(-l0)))))/(1-sum(nk/(exp(l0)*(1-nk)+nk)))
        #print(l1)
      }
    }        
  }
  mle_lam <- l1                                         #MLE of lambda
  mle_psi <- l1/(1-exp(-l1))                            #MLE of psi
  
  pk <- -1/l1*log(1-nk*(1-exp(-l1)))   
  ml <- (-N)*log(exp(l1)-1)+sum(Nk*log(exp(l1*pk)-1))	  #maximum log-likelihood 
  mle_p <- array(0,length(sel))  
  mle_p[sel>0] <- pk                                    #MLE of lineage frequencies
  out <- list(ml, mle_lam, mle_psi, mle_p)
  out	
}
