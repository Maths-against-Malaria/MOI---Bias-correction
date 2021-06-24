#' Derives the BCMLE of MOI parameter and frequency spectra
#'
#' @description derives the bias-corrected maximum-likelihood estimate (BCMLE)
#'   of the MOI parameter (Poisson parameter) and the lineage (allele)
#'   frequencies.
#'
#' @param N integer; sample size
#' @param Nk integer vector; number of lineage prevalence counts in a dataset.
#'   for a simulated data this is simply derived as \code{colSums(dataset)}. To
#'   derive the MLE and lineage prevalence counts for a real dataset please
#'   refer to the package \link[MLMOI]{moimle}.
#'
#' @return list; 
#' 1. bcmle_lam...the BCMLE of MOI parameter lambda 
#' 2. bcmle_psi...the BCMLE of mean MOI psi
#' 3. bcmle_p...the BCMLE of lineage frequencies
#'
#' @export
#'
#' @examples
#' \donotrun{
#' m <- cpoiss(2, 150) #lambda = 2, N = 150
#' p <- c(0.6,0.4) #lineage frequencies
#' dataset <- mnom(m, p)
#' Nk <- colSums(dataset)
#' BCMLE(150, Nk)
#' }
BCMLE <- function(N, Nk){
    mle <- MLE(N,Nk)
    mle_lam <- mle[[2]]
    mle_p <- mle[[4]]
    
    bias <- second_order_bias(N, mle_lam, mle_p)
    bias_lam <- bias[[1]]
    bias_p <- bias[[2]]
    
    bcmle_lam <- mle_lam - bias_lam                 #bias-corrected MLE of lambda
    bcmle_psi <- bcmle_lam/(1 - exp(-bcmle_lam))    #bias-corrected MLE of psi
    
    bcmle_p <- mle_p - bias_p                       #bias-corrected MLE of lambda lineage frequencies
    
    out <- list(bcmle_lam, bcmle_psi, bcmle_p)
    out
}
