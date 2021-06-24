#' Derives the HBCMLE3 of MOI parameter and frequency spectra
#'
#' @description derives the 3rd version of heuristically bias-corrected
#'   maximum-likelihood estimate (HBCMLE3) of the MOI parameter (Poisson
#'   parameter) and the lineage (allele) frequencies.
#'
#' @param N integer; sample size
#' @param Nk integer vector; number of lineage prevalence counts in a dataset.
#'   for a simulated data this is simply derived as \code{colSums(dataset)}. To
#'   derive the MLE and lineage prevalence counts for a real dataset please
#'   refer to the package \link[MLMOI]{moimle}.
#'
#' @return list; 
#' 1. hbcmle_3_lam...the HBCMLE3 of MOI parameter lambda 
#' 2. hbcmle_3_psi...the HBCMLE3 of mean MOI psi 
#' 3. hbcmle_3_p...the HBCMLE3 of lineage frequencies
#'
#' @export
#'
#' @examples
#' \donotrun{
#' m <- cpoiss(2, 150) #lambda = 2, N = 150
#' p <- c(0.6,0.4) #lineage frequencies
#' dataset <- mnom(m, p)
#' Nk <- colSums(dataset)
#' HBCMLE3(150, Nk)
#' }
HBCMLE3 <- function(N, Nk){
    mle <- MLE(N, Nk)
    mle_lam <- mle[[2]]
    mle_p <- mle[[4]]
    
    bcmle <- BCMLE(N, Nk)
    bcmle_lam <- bcmle[[1]]
    bcmle_p <- bcmle[[3]]
    
    bias <- second_order_bias(N, mle_lam, mle_p)                   #second-order bias evaluated at the MLE
    bias_lam <- bias[[1]]
    bias_p <- bias[[2]]
    
    p_pathologic <- prob_pathological(N, bcmle_lam, bcmle_p)    #probability of pathological data evaluated at the BCMLE
    p_regular <- 1 - p_pathologic                               #probability of regular data evaluated at the BCMLE
    
    hbcmle_3_lam <- p_regular*mle_lam - bias_lam                #HBCMLE2 of lambda
    hbcmle_3_psi <- hbcmle_3_lam/(1 - exp(-hbcmle_3_lam))       #HBCMLE2 of psi
    hbcmle_3_p <- bcmle_p                                       #HBCMLE2 of lineage frequencies
    
    out <- list(hbcmle_3_lam, hbcmle_3_psi, hbcmle_3_p)
    out
}
