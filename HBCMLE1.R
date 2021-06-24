#' Derives the HBCMLE1 of MOI parameter and frequency spectra
#'
#' @description derives the 1st version of heuristically bias-corrected
#'   maximum-likelihood estimate (HBCMLE1) of the MOI parameter (Poisson
#'   parameter) and the lineage (allele) frequencies.
#'
#' @param N integer; sample size
#' @param Nk integer vector; number of lineage prevalence counts in a dataset.
#'   for a simulated data this is simply derived as \code{colSums(dataset)}. To
#'   derive the MLE and lineage prevalence counts for a real dataset please
#'   refer to the package \link[MLMOI]{moimle}.
#'
#' @return list; 
#' 1. hbcmle_1_lam...the HBCMLE1 of MOI parameter lambda 
#' 2. hbcmle_1_psi...the HBCMLE1 of mean MOI psi 
#' 3. hbcmle_1_p...the HBCMLE1 of lineage frequencies
#'
#' @export
#'
#' @examples
#' \donotrun{
#' m <- cpoiss(2, 150) #lambda = 2, N = 150
#' p <- c(0.6,0.4) #lineage frequencies
#' dataset <- mnom(m, p)
#' Nk <- colSums(dataset)
#' HBCMLE1(150, Nk)
#' }
HBCMLE1 <- function(N, Nk){
    mle <- MLE(N, Nk)
    mle_lam <- mle[[2]]
    mle_p <- mle[[4]]
    
    bcmle <- BCMLE(N, Nk)
    bcmle_lam <- bcmle[[1]]
    bcmle_p <- bcmle[[3]]
    
    p_pathologic <- prob_pathological(N, mle_lam, mle_p)    #probability of pathological data evaluated at the MLE
    p_regular <- 1 - p_pathologic                           #probability of regular data evaluated at the MLE
    
    hbcmle_1_lam <- p_regular*bcmle_lam                     #HBCMLE1 of lambda
    hbcmle_1_psi <- hbcmle_1_lam/(1 - exp(-hbcmle_1_lam))   #HBCMLE1 of psi
    hbcmle_1_p <- bcmle_p                                   #HBCMLE1 of lineage frequencies
    
    out <- list(hbcmle_1_lam, hbcmle_1_psi, hbcmle_1_p)
    out
}
