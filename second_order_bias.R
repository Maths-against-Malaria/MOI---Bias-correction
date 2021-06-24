#' Derives the approximated second-order bias of the MLE
#' 
#' @description this function in fact derives the second-order bias - O(N^2)
#' in the paper this was specified by the function B evaluated at the 
#' estimates (equation 18). 
#'
#' @param N integer; sample size 
#' @param lambda numeric; MOI parameter
#' @param p numeric vector; lineage frequencies
#'
#' @return list
#' 1. second-order bias of the MLE of MOI parameter lambda
#' 2. second-order bias of the MLE of lineage frequencies
#' 
#' @export
#'
#' @examples
#' \donotrun{
#' second_order_bias (150, 1.23, c(0.64, 0.36))
#' }
#' 
second_order_bias <- function(N, lambda, p){
    lep <- lambda*p
    dk <- exp(lep)-1
    d <- sum(dk)
    d0 <- 1/(exp(lambda)-1)
    x <- (1- d*d0) 
    y <- N*(d0 + 1)
    den <- y*x
    nom<- (d0 + 1/2)*d - d0*((d^2) - sum(dk^2))/(2*x)

    bias_lam <- nom/den                     #second-order bias of the lambda estimate
    
    nomp <- (dk - p*d)*(d0 + 0.5 - (1/lambda)) + d0*(dk^2)/2 + d0*(d*(p*d - dk) + (d0*dk - p)*(sum(dk^2)))/(2*x)
    
    bias_p <- nomp/(den*lambda)            #second-order bias of lineage frequency estimates
    
    out <- list(bias_lam, bias_p)
    out
}
