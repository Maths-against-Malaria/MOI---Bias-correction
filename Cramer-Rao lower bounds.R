#' Cramer-Rao lower bounds (CRLB) of the model parameters
#'
#' @param N integer; sample size
#' @param lambda numeric; MOI parameter
#' @param p numeric vector; lineage frequencies
#'
#' @return Cramer-Rao lower bounds (CRLB) of the model parameters for the
#'   specified model parameters and sample size
#' @export
#'
#' @examples
crlb <- function(N, lambda, p){
    p <- sort(p,decreasing = T)
    dk <- exp(lambda*p)-1
    d <- sum(dk)
    d0 <- 1/(exp(lambda)-1)
    x <- (1 - d*d0)*N 
    crlam11 <- d/(x*(d0 + 1))
    cr11 <- (d0 + 1)*((1 - lambda*d0)^2)*d/x
    crkk <- (dk/N + (d*(p^2) - 2*p*dk + d0*(dk^2))/x)/((d0 + 1)*(lambda^2))
    c(cr11,crkk)
}
