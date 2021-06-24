#' Generates negative binomial distributed random numbers
#'
#' @description This function generates N random numbers from Negative binomial
#'   distribution with parameters (success, p)
#'
#' @param N integer; sample size
#' @param success numeric; the number of successes  
#' @param p numeric; the probability of success
#'
#' @return vector; MOI for N samples
#' @export
#'
#' @examples
#' \donotrun{
#' cnegb(150, 2, 0.8)
#' }
cnegb <- function(N, success, p){
    m <- 100 # to accelerate computation it is assumed that m<100 is generically drawn
    out <- rep(0,N)
    x <- runif(N,min=0,max=1)
    p0 <- pnbinom(0,size = success, prob =  p)
    nc <- 1/(1 - p0)
    pvec <- (pnbinom(1:m,size = success, prob = p) - p0)*nc
    pvec <- c(pvec,1)
    for (i in 1:N){
        k <- 1
        while(x[i] > pvec[k]){
            k <- k+1
        }
        if(k==m){ # if a m>=100 is drawn this is executed
            k <- k+1
            a <- dnbinom(k, size = success, prob = p)*nc
            b <- pvec[m]+a
            while(x[i]>b){
                k <- k+1
                a <- a*(1-p)*(k+success-1)/k
                b <- b+a
            }
        }
        out[i] <- k
    }
    out
}
