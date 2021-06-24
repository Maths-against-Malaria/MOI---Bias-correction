#' Generates multinomially distributed vectors
#'
#' @description function return N multinomial vectors, where the k-th vector is
#'   distributed as Mult(m_k;p_1,...,p_n) where m=(m_1,...,m_N),
#'   p=(p_1,...,p_n)...frequency distribution
#'
#' @param m vector; MOI (#super-infections) of N samples
#' @param p vector; lineage frequency distribution
#'
#' @return matrix; dataset of size N x n where N is the sample size and n is the
#'   number of lineages. Each row corresponds to a sample which is the vector of
#'   the super-infections corresponding to that sample.
#'
#' @export
#' \donotrun{
#' m <- cpoiss(2, 150)
#' p <- c(0.6,0.4) 
#' mnom(m, p)
#' }
mnom <- function(m, p) { 
  N <-length(m)
  out<-array(0,dim=c(N,length(p)))
  for(k in 1:N){
    out[k,]=rmultinom(1,m[k],p)
  }
  out
}
