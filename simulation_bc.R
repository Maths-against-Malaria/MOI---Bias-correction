#' Simulation with nested approach
#'
#' @param S integer; number of simulation steps
#' @param ssize vector; sample sizes, e.g., c(400, 200, 100, 50)
#' @param linfreq vector; lineage-frequency distribution
#' @param path path to the file
#' @param lambda numeric; MOI parameter
#'
#' @return the simulation results are stored in a txt file specified by path
#' @export
#'
simu_CP <- function(S, ssize, linfreq, path, lambda){
    n <- length(p)
    p <- as.vector(round(linfreq, digits = 2))
    sz <- length(ssize)
    ssize <- sort(ssize, decreasing = T)
    NN <- ssize[1]
    simfinal <- rep(list(matrix(NA, S, 21 + 5*n)), sz)
    while (sim < S){
        M <- mnom(NN, lambda, p, n)
        NkN <- colSums(M)
        if (sum(NkN) <= NN || max(NkN) == NN){
            
        }else{
            sim <- sim + 1
            mleN <- MLE(NkN, NN)
            bcmleN <- BCMLE(NkN, NN)
            hbcmle1N <- HBCMLE1(NkN, NN)
            hbcmle2N <- HBCMLE2(NkN, NN)
            hbcmle3N <- HBCMLE3(NkN, NN)
            
            est_lambda_psi_N <- c(mleN[[1]], mleN[[2]],
                                bcmleN[[1]], bcmleN[[2]],
                                hbcmle1N[[1]], hbcmle1N[[2]], 
                                hbcmle2N[[1]], hbcmle2N[[2]],
                                hbcmle3N[[1]], hbcmle3N[[2]])
            
            mlepn <- mleN[[3]]
            bcmlepn <- bcmleN[[3]] 
            bcmleqpn <- hbcmle1N[[3]] 
            bcmleqqpn <- hbcmle2N[[3]] 
            bcmleqqqpn <- hbcmle3N[[3]] 
            
            eun <- sqrt(sum((p - mlepn)^2)) ##Euclidean distance between true p and mle of p
            eubcn <- sqrt(sum((p - bcmlepn)^2)) ##Euclidean distance between true p and bias corrected 
            eubcqn <- sqrt(sum((p - bcmleqpn)^2)) ##Euclidean distance between true p and bias corrected
            eubcqqn <- sqrt(sum((p - bcmleqqpn)^2)) ##Euclidean distance between true p and bias corrected
            eubcqqqn <- sqrt(sum((p - bcmleqqqpn)^2)) ##Euclidean distance between true p and bias corrected
            
            pn_nonzero <- p[mlepn>0]
            mlepn_nonzero <- mlepn[mlepn>0]
            bcmlepn_nonzero <- bcmlepn[bcmlepn>0]
            bcmleqpn_nonzero <- bcmleqpn[bcmleqpn>0]
            bcmleqqpn_nonzero <- bcmleqqpn[bcmleqqpn>0]
            bcmleqqqpn_nonzero <- bcmleqqqpn[bcmleqqqpn>0]
            
            kln <- sum(mlepn_nonzero*log(mlepn_nonzero/pn_nonzero)) ##Kullback-Leibler between true p and mle of p
            klbcn <- sum(bcmlepn_nonzero*log(bcmlepn_nonzero/pn_nonzero)) ##Kullback-Leibler between true p and bias corrected 
            klbcqn <- sum(bcmleqpn_nonzero*log(bcmleqpn_nonzero/pn_nonzero)) ##Kullback-Leibler between true p and bias corrected 
            klbcqqn <- sum(bcmleqqpn_nonzero*log(bcmleqqpn_nonzero/pn_nonzero)) ##Kullback-Leibler between true p and bias corrected 
            klbcqqqn <- sum(bcmleqqqpn_nonzero*log(bcmleqqqpn_nonzero/pn_nonzero)) ##Kullback-Leibler between true p and bias corrected 
            
            simfinal[[1]][sim,] <- c(est_lambda_psi_N,
                                     mlepn, bcmlepn, bcmleqpn, bcmleqqpn, bcmleqqqpn,
                                     eun, eubcn, eubcqn, eubcqqn, eubcqqqn,
                                     kln, klbcn, klbcqn, klbcqqn, klbcqqqn)
            k <- 1
            for (N in ssize[-1]) {
                k <- k + 1
                Nk <- colSums(M[1:N,])
                if (sum(Nk) <= N || max(Nk) == N){
                    Nk <- NA
                    newdata <- innersamplegenerator_CP(Nk,N,lambda,p,n)
                    M <- newdata[[1]]
                    Nk <- newdata[[2]]
                    
                    mle <- MLE(Nk, N)
                    bcmle <- BCMLE(Nk, N)
                    hbcmle1 <- HBCMLE1(Nk, N)
                    hbcmle2 <- HBCMLE2(Nk, N)
                    hbcmle3 <- HBCMLE3(Nk, N)
                    
                    est_lambda_psi <- c(mle[[1]], mle[[2]],
                                        bcmle[[1]], bcmle[[2]],
                                        hbcmle1[[1]], hbcmle1[[2]], 
                                        hbcmle2[[1]], hbcmle2[[2]],
                                        hbcmle3[[1]], hbcmle3[[2]])
                    
                    mlep <- mle[[3]]
                    bcmlep <- bcmle[[3]] 
                    bcmleqp <- hbcmle1[[3]] 
                    bcmleqqp <- hbcmle2[[3]] 
                    bcmleqqqp <- hbcmle3[[3]] 
                    
                    eu <- sqrt(sum((p - mlep)^2)) ##Euclidean distance between true p and mle of p
                    eubc <- sqrt(sum((p - bcmlep)^2)) ##Euclidean distance between true p and bias corrected 
                    eubcq <- sqrt(sum((p - bcmleqp)^2)) ##Euclidean distance between true p and bias corrected
                    eubcqq <- sqrt(sum((p - bcmleqqp)^2)) ##Euclidean distance between true p and bias corrected
                    eubcqqq <- sqrt(sum((p - bcmleqqqp)^2)) ##Euclidean distance between true p and bias corrected
                    
                    p_nonzero <- p[mlep>0]
                    mlep_nonzero <- mlep[mlep>0]
                    bcmlep_nonzero <- bcmlep[bcmlep>0]
                    bcmleqp_nonzero <- bcmleqp[bcmleqp>0]
                    bcmleqqp_nonzero <- bcmleqqp[bcmleqqp>0]
                    bcmleqqqp_nonzero <- bcmleqqqp[bcmleqqqp>0]
                    
                    
                    kl <- sum(mlep_nonzero*log(mlep_nonzero/p_nonzero)) ##Kullback-Leibler between true p and mle of p
                    klbc <- sum(bcmlep_nonzero*log(bcmlep_nonzero/p_nonzero)) ##Kullback-Leibler between true p and bias corrected 
                    klbcq <- sum(bcmleqp_nonzero*log(bcmleqp_nonzero/p_nonzero)) ##Kullback-Leibler between true p and bias corrected 
                    klbcqq <- sum(bcmleqqp_nonzero*log(bcmleqqp_nonzero/p_nonzero)) ##Kullback-Leibler between true p and bias corrected 
                    klbcqqq <- sum(bcmleqqqp_nonzero*log(bcmleqqqp_nonzero/p_nonzero)) ##Kullback-Leibler between true p and bias corrected 
                    
                    simfinal[[k]][sim,] <- c(est_lambda_psi,
                                             mlep, bcmlep, bcmleqp,bcmleqqp,bcmleqqqp,
                                             eu, eubc, eubcq, eubcqq, eubcqqq, 
                                             kl, klbc, klbcq, klbcqq, klbcqqq)
                    
                }else{
                    mle <- MLE(Nk, N)
                    bcmle <- BCMLE(Nk, N)
                    hbcmle1 <- HBCMLE1(Nk, N)
                    hbcmle2 <- HBCMLE2(Nk, N)
                    hbcmle3 <- HBCMLE3(Nk, N)
                    
                    est_lambda_psi <- c(mle[[1]], mle[[2]],
                                        bcmle[[1]], bcmle[[2]],
                                        hbcmle1[[1]], hbcmle1[[2]], 
                                        hbcmle2[[1]], hbcmle2[[2]],
                                        hbcmle3[[1]], hbcmle3[[2]])
                    
                    mlep <- mle[[3]]
                    bcmlep <- bcmle[[3]] 
                    bcmleqp <- hbcmle1[[3]] 
                    bcmleqqp <- hbcmle2[[3]] 
                    bcmleqqqp <- hbcmle3[[3]] 
                    
                    eu <- sqrt(sum((p - mlep)^2)) ##Euclidean distance between true p and mle of p
                    eubc <- sqrt(sum((p - bcmlep)^2)) ##Euclidean distance between true p and bias corrected 
                    eubcq <- sqrt(sum((p - bcmleqp)^2)) ##Euclidean distance between true p and bias corrected
                    eubcqq <- sqrt(sum((p - bcmleqqp)^2)) ##Euclidean distance between true p and bias corrected
                    eubcqqq <- sqrt(sum((p - bcmleqqqp)^2)) ##Euclidean distance between true p and bias corrected
                    
                    p_nonzero <- p[mlep>0]
                    mlep_nonzero <- mlep[mlep>0]
                    bcmlep_nonzero <- bcmlep[bcmlep>0]
                    bcmleqp_nonzero <- bcmleqp[bcmleqp>0]
                    bcmleqqp_nonzero <- bcmleqqp[bcmleqqp>0]
                    bcmleqqqp_nonzero <- bcmleqqqp[bcmleqqqp>0]
                    
                    kl <- sum(mlep_nonzero*log(mlep_nonzero/p_nonzero)) ##Kullback-Leibler between true p and mle of p
                    klbc <- sum(bcmlep_nonzero*log(bcmlep_nonzero/p_nonzero)) ##Kullback-Leibler between true p and bias corrected 
                    klbcq <- sum(bcmleqp_nonzero*log(bcmleqp_nonzero/p_nonzero)) ##Kullback-Leibler between true p and bias corrected 
                    klbcqq <- sum(bcmleqqp_nonzero*log(bcmleqqp_nonzero/p_nonzero)) ##Kullback-Leibler between true p and bias corrected 
                    klbcqqq <- sum(bcmleqqqp_nonzero*log(bcmleqqqp_nonzero/p_nonzero)) ##Kullback-Leibler between true p and bias corrected 
                    
                    simfinal[[k]][sim,] <- c(est_lambda_psi, mlep, bcmlep, bcmleqp,bcmleqqp,bcmleqqqp,
                                             eu, eubc, eubcq, eubcqq, eubcqqq, kl, klbc, klbcq, klbcqq, klbcqqq)
                }
            }
        }
    }
    for(d in 1:sz){
        B <- c(lambda, apply(simfinal[[d]],2,mean),apply(simfinal[[d]],2,var), ssize[d], p)
        write.table(t(B),paste(path,"/data-n",toString(n), "-maxfreq", round(p[1],digits = 2), ".txt",sep=""),
                    append=TRUE, sep=" ", col.names=FALSE, row.names=FALSE)
    }
}


######################################################################


#' For nested simulation (above), generates random dataset of size N
#'
#' @param Nk integer vector; number of lineage prevalence counts in a dataset.
#'   for a simulated data this is simply derived as \code{colSums(dataset)}. To derive
#'   the MLE and lineage prevalence counts for a real dataset please refer to
#'   the package \link[MLMOI]{moimle}.
#' @param N integer; sample size
#' @param lambda numeric; MOI parameter
#' @param p vector; lineage-frequency distribution
#' @param n integer; number of lineages
#'
#' @return new random of size N and Nk's 
#' @export
#'
innersamplegenerator_CP <- function(Nk,N,lambda,p,n) {
    while (is.na(Nk[1]) == TRUE){
        MN <- mnom(N, lambda, p, n)
        Nk <- colSums(MN)
        if (sum(Nk) <= N || max(Nk) == N){
            Nk <- NA
        }else{
            M <- MN
        }
    } 
    list(M,Nk)
}

