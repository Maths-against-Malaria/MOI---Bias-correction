######################################################################
#------------------------Estimating parameters------------------------
######################################################################

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

######################################################################

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

######################################################################

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

######################################################################

HBCMLE2 <- function(N, Nk){
    mle <- MLE(N, Nk)
    mle_lam <- mle[[2]]
    mle_p <- mle[[4]]
    
    bcmle <- BCMLE(N, Nk)
    bcmle_lam <- bcmle[[1]]
    bcmle_p <- bcmle[[3]]
    
    bias <- second_order_bias(N, bcmle_lam, bcmle_p)               #second-order bias evaluated at the BCMLE
    bias_lam <- bias[[1]]
    bias_p <- bias[[2]]
    
    p_pathologic <- prob_pathological(N, bcmle_lam, bcmle_p)    #probability of pathological data evaluated at the BCMLE
    p_regular <- 1 - p_pathologic                               #probability of regular data evaluated at the BCMLE
    
    hbcmle_2_lam <- p_regular*(mle_lam - bias_lam)              #HBCMLE2 of lambda
    hbcmle_2_psi <- hbcmle_2_lam/(1 - exp(-hbcmle_2_lam))       #HBCMLE2 of psi
    hbcmle_2_p <- mle_p - bias_p                                #HBCMLE2 of lineage frequencies
    
    out <- list(hbcmle_2_lam, hbcmle_2_psi, hbcmle_2_p)
    out
}

######################################################################

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

######################################################################

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

######################################################################

prob_pathological <- function(N, lambda, p){
    lep <- lambda*p
    dk <- exp(lep) - 1
    d <- sum(dk)
    d0 <- 1/(exp(lambda) - 1)
    x <- (1 - d*d0) 
    y <- N*(d0 + 1)
    den <- y*x
    nom<- (d0 + 1/2)*d - d0*((d^2) - sum(dk^2))/(2*x)
    
    q1 <- (d*d0)^N 
    q2 <- sum((dk*d0)^N)
    q3 <- (1 - prod(1 - (1 - exp(-lep))^N))/(1 - exp(-lambda))^N
    
    q1[is.nan(q1)==T] <- 0
    q2[is.nan(q2)==T] <- 0
    q3[is.nan(q3)==T] <- 0
    
    p_pathologic <- q3 + q1 - q2
    p_pathologic
}

######################################################################

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

######################################################################
#----------------------------Generating data--------------------------
######################################################################

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

######################################################################

mnom <- function(m, p) { 
    N <-length(m)
    out<-array(0,dim=c(N,length(p)))
    for(k in 1:N){
        out[k,]=rmultinom(1,m[k],p)
    }
    out
}

######################################################################

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

######################################################################