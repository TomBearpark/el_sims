################################################################################
################################################################################

#################
## Example (How to use get.aVar)
#
# Generate the data JUST ONCE at the beginning
#
# nsim    <- 100000
# data    <- gen_data(n = nsim)
# Xdraw   <- as.matrix(data$df[,-1])
# Xtilde  <- cross_mom(Xdraw)
#
#
# just loop over betas the following mf
# aux <- get.aVar(Xdraw,Xtilde,data$model.specs$beta)
#
#
#################


# Main Function
get.aVar <- function(X.draw,X.tilde.draw,beta) {
  # Asymptotic Variance of MLE
  message("------mle-----------")
  
  mle <- avar.mle(X = X.draw, beta = beta)
  message("------mom-----------")
  # Asymptotic Variance of MoM
  mom <- avar.mom(X = X.draw, beta = beta)
  # Asymptotic Variance of GMM/CUE/EL
  message("------gmm-----------")
  gmm <- avar.gmm(X = X.draw, XX = X.tilde.draw, beta = beta)
  
  return(list(aVar.mle = mle,
              aVar.mom = mom,
              aVar.gmm = gmm))
}


################################################################################
################################################################################
## Auxiliary functions

avar.mle <-  function(X, beta) {
  sims <- nrow(X)
  k    <- length(beta)
  
  store <- matrix(0, k, k)
  for (s in seq_len(sims)) {
    xi <- X[s,]
    xb <- sum(xi*beta)
    
    if(pnorm(xb) %in% c(0,1)) next
    
    varcov <- (dnorm(xb)^2/(pnorm(xb)*(1-pnorm(xb))))*xi%*%t(xi)
    
    store <- varcov + store    
  }  
  
  avar <- solve(store/sims)
  
  return(diag(avar))
}

avar.mom <- function(X,beta) {
  sims <- nrow(X)
  k <- length(beta)

  G.store <- matrix(0, k, k)
  V.store <- matrix(0, k, k)
  for (s in seq_len(sims)) {
    xi  <- X[s,]
    xb  <- sum(xi*beta)
    Phi <- pnorm(xb)
    xx  <- xi%*%t(xi)
    if(dnorm(xb) %in% c(0,1)) next
    G <- dnorm(xb)*xx
    V <- Phi*(1-Phi)*xx
    
    G.store <- G.store + G
    V.store <- V.store + V
  }

  G <- G.store/sims; V <- V.store/sims 
  
  Ginv <- solve(G)
  
  avar <- Ginv %*% V %*% t(Ginv)
  
  return(diag(avar))
}

cross_mom <- function(X) {
  colX <- ncol(X)
  
  G <- apply(X,2,function(x) X*x) 
  G <- matrix(G, nrow = nrow(X), ncol = colX^2)
  
  selcol <- c()
  for (k in seq_len(colX)) {
    trues  <- rep(T,colX-k+1)
    falses <- rep(F,k-1)
    selcol <- append(selcol,c(falses,trues))
  }
  
  return(G[,selcol])
}


avar.gmm <- function(X,XX,beta) {
  sims <- nrow(X)
  k    <- length(beta)
  ktil <- k*(k+1)/2
  
  G.store <- matrix(0, ktil, k)
  V.store <- matrix(0, ktil, ktil)
  for (s in seq_len(sims)) {
    xi    <- X[s,]
    xitil <- XX[s,]
    xb    <- sum(xi*beta)
    Phi <- pnorm(xb)
    if(pnorm(xb) %in% c(0,1)) next
    xx  <- xitil%*%t(xi)
    
    G <- dnorm(xb)*xx
    V <- Phi*(1-Phi)*xitil%*%t(xitil)
    
    G.store <- G.store + G
    V.store <- V.store + V
  }
  
  G <- G.store/sims; V <- V.store/sims 
  
  avar <- solve(t(G) %*% solve(V) %*% G)
  
  return(diag(avar))
}
