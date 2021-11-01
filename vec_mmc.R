

moments_mom_alt <- function(theta, data) {
  X <- as.matrix(data[,-1])
  y <- data$y
  
  eps <- y - pnorm(as.matrix(X) %*% theta)
  G <- apply(X,2,function(x) x*eps) 
  
  G
}


moments_two_step <- function(theta, data) {
  X <- as.matrix(data[,-1])
  y <- data$y
  
  colX <- ncol(X)
  
  eps <- (y - pnorm(as.matrix(X) %*% theta))[,1,drop=T]
  G <- apply(X,2,function(x) X*x*eps) 
  G <- matrix(G, nrow = nrow(X), ncol = colX^2)
  
  selcol <- c()
  for (k in seq_len(colX)) {
    trues  <- rep(T,colX-k+1)
    falses <- rep(F,k-1)
    selcol <- append(selcol,c(falses,trues))
  }
  
  G[,selcol]
}