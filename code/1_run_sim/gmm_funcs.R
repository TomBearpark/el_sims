#   ____________________________________________________________________________
#   GMM functions                                                           ####
#   For details of the gmm package, see
###    https://cran.csiro.au/web/packages/gmm/vignettes/gmm_with_R.pdf

require(dplyr)
require(purrr)

# Helper function for getting all combinations of a list, returns a handy tibble
gen_combinations <- function(x_list) {
  n <- length(x_list)
  ll <- tibble(v1 = NA, v2 = NA, .rows = choose(n, 2) + n)
  count <- 1
  for (ii in 1:n) {
    for (jj in 1:n) {
      if (ii > jj) next
      ll$v1[count] <- x_list[ii]
      ll$v2[count] <- x_list[jj]
      count <- count + 1
    }
  }
  ll
}

# Generate vector of moments to minimize in GMM procedure for all 55 moms.

moments_two_step <- function(theta, data) {
  X <- as.matrix(data[,-1])
  y <- data$y
  
  colX <- ncol(X)
  
  eps_hat <- (y - pnorm(as.matrix(X) %*% theta))[,1,drop=T]
  G <- apply(X,2,function(x) X*x*eps_hat) 
  G <- matrix(G, nrow = nrow(X), ncol = colX^2)
  
  selcol <- c()
  for (k in seq_len(colX)) {
    trues  <- rep(T,colX-k+1)
    falses <- rep(F,k-1)
    selcol <- append(selcol,c(falses,trues))
  }
  
  G[,selcol]
}

D_two_step <- function(theta, data){
 stop("not implemented yet") 
}

# Generate moments for method of moments estimator 
moments_mom <- function(theta, data) {
  X <- as.matrix(data[,-1])
  y <- data$y
  
  eps_hat <- y - pnorm(as.matrix(X) %*% theta)
  G <- apply(X,2,function(x) x*eps_hat) 
  
  G
}

# Derivative for the method of moments estimator
D_mom <- function(theta, data){
  X <- as.matrix(data[,-1])
  - dnorm(theta %*% X ) 
}


## Wrapper function for running the gmms
est.GMM <- function(data, type = "twoStep", init = NULL) {
  stopifnot(type %in% c("mom", "twoStep", "cue", "EL"))

  # Initialise with linear regression
  if(is.null(init)){
    init <- lm(y ~ . - 1, data)$coefficients %>% as.matrix()  
  }
  
  # Start timer
  tic()
  iter.max <- 100000
  # Get the right moment function
  if (type == "mom") {
  
    moments <- moments_mom
    est <- gmm(g = moments, x = data, t0 = init, wmatrix = "ident", 
               optfct="nlminb", control=list(iter.max=iter.max))  # this is exaclty identified, should be the same without W
    converge <- est$algoInfo$convergence
 
  } else if (type %in% c("twoStep")) {
    
    moments <- moments_two_step
    est <- gmm(g = moments, x = data, t0 = init, type = "twoStep", 
               optfct="nlminb", control=list(iter.max=iter.max)) 
    converge <- est$algoInfo$convergence
  
  } else if (type == "cue"){
    
    moments <- moments_two_step
    est <- gel(g = moments, x = data, tet0 = init, type = "CUE", 
               optfct="nlminb", control=list(iter.max=iter.max)) 
    converge <- est$conv_par
  
  } else if (type == "EL"){
    moments <- moments_two_step
    est <- gel(g = moments, x = data, tet0 = init, type = "EL", 
               optfct="nlminb", control=list(iter.max=iter.max)) 
    converge <- est$conv_par
  }
  # End the timer
  tt <- toc()
  
  list(beta.hat = est$coefficients, converge = converge, 
       time = tt$toc - tt$tic)
}
