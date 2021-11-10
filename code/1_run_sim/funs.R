gen_sigma <- function(k, rho){
  sigma <- diag(nrow = k)
  for(i in 1:k) for(j in 1:k) sigma[i,j] <- rho^(abs(i - j))
  sigma
}

gen_sigma_less_fancy <- function(k, rho){
  sigma <- matrix(rho, nrow = k, ncol = k)
  diag(sigma) <- 1
  sigma
}

gen_data <- function(k = 10, n = 1000, constant = TRUE, 
                     X.mu = NULL, X.sigma = "I", rho = 0, blim = 1, fix_beta = FALSE){
  
  if(X.sigma == "I") {
    message ("using identity covariance matrix for X")
    X.sigma <- diag(k)
  } else if(X.sigma == "decay") {
    if(rho == 0) stop("need to specify rho")
    X.sigma <- gen_sigma(k = k, rho = rho)
  }else if(X.sigma == "diagish"){
    if(rho == 0) stop("need to specify rho")
    X.sigma <- gen_sigma_less_fancy(k = k, rho = rho)
  }else{
    stop("that structure on X.sigma is not currently implemented!")
  }
  
  # If not specified the true betas are drawn from a uniform on [-1,1]
  if (fix_beta == FALSE) {
    beta <- runif(k, min = -blim, max = blim) 
  }else{
    beta <- rep(1, k)
  }
 	beta <- matrix(beta) 
  # and centered on 0
  if (is.null(X.mu)) X.mu <- rep(0,k)
  
  # Fill in random draws for the covariates 
  X <- MASS::mvrnorm(n = n, mu = X.mu, Sigma = X.sigma)

  # add constant if requested
  if (constant == TRUE) X[,1] <- 1
  
  # Draw a error term: normal due to probit
  u <- rnorm(n)
  
  # Generate binary Y values
  Y <- 1*(X %*% beta >= u) 
  
  # Calculate the population residuals 
  eps <- Y - dnorm(X %*% beta)
  
  # Clean up outputs 
  out <- data.frame(Y, X) 
  names(out) <- c("y", paste0("x", seq(1, k)))
  
  model.specs <- list('beta' = beta,
                      'regressors' = k,
                      'sample' = n,
                      'constant' = constant,
                      'X.mu' = X.mu,
                      'X.sigma' = X.sigma, 
                      'rho' = rho)
  
  return(list(df = out, u = u, model.specs = model.specs, eps = eps))
}


est.ML <- function(data) {
  
  tic()
  est <- glm(formula = y ~ . - 1, 
                family = binomial(link = "probit"), 
                data = data)
  tt <- toc()
  
  # I'm coding it as a list to allow us considering SEs later
  return(list(beta.hat = est$coefficients, converge =1- 1*est$converged, 
              time = tt$toc - tt$tic))
}

get_b_tag <- function(blim){
  if(!blim == 1){
    btag <- paste0("_b", str_replace(blim, "[.]", "_"))
  }else{
    btag <- ""
  }
  btag
}

get_var_tag <- function(X.sigma, rho){
  if(X.sigma == "I"){
    var_tag <- ""
  }else if(X.sigma == "decay"){
    var_tag <- paste0("_decay_rho", str_replace(rho, "[.]", "_"))
  }else if(X.sigma == "diagish"){
    var_tag <- paste0("_diagish_rho", str_replace(rho, "[.]", "_"))
  }else{
    stop("Didn't implement this X.sigma, cheeky!")
  }
  var_tag
}


get_beta_tag <- function(fix_beta){
  if (fix_beta) {
    beta_tag <- "_beta1" 
  }else {
    beta_tag <- ""
  }
  beta_tag
}
