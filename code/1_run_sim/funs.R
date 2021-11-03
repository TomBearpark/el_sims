gen_sigma <- function(k, rho){
  sigma <- diag(nrow = k)
  for(i in 1:k) for(j in 1:k) sigma[i,j] <- rho^(abs(i - j))
  sigma
}

gen_data <- function(beta = NULL, k = 10, n = 1000, constant = TRUE, 
                     X.mu = NULL, X.sigma = "I", rho = 0){
  
  if(X.sigma == "I") {
    message ("using identity covariance matrix for X")
  } else if(X.sigma == "decay") {
    if(rho == 0) stop("need to specify rho")
  }else{
    stop("that structure on X.sigma is not currently implemented!")
  }
  
  # If not specified the true betas are drawn from a uniform on [-1,1]
  if (is.null(beta)) beta <- runif(k, min = -1, max = 1) 
  
  # If not specified covariates are independent
  if (X.sigma == "I")  {
    X.sigma <- diag(k)
  }else{
    X.sigma <- gen_sigma(k = k, rho = rho)
  }
  
  # and centered on 0
  if (is.null(X.mu)) X.mu <- rep(0,k)
  
  # Fill in random draws for the covariates 
  X <- MASS::mvrnorm(n = n, mu = X.mu, Sigma = X.sigma)

  # add constant if requested
  if (constant == TRUE) X[,1] <- 1
  
  # Draw a error term: normal due to probit
  epsilon <- rnorm(n)
  
  # Generate binary Y values
  Y <- 1*(X %*% beta >= epsilon) 
  
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
  
  return(list(df = out, eps = epsilon, model.specs = model.specs))
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

write_fail <- function(data.obj){
  data.obj
}

get_var_tag <- function(X.sigma, rho){
  if(X.sigma == "I"){
    var_tag <- ""
  }else if(X.sigma == "decay"){
    var_tag <- paste0("_decay_rho", str_replace(rho, "[.]", "_"))
  }
  var_tag
}
