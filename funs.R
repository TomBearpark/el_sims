gen_data <- function(beta = NULL, k = 10, n = 1000, constant = TRUE, 
                     X.mu = NULL, X.sigma = NULL){
  
  # If not specified the true betas are drawn from a uniform on [-1,1]
  if (is.null(beta)) beta <- runif(k, min = -1, max = 1) 
  
  # If not specified covariates are independent
  if (is.null(X.sigma))  X.sigma <- diag(k)
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
                      'X.sigma' = X.sigma)
  
  return(list(df = out, eps = epsilon, model.specs = model.specs))
}


est.ML <- function(data) {
  est <- glm(formula = y ~ . - 1, 
                family = binomial(link = "probit"), 
                data = data)
  
  # I'm coding it as a list to allow us considering SEs later
  results <- list(beta.hat = est$coefficients)
  return(results)
}


est.ALL <- function() {
  ## ML probit model
  
  ## Method of moments estimator
  
  ## 1-step GMM with I
  
  ## Two-step efficient GMM
  
  ## CUE
  
  ## EL
}

