n <- 100
k <- 10
beta <- matrix(rep(1, k + 1), ncol = 1)

gen_data <- function(k, n, beta){
  
  # Initialise covariate matrix, always including a constant 
  X <- matrix(nrow = n, ncol = k + 1)
  X[,1] <- 1
  
  # Fill in random draws for the covariates
  for(ii in 2:(k + 1)) X[,ii] <- rnorm(n)
  
  # Draw a error term: normal due to probit
  epsilon <- rnorm(n)
  
  # Generate binary Y values
  Y <- 1*(X %*% beta >= epsilon) 
  
  # Clean up outputs 
  out <- data.frame(Y, X) 
  names(out) <- c("y", paste0("x", seq(0, k)))
  out
}

# df <- gen_data(k, n, beta)
# sum(df$y) / length(df$y)

