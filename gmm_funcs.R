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

  # Format
  X <- data %>% dplyr::select(starts_with("x"))
  y <- data$y

  # Get a list of all the moment conditions we want to include
  moms <- gen_combinations(names(X))
  # Pre compute epsilons
  eps <- y - pnorm(as.matrix(X) %*% theta)

  # Get a vector of the sample moments for each moment condition
  G <- map_dfc(
    1:dim(moms)[1],
    function(ii) {
      (X[moms$v1[ii]] * X[moms$v2[ii]] * eps)
    }
  )

  # Return n*q matrix, q is number of moments
  as.matrix(G)
}

# Generate vector of moments to minimize in GMM procedure for all 55 moms.
moments_mom <- function(theta, data) {

  # Format
  X <- data %>% dplyr::select(starts_with("x"))
  y <- data$y

  # Pre compute epsilons
  eps <- y - pnorm(as.matrix(X) %*% theta)

  # Get a vector of the sample moments for each moment condition
  G <- map_dfc(
    names(X),
    function(ii) {
      (X[ii] * eps)
    }
  )

  # Return n*k matrix, k is number of moments/regressors
  as.matrix(G)
}



est.GMM <- function(data, type = "twoStep") {
  stopifnot(type %in% c("mom", "twoStep", "cue"))

  # Initialise with linear regression
  init <- lm(y ~ . - 1, data)$coefficients %>% as.matrix()

  # Get the right moment function
  if (type == "mom") {
    moments <- moments_mom
    est <- gmm(g = moments, x = data, init, wmatrix = "ident")
  } else if (type %in% c("cue", "twoStep")) {
    moments <- moments_two_step
    est <- gmm(g = moments, x = data, init, type = type)
  } 

  list(beta.hat = est$coefficients)
}
