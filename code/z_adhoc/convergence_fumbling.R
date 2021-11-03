###############################################################################
## Load stuff
if(!require(pacman)) install.packages('pacman')
pacman::p_load('MASS', 'tidyverse', 'gmm', 'furrr', 'tictoc')
theme_set(theme_bw())

## TO do - set an out_dir outside of github so we don't save huge files there 

# Set directories
user <- Sys.getenv("USER")
if( user == "tombearpark"){
  dir <- file.path('/Users/tombearpark/Documents/GitHub/el_sims/')
}else if (user == "bearpark"){
  # For the server
  dir <- file.path("/home/bearpark/el_sims/")
} else{
  dir <- 'D:\\Dropbox\\Università\\PhD\\II Year\\Fall\\ECO519 - Non-linear Econometrics\\psets\\el_sims\\code\\'
} 

fig_loc <- file.path(dir, "fig/")
tab_loc <- file.path(dir, "out/")
dir.create(fig_loc, showWarnings = FALSE); dir.create(tab_loc, showWarnings = FALSE)

# Load auxiliary funs
source(file.path(dir, "code/1_run_sim/funs.R"))
source(file.path(dir, "code/1_run_sim/gmm_funcs.R"))

###############################################################################
set.seed(1)

D_mom <- function(theta, data){
  X <- as.matrix(data[,-1])
  - t(X) %*% diag(dnorm(X%*%theta)[,1,drop=T])
}

D_mom_tom <- function(theta, data){
  X <- as.matrix(data[,-1])
  y <- data$y
  k <- dim(X)[2]
  eps_hat <- y - pnorm(as.matrix(X) %*% theta)
  D <- matrix(ncol = k, nrow = k)
  for(k1 in 1:k){
    for(k2 in 1:k){
      D[k1, k2] <- mean(X[,k1] * X[,k2] * eps_hat)
    }
  }
  D
}

k <- 10
data.obj <- gen_data(k = k, n = 1000, X.sigma = "I", rho = 0)
data <- data.obj$df
beta <- data.obj$model.specs$beta

D_mom(theta = beta, data = data)
D_mom_tom(theta = beta, data = data)


init <- lm(y ~ . - 1, data)$coefficients %>% as.matrix()  

## See how convergence is looking... 

## mom
gmm(g = moments_mom, x = data, t0 = init, wmatrix = "ident", 
           gradv = D_mom) 

gmm(g = moments_mom, x = data, t0 = init, wmatrix = "ident", 
           gradv = D_mom_tom) 

gmm(g = moments_mom, x = data, t0 = init, wmatrix = "ident", 
    optfct="nlminb", control=list(iter.max=10000)) 

gmm(g = moments_mom, x = data, t0 = init, wmatrix = "ident", 
     control=list(maxit=1000), method = "BFGS")

gmm(g = moments_mom, x = data, t0 = init, wmatrix = "ident", 
    control=list(maxit=100000))

## two step
gmm(g = moments_two_step, x = data, t0 = init, type = "twoStep")
gmm(g = moments_two_step, x = data, t0 = init, type = "twoStep", 
    control=list(maxit=100000))
    optfct="nlminb", control=list(iter.max=10000), 
    lower = rep(-1, k), upper = rep(1, k)) 
beta

## CUE
gel(g = moments_two_step, x = data, tet0 = init, type = "CUE")
gel(g = moments_two_step, x = data, tet0 = init, type = "CUE", 
    optfct="nlminb", control=list(iter.max=10000))

## EL
gel(g = moments_two_step, x = data, tet0 = init, type = "EL")
gel(g = moments_two_step, x = data, tet0 = init, type = "EL", 
    optfct="nlminb", control=list(iter.max=10000))


