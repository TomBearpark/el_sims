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
  - t(X) %*% diag(dnorm(X%*%theta)[,1,drop=T]) %*% X
}

data.obj <- gen_data(k = 10, n = 1000, X.sigma = "I", rho = 0)
data <- data.obj$df
beta <- data.obj$model.specs$beta


init <- lm(y ~ . - 1, data)$coefficients %>% as.matrix()  
est <- gmm(g = moments_mom, x = data, t0 = init, wmatrix = "ident", 
           itermax=1000000, gradv = D_mom) # 






