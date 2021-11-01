## A little bit of timing
pacman::p_load('microbenchmark','gmm', 'tidyverse')
user <- Sys.getenv("USER")
if( user == "tombearpark"){
  dir <- file.path('/Users/tombearpark/Documents/GitHub/el_sims/')
} else {
  dir <- 'D:\\Dropbox\\Università\\PhD\\II Year\\Fall\\ECO519 - Non-linear Econometrics\\psets\\el_sims\\'
} 

source(file.path(dir,"funs.R"))
source(file.path(dir,"gmm_funcs.R"))
source(file.path(dir,"vec_mmc.R"))


## Load data and initialise with OLS
n <- 10000
k <- 5
data.obj <- gen_data(k = k, n = n)
data <- data.obj$df
beta <- data.obj$model.specs$beta
init <- lm(y ~ . - 1, data)$coefficients %>% as.matrix()

## Compare for mom
# find that filippo's  version is abotu 30% faster, so swap that out
est <- gmm(g = moments_mom, x = data, t0 = init, wmatrix = "ident")
est2 <- gmm(g = moments_mom_alt, x = data, t0 = init, wmatrix = "ident")
m1 <- microbenchmark(
  est <- gmm(g = moments_mom, x = data, t0 = init, wmatrix = "ident"),
  est2 <- gmm(g = moments_mom_alt, x = data, t0 = init, wmatrix = "ident"), 
  times = 10
); m1

## Compare for full gmm
est  <- gmm(g = moments_two_step, x = data, t0 = init, type = "twoStep")
est2 <- gmm(g = moments_two_step_alt, x = data, t0 = init, type = "twoStep")
m2 <- microbenchmark(
  est  <- gmm(g = moments_two_step, x = data, t0 = init, type = "twoStep"),
  est2 <- gmm(g = moments_two_step_alt, x = data, t0 = init, type = "twoStep"), 
  times = 10
); m2
