###############################################################################
## Load stuff
if(!require(pacman)) install.packages('pacman')
pacman::p_load('MASS', 'tidyverse', 'gmm', 'doParallel','doRNG', 'furrr')
theme_set(theme_bw())

# Set directories
user <- Sys.getenv("USER")
if( user == "tombearpark"){
  dir <- file.path('/Users/tombearpark/Documents/GitHub/el_sims/')
} else {
  dir <- 'D:\\Dropbox\\Università\\PhD\\II Year\\Fall\\ECO519 - Non-linear Econometrics\\psets\\el_sims\\'
} 

fig_loc <- file.path(dir, "fig/")
tab_loc <- file.path(dir, "tab/")

# Load auxiliary funs
source(file.path(dir,"funs.R"))
source(file.path(dir,"gmm_funcs.R"))

###############################################################################
set.seed(8894)

n <- 1000
k <- 10

run_sim <- function(i, k, n){
  
  data.obj <- gen_data(k = k, n = n)
  data.df <- data.obj$df
  beta <- data.obj$model.specs$beta
  
  ## Maximum likelihood
  ml   <- est.ML(data.df)
  
  ## Method of moments
  mom  <- est.GMM(data.df, type = "mom")
  
  ## Twostep GMM
  gmm2 <- est.GMM(data.df, type = "twoStep")
  
  ## CUE
  cue  <- est.GMM(data.df, type = "cue", init = beta)
  
  ## EL
  el  <- est.GMM(data.df, type = "EL", init = beta)
  
  ## Store results
  tibble(i = i, 
         beta = beta, 
         ml = ml$beta.hat, 
         mom = mom$beta.hat[,1], 
         gmm = gmm2$beta.hat[,1], 
         cue = cue$beta.hat,
         el = el$beta.hat) 
}
plan(multisession, workers = 8)

results <- future_map_dfr(1:1000, run_sim, k = k, n = n)
