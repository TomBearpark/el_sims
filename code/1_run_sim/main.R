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
tab_loc <- file.path(dir, "tab/")
dir.create(fig_loc, showWarnings = FALSE); dir.create(tab_loc, showWarnings = FALSE)

# Load auxiliary funs
source(file.path(dir, "code/1_run_sim/funs.R"))
source(file.path(dir, "code/1_run_sim/gmm_funcs.R"))

###############################################################################
seed <- 8894; set.seed(seed)

n <- 1000
k <- 10

run_sim <- function(i, k, n, X.sigma = "I", rho = NULL){
  print(i); tic()
  
  data.obj <- gen_data(k = k, n = n, X.sigma = X.sigma, rho = rho)
  data.df <- data.obj$df
  beta <- data.obj$model.specs$beta
  
  ## Maximum likelihood
  ml   <- est.ML(data.df)
  
  ## Method of moments
  mom  <- est.GMM(data.df, type = "mom")
  
  ## Twostep GMM
  gmm2 <- est.GMM(data.df, type = "twoStep")
  
  ## CUE
  cue  <- est.GMM(data.df, type = "cue")
  
  ## EL
  el  <- est.GMM(data.df, type = "EL")
  
  tt <- toc()
  
  ## Store results
  tibble(i = i, 
         var = names(data.obj$df)[-1], 
         beta = beta, 
         ml = ml$beta.hat, 
         mom = mom$beta.hat[,1], 
         gmm = gmm2$beta.hat[,1], 
         cue = cue$beta.hat[,1],
         el = el$beta.hat[,1], 
         time = tt$toc - tt$tic
         ) 
}

run_study <- function(n, k, X.sigma, rho = 0, 
                      ncores = 50, ndraws = 1000, seed = 8894){

  plan(multisession, workers = ncores)
  results <- future_map_dfr(1:ndraws, run_sim, k = k, n = n, 
                            X.sigma = X.sigma, rho = rho, 
                          .options = furrr_options(seed = seed), 
                          .progress = TRUE)
  
  if(X.sigma == "I"){
    var_tag <- ""
  }else if(X.sigma == "decay"){
    var_tag <- paste0("_decay_rho", str_replace(rho, "[.]", "_"))
  }
  write_csv(results, 
            file = file.path(tab_loc, 
                      paste0("sim", ndraws, "_k", k, "_n", n, var_tag ,".csv")))
}

# # Naive first step
# run_study(n = 1000,  k = 10, X.sigma = "I", ncores = 50)
run_study(n = 10000, k = 10, X.sigma = "I", ncores = 50)

# run_study(n = 1000, k = 10, X.sigma = "decay", rho = 0.5, ndraws = 1000, ncores = 50)

# results %>%
#   mutate(across(ml:el, ~.x-beta)) %>%
#   pivot_longer(cols = ml:el) %>%
#   ggplot() +
#   geom_density(aes(x = value, color = var)) +
#   facet_wrap(vars(name))
