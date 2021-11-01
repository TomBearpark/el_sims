###############################################################################
## Load stuff
if(!require(pacman)) install.packages('pacman')
pacman::p_load('MASS', 'tidyverse', 'gmm', 'doParallel','doRNG', 'furrr', 'tictoc')
theme_set(theme_bw())

## TO do - set an out_dir outside of github so we don't save huge files there 

# Set directories
user <- Sys.getenv("USER")
if( user == "tombearpark"){
  code_dir <- file.path('/Users/tombearpark/Documents/GitHub/el_sims/')
}else if (user == "bearpark"){
  # For the server
  code_dir <- file.path("/home/bearpark/el_sims/")
} else{
  code_dir <- 'D:\\Dropbox\\Università\\PhD\\II Year\\Fall\\ECO519 - Non-linear Econometrics\\psets\\el_sims\\'
} 

fig_loc <- file.path(code_dir, "fig/")
tab_loc <- file.path(code_dir, "tab/")
dir.create(fig_loc, showWarnings = FALSE); dir.create(tab_loc, showWarnings = FALSE)

# Load auxiliary funs
source(file.path(code_dir,"funs.R"))
source(file.path(code_dir,"gmm_funcs.R"))

###############################################################################
seed <- 8894; set.seed(seed)

n <- 1000
k <- 10

run_sim <- function(i, k, n){
  print(i); tic()
  
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

# run_sim(i = 1, k = k, n = n)


plan(multisession, workers = 8)
results <- future_map_dfr(1:8, run_sim, k = k, n = n, 
                          .options = furrr_options(seed = seed), 
                          .progress = TRUE)

write_csv(results, file = file.path(tab_loc, 
                                    paste0("sim_k", k, "_n", n, ".csv")))

# results %>% 
#   mutate(across(ml:el, ~.x-beta)) %>% 
#   pivot_longer(cols = ml:el) %>% 
#   ggplot() + 
#   geom_density(aes(x = value)) + 
#   facet_wrap(~name, scales="free")
