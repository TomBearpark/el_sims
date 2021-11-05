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
seed <- 8894; set.seed(seed)

run_sim <- function(i, k, n, X.sigma = "I", rho = NULL, blim, 
                    write_fail = FALSE, fix_beta = FALSE){
  print(i)
  
  data.obj <- gen_data(k = k, n = n, X.sigma = X.sigma, rho = rho, 
                       blim = blim, fix_beta = fix_beta)
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
  
  ## Store results
  
  times <- tibble(i = i,          
                  ml = ml$time, 
                  mom = mom$time, 
                  gmm = gmm2$time, 
                  cue = cue$time,
                  el = el$time)
  
  fails <- tibble(i = i,          
                  ml = ml$converge, 
                  mom = mom$converge, 
                  gmm = gmm2$converge, 
                  cue = cue$converge,
                  el = el$converge)
  
  coefs <- tibble(i = i, 
         var = names(data.obj$df)[-1], 
         beta = beta, 
         ml = ml$beta.hat, 
         mom = mom$beta.hat, 
         gmm = gmm2$beta.hat, 
         cue = cue$beta.hat,
         el = el$beta.hat) 
  
  # if(sum(fails[,-1])&write_fail){
  #   write_csv(data.obj$df, file = paste0(tab_loc, "failed_data/data_obj",i,".csv"))
  #   write_csv(tibble(beta), file = paste0(tab_loc, "failed_data/beta",i,".csv"))
  # }
  
  return(list(times = times, fails = fails, coefs = coefs))
}

run_study <- function(n, k, X.sigma, rho = 0, 
                      ncores = 50, ndraws = 1000, seed = 8894, blim = 1, 
                      fix_beta = FALSE){
  
  # Tag for file names
  var_tag <- get_var_tag(X.sigma = X.sigma, rho = rho)
  btag <- get_b_tag(blim)
  beta_tag <- get_beta_tag(fix_beta)
  
  file <- paste0("sim", ndraws, "_k", k, "_n", n, var_tag, btag, beta_tag, ".csv")
  
  
  # Set up parallel compute
  message(file)
  plan(multisession, workers = ncores)
  results <- future_map(1:ndraws, run_sim, k = k, n = n, 
                            X.sigma = X.sigma, rho = rho, blim = blim, 
                            fix_beta = fix_beta, 
                          .options = furrr_options(seed = seed), 
                          .progress = TRUE)
  
  coefs <- map_dfr(1:ndraws, function(x) results[[x]]$coefs)
  times <- map_dfr(1:ndraws, function(x) results[[x]]$times)
  converge <- map_dfr(1:ndraws, function(x) results[[x]]$fails)
  
  write_csv(coefs, file = file.path(tab_loc, "coefs/", file))
  write_csv(times, file = file.path(tab_loc, "times/", file))
  write_csv(converge, file = file.path(tab_loc, "converge/", file))
}


# Local testing... 

# n <- 1000
# for(k in c(5, 10)){
#   run_study(n = n, k = k, X.sigma = "diagish", 
#             ndraws = 1000, ncores = 50, rho = 0.5, fix_beta = TRUE)
# }



# Run stuff on the server... 
full_run <- FALSE

if(full_run){
  ndraws <- 5000
  ncores <- 50
  n      <- 1000
  
  for(k in c(2, 5, 10)){
    run_study(n = n, k = k, X.sigma = "I", ndraws = ndraws, ncores = ncores)
    for(rho in c(0.5, 0.9)){
      run_study(n = n, k = k, X.sigma = "decay",   rho = rho, ndraws = ndraws, ncores = ncores)  
      run_study(n = n, k = k, X.sigma = "diagish", rho = rho, ndraws = ndraws, ncores = ncores)
    }
  }
}



