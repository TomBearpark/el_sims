###########################################################################
### Set up and paths 

## Load stuff
if(!require(pacman)) install.packages('pacman')
pacman::p_load('tidyverse')
theme_set(theme_bw())

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

tab_loc <- file.path(code_dir, "out/")
out     <- file.path(code_dir, "results/")

fig_out <- file.path(out, "fig/")
tab_out <- file.path(out, "tab/")
  
# Source some functions 
source(file.path(code_dir, "code/1_run_sim/funs.R"))

###########################################################################
### 1. Functions 

load_data <- function(k, n, X.sigma, rho = 0, tab_loc, blim = 1, 
                      ndraws = 5000, fix_beta = FALSE){
  
  var_tag  <- get_var_tag(X.sigma = X.sigma, rho = rho)
  btag     <- get_b_tag(blim)
  beta_tag <- get_beta_tag(fix_beta)
  
  file <- paste0("sim", ndraws, "_k", k, "_n", n, var_tag, btag, beta_tag, ".csv")
  
  df      <- read_csv(file.path(tab_loc, "coefs/", file))
  conv    <- read_csv(file.path(tab_loc, "converge/", file))
  times   <- read_csv(file.path(tab_loc, "times/", file))
  
  df %>% 
    mutate(across(ml:el, ~(.x-beta))) %>%
    pivot_longer(cols = ml:el, values_to = "error", names_to = "estimator") %>% 
    left_join(pivot_longer(conv, ml:el, names_to = "estimator", 
                           values_to = "converge")) %>% 
    left_join(pivot_longer(times, ml:el, names_to = "estimator", 
                           values_to = "time")) %>% 
    mutate(k = !!k, n = !!n, x_var = !!var_tag, blim = btag, rho = rho) %>% 
    mutate(estimator = factor(estimator, levels = c("ml", "mom", "gmm", "cue", "el"))) %>% 
    mutate(covar = factor(var, levels = paste0("x", 1:k)))
}

plot_results <- function(df, var_name, scales = "free", xlim = TRUE){
  p <- ggplot(df) +
    geom_density(aes(x = error, color = covar, group = covar)) +
    facet_wrap(vars(estimator, k), scales = scales, nrow = 5) + 
    ggtitle(paste0("Estimated coefs minus the real values, X.sigma = ", var_name))
  if(xlim) p <- p + xlim(c(-1, 1))
  p
}

# mean-bias, median-bias, standard deviation, 
# median absolute di§erence from the median, 
# root mean squared error, and median absolute error.

rmse_tab <- function(df){
  df  %>% 
    mutate(sq_error = error^2) %>% 
    group_by(estimator, k, x_var) %>% 
    summarize(rmse = sqrt(mean(sq_error)), 
              mean_bias = mean(error), 
              median_bias = median(error), 
              sd = sd(error), 
              mae = median(abs(error)))
}

load_all_k <- function(k_list, n, X.sigma, rho, tab_loc){
  map_dfr(k_list, load_data, n = n, X.sigma = X.sigma, rho = rho, tab_loc = tab_loc) 
}

get_opts <- function(k_list, X.sigma_list, rho_list){
  expand_grid(k = k_list, X.sigma = X.sigma_list, rho = rho_list) %>% 
    mutate(rho = ifelse(X.sigma == "I", "0", rho)) %>% distinct()
}

load_all <- function(k_list = c(2,5,10), X.sigma_list = c("I", "diagish", "decay"), 
                       rho_list = c(.5, .9), n_draws = 5000, tab_loc){
  opts <- 
  pmap_dfr(opts, load_data, n =1000, tab_loc = tab_loc, ndraws = 5000)
}

sim_name <- function(X.sigma, rho){
  if(X.sigma == "I") return("I")
  else return(paste0(X.sigma, " ", rho))
}

bar_plot <- function(df, plot_var){
  ggplot(df) + 
    geom_col(aes(x = estimator, y = .data[[plot_var]])) + 
    facet_wrap(vars(k, x_var), scales = "free", nrow = 3)
}



###########################################################################
### 2.1 Error density plots

X.sigma_list = c("I", "diagish", "decay")
k_list = c(2,5,10)
rho_list = c(.5, .9)

opts <- get_opts(k_list = k_list, X.sigma_list = X.sigma_list, rho_list = rho_list)

## Load all the data
df <- load_all(tab_loc = tab_loc)


## Error density plots 
dir_out <- file.path(fig_out, "error_density/")
dir.create(dir_out, showWarnings = FALSE)

main <- opts %>% select(-k) %>% distinct()
for(ii in 1:dim(main)[1]){
  print(ii)
  val <- opts[ii,]
  X.sigma <- val$X.sigma; rho <- val$rho
  
  plot_df <- df %>% 
    filter(X.sigma == !!X.sigma, rho == !!rho)  
  
  out_name <- sim_name(X.sigma, rho)
  
  p <- plot_df %>% 
    plot_results(var_name = out_name)
  ggsave(p, file = paste0(dir_out, out_name, ".png"), height = 7, width = 12)
  
  p2 <- plot_df %>% 
    filter(converge == 0) %>% 
    plot_results(var_name = out_name, xlim = FALSE)
  ggsave(p2, file = paste0(dir_out, out_name, "free.png"), height = 7, width = 12)
}

###########################################################################
### 2.2 Bar plot comparisons

dir_out2 <- file.path(fig_out, "bar_comparison/")
dir.create(dir_out2, showWarnings = FALSE)

res <- df %>% rmse_tab()
comp_vars <- names(res)[4:length(names(res))]


for (cc in comp_vars){
  print(cc)
  
  dir_outcc <- file.path(fig_out, "bar_comparison/", cc, "/")
  dir.create(dir_outcc, showWarnings = FALSE)
  
  df %>% rmse_tab() %>% 
    bar_plot(plot_var = cc) %>% 
    ggsave(file = paste0(dir_outcc, "raw.png"), height = 7, width = 12)
  
  df %>% rmse_tab() %>% filter(estimator != "cue") %>% 
    bar_plot(plot_var = cc) %>% 
    ggsave(file = paste0(dir_outcc, "no_cue.png"), height = 7, width = 12)
  
  df %>% filter(between(error, -1, 1)) %>% rmse_tab()  %>% 
    bar_plot(plot_var = cc) %>% 
    ggsave(file = paste0(dir_outcc, "limit1_1.png"), height = 7, width = 12)
  
  df %>% filter(converge == 0) %>% rmse_tab()   %>% 
    bar_plot(plot_var = cc) %>% 
    ggsave(file = paste0(dir_outcc, "converged.png"), height = 7, width = 12)
  
  df %>% filter(converge == 0, estimator != "cue") %>% rmse_tab()   %>% 
    bar_plot(plot_var = cc) %>% 
    ggsave(file = paste0(dir_outcc, "no_cue_converged.png"), height = 7, width = 12)
}

df %>% 
  rmse_tab() %>%
  write_csv(file = paste0(dir_out2, "raw.csv"))


## Lets find a bad draw, check it out a bit... 
df %>% 
  filter(converge == 0) %>% 
  filter(abs(error) > 100000)  %>% 
  write_csv(file = paste0(tab_loc, "terrible_cues.csv"))
