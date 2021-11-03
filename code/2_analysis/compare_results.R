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

source(file.path(code_dir, "code/1_run_sim/funs.R"))

tab_loc <- file.path(code_dir, "out/")

### Load the data 
k <- 10
n <- 1000
ndraws <- 1000


load_data <- function(k, n, X.sigma, rho = 0, tab_loc){
  
  var_tag <- get_var_tag(X.sigma = X.sigma, rho = rho)
  file    <- paste0("sim", ndraws, "_k", k, "_n", n, var_tag ,".csv")
  
  df      <- read_csv(file.path(tab_loc, "coefs/", file))
  conv    <- read_csv(file.path(tab_loc, "converge/", file))
  times   <- read_csv(file.path(tab_loc, "times/", file))
  
  df %>% 
    mutate(across(ml:el, ~.x-beta)) %>% 
    pivot_longer(cols = ml:el, values_to = "bias", names_to = "estimator") %>% 
    left_join(pivot_longer(conv, ml:el, names_to = "estimator", 
                           values_to = "converge")) %>% 
    left_join(pivot_longer(times, ml:el, names_to = "estimator", 
                           values_to = "time")) %>% 
    mutate(k = !!k, n = !!n, x_var = !!var_tag)
}


### Make some plots... 
df <- map_dfr(c(5, 10, 12), load_data, n = 1000, X.sigma = "I", tab_loc = tab_loc)
df %>%
  ggplot() +
  geom_density(aes(x = bias, color = var)) +
  xlim(c(-1, 1)) +
  facet_wrap(vars(estimator, k), scales = "free", nrow = 5) + 
  ggtitle("Estimated coefs mins the real values, X.sigma = I")


df <- map_dfr(c(5, 10, 12), load_data, n = 1000, X.sigma = "decay",rho = 0.9, tab_loc = tab_loc)
df %>%
  ggplot() +
  geom_density(aes(x = bias, color = var)) +
  xlim(c(-1, 1)) +
  facet_wrap(vars(estimator, k), scales = "free", nrow = 5) + 
  ggtitle("Estimated coefs mins the real values, X.sigma = decay 0.9")


df <- map_dfr(c(5, 10, 12), load_data, n = 1000, X.sigma = "decay",rho = 0.5, tab_loc = tab_loc)
df %>%
  ggplot() +
  geom_density(aes(x = bias, color = var)) +
  xlim(c(-1, 1)) +
  facet_wrap(vars(estimator, k), scales = "free", nrow = 5) + 
  ggtitle("Estimated coefs mins the real values, X.sigma = decay 0.5")


df %>% 
  group_by(converge, estimator) %>% 
  tally()


# Variance of 
df %>% 
  pivot_longer(cols = ml:el) %>%
  group_by(var, name) %>% 
  summarize(sd = sd(value))

# Run time for a full set of estimators  
df %>% ggplot() + 
  geom_density(aes(x = time))
