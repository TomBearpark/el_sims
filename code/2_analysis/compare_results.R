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

tab_loc <- file.path(code_dir, "tab/")

### Load the data 
k <- 10
n <- 1000
# df <- read_csv(file.path(tab_loc, paste0("sim_k", k, "_n", n, ".csv")))

df <- read_csv(file.path(tab_loc, "sim1000_k10_n1000.csv")) %>% 
  mutate(rho = 0)
df1 <- read_csv(file.path(tab_loc, "sim1000_k10_n1000_decay_rho0_5.csv")) %>% 
  mutate(rho = 0.5)

### Make some plots... 
df1 %>%
  mutate(across(ml:el, ~.x-beta)) %>%
  pivot_longer(cols = ml:el) %>%
  ggplot() +
  geom_density(aes(x = value, color = var)) +
  xlim(c(-1, 1)) +
  facet_wrap(vars(name), scales = "free") + 
  ggtitle("Estimated coefs mins the real values")

# Variance of 
df %>% 
  pivot_longer(cols = ml:el) %>%
  group_by(var, name) %>% 
  summarize(sd = sd(value))

# Run time for a full set of estimators  
df %>% ggplot() + 
  geom_density(aes(x = time))
