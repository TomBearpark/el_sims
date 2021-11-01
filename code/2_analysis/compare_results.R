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
df <- read_csv(file.path(tab_loc, paste0("sim_k", k, "_n", n, ".csv")))

### Make some plots... 
df %>%
  mutate(across(ml:el, ~.x-beta)) %>%
  pivot_longer(cols = ml:el) %>%
  ggplot() +
  geom_density(aes(x = value, color = var)) +
  facet_wrap(vars(name)) + 
  ggtitle("Estimated coefs mins the real values")

df %>% 
  pivot_longer(cols = ml:el) %>%
  group_by(var, name) %>% 
  summarize(sd = sd(value))
  
df %>% ggplot() + 
  geom_density(aes(x = time))