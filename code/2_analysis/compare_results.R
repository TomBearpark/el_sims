###########################################################################
### Set up and paths 

## Load stuff
if(!require(pacman)) install.packages('pacman')
pacman::p_load('tidyverse', 'latex2exp', 'stargazer', 'scales')
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
source(file.path(code_dir, "code/1_run_sim/get_Avar.R"))

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
    mutate(covar = factor(var, levels = paste0("x", 1:!!k)))
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

load_all_k <- function(k_list, n, X.sigma, rho, tab_loc, fix_beta = FALSE){
  map_dfr(k_list, load_data, n = n, X.sigma = X.sigma, rho = rho, tab_loc = tab_loc, 
          fix_beta = fix_beta) 
}

get_opts <- function(k_list, X.sigma_list, rho_list){
  expand_grid(k = k_list, X.sigma = X.sigma_list, rho = rho_list) %>% 
    mutate(rho = ifelse(X.sigma == "I", "0", rho)) %>% distinct()
}

load_all <- function(k_list = c(2,5,10), X.sigma_list = c("I", "diagish", "decay"), 
                       rho_list = c(.5, .9), n_draws = 5000, tab_loc, 
                     fix_beta = FALSE){
  opts <- get_opts(k_list = k_list, X.sigma_list = X.sigma_list, rho_list = rho_list)
  pmap_dfr(opts, load_data, n =1000, tab_loc = tab_loc, ndraws = n_draws, 
           fix_beta = fix_beta)
}

sim_name <- function(X.sigma, rho){
  if(X.sigma == "I") return("I")
  else return(paste0(X.sigma, " ", rho))
}

add_labs <- function(df){
  df$lab2 <- factor(df$x_var, levels = c("","_diagish_rho0_5", "_diagish_rho0_9", 
                                         "_decay_rho0_5", "_decay_rho0_9"), 
                    labels = c(TeX("$I_{k-1}$"),
                               TeX("$\\Sigma_{0.5}$"), TeX("$\\Sigma_{0.9}$"), 
                               TeX("$\\Xi_{0.5}$"), TeX("$\\Xi_{0.9}$")))
  
  df$lab1 <- factor(df$k, levels = c(2, 5, 10), 
                    labels = c(TeX("$k=2$"), TeX("$k=5$"), TeX("$k=10$")))

  df
}

bar_plot <- function(df, plot_var){
  
  df <- add_labs(df)
  nrow <- length(unique(df$lab1))
  ggplot(df) + 
    geom_col(aes(x = estimator, y = .data[[plot_var]]), fill = "#FC4E07", 
             alpha = .8, color = "black", width = 0.5) + 
    facet_wrap(vars(lab1, lab2), scales = "free", labeller = "label_parsed", 
               nrow = 3) + 
    xlab("")
}


add_table_labs <- function(df){
  df$`X distribution` <- case_when(
    df$x_var == "" ~ "I", 
    df$x_var == "_diagish_rho0_5" ~ '$\\Sigma_{0.5}$', 
    df$x_var == "_diagish_rho0_9" ~ '$\\Sigma_{0.9}$', 
    df$x_var == "_decay_rho0_5" ~   '$\\Xi_{0.5}$', 
    df$x_var == "_decay_rho0_9" ~   '$\\Xi_{0.9}$'
  )
  df
}



###########################################################################
### 2.1 Error density plots

X.sigma_list <- c("I", "diagish", "decay")
k_list       <- c( 2,  5,  10)
rho_list     <- c(.5, .9)

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

# Add a tag for if somethign didn't converge 
df <- df %>% group_by(i) %>% mutate(fail = sum(converge))  %>% ungroup()


cc <- comp_vars[1]
r <- df %>% rmse_tab()
r %>% filter(estimator != "cue") %>% bar_plot("rmse")


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
  
  df %>% filter(fail == 0) %>% rmse_tab()   %>% 
    bar_plot(plot_var = cc) %>% 
    ggsave(file = paste0(dir_outcc, "converged.png"), height = 7, width = 12)
  
  df %>% filter(fail == 0, estimator != "cue") %>% rmse_tab()   %>% 
    bar_plot(plot_var = cc) %>% 
    ggsave(file = paste0(dir_outcc, "no_cue_converged.png"), height = 7, width = 12)
}

df %>% 
  rmse_tab() %>%
  write_csv(file = paste0(dir_out2, "raw.csv"))


## Output tex tables
for(k in k_list){
  df %>% 
    rmse_tab() %>% add_table_labs() %>% 
    select(-x_var) %>% 
    relocate(estimator, `X distribution`) %>% 
    filter(k == !!k) %>% 
    mutate(across(rmse:mae, round, 5)) %>% 
    data.frame() %>%
    stargazer(summary = FALSE, rownames = FALSE) %>% 
    capture.output() %>% 
    writeLines(paste0(tab_out, "tab_k", k,".tex"))
}

## compute times
df %>% 
  group_by(estimator) %>% summarize(mean_time = mean(time), 
                                    median_time = median(time),
                                    max_time = max(time)
                                    ) %>% 
  data.frame() %>% 
  stargazer(summary = FALSE, rownames = FALSE)

df %>% 
  ggplot() + 
  geom_density(aes(x = log(time), color = estimator), bw = 1) 

## Lets find a bad draw, check it out a bit... 
df %>% 
  filter(converge == 0) %>% 
  filter(abs(error) > 100000)  %>% 
  write_csv(file = paste0(tab_loc, "terrible_cues.csv"))

###########################################################################
## 2.3 Fixed beta results 


df2 <- load_all(tab_loc = tab_loc, fix_beta = TRUE, k_list = c(2, 5))
dir_out3 <- file.path(fig_out, "bar_comparison/fix_beta/")
dir.create(dir_out3, showWarnings = FALSE)

cc <- "rmse"

dir_outcc <- file.path(dir_out3, cc, "/")
dir.create(dir_outcc, showWarnings = FALSE)

df2 %>% rmse_tab() %>% 
  bar_plot(plot_var = cc) %>% 
  ggsave(file = paste0(dir_outcc, "raw.png"), height = 7, width = 12)

df2 %>% rmse_tab() %>% filter(estimator != "cue") %>% 
  bar_plot(plot_var = cc) %>% 
  ggsave(file = paste0(dir_outcc, "no_cue.png"), height = 7, width = 12)


## Get the asymtotic variances 

gen_ratio_df <- function(X.sigma, rho, k = 5, filter_cue = TRUE){

  message(paste0(X.sigma, "------", rho))
  nsim    <- 100000
  data    <- gen_data(n = nsim, k = k, fix_beta = TRUE, 
                      X.sigma = X.sigma, rho = rho)
  Xdraw   <- as.matrix(data$df[,-1])
  Xtilde  <- cross_mom(Xdraw)
  # just loop over betas the following mf
  aux <- get.aVar(Xdraw,Xtilde,data$model.specs$beta)
  
  message("------formatting data-----------")
  
  avars <- aux %>% 
    as_tibble() %>% 
    mutate(var = paste0("x", 1:k)) %>% 
    mutate(aVar.cue = aVar.gmm, aVar.el = aVar.gmm) %>% 
    pivot_longer(cols = -var, values_to = 'avar') %>% 
    mutate(estimator  = substr(name, 6, str_length(name))) %>% 
    mutate(estimator  = ifelse(estimator == "mle", "ml", estimator)) %>% 
    select(-name)
    
  TT <- 1000
  
  plot_df <- load_data(k = k, n = 1000, X.sigma = 
                          X.sigma, rho = rho, tab_loc = tab_loc, 
                        fix_beta = TRUE) %>% 
      filter(k == !!k) %>% 
      group_by(k, estimator, var) %>% 
      summarize(var_mc = sd(sqrt(TT) * error)^2) %>% 
      left_join(avars) %>% 
      mutate(ratio = var_mc / avar)
  
  if(filter_cue)  {
    plot_df <- plot_df %>% filter(estimator != "cue")
    levs <- c("ml", "mom", "gmm", "el")
  }else{
    levs <- c("ml", "mom", "gmm", "cue", "el")
  }

  plot_df %>% 
    mutate(x_var = get_var_tag(X.sigma = X.sigma, rho = rho)) %>% 
    add_table_labs()
    
}

main <- main %>% mutate(rho = as.numeric(rho))

plot_df <- 
  pmap_dfr(main, gen_ratio_df, k = 2, filter_cue = FALSE) %>% 
  bind_rows(
    pmap_dfr(main, gen_ratio_df, k = 5, filter_cue = FALSE)
  )


pp <- plot_df %>% 
  filter(estimator != "cue", var == "x2") %>% 
  mutate(estimator = factor(estimator, levels = c("ml", "mom", "gmm", "el"))) %>% 
  mutate(ratio  = (ratio - 1)*100) %>% 
  mutate(lab = factor( x_var, levels = c("","_diagish_rho0_5", "_diagish_rho0_9", 
                                     "_decay_rho0_5", "_decay_rho0_9"), 
                labels = c(TeX("$I_{k-1}$"),
                           TeX("$\\Sigma_{0.5}$"), TeX("$\\Sigma_{0.9}$"), 
                           TeX("$\\Xi_{0.5}$"), TeX("$\\Xi_{0.9}$")))) %>% 
  mutate(k = factor(k, levels = c(2, 5), 
         labels = c(TeX("$k=2$"), TeX("$k=5$")))) %>% 
  ggplot() + 
  geom_col(aes(x = estimator, y = ratio, fill = "#FC4E07"), 
                color = "black", position = "dodge", alpha = .8, width = 0.5)+ 
  theme(legend.position = "none") + 
  scale_y_continuous(labels = scales::percent_format(scale = 1, accuracy = 1L), 
                     limits = c(0, NA)) + 
  coord_flip() + xlab("") + ylab("") + 
  facet_wrap(vars(k, lab), scales="free", nrow = 2, labeller = "label_parsed")

pp
avar_path <- file.path(fig_out, "avar_comparison.png")
ggsave(pp, file = avar_path, height = 5, width = 12)







## Scraps...
# gen_ratio_plot <- function(X.sigma, rho, k = 5, filter_cue = TRUE){
#   
#   nsim    <- 100000
#   data    <- gen_data(n = nsim, k = k, fix_beta = TRUE, 
#                       X.sigma = X.sigma, rho = rho)
#   Xdraw   <- as.matrix(data$df[,-1])
#   Xtilde  <- cross_mom(Xdraw)
#   # just loop over betas the following mf
#   aux <- get.aVar(Xdraw,Xtilde,data$model.specs$beta)
#   
#   avars <- aux %>% 
#     as_tibble() %>% 
#     mutate(var = paste0("x", 1:k)) %>% 
#     mutate(aVar.cue = aVar.gmm, aVar.el = aVar.gmm) %>% 
#     pivot_longer(cols = -var, values_to = 'avar') %>% 
#     mutate(estimator  = substr(name, 6, str_length(name))) %>% 
#     mutate(estimator  = ifelse(estimator == "mle", "ml", estimator)) %>% 
#     select(-name)
#   
#   TT <- 1000
#   
#   plot_df <- load_data(k = k, n = 1000, X.sigma = 
#                          X.sigma, rho = rho, tab_loc = tab_loc, 
#                        fix_beta = TRUE) %>% 
#     filter(k == !!k) %>% 
#     group_by(k, estimator, var) %>% 
#     summarize(var_mc = sd(sqrt(TT) * error)^2) %>% 
#     left_join(avars) %>% 
#     mutate(ratio = var_mc / avar)
#   
#   if(filter_cue)  {
#     plot_df <- plot_df %>% filter(estimator != "cue")
#     levs <- c("ml", "mom", "gmm", "el")
#   }else{
#     levs <- c("ml", "mom", "gmm", "cue", "el")
#   }
#   
#   plot_df %>% 
#     mutate(estimator = factor(estimator, levels = levs)) %>% 
#     mutate(ratio  = (ratio - 1)*100) %>% 
#     ggplot() + 
#     geom_col(aes(x = estimator, y = ratio, fill = var), 
#              color = "black", position = "dodge", alpha = .7) + 
#     ylab("Ratio of aVar and estimated var")
# }
# 
# gen_ratio_plot