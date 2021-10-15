###############################################################################
## Load stuff
if(!require(pacman)) install.packages('pacman')
pacman::p_load('MASS', 'tidyverse', 'gmm')
theme_set(theme_bw())

# Set directories
user <- Sys.getenv("USER")
if( user == "tombearpark"){
  dir <- '/Users/tombearpark/Documents/GitHub/metricsfolks/ECO513/psets/ps2/'
} else {
  dir <- 'D:\\Dropbox\\Università\\PhD\\II Year\\Fall\\ECO519 - Non-linear Econometrics\\psets\\el_sims\\'
} 

fig_loc <- paste0(dir, "fig/")
tab_loc <- paste0(dir, "tab/")

# Load auxiliary funs
source(paste0(dir,"funs.R"))

###############################################################################
set.seed(8894)
n <- 10000
k <- 10
data.obj <- gen_data(k = k, n = n)


est.ML(data.obj$df)
data.obj$model.specs$beta

