###############################################################################
## Load stuff
if(!require(pacman)) install.packages('pacman')
pacman::p_load('MASS', 'tidyverse', 'gmm')
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

n <- 5000
k <- 5

data.obj <- gen_data(k = k, n = n)
data.df <- data.obj$df
beta <- data.obj$model.specs$beta

## Maximum likelihood
ml <- est.ML(data.df)

## Method of moments
mom <- est.GMM(data.df, type = "mom")

## Twostep GMM
gmm2 <- est.GMM(data.df, type = "twoStep")

## CUE
cue <- est.GMM(data.df, type = "cue")


