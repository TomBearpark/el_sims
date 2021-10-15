###############################################################################
user <- Sys.getenv("USER")
if( user == "tombearpark"){
  dir <- '/Users/tombearpark/Documents/GitHub/metricsfolks/ECO513/psets/ps2/'
} else {
  dir <- ''
} 

fig_loc <- paste0(dir, "fig/")
tab_loc <- paste0(dir, "tab/")


if(!require(pacman)) install.packages('pacman')
pacman::p_load('MASS', 'tidyverse')
theme_set(theme_bw())






library(mvtnorm)

