##################################################################################
## Organizing BayPass Results
## Author Daniel Anstett
## 
##
## Last Modified June 20, 2021
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)
setwd("/Users/daniel_anstett/Dropbox/AM_Workshop/12_baypass_results")
baye_1 <- read.table("cov_model_rand_1_summary_betai_reg.out",sep="\t")
