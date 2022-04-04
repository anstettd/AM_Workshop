##################################################################################
## Climate associated SNP distribution across chromosomes
## Author Daniel Anstett
## 
## 
## Last Modified Marc 31, 2021
###################################################################################

#Import libraries
library(tidyverse)
library(Hmisc)

#Import data
climate <- read_csv("Donor_selection/Data/climate_pop.csv")

climate_corr <- as.matrix(climate[,6:14])
climate_all <- rcorr(climate_corr)
climate_cor <- as.data.frame(climate_all$r)
write_csv(climate_cor, "Genomics_scripts/Data/climate_corr.csv")



