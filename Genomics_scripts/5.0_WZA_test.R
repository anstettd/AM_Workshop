##################################################################################
## Lostruct Test
## Author Daniel Anstett
## 
## 
## Last Modified August 2, 2022
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)
library(Kendall)

#Import Climate Data
climate <- read_csv("Genomics_scripts/Data/env_baseline.csv",col_names=FALSE)
climate2 <- filter(climate[1:2,])
climate5 <- filter(climate[5,])
climate_wza <- rbind(climate2,climate5)

#Import genotype data
all_data <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/WZA_vignette/1_0.5_192.alleleFreqs.csv", header = F)