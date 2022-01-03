##################################################################################
## Gradient forest for full SNP dataset
## Author Daniel Anstett
## 
## 
## Last Modified Dec 13, 2021
###################################################################################
##Libraries
library(tidyverse)
library(gradientForest)

##Import Data
snp_clim_bf10 <- read_csv("Genomics_scripts/Data/snp_clim_bay.csv")
snp_clim_bf10NA <- snp_clim_bf10 %>%
  select_if(~ !any(is.na(.)))
#snp_clim_ful <- read_csv("Genomics_scripts/Data/snp_clim_full.csv") # full data


## Generate specific dataframes for GF model
env_site <- snp_clim_bf10 %>% select(MAT,MAP,PAS,EXT,CMD,Tave_wt,Tave_sm,PPT_wt,PPT_sm)
test_snp <- snp_clim_bf10[,15:114] %>%
  select_if(~ !any(is.na(.)))

env_site<-as.data.frame(env_site)
test_snp<-as.matrix(test_snp)
rownames(test_snp)<-as.character(c(1:dim(test_snp)[1]))



# Gradient Forest Model
gf <- gradientForest(cbind(env_site, test_snp),
                      predictor.vars = colnames(env_site), response.vars = colnames(test_snp),
                      ntree = 500, transform = NULL, compact = T,
                      nbin = 201, maxLevel = 5 , corr.threshold = 0.5)

