##################################################################################
## Gradient forest for full SNP dataset
## Author Daniel Anstett
## 
## 
## Last Modified Dec 13, 2021
###################################################################################

# Clear environment
rm(list = ls())

# Get this package retrieving function
## This function will automatically load packages that you already have
## and will install packages you don't yet have then load them
ipak <- function(pkg){
  # Function written by Dr. Evan Fricke
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = T)
  sapply(pkg, require, character.only = T)
}

# Define the packages that the script needs
myPackages <- c("randomForest", "extendedForest", "gradientForest")

# Load the packages
ipak(myPackages)

# Install packages
install.packages("extendedForest", repos="http://R-Forge.R-project.org")
install.packages("gradientForest", repos="http://R-Forge.R-project.org")

###################################################################################
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
test_snp<-as.data.frame(test_snp)


df_in_1<-data.frame(env_site, test_snp)
pred<-colnames(env_site)
resp<-colnames(test_snp)

#resp<-as.factor(resp)
#resp<-droplevels(resp)

# Gradient Forest Model
gf <- gradientForest(df_in_1,
                      predictor.vars = pred, response.vars = resp,
                      ntree = 500, transform = NULL, compact = T,
                      nbin = 201 , corr.threshold = 0.5)


