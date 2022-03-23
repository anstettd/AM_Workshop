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
myPackages <- c("randomForest", "extendedForest", "gradientForest","tidyverse")

# Load the packages
ipak(myPackages)

# Install packages
#install.packages("extendedForest", repos="http://R-Forge.R-project.org")
#install.packages("gradientForest", repos="http://R-Forge.R-project.org")

#library(randomForest)
#library(gradientForest)
#library(tidyverse)

###################################################################################
##Import Data
snp_clim_bf10 <- read_csv("Genomics_scripts/Data/snp_clim_bay.csv")
#snp_clim_bf10NA <- snp_clim_bf10 %>%
#  select_if(~ !any(is.na(.)))
####snp_clim_ful <- read_csv("Genomics_scripts/Data/snp_clim_full.csv") # full data


## Generate specific dataframes for GF model
env_site <- snp_clim_bf10 %>% select(MAT,MAP,PAS,EXT,CMD,Tave_wt,Tave_sm,PPT_wt,PPT_sm)
#test_snp <- snp_clim_bf10[,15:114] %>%
#  select_if(~ !any(is.na(.)))
#write_csv(test_snp,"Genomics_scripts/Data/test_snp.csv")

#Import filtered test snp datatset
test_snp <- read.csv("Genomics_scripts/Data/test_snp.csv")


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






##################################################################################################################
#Basic Plots

#Importance Plot
plot(gf, plot.type = "O")

most_important <- names(importance(gf))[1:9]

##Split Density Plot
#The second plot is the splits density plot (plot.type="S"), which shows binned 
#split importance and location on each gradient (spikes), kernel density of splits 
#(black lines), of observations(red lines) and of splits standardised by 
#observations density (blue lines). Each distribution integrates to predictor 
#importance. These show where important changes in the abundance of
#multiple species are occurring along the gradient; they indicate a composition change rate
plot(gf, plot.type = "S", imp.vars = most_important,leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
     cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))

##Cumulaive Plot
#for each species shows cumulative importance distributions of splits improvement scaled by
#R2 weighted importance, and standardised by density of observations. These show cumulative
#change in abundance of individual species, where changes occur on the gradient, and the species
#changing most on each gradient.
plot(gf, plot.type = "C", imp.vars = most_important,show.overall = F, legend = T, leg.posn = "topleft",
       leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4,cex.axis = 0.6, line.ylab = 0.9, 
       par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0,0.3, 0, 0)))


plot(gf, plot.type = "P", show.names = F, horizontal = F, cex.axis = 1, cex.labels = 0.7, line = 2.5)
















