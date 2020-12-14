###################################################################################
### SCRIPT PURPOSE: Select climatic predictor variables for SDM models
# Modified from Angert et al. 2018, American Naturalist
# Author: Amy Angert
# last update:  13 Dec 2020

### PROGRAM FUNCTIONS: 
  # Inspect distribution of each bioclimatic variable
  # Examine correlations among predictor varaibles
  # Estimate variable importance using a simple GLM

###################################################################################

################################################################################# 
### LOAD LIBRARIES AND PREPARE INPUTS

## CLEAR WORKSPACE
rm(list = ls(all.names = TRUE))

## LIBRARIES
library(tidyverse)
library(GGally) # for visualizing correlation matrix

## INPUTS
dat <- read_csv("SDM/data_files/biovars.csv")

# Precipitation variables (bio12-bio19) should be log-transformed
dat <- dat %>% 
  mutate(bio12=log(bio12+0.5),
            bio13=log(bio13+0.5),
            bio14=log(bio14+0.5),
            bio15=log(bio15+0.5),
            bio16=log(bio16+0.5),
            bio17=log(bio17+0.5),
            bio18=log(bio18+0.5),
            bio19=log(bio19+0.5))

################################################################################


################################################################################
### INSPECT DISTRIBUTIONS AND COLLINEARITY

# pull bioclim variables out into their own dataframe
cor.dat <- dat %>% 
  dplyr::select(bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10, 
         bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19)

# scatterplot matrix with univariate distribution of each variable along diagonal
cor.mat <- ggpairs(cor.dat)
cor.mat

################################################################################


################################################################################
### EXAMINE UNIVARIATE IMPORTANCE OF EACH CANDIDATE PREDICTOR
## This function  is borrowed and modified from Niklaus E. Zimmermann
## It evaluates the predictive power of individual predictor variables by running a logistic regression of presence/absence against linear and quadratic terms of each variable
## Metrics of predictive power = adjusted deviance, AIC

## Variable Importance (VIP) function
# Inputs:
  # tr.spp = data frame with pres/abs
  # tr.var = data frame with candidate predictor variables
  # pres = column number with pres/abs
  # pf = column number where candidate predictors start
  # pl = column number where candidate predictors end
varimp.glm = function(tr.spp, tr.var, pres, pf, pl) {
  tmp.mat = matrix(ncol=2, nrow=(pl-pf+1))
  for (i in pf:pl) {
    tmp = glm(tr.spp[,pres] ~ tr.var[,i] + I((tr.var[,i])^2), na.action=na.omit, family=binomial)
    tmp.mat[(i-pf+1),1] = tmp$aic
    tmp.mat[(i-pf+1),2] = (1-(tmp$deviance/tmp$null.deviance))
    }
    return(tmp.mat)
  }

## Set up inputs for VIP function
tr.vip <- dat %>% # keep only Pres/Abs & bioclim predictors
  dplyr::select(presabs, 
         bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10, 
         bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19)
tr.vip <- as.data.frame(tr.vip)
pres=1                     # column number with presence:absence
v.start=2                  # column number where predictor variables start
v.stop=ncol(tr.vip)        # column number where predictor variables end
# v.num=v.stop-1             # number of predictor variables

## RUN VIP function
dev.fit = varimp.glm(tr.vip, tr.vip, pres, v.start, v.stop)

# spruce up labels
x.labs=as.data.frame(names(tr.vip[2:v.stop]))
dev.fit =cbind(as.data.frame(dev.fit), x.labs) # output matrix; col=1 AIC, col=2 Adj deviance
names(dev.fit) = c("AIC", "AdjDev", "Var")

# sort by descending order of explanatory power
dev.fit = dev.fit[order(dev.fit$AIC),]
dev.fit 

# basic barplot of explanatory power
ggplot(dev.fit, aes(x=Var, y=AdjDev)) + geom_bar(stat="identity")

################################################################################


################################################################################
### VARIABLE SELECTION BASED ON IMPORTANCE AND COLLINEARITY

## (NOTE: variable selection is altered from Angert et al. 2018, due to addition of occupancy dataset to training input &/or use of 1961-1990 normals)

# Single best variable = bio15
# bio15 has no |r|>0.7 with other bioclim
# 2nd best variable = bio10
# If bio10 is included, then bio1, bio5, bio8, and bio9 should be excluded (|r|>0.7)
# Next best variable of those remaining = bio14
# If bio14 is included, then bio17 and bio18 should be excluded (|r|>0.7)
# Next best bariable = bio12
# If bio12 is included, then bio13, bio16, and bio19 should be excluded (|r|>0.7)
# Next best variable = bio11
# If bio11 is included, then bio1, bio6, and bio8 should be excluded (|r|>0.7)
# Next best variable = bio4
# If bio4 is included, then bio7 should be excluded (|r|>0.7)
# Next best of remaining variables = bio3
# bio3 has no |r|>0.7 with other bioclim
# Choose last remaining = bio2

## Final list of variables to include:
	#bio15, bio10, bio14, bio12, bio11, bio4, bio3, bio2

## Save file with selected predictors
dat.input <- dat %>% 
  select(set, lat, long, el, presabs, bio15, bio10, bio14, bio12, bio11, bio4, bio3, bio2)
write_csv(dat.input, "SDM/data_files/sdm_input.csv")

################################################################################
