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
library(tidyverse) # for gpplot and data manipulations
library(GGally) # for visualizing correlation matrix

## INPUTS
dat <- read_csv("SDM/data_files/biovars.csv")

################################################################################


################################################################################
### INSPECT DISTRIBUTIONS AND COLLINEARITY

# pull bioclim variables out into their own dataframe
cor.dat <- dat %>% 
  dplyr::select(bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10, 
         bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19)

# scatterplot matrix with univariate distribution of each variable along diagonal
cor.mat <- ggpairs(cor.dat)
cor.mat # view this in conjunction with variable importance, below

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
#ggplot(dev.fit, aes(x=Var, y=AdjDev)) + geom_bar(stat="identity")

################################################################################


################################################################################
### VARIABLE SELECTION BASED ON IMPORTANCE AND COLLINEARITY

## (NOTE: variable selection is altered from Angert et al. 2018, due to addition of occupancy dataset to training input &/or use of 1961-1990 normals)

# Single best variable = bio10. INCLUDE
# bio10 has |r|>0.7 with bio1, bio5, bio8, bio9, bio11
# If bio10 is included, then bio1, bio5, bio8, bio9, and bio11 should be excluded (|r|>0.7)

# 2nd best variable = bio1, but that is collinear with bio10. EXCLUDE

# 3rd best variable = bio5, but that is collinear with bio10. EXCLUDE

# 4th best variable = bio14. INCLUDE
# bio14 has |r|>0.7 with bio1, bio17, bio18
# If bio14 is included, then bio1, bio17 and bio18 should be excluded (|r|>0.7)

# 5th best variable = bio11, but that is collinear with bio10. EXCLUDE

# 6th best variable = bio15. INCLUDE
# bio15 is not too collinear with any other variable.

# 7th best variable = bio9, but that is collinear with bio10. EXCLUDE

# 8th best variable = bio8, but that is collinear with bio10. EXCLUDE

# 9th best variable = bio17, but that is collinear with bio14. EXCLUDE

# 10th best variable = bio18, but that is collinear with bio14. EXCLUDE.

# 11th best variable = bio12. INCLUDE
# bio12 has |r|>0.7 with bio13, bio16, bio19
# If bio12 is included, then bio13, bio16 and bio19 should be excluded (|r|>0.7)

# 12th best variable = bio6. INCLUDE
# bio6 has |r|>0.7 with bio1, bio8, bio11 but these have already been excluded

# 13th best variable = bio16, but this is collinear with bio12. EXCLUDE

# 14th best variable = bio13, but this is collinear with bio12. EXCLUDE

# 15th best variable = bio19, but this is collinear with bio12. EXCLUDE

# 16th best variable = bio3. INCLUDE
# bio3 has |r|>0.7 with bio4
# If bio3 is included, then bio4 should be excluded (|r|>0.7)

# 17th best variable = bio2
# bio2 has |r|>0.7 with bio7
# If bio2 is included, then bio7 should be excluded (|r|>0.7)

# 18th best variable = bio4, but this is collinear with bio3. EXCLUDE

# last variable = bio7, but this is collinear with bio2

## Final list of variables to INCLUDE: bio10, bio14, bio15, bio12, bio6, bio3 bio2

## Final list of variables to EXCLUDE: bio1, bio5, bio8, bio9, bio11, bio17, bio18, bio13, bio16, bio19, bio4, bio7

## Save file with selected predictors
dat.input <- dat %>% 
  dplyr::select(Master.ID, Latitude, Longitude, Elevation, presabs, bio10, bio14, bio15, bio12, bio6, bio3, bio2)
write_csv(dat.input, "SDM/data_files/sdm_input.csv")

################################################################################
