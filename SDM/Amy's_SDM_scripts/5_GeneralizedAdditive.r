###################################################################################
### SCRIPT PURPOSE: SDM based on generalized additive model of presence/absence ~ bioclimatic predictors
# Modified from Angert et al. 2018, American Naturalist
# Author: Amy Angert
# last update:  14 Dec 2020

## OVERALL WORKFLOW:
#	Build full and reduced model logistic GLMs
#	Perform resubstitution and cross-validation to assess accuracy

###################################################################################

### LOAD LIBRARIES AND PREPARE INPUTS

## LIBRARIES
library(gam) # for GAM models
library(PresenceAbsence) # for accuracy stats
library(DAAG) # for cross-validation sampling

## INPUTS
# read in file
dat <- read_csv('SDM/data_files/sdm_input.csv')

# slim dataframe to conform to structure required by mod.form function
# (see script modforms.R)
dat.input <- dat %>% 
  select(presabs, bio15, bio10, bio14, bio12, bio11, bio4, bio3, bio2)

################################################################################


################################################################################
### BUILD FULL AND REDUCED GAM MODELS 

## Source code file with function for generating lists of model terms
source("SDM/R_Code/modforms.R")

## GAM model 1: all smoothers using default df (smoothers=4)
mod1.GAM <- gam(mod.form.4(dat.input, 1, 2), family=binomial, data=dat.input)
#mod1.fit <- 100*(1-mod1.GAM$deviance/mod1.GAM$null.deviance) 

## GAM model 2: stepwise with different scopes
# scope list gives range of possible smoothing degrees for each variable
# Daniel: Changed "step.gam" to "step.Gam" but Amy changed back
# Amy changing to mod1 as basis and deleting old basis mod2 (smoothers=3)
mod2.GAM <- step.gam(mod1.GAM, scope=list( 
		"bio2"	=~1+ bio2	+s(bio2,2) +s(bio2,3)	+s(bio2,4),
		"bio3"	=~1+ bio3 	+s(bio3,2) 	+s(bio3,3) 	+s(bio3,4),
		"bio4"	=~1+ bio4 	+s(bio4,2) 	+s(bio4,3) 	+s(bio4,4),
		"bio10"	=~1+ bio10 	+s(bio10,2) +s(bio10,3) +s(bio10,4),
		"bio11"	=~1+ bio11 	+s(bio11,2) +s(bio11,3) +s(bio11,4),
		"bio12"	=~1+ bio12 	+s(bio12,2) +s(bio12,3)	+s(bio12,4),
		"bio14"	=~1+ bio14 	+s(bio14,2) +s(bio14,3)	+s(bio14,4),
		"bio15"	=~1+ bio15 	+s(bio15,2) +s(bio15,3)	+s(bio15,4)),
		trace=F)
#mod4.fit <- 100*(1-mod4.GAM$deviance/mod4.GAM$null.deviance)  
save(mod2.GAM, file="SDM/Output/GAM.mod2.Rda")

################################################################################


################################################################################
######## {If necessary} LOAD SAVED MODEL OBJECTS AND PREDICTIONS

mod = get(load("SDM/Output/GAM.mod2.Rda"))
pred = predict(mod, type="response")

################################################################################


################################################################################
### ACCURACY CALCULATIONS

source("SDM/R_code/accuracy.R")

modl="mod2.GAM" # label to keep track of which model

## METHOD 1: resusbtitution
## i.e., use full input dataframe as testing dataframes for accuracy calculation 
accs <- accuracy(dat.input, pred, modl)
accs$thresh = c("SensSpec", "Kappa")
accs$model = "GAM.mod2"
save(accs, file="SDM/Output/GAM.mod2.accs.Rda")

## METHOD 2: k-fold cross-validation
## i.e., hold out 1 of k folds as testing dataframe for accuracy calculation
## (results will vary slightly if re-run due to random division of data set into folds)
cv.accs <- cv.accuracy(mod, dat.input, modl)
cv.accs$thresh = c("SensSpec", "Kappa")
cv.accs$model = "GAM.mod2"
save(cv.accs, file="SDM/Output/GAM.mod2.cvaccs.Rda")

################################################################################
