###################################################################################
### SCRIPT PURPOSE: SDM based on random forest models of presence/absence ~ bioclimatic predictors
# Modified from Angert et al. 2018, American Naturalist
# Author: Amy Angert
# last update:  14 Dec 2020

## OVERALL WORKFLOW:
#	Build RF model
#	Perform resubstitution and cross-validation to assess accuracy

###################################################################################
### LOAD LIBRARIES AND PREPARE INPUTS

## CLEAR WORKSPACE
rm(list = ls(all.names = TRUE))

## LIBRARIES
library(randomForest) # for RF models
library(PresenceAbsence) # for accuracy stats
library(DAAG) # for cross-validation sampling

## INPUTS
# read in file
dat <- read_csv('SDM/data_files/sdm_input.csv')

# slim dataframe to conform to structure required by mod.form function
# (see script modforms.R)
dat.input <- dat %>% 
  select(presabs, bio10, bio14, bio15, bio12, bio6, bio3, bio2)

################################################################################


################################################################################
### RANDOM FORESTS MODEL

source("SDM/Amy's_SDM_scripts/modforms.R")

mod1.RF <- randomForest(mod.form.lin(dat.input, 1, 2), importance=T, keep.forest=T, data=dat.input)

save(mod1.RF, file="RF.mod1.Rda")	

################################################################################


################################################################################
### {If necessary} LOAD SAVED MODEL AND CALCULATE PREDICTIONS

mod = get(load("RF.mod1.Rda"))
pred = predict(mod1.RF, type="prob")[,2]

################################################################################


################################################################################
### ACCURACY CALCULATIONS, MODEL=RF

source("SDM/Amy's_SDM_scripts/accuracy.R")

modl = "mod1.RF" # add var to keep track of model, needed for accuracy function below

# due to nature of model, the only accuracy calculations are inherently cross-validated
acc <- accuracy(dat.input, pred, modl)
acc$thresh = c("SensSpec", "Kappa")
acc$model = "RF.mod1"
save(acc, file="SDM/Output/RF.mod1.accs.Rda")

################################################################################
