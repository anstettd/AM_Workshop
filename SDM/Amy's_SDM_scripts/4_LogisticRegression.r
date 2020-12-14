################################################################################# 
### SCRIPT PURPOSE: SDM based on generalized linear model (i.e., logistic regression of presence/absence ~ bioclimatic predictors)
# Modified from Angert et al. 2018, American Naturalist
# Author: Amy Angert
# last update:  13 Dec 2020

## OVERALL WORKFLOW:
  #	Build full and reduced model logistic GLMs
  #	Perform resubstitution and cross-validation to assess accuracy

###################################################################################


################################################################################# 
### LOAD LIBRARIES AND PREPARE INPUTS

## LIBRARIES
library(PresenceAbsence) # for accuracy stats
library(DAAG) # for cross-validation sampling

## INPUTS
# read in file
dat <- read_csv('SDM/data_files/sdm_input.csv')

# if necessary, filter rows to appropriate ratio of pseudoabsences:presences 
# (not necessary here; all rows are used in GLM)

# slim dataframe to conform to structure required by mod.form function
# (see script modforms.R)
dat.input <- dat %>% 
  select(presabs, bio15, bio10, bio14, bio12, bio11, bio4, bio3, bio2)

################################################################################
### BUILD FULL AND REDUCED GLM MODELS 

## Source code file with function for generating lists of model terms
source("SDM/Amy's_SDM_scripts/modforms.R")

## Build initial model with all variables: mod1.LR
mod1.GLM <- glm(mod.form.quad(dat.input, 1, 2), family=binomial, data=dat.input)
#mod1.fit=100*(1-mod1.GLM$deviance/mod1.GLM$null.deviance)
mod1.pred <- predict(mod1.GLM, type="response") # model prediction
#summary(mod1.GLM) # full model summary stats
save(mod1.GLM, file="SDM/Output/GLM.mod1.Rda") # save model object

## Build most parsimonious model using backwards variable reduction: mod2.LR
mod2.GLM <- step(mod1.GLM, trace=F) # backwards stepwise variable reduction
#mod2.fit=100*(1-mod2.GLM$deviance/mod2.GLM$null.deviance) # model fit
mod2.pred <- predict(mod2.GLM, type="response") # model prediction
#summary(mod2.GLM)  # reduced model summary	
save(mod2.GLM, file="SDM/Output/GLM.mod2.Rda") # save model object

################################################################################


################################################################################
### {IF NECESSARY} LOAD FINAL REDUCED MODELS AND THEIR PREDICTIONS 

mod <- get(load("SDM/Output/GLM.mod2.Rda"))
pred <- predict(mod, type="response") # model prediction

################################################################################


################################################################################
### ACCURACY CALCULATIONS

source("SDM/R_code/accuracy.R")

modl="mod2.GLM"  # label for the model, tracking purposes only

## METHOD 1: resusbtitution
## i.e., use full input dataframe as testing dataframes for accuracy calculation 
acc <- accuracy(dat.input, pred, modl)
acc$thresh = c("SensSpec", "Kappa")
acc$model = "GLM.mod2"
save(accs, file="SDM/Output/GLM.mod2.accs.Rda")

## METHOD 2: k-fold cross-validation
## i.e., hold out 1 of k folds as testing dataframe for accuracy calculation
## (results will vary slightly if re-run due to random division of data set into folds)

cv.accs <- cv.accuracy(mod, dat, modl)
cv.accs = c("SensSpec", "Kappa")
cv.accs$model = "GLM.mod2"
save(cv.accs, file="SDM/Output/GLM.mod2.cvaccs.Rda")

################################################################################
