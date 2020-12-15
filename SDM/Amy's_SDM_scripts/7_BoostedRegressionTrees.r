###################################################################################
### SCRIPT PURPOSE: SDM based on random forest models of presence/absence ~ bioclimatic predictors
# Modified from Angert et al. 2018, American Naturalist
# Author: Amy Angert
# last update:  14 Dec 2020

## OVERALL WORKFLOW:
#	Build boosted regression tree model of presence/absence ~ bioclimatic predictor variables
#	Perform resubstitution and cross-validation to assess accuracy

###################################################################################
### LOAD LIBRARIES AND PREPARE INPUTS

## CLEAR WORKSPACE
rm(list = ls(all.names = TRUE))

## LIBRARIES
library(gbm) # for BRT models
library(dismo) # for gbm.step function (per Elith et al. 2008, J Animal Ecol 77:802-813)
library(PresenceAbsence) # for accuracy stats

## INPUTS
# read in file
dat <- read_csv('SDM/data_files/sdm_input.csv')

# slim dataframe to conform to structure required by mod.form function
# (see script modforms.R)
dat.input <- dat %>% 
  dplyr::select(presabs, bio10, bio15, bio11, bio17, bio12, bio3, bio2)
dat.input <- as.data.frame(dat.input)

################################################################################


################################################################################
### BOOSTED REGRESSION TREE MODELS

## Notes about model optimization 
# Goal 1 = lowest deviance (i.e., low residual variance)
# Goal 2 = fewest trees (i.e., parsimony)
# These two goals might conflict, in that there is often trade-off between a less parsimonious model (high n trees) with low deviance (more variance explained) and a more parsimonious model (few trees) with higher deviance (less variance explained)
# Below we will test a range of tree complexities and learning rates to "optimize" (really, tinker) with these outcomes

# Set up variables for gbm function
resp <- paste("as.factor(",colnames(dat.input[1]),")",sep="") # assign response to column number
n.col <- ncol(dat.input) # number of columns
preds <- 2:n.col # assign predictors to column numbers
	
# Basic BRT model with default parameter settings
mod1.BRT <- gbm.step(data=dat.input, gbm.x=preds, gbm.y=1, family="bernoulli", tree.complexity=2, learning.rate=0.01, bag.fraction=0.75, n.folds=5, n.trees=50, plot.main=TRUE, keep.fold.fit=TRUE, step.size=15) 
#save(mod1.BRT, file="SDM/Output/BRT.mod1.Rda")
 		
# Learning rate adjusted up 
mod2.BRT <- gbm.step(data=dat.input, gbm.x=2:8, gbm.y=1, family="bernoulli", tree.complexity=2, learning.rate=0.001, bag.fraction=0.75, n.folds=5, plot.main=TRUE, keep.fold.fit=TRUE, step.size=15)
#save(mod2.BRT, file="SDM/Output/BRT.mod2.Rda")

## Learning rate adjusted down
mod3.BRT <- gbm.step(data=dat.input, gbm.x=preds, gbm.y=1, family="bernoulli", tree.complexity=2, learning.rate=0.1, bag.fraction=0.75, n.folds=5, plot.main=TRUE, keep.fold.fit=TRUE, step.size=15)
save(mod3.BRT, file="SDM/Output/BRT.mod3.Rda")

## Tree complexity adjusted up (this model fails)
mod4.BRT <- gbm.step(data=dat.input, gbm.x=preds, gbm.y=1, family="bernoulli", tree.complexity=3, learning.rate=0.1, bag.fraction=0.75, n.folds=5, plot.main=TRUE, keep.fold.fit=TRUE)#, step.size=5) 
#save(mod4.BRT, file="SDM/Output/BRT.mod4.Rda")

## Tree complexity and learning rate both adjusted up
mod5.BRT <- gbm.step(data=dat.input, gbm.x=preds, gbm.y=1, family="bernoulli", tree.complexity=3, learning.rate=0.01, bag.fraction=0.75, n.folds=5, plot.main=TRUE, keep.fold.fit=TRUE, step.size=5)
#save(mod5.BRT, file="SDM/Output/BRT.mod5.Rda")
	
## Learning rate and tree complexity both adjusted down
mod6.BRT <- gbm.step(data=dat.input, gbm.x=preds, gbm.y=1, family="bernoulli", tree.complexity=1, learning.rate=0.1, bag.fraction=0.75, n.folds=5, plot.main=TRUE, keep.fold.fit=TRUE, step.size=5)
#save(mod6.BRT, file="SDM/Output/BRT.mod6.Rda")

## examine BRT model performance
# this pulls n trees and deviance from each model's output into a data frame
perf <- matrix(ncol=2, nrow=6)
for (j in 1:6) {
	mod = get(paste("mod",j,".BRT", sep=""))
	if((is.numeric(mod$cv.statistics[[1]]))==T){
	  perf[j,1] = round(mod$cv.statistics[[1]],3)
	  perf[j,2] = mod$n.trees
				} else {
				  perf[j,1] = NA # failed model 
				}
  }
perf <- as.data.frame(perf)
row.names(perf) <- c("mod1", "mod2","mod3", "mod4", "mod5", "mod6")
colnames(perf) <- c("dev", "ntrees")
perf

ggplot(perf, aes(x=ntrees, y=dev)) + geom_point()
# Model 3 is the winner - it has the lowest deviance and the second fewest trees

################################################################################


################################################################################
### LOAD FINAL MODEL AND PREDICTIONS 

mod = get(load(paste("SDM/Output/BRT.mod3.Rda", sep="")))
pred = mod$fitted

################################################################################
### ACCURACY CALCULATIONS, MODEL=BRT

source("SDM/Amy's_SDM_scripts/accuracy.R")

modl="BRT.mod3" # required label for accuracy function              

## be sure to adjust model #
acc <- accuracy.brt(dat.input, pred, modl)
acc$thresh = c("SensSpec", "Kappa")
acc$model = "BRT.mod3"
save(acc, file="SDM/Output/BRT.mod3.accs.Rda")

################################################################################
