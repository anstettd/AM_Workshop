###################################################################################
### SCRIPT PURPOSE: SDM based on random forest models of presence/absence ~ bioclimatic predictors
# Modified from Angert et al. 2018, American Naturalist
# Author: Amy Angert
# last update:  14 Dec 2020

## OVERALL WORKFLOW:
#	Build Random Forests model of presence/absence ~ bioclimatic predictor variables
#	Perform resubstitution and cross-validation to tune parameters and assess accuracy

###################################################################################
### LOAD LIBRARIES AND PREPARE INPUTS

## CLEAR WORKSPACE
rm(list = ls(all.names = TRUE))

## LIBRARIES
library(randomForest) # for RF models
library(PresenceAbsence) # for accuracy stats

## INPUTS
# read in file
dat <- read_csv('SDM/data_files/sdm_input.csv')

# slim dataframe to conform to structure required by mod.form function
# (see script modforms.R)
dat.input <- dat %>% 
  dplyr::select(presabs, bio10, bio15, bio11, bio17, bio12, bio3, bio2)

################################################################################


################################################################################
### TEST RANGE OF PARAMETER VALUES FOR OPTIMIZING RF MODEL

source("SDM/Amy's_SDM_scripts/modforms.R")

# Number of trees
trees <- c("50", "500","1000", "1500") # trees 

# Portion of data split for cross validation
cv_folds <- c(10, 5, 3) 

# Number of variables randomly sampled as candidates at each split
my_mtry <- c(2,3,4) 

for(l in 1:length(my_mtry)){
  this_mtry <- my_mtry[l]
  for(m in 1:length(cv_folds)){
    this_cv_fold <- cv_folds[m]
    for(j in 1:length(trees)){
      this_tree <- trees[j]; this_tree <- as.numeric(this_tree)	
        # run model
        mod1.RF <- randomForest(mod.form.lin(dat.input, 1, 2), mtry=this_mtry, 
                               importance=T, keep.forest=T, data=dat.input, 
                               ntree=this_tree, cv.fold=this_cv_fold)           
        # predictions from model
        mod1.pred = predict(mod1.RF, type="prob")[,2] 
        # assign informative names 			
        assign(paste("RF_cv", this_cv_fold,
                     "_ntree",this_tree,
                     "_mtry", this_mtry, sep=""), mod1.RF)
        assign(paste("RF_cv",this_cv_fold,
                     "_ntree",this_tree,
                     "_mtry", this_mtry,".pred", sep=""), mod1.pred)
      }
    }
  }

source("SDM/Amy's_SDM_scripts/accuracy.R")
accs = c()
for(l in 1:length(my_mtry)){
  this_mtry <- my_mtry[l]
  for(m in 1:length(cv_folds)){
    this_cv_fold = cv_folds[m]
    for(j in 1:length(trees)){
      this_tree = trees[j]; this_tree <- as.numeric(this_tree)	
        pred = get(paste("RF_cv",this_cv_fold,
                         "_ntree",this_tree,
                         "_mtry", this_mtry, ".pred", sep=""))
        # add var to keep track of model
        modl = paste("RF_cv",this_cv_fold,
                     "_ntree",this_tree,"_mtry", this_mtry, sep="")      
        temp = accuracy(dat.input, pred, modl)	
        temp$thresh = c("SensSpec")
        temp$model = modl
        temp$mtry <- my_mtry[l]
        temp$cvfold <- cv_folds[m]
        temp$ntree <- trees[j]
        accs = rbind(accs, temp)
      }
    }
  }
accs[order(accs$AUC), ]

################################################################################
### RANDOM FORESTS MODEL

source("SDM/Amy's_SDM_scripts/modforms.R")

# default parameter settings
mod1.RF <- randomForest(mod.form.lin(dat.input, 1, 2), importance=T, keep.forest=T, data=dat.input)
save(mod1.RF, file="RF.mod1.Rda")	

# tuned as per above
mod2.RF <- randomForest(mod.form.lin(dat.input, 1, 2), mtry=3, ntree=1500, cv.fold=3, importance=T, keep.forest=T, data=dat.input)
save(mod2.RF, file="RF.mod2.Rda")	

################################################################################


################################################################################
### {If necessary} LOAD SAVED MODEL AND CALCULATE PREDICTIONS

mod = get(load("RF.mod2.Rda"))
pred = predict(mod2.RF, type="prob")[,2]

################################################################################


################################################################################
### ACCURACY CALCULATIONS, MODEL=RF

source("SDM/Amy's_SDM_scripts/accuracy.R")

modl = "mod2.RF" # add var to keep track of model, needed for accuracy function below

# due to nature of model, the accuracy calculations are inherently cross-validated
acc <- accuracy(dat.input, pred, modl)
acc$thresh = c("SensSpec", "Kappa")
acc$model = "RF.mod2"
save(acc, file="SDM/Output/RF.mod2.accs.Rda")

################################################################################
