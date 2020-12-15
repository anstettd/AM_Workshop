#################################################################################
### SCRIPT PURPOSE: create and evaluate a maximum entropy ("MaxEnt") species distribution model 
# Modified from Angert et al. 2018, American Naturalist
# Author: Amy Angert
# last update:  15 Dec 2020

## OVERALL WORKFLOW:
#	Build Maxent model of presence/(pseudo)absence ~ bioclimatic predictor variables
# Optimize Maxent parameter settings
#	Perform resubstitution and cross-validation to assess accuracy

## IMPORTANT NOTES ABOUT RUNNING MAXENT
# Maxent is an open-source java program that should be downloaded here:
# https://biodiversityinformatics.amnh.org/open_source/maxent/
# After download, the maxent.jar file should be moved into the dismo package's java library within R
# On a mac OS, the path is MacHD > Library > Frameworks > R.Framework > Resources > Library > dismo > java
# On a PC, the path is ~/R/win-library/[version]/dismo/java
# You will also need to have Java installed: www.java.com

##################################################################################### LOAD LIBRARIES AND PREPARE INPUTS

## CLEAR WORKSPACE
rm(list = ls(all.names = TRUE))

## LIBRARIES
library(dismo) # for maxent model syntax
library(rJava) # for maxent program
library(ENMeval) # for optimizing model parameters
library(raster) # for making/reading bioclim rasters for ENMeval
library(PresenceAbsence) # for accuracy stats

## INPUTS
# read in file
dat <- read_csv('SDM/data_files/sdm_input.csv')

# slim dataframe to conform to structure required by mod.form function
# (see script modforms.R)
dat.input <- dat %>% 
  dplyr::select(presabs, bio10, bio14, bio15, bio12, bio6, bio3, bio2)

#############################################################
### MAXENT MODEL

# Temporarily give MaxEnt more memory, if needed
#options(java.parameters = "-Xmx1g" )

# model with default parameter settings
mod1.MAX <- maxent(dat.input[,c(2:8)], dat.input$presabs)

# run ENMeval to optimize parameter settings


# model with parameters adjusted based on ENMeval
mod2.MAX <- maxent(dat.input[,c(2:8)], dat.input$presabs, args=c("betamultiplier=2.5", "threshold=FALSE")) 
	
save(mod2.MAX, file="SDM/Output/MAX.mod2.Rda")

##################################################################




##################################################################
### LOAD FINAL MODEL AND ITS PREDICTIONS 

mod <- get(load("SDM/Output/MAX.mod2.Rda"))

##################################################################


##################################################################
### CALCULATE ACCURACY, MODEL = MAXENT

source("SDM/Amy's_SDM_scripts/accuracy.R")
modl = "mod2.MAX" # add var to keep track of model

# Resubstitution
acc <- accuracy.max(dat.input, mod, modl)
acc$thresh = c("SensSpec", "Kappa")
acc$model = "MAX.mod2"
save(acc, file="SDM/Output/MAX.mod2.accs.Rda")
	
# Cross-validation
x.fold=5
n.col=8
cv.acc <- cv.accuracy.max(dat.input, mod, modl, x.fold, n.col)
cv.acc$thresh = c("SensSpec", "Kappa")
cv.acc$model = "MAX.mod1"
save(cv.acc, file="SDM/Output/MAX.mod2.cvaccs.Rda")

################################################################################




################################################################################
######## START SOME EXTRA STUFF

## compare model variable importance
#setwd(path.obj)
library(dismo)
par(ask=T)
par(mfrow=c(1,1))
for (i in 1:10) {
	mod = get(load(paste("MAX.mod1.",i,".Rda", sep="")))   
	plot(mod) 
	}

## compare variable responses
#setwd(path.obj)
library(dismo)
par(ask=T)
par(mfrow=c(2,4))
for (i in 1:10) {
	mod = get(load(paste("MAX.mod1.",i,".Rda", sep="")))   
	response(mod)
	}
	
######## END SOME EXTRA STUFF
################################################################################
