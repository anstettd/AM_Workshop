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
# read in file of presences and pseudoabsences and their associated bioclim values
dat <- read_csv('SDM/data_files/sdm_input.csv')

# slim dataframe to conform to structure required by mod.form function
# (see script modforms.R)
# this is unprojected
dat.input <- dat %>% 
  dplyr::select(presabs, bio10, bio15, bio11, bio17, bio12, bio3, bio2)

# read in rasters of bioclim variables for us in ENMeval
# these are in LCC projection
preds.lcc <- stack("SDM/data_files/bio2.grd", 
               "SDM/data_files/bio3.grd", 
               "SDM/data_files/bio10.grd", 
               "SDM/data_files/bio11.grd", 
               "SDM/data_files/bio12.grd", 
               "SDM/data_files/bio15.grd", 
               "SDM/data_files/bio17.grd")

# unproject rasters to match pres/abs
prj.lcc <- "+proj=lcc +lon_0=-95 +lat_1=49 +lat_2=77 +type=crs"  #specify current projection 
proj4string(preds.lcc) <- CRS(prj.lcc) #define current projection of rasters
prj.wgs = "+proj=longlat + type=crs" #specify unprojected coordinate system
preds = projectRaster(preds.lcc, crs=CRS(prj.wgs)) #transform raster projection

#################################################################################


#################################################################################
### RUN ENMeval to optimize parameter settings
## see here for user guide: https://cran.r-project.org/web/packages/ENMeval/vignettes/ENMeval-vignette.html#eval

dat <- as.data.frame(dat)
dat.pres <- dat %>% filter(presabs==1) %>% dplyr::select(Longitude, Latitude) 
dat.abs <- dat %>% filter(presabs==0) %>% dplyr::select(Longitude, Latitude) 

tuna = ENMevaluate(dat.pres, preds, bg.coords=dat.abs, method="randomkfold", kfold=5, algorithm='maxent.jar')

# see table of evaluation metrics
tuna@results

# find model with lowest AICc
bestmod <- which(tuna@results$delta.AICc == 0) 
tuna@results[bestmod,]
# this identifies the model with linear and quadratic features (fc="LQ") and a low regularization multiplier (rm=0.5) settings

#############################################################
### RUN MAXENT MODEL

# Temporarily give MaxEnt more memory, if needed
#options(java.parameters = "-Xmx1g" )

# model with default parameter settings
mod1.MAX <- maxent(dat.input[,c(2:8)], dat.input$presabs)

# model with parameters adjusted based on ENMeval (LQ_0.5)
mod2.MAX <- maxent(dat.input[,c(2:8)], dat.input$presabs, args=c(RMvalues=c(0.5), fc=c("LQ")))
	
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
