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
  dplyr::select(presabs, bio10, bio15, bio11, bio17, bio12, bio3, bio2)

#################################################################################
### MAKE RASTERS of point-associated bioclim for use in ENMeval

# Order the dataframe so when it is rasterized presence records have their original values
dat <- dat[order(rev(dat$presabs),decreasing=TRUE), ]

# Convert to spatial points
coordinates(dat) = ~Longitude + Latitude 
prj.wgs = "+proj=longlat +ellps=WGS84"
proj4string(dat) = CRS(prj.wgs)

# Load existing raster (90m DEM from hydrosheds) 
dem <- raster(file.choose()) # click on to .bil file within na_dem_30s_bil folder

# rasterize thinned pseudo-absence records for Maxent Background
vars <- c("bio2", "bio3", "bio10", "bio11", "bio12", "bio15", "bio17")

for(i in 1:length(vars)){
  dat_rast <- rasterize(dat, dem, field=vars[i], fun='first', na.rm=TRUE, background=NA) # consider changing function to mean
  assign(vars[i], dat_rast)
}

# stack rasters
preds <- stack(bio2, bio3, bio10, bio11, bio12, bio15, bio17)
names(preds) <- vars

#################################################################################
### RUN ENMeval to optimize parameter settings

dat <- as.data.frame(dat)
dat.pres <- dat %>% filter(presabs==1) %>% dplyr::select(Longitude, Latitude) 
dat.abs <- dat %>% filter(presabs==0) %>% dplyr::select(Longitude, Latitude) 

tuna = ENMevaluate(dat.pres, preds, bg.coords=dat.abs, method="randomkfold", kfold=5, algorithm='maxent.jar')

# see table of evaluation metrics
tuna@results

# find model with lowest AICc
which(tuna@results$delta.AICc == 0)

#############################################################
### RUN MAXENT MODEL

# Temporarily give MaxEnt more memory, if needed
#options(java.parameters = "-Xmx1g" )

# model with default parameter settings
mod1.MAX <- maxent(dat.input[,c(2:8)], dat.input$presabs)

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
