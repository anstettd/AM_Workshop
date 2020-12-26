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
# Read in file of presences and pseudoabsences and their associated bioclim values
dat <- read_csv('SDM/data_files/sdm_input.csv')

# Slim dataframe to conform to structure required by mod.form function
# (see script modforms.R)
# This is unprojected
dat.input <- dat %>% 
  dplyr::select(presabs, bio10, bio15, bio11, bio17, bio12, bio3, bio2) %>% 
  as.data.frame()

# Read in rasters of bioclim variables for us in ENMeval
# These are in LCC projection
preds.lcc <- stack("SDM/data_files/bio2.grd", 
               "SDM/data_files/bio3.grd", 
               "SDM/data_files/bio10.grd", 
               "SDM/data_files/bio11.grd", 
               "SDM/data_files/bio12.grd", 
               "SDM/data_files/bio15.grd", 
               "SDM/data_files/bio17.grd")

# Unproject rasters to match pres/abs
prj.lcc <- "+proj=lcc +lon_0=-95 +lat_1=49 +lat_2=77 +type=crs"  #specify current projection 
proj4string(preds.lcc) <- CRS(prj.lcc) #define current projection of rasters
prj.wgs = "+proj=longlat + type=crs" #specify unprojected coordinate system
preds = projectRaster(preds.lcc, crs=CRS(prj.wgs)) #transform raster projection

#################################################################################


#################################################################################
### RUN ENMeval to optimize parameter settings
## See here for user guide: https://cran.r-project.org/web/packages/ENMeval/vignettes/ENMeval-vignette.html#eval

dat.pres <- dat %>% filter(presabs==1) %>% dplyr::select(Longitude, Latitude) %>% as.data.frame()
dat.abs <- dat %>% filter(presabs==0) %>% dplyr::select(Longitude, Latitude) %>% as.data.frame()

tuna = ENMevaluate(dat.pres, preds, bg.coords=dat.abs, method="randomkfold", kfold=5, algorithm='maxent.jar')

# See table of evaluation metrics
tuna@results

# Find model with lowest AICc
bestmod <- which(tuna@results$delta.AICc == 0) 
tuna@results[bestmod,]
# This identifies the model with linear and quadratic features only (no hinge, threshold, or product features) and a low regularization multiplier of 0.5

# Save the best model 
mod.MAX <- tuna@models[[bestmod]]
save(mod.MAX, file="SDM/Output/MAX.mod.Rda")


