################################################################################# 
### SCRIPT PURPOSE: Prepare climatic predictor variables associated with each presence and pseudoabsence point coordinate
# Modified from Angert et al. 2018, American Naturalist
# Author: Amy Angert
# last update:  13 Dec 2020

## OVERALL WORKFLOW:
# !!! This assumes you have run climateNA external to this script and have its output.
# Read in climateNA output file with monthly temperature and precipitation records for 1961-1990
# Calculate bioclimatic variables from monthly temperature and precipitation

################################################################################# 


################################################################################# 
### LOAD LIBRARIES AND PREPARE INPUTS

## CLEAR WORKSPACE
rm(list = ls(all.names = TRUE))

## LIBRARIES
library(tidyverse) # for data manipulation
library(dismo) # for biovars function

## INPUTS
clim <- read_csv("SDM/data_files/points_Normal_1961_1990MSY.csv")

# Ensure no missing values
complete <- clim[complete.cases(clim), ] # should be same dim as clim if there are no missing values

################################################################################# 


################################################################################# 
### CALCULATE DERIVED BIOCLIMATIC VARIABLES
## 19 bioclimatic variables from Worldclim (https://www.worldclim.org/data/bioclim.html)

# Requires matrices of Tmin, Tmax, and Precip
tmin <- clim %>% 
  dplyr::select(Tmin01, Tmin02, Tmin03, Tmin04, Tmin05, Tmin06, Tmin07, Tmin08, Tmin09, Tmin10, Tmin11, Tmin12) %>% 
  as.matrix()

tmax <- clim %>% 
  dplyr::select(Tmax01, Tmax02, Tmax03, Tmax04, Tmax05, Tmax06, Tmax07, Tmax08, Tmax09, Tmax10, Tmax11, Tmax12) %>% 
  as.matrix()

prec <- clim %>% 
  dplyr::select(PPT01, PPT02, PPT03, PPT04, PPT05, PPT06, PPT07, PPT08, PPT09, PPT10, PPT11, PPT12) %>% 
  as.matrix()

# Calculate bio1-bio19
bio <- biovars(prec, tmin, tmax)
bio <- data.frame(bio)

# Bind back to pres/abs dataset
bioclim <- cbind(clim, bio) %>% 
  dplyr::select(Master.ID=ID1, presabs=ID2, Latitude, Longitude, Elevation, 
                bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10, 
                bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19) %>% 
         mutate(bio12=log(bio12+0.5), 
                bio13=log(bio13+0.5), 
                bio14=log(bio14+0.5), 
                bio15=log(bio15+0.5), 
                bio16=log(bio16+0.5), 
                bio17=log(bio17+0.5), 
                bio18=log(bio18+0.5), 
                bio19=log(bio19+0.5))
write_csv(bioclim, "SDM/data_files/biovars.csv")

#################################################################################
