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

## LIBRARIES
library(tidyverse) # for data manipulation
library(dismo) # for biovars function

## INPUTS
clim <- read_csv("SDM/data_files/climatena.csv")

# Ensure no missing values
missing <- dat1[complete.cases(clim), ] # should be dim 0 if there are no missing values

################################################################################# 


################################################################################# 
### CALCULATE DERIVED BIOCLIMATIC VARIABLES

## Requires matrices of Tmin, Tmax, and Precip
tmin <- clim %>% 
  select(Tmin01, Tmin02, Tmin03, Tmin04, Tmin05, Tmin06, Tmin07, Tmin08, Tmin09, Tmin10, Tmin11, Tmin12) %>% 
  as.matrix()

tmax <- clim %>% 
  select(Tmax01, Tmax02, Tmax03, Tmax04, Tmax05, Tmax06, Tmax07, Tmax08, Tmax09, Tmax10, Tmax11, Tmax12) %>% 
  as.matrix()

prec <- clim %>% 
  select(PPT01, PPT02, PPT03, PPT04, PPT05, PPT06, PPT07, PPT08, PPT09, PPT10, PPT11, PPT12) %>% 
  as.matrix()

## Calculate bio1-bio19
bio <- biovars(prec, tmin, tmax)
bio <- data.frame(bio)

# Bind back to pres/abs dataset
bioclim <- cbind(clim, bio)
write_csv(bioclim, "SDM/data_files/")

#################################################################################