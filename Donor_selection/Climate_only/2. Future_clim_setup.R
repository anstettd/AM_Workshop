##############################################################################
### SCRIPT PURPOSE: Transform future climate at select sites to bioclim

# Modified from Angert et al. 2018, American Naturalist
# Author: Daniel Anstett & Amy Angert
# last update:  5 Jan 2021

## OVERALL WORKFLOW:
# !!! For future climate data, this assumes you have run climateNA (http://climatena.ca/) external to this script and have its output.
# Calculate bioclimatic variables from monthly temperature and precipitation
##############################################################################

## CLEAR WORKSPACE
rm(list = ls(all.names = TRUE))

## LIBRARIES
library(tidyverse) # for data manipulation
library(dismo) # for biovars function
library(rgdal) # for transforming projections

## INPUTS
## (just points files for now; rasters loaded below)
gcc.clim <- read_csv("Donor_selection/data/timeseries_lat_4GCMs_AMT.csv")

# Ensure no missing values
complete <- gcc.clim[complete.cases(gcc.clim), ] # should be same dim as clim if there are no missing values

##############################################################################


##############################################################################
### CALCULATE DERIVED BIOCLIMATIC VARIABLES
## 19 bioclimatic variables from Worldclim (https://www.worldclim.org/data/bioclim.html)

# Requires matrices of Tmin, Tmax, and Precip
tmin <- gcc.clim %>% 
  dplyr::select(Tmin01, Tmin02, Tmin03, Tmin04, Tmin05, Tmin06, Tmin07, Tmin08, Tmin09, Tmin10, Tmin11, Tmin12) %>% 
  as.matrix()

tmax <- gcc.clim %>% 
  dplyr::select(Tmax01, Tmax02, Tmax03, Tmax04, Tmax05, Tmax06, Tmax07, Tmax08, Tmax09, Tmax10, Tmax11, Tmax12) %>% 
  as.matrix()

prec <- gcc.clim %>% 
  dplyr::select(PPT01, PPT02, PPT03, PPT04, PPT05, PPT06, PPT07, PPT08, PPT09, PPT10, PPT11, PPT12) %>% 
  as.matrix()

# Calculate bio1-bio19
bio <- biovars(prec, tmin, tmax)
bio <- data.frame(bio)

# Bind back to pres/abs dataset
bioclim.f <- cbind(gcc.clim, bio) %>% 
  dplyr::select(GCM, ID1, ID2, Latitude, Longitude, Elevation, 
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
write_csv(bioclim.f, "Donor_selection/data/gcc.clim.csv")

##############################################################################

