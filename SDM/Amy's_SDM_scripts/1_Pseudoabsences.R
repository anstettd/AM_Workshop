################################################################################# 
### SCRIPT PURPOSE: Generate pseudoabsence records 
  # Modified from Angert et al. 2018, American Naturalist
  # Author: Amy Angert
  # last update:  13 Dec 2020

## Simplifications and changes made for AM workshop
  # Converted to tidyverse syntax (unless manipulating spatial files)
  # Using relative paths within R project instead of setwd commands
  # Chose 1 thinned presence input file so that everything is not done across i=1:10 replicates (using i=6 from Am Nat paper, chosen at random)
  # Building model with full presence dataset (thinned herbarium records + occupancy points) instead of reserving occupancy points for testing
  # Exporting one master set of points instead of separate sets for each ratio of pseudabsences:presences

## Still to do
  # Add missing code to draw background points in the first place

###############################################################################



############################################################################### 
### LOAD LIBRARIES AND INPUTS

## CLEAR WORKSPACE
rm(list = ls(all.names = TRUE))

## LIBRARIES
library(tidyverse)
library(sp)
library(raster)
library(rgeos)

## INPUTS
# Read in presence records
# (these are from herbarium collections, thinned to reduce oversampled areas, and randomized occupancy surveys; see script 0)
pres <- read_csv("SDM/data_files/presences.csv")

# convert to spatial points file
coordinates(pres) <- ~Longitude+Latitude
projection(pres) <- CRS('+proj=longlat')

################################################################################


################################################################################
### DRAW 20K BACKGROUND POINTS
    
## FIGURE OUT WHERE THIS CODE IS...    
    
# For now, using pre-existing background points
# ***THIS NEEDS TO CHANGE***
pseudo <- read_csv("SDM/data_files/all.records.aug.31.csv") %>% 
      filter(DATASET=="herb" & PRESABS==0) %>% 
      dplyr::select(MASTER.ID, DATASET, PRESABS, YEAR, Latitude, Longitude, Elevation)

# Convert to spatial points files
coordinates(pseudo) = ~Longitude + Latitude
projection(pseudo) = CRS('+proj=longlat')

################################################################################
    

################################################################################
### SUBSAMPLE FROM BACKGROUND POINTS TO MAKE PSEUDOABSENCES 
## RATIONALE: Pseudoabsences should fall in sampled regions yet not in occupied habitat. Here, the zerodist2 function is used to sample background points from within donut-shaped areas surrounding each presence record, so that resulting pseudoabsences are >0.6 km but <80 km from a presence.

# Keep any background points within 80 km of a presence record
ok = zerodist2(pseudo, pres, 80)
psu1 = pseudo[unique(ok[,1]),] # not using tidyverse here - doesn't work on spatial data frames

# Drop any background points that are closer then 0.6 km from a presence record
drop = zerodist2(pseudo, pres, 0.600)
psu2 = psu1[-drop[,1],]
    
# Randomly select n pseudoabsences, where n= number of presences (i.e., 1:1 ratio abs:pres)
# (Angert et al. 2018 showed this matches natural prevalence and also reduces overfitting)
npsu = length(pres) 
samp = sample(1:nrow(psu2), npsu, replace=FALSE) 	
psu2.use = psu2[samp,]

################################################################################


################################################################################
## APPEND PRESENCES & PSEUDOABSENCES

# Remove spatial class
pres <- as.data.frame(pres)
psu2.use = as.data.frame(psu2.use)

# Bind rows & format as input for ClimateNA
# (required column names = ID1, ID2, lat, long, el; must have exactly and only these)
points <- bind_rows(pres, psu2.use) %>% 
  dplyr::select(ID1=MASTER.ID, ID2=PRESABS, lat=Latitude, long=Longitude, el=Elevation)

# Write to file
write_csv(points, "SDM/data_files/points.csv")

################################################################################

# These files should be imported into the ClimateNA graphical user interface
# In ClimateNA, use the time series option to get 30-year average of monthly temperature and precipitation for 1961-1990
# Save ClimateNA output with default name: points_Normal_1961_1990MSY.csv
