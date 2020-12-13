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

## RATIONALE FOR STEP 1: Pseudoabsences should fall in sampled regions yet not in occupied habitat. Here, the zerodist2 function is used to sample background points from within donut-shaped areas surrounding each presence record, so that resulting pseudoabsences are >0.6 km but <80 km from a presence.

# Keep any background points within 80 km of a presence record
ok = zerodist2(pseudo, pres, 80)
psu1 = pseudo[unique(ok[,1]),] # not using tidyverse here - doesn't work on spatial data frames

# Drop any background points that are closer then 0.6 km from a presence record
drop = zerodist2(pseudo, pres, 0.600)
psu2 = psu1[-drop[,1],]
    
## RATIONALE FOR STEP 2: Different ENM algorithms are optimized for different ratios of presence records relative to pseudoabsences.
    
# 10:1 ratio abs:pres for GLM and GAM algorithms
npsu.glm = 10*length(pres) 
samp.glm = sample(1:nrow(psu2), npsu.glm, replace=FALSE) #subsample from anything meeting distance criteria	
psu2.glm = psu2[samp.glm,]
psu2.glm$set = "glm_gam"
psu2.glm <- as.data.frame(psu2.glm)

# 4:1 ratio abs:pres for RF and MAX algorithms
npsu.rf = 4*length(pres) 
samp.rf = sample(1:nrow(psu2.glm), npsu.rf, replace=FALSE) #subsample from those retained for GLM/GAM	
psu2.rf = psu2.glm[samp.rf,]
psu2.rf$set = "rf_max"

# 1:1 ratio abs:pres for BRT algorithm
npsu.brt = length(pres) 
samp.brt = sample(1:nrow(psu2.rf), npsu.brt, replace=FALSE) #subsample from those retained for RF/MAX	
psu2.brt = psu2.rf[samp.brt,]
psu2.brt$set = "brt"

## MERGE into one set that can be filtered for appropriate input to each algorithm
# What we want is brt + (rf-brt) + (glm-rf)
psu2.rf_brt <- anti_join(psu2.rf, psu2.brt, by="MASTER.ID") 
psu2.glm_rf <- anti_join(psu2.glm, psu2.rf, by="MASTER.ID") 

psu.all <- bind_rows(psu2.brt, psu2.rf_brt, psu2.glm_rf)

# check numbers
psu.all %>% count(set) #should have 303 brt, 1212-303=909 rf/max, 3030-1212=1818 glm/gam
psu.all %>% n_distinct("MASTER.ID") #should each be unique, n=3030

################################################################################


################################################################################
## APPEND PRESENCES & PSEUDOABSENCES

# Remove spatial class
pres <- as.data.frame(pres)
pres$set <- "all"

# Bind rows & format as input for ClimateNA
# (required column names = ID1, ID2, lat, long, el; must have exactly and only these)
points <- bind_rows(pres, psu.all) %>% 
  dplyr::select(ID1=set, ID2=PRESABS, lat=Latitude, long=Longitude, el=Elevation)

# Write to file
write_csv(points, "SDM/data_files/points.csv")

################################################################################

# These files should be imported into the ClimateNA graphical user interface
# In ClimateNA, use the time series option to get 30-year average of monthly temperature and precipitation for 1961-1990
# Save ClimateNA output as climatena.csv
