################################################################################# 
### SCRIPT PURPOSE: Generate pseudoabsence records 
  # Based on methods used in Angert et al. 2018, Am Nat
  # Author: Amy Angert
  # last update:  11 Dec 2020

## Simplifications and changes made for AM workshop
  # Converted to tidyverse syntax (unless manipulating spatial files)
  # Using relative paths within R project instead of setwd commands
  # Chose 1 thinned presence input file so that everything is not done across i=1:10 replicates (using i=6 from Am Nat paper, chosen at random)
  # Building model with full presence dataset (thinned herbarium records + occupancy points) instead of reserving occupancy points for testing

## Still to do
  # Add missing code to draw background points in the first place
  # Add instructions for downloading climate records
  # Delete merges below where climate variables are pulled in from another pre-existing file

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
samp.glm = sample(1:nrow(psu2), npsu.glm, replace=FALSE) 	
psu2.glm = psu2[samp.glm,]
assign(paste("pseu","a",sep=""), psu2.glm)

# 4:1 ratio abs:pres for RF and MAX algorithms
npsu.rf = 4*length(pres) 
samp.rf = sample(1:nrow(psu2), npsu.rf, replace=FALSE) 	
psu2.rf = psu2[samp.rf,]
assign(paste("pseu","b",sep=""), psu2.rf)

# 1:1 ratio abs:pres for BRT algorithm
npsu.brt = length(pres) 
samp.brt = sample(1:nrow(psu2), npsu.brt, replace=FALSE) 	
psu2.brt = psu2[samp.brt,]
assign(paste("pseu","c",sep=""), psu2.brt)

################################################################################


################################################################################
## APPEND PRESENCES & PSEUDOABSENCES

# Remove spatial class
pres <- as.data.frame(pres)
pseua <- as.data.frame(pseua)
pseub <- as.data.frame(pseub)
pseuc <- as.data.frame(pseuc)

# Bind rows & format as input for ClimateNA
# (required columns = ID1, ID2, lat, long, el)
points.glmgam <- bind_rows(pres, pseua) %>% 
  mutate(ID2=MASTER.ID) %>% 
  dplyr::select(ID1=MASTER.ID, ID2, lat=Latitude, long=Longitude, el=Elevation)
points.rfmax <- bind_rows(pres, pseub) %>% 
  mutate(ID2=MASTER.ID) %>% 
  dplyr::select(ID1=MASTER.ID, ID2, lat=Latitude, long=Longitude, el=Elevation)
points.brt <- bind_rows(pres, pseuc) %>% 
  mutate(ID2=MASTER.ID) %>% 
  dplyr::select(ID1=MASTER.ID, ID2, lat=Latitude, long=Longitude, el=Elevation)

# Write to files
write_csv(points.glmgam, "SDM/data_files/points_glmgam.csv")
write_csv(points.rfmax, "SDM/data_files/points_rfmax.csv")
write_csv(points.brt, "SDM/data_files/points_brt.csv")

################################################################################

# These files should be imported to ClimateNA graphical user interface
# Download monthly temperature and precipitation records for each point


	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	