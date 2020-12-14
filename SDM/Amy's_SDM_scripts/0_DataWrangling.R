################################################################################# 
### SCRIPT PURPOSE: Cull and reformat input files 
# Reworking files used in Angert et al. 2018, Am Nat so that starting inputs are representative for AM workshop
# Author:  Amy Angert
# last update:  12 Dec 2020

############################################################################### 


############################################################################### 
### LOAD LIBRARIES AND INPUTS

## CLEAR WORKSPACE
rm(list = ls(all.names = TRUE))

## LIBRARIES
library(tidyverse)

## INPUTS
# Master file
# (this file contains rows for all cleaned herbarium records, all occupancy points, and 20K background points and columns of climate variables)
all <- read_csv("SDM/data_files/all.records.aug.31.csv")

# Thinned herbarium records
# (this file has herbarium records that have been thinned to reduce oversampling (per folder herb_thinning), but it is missing elevation and year)
# it is replicate 6 of the 10 used in Am Nat 2018
thin <- read_csv("SDM/data_files/herbarium_thinned.csv")

############################################################################### 

############################################################################### 
## WRANGLE UNUSED AND MISSING COLUMNS

# Remove climate columns from all
all <- all %>% dplyr::select(MASTER.ID, DATASET, PRESABS, YEAR, Latitude, Longitude, Elevation)

# Pull out occupancy presence records
occ <- all %>% filter(DATASET=="occ" & PRESABS==1)

# Pull out herbarium presence records
# (this file has not been thinned, but it has useful covariates)
herb <- all %>% filter(DATASET=="herb" & PRESABS==1)

# Retain rows of herb that match those selected during thinning
herb.thin <- left_join(thin, herb, by="MASTER.ID") %>% dplyr::select(MASTER.ID, DATASET, PRESABS, YEAR, Latitude, Longitude, Elevation)

# Append occupancy and thinned herbarium records for complete set of usable presences
# (this will become our input for script 1)
presences <- bind_rows(herb.thin, occ)
write_csv(presences, "SDM/data_files/presences.csv")

############################################################################### 