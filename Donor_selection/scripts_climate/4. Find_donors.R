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
library(rgdal) # for transforming projections
library(raster) # for handeling rasters

## INPUTS

#Data
gcc.clim <- read_csv("Donor_selection/data/gcc.clim.csv")
# slim dataframe to just bio10 (average temp of warmest quarter)
gcc.clim <- gcc.clim %>% dplyr::select(GCM, ID1, ID2, Latitude, Longitude, Elevation, 
                         bio10)

#Rasters
bio10.6190<-raster("Donor_selection/data/bio10.grd") #bring in 1961 to 1990 bio10 raster
SDM<-raster("SDM/Output/UnweightedEnsemble_Unprojected.grd") #bring in ensemble SDM raster







