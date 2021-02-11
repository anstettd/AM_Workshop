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
library(raster) # for handeling rasters
library(cowplot)

## INPUTS

#Data
gcc.clim <- read_csv("Donor_selection/data/gcc.clim.csv")
# slim dataframe to just bio10 (average temp of warmest quarter)
gcc.clim <- gcc.clim %>% dplyr::select(GCM, ID1, ID2, Latitude, Longitude, Elevation, bio10)



gcc.clim.S15<-gcc.clim %>% filter(ID1=="S15")

#Climate Raster
bio10<-raster("Donor_selection/data/bio10.clip.grd") #bring in 1961 to 1990 bio10 raster
plot(bio10, main="Bio 10 Across Range")

bio10.S15.rcp45.2055 <- bio10
bio10.S15.rcp45.2055[bio10<(gcc.clim.S15$bio10[1]-0.5)] <- NA
bio10.S15.rcp45.2055[bio10>(gcc.clim.S15$bio10[1]+0.5)] <- NA
plot(bio10.S15.rcp45.2055, main="2055 at RCP 4.5")

bio10.S15.rcp45.2085 <- bio10
bio10.S15.rcp45.2085[bio10<(gcc.clim.S15$bio10[2]-0.5)] <- NA
bio10.S15.rcp45.2085[bio10>(gcc.clim.S15$bio10[2]+0.5)] <- NA
plot(bio10.S15.rcp45.2085, main="2085 at RCP 4.5")

bio10.S15.rcp85.2055 <- bio10
bio10.S15.rcp85.2055[bio10<(gcc.clim.S15$bio10[3]-0.5)] <- NA
bio10.S15.rcp85.2055[bio10>(gcc.clim.S15$bio10[3]+0.5)] <- NA
plot(bio10.S15.rcp85.2055, main="2055 at RCP 8.5")

bio10.S15.rcp85.2085 <- bio10
bio10.S15.rcp85.2085[bio10<(gcc.clim.S15$bio10[4]-0.5)] <- NA
bio10.S15.rcp85.2085[bio10>(gcc.clim.S15$bio10[4]+0.5)] <- NA
plot(bio10.S15.rcp85.2085, main="2085 at RCP 8.5")









#Sample without referencing
bio10.S15.rcp45.2055 <- bio10
bio10.S15.rcp45.2055[bio10<(20.35-1)] <- NA
bio10.S15.rcp45.2055[bio10>(20.35+1)] <- NA
plot(bio10.S15.rcp45.2055)

bio10.S15.rcp45.2055 <- bio10
bio10.S15.rcp45.2055[bio10<(20.35-0.5)] <- NA
bio10.S15.rcp45.2055[bio10>(20.35+0.5)] <- NA
plot(bio10.S15.rcp45.2055)










