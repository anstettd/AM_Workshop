##############################################################################
### SCRIPT PURPOSE: Clip climate ratser with range extent

# Modified from Angert et al. 2018, American Naturalist
# Author: Daniel Anstett
# last update:  3 Feb 2021

## OVERALL WORKFLOW:
# Assumes you have bioclim raster and shape file for range extent
# Produces new rater that is contrained by range extent shape file
##############################################################################

#Import libraries
library(sf)
library(raster)
library(tidyverse)
library(tmap)


#bring in raster
nlayers(stack("Donor_selection/data/bio10.grd")) #has just one layer
bio10<-raster("Donor_selection/data/bio10.grd") #bring in 1961 to 1990 bio10 raster
st_crs(bio10) #in WGS1984
states<-ne_states(country=c("canada","united states of america"), #bring in NA Map
  returnclass= "sf")

tm_shape(states)+
  tm_raster(bio10)


#bring in range map


#confirm projection








