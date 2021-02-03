##############################################################################
### SCRIPT PURPOSE: Clip climate ratser with range extent

# Modified from Angert et al. 2018, American Naturalist
# Author: Daniel Anstett
# last update:  3 Feb 2021

## OVERALL WORKFLOW:
# Assumes you have bioclim raster and shape file for range extent
# Produces new rater that is contrained by range extent shape file
##############################################################################

library(sf)
library(tmap)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)

#bring in raster
bio10<-raster()

#bring in range map


#confirm projection








calo <- states %>%
  filter(name_en=="Oregon" |
           name_en=="California")


#plot selected states, by wes bbox
tm_shape(calo) +
  tm_borders()
