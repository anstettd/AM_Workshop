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

library(raster)
library(rgdal)
library(sf)
library(tmap)
library(tidyverse)
library(rnaturalearth)

#Transform M. cardinalis distribution raster into a shape file
#sdm <- raster("SDM/Output/ThresholdedEnsemble_Unprojected.grd") #bring in ensemble SDM raster
#card <- rasterToPolygons(sdm)
#raster::shapefile(card,"SDM/Output/c_range.shp")

# Import M.cardinalis ensamble range extent as sf polygon
c_range <- st_read("SDM/Output/c_range.shp")

##visualize not working
#Setup legal boundries
states<-ne_states(country=c("united states of america"), returnclass= "sf")
west <- states %>% filter(name_en=="Oregon" | name_en=="California" | name_en=="Nevada")
#plot using tmap
tm_shape(west)+
tm_borders()+
tm_shape(c_range)+
  tm_polygons()
  
#Import climate variable
lcc10<-raster("Donor_selection/data/bio10.grd") #bring in 1961 to 1990 bio10 raster
#re-project bio10 raster into WGS 1984 (EPSG 4326)
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
bio10 <- projectRaster(lcc10, crs=EPSG4326)

#CRS match
crs(c_range) 
crs(bio10)

#Clip
bio10.clip <- raster::crop(bio10, c_range)




