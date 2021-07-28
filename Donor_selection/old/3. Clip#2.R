##############################################################################
### SCRIPT PURPOSE: Clip climate ratser with range extent

# Modified from Angert et al. 2018, American Naturalist
# Author: Daniel Anstett
# last update:  9 Feb 2021

## OVERALL WORKFLOW:
# Assumes you have bioclim raster and shape file for range extent
# Produces new rater that is contrained by range extent shape file
##############################################################################

#Import libraries
library(rgdal)
library(sf)
library(tmap)
library(tidyverse)
library(rnaturalearth)
library(raster)

#Transform M. cardinalis distribution raster into a shape file
sdm <- raster("SDM/Output/ThresholdedEnsemble_Unprojected.grd") #bring in ensemble SDM raster
card <- rasterToPolygons(sdm, dissolve=TRUE)
raster::shapefile(card,"SDM/Output/c_range.shp")

# Import M.cardinalis ensamble range extent as sf polygon
c_range <- st_read("SDM/Output/c_range.shp")

#Import climate variable
lcc10<-raster("Donor_selection/data/bio_clim/bio10.grd") #bring in 1961 to 1990 bio10 raster
#re-project bio10 raster into WGS 1984 (EPSG 4326)
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
bio10 <- projectRaster(lcc10, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
c_range <- st_transform(c_range, crs = 4326) # reproject to WGS 1984 (EPSG 4326)


#CRS match
crs(c_range) 
crs(bio10)

#Clip
bio10.clip <- raster::crop(bio10, extent(c_range))
bio10.mask<- mask(bio10.clip,c_range)

#writeRaster(bio10.mask, file="Donor_selection/data/bio10.clip.grd")


  
  
  
  
  
  
  
  
  
  
  