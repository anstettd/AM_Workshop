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
library(stars)
library(tidyverse)
library(tmap)
library(rnaturalearth)


#bring in climate raster
nlayers(stack("Donor_selection/data/bio10.grd")) #has just one layer
bio10<-raster("Donor_selection/data/bio10.grd") #bring in 1961 to 1990 bio10 raster

#bring in range map
nlayers(stack("SDM/Output/ThresholdedEnsemble_Unprojected.grd")) #has just one layer
sdm<-raster("SDM/Output/ThresholdedEnsemble_Unprojected.grd") #bring in ensemble SDM raster
#sdm_s<-raster("SDM/Output/UnweightedEnsemble_LCCProjection.grd") #bring in ensemble SDM raster

#covert raster to sf file












#convert raster to polygon
sdm.poly<-rasterToPolygons(sdm, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE)
bio10.poly<-rasterToPolygons(bio10, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE)

#plot
states<-ne_states(country=c("canada","united states of america"), #bring in NA Map
                  returnclass= "sf")
tm_shape(states)+ #not working right now
  m_polygon(bio10)

#confirm projection
st_crs(bio10.poly) == st_crs(sdm.poly) #not in the same coordinate system
st_crs(sdm.poly) 
st_crs(bio10.poly)

#re-project SDM
sdm.pt<- st_transform(sdm.poly,st_crs(bio10.poly))










#confirm projection
crs(sdm) #in WGS1984
crs(bio10) #in WGS1984
crs(bio10) <- 



st_crs(bio10) == st_crs(sdm)
st_crs(bio10) == st_crs(sdm_s)



crs(smd_s)







#import rasters
bio10 <- read_stars("Donor_selection/data/bio10.grd")
SDM <- read_stars("SDM/Output/ThresholdedEnsemble_Unprojected.grd")









#To Do
# get SDM into .shp file
# get to files into same coordinate system
# see https://www.earthdatascience.org/courses/earth-analytics/lidar-raster-data-r/crop-raster-data-in-r/
# https://www.nceas.ucsb.edu/sites/default/files/2020-04/OverviewCoordinateReferenceSystems.pdf




#sdm_s<-raster("SDM/Output/UnweightedEnsemble_LCCProjection.grd") #bring in ensemble SDM raster LCC





