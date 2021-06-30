##############################################################################
### SCRIPT PURPOSE: Clip climate ratser with range extent

# Modified from Angert et al. 2018, American Naturalist
# Author: Daniel Anstett
# last update:  June 23 2021

## OVERALL WORKFLOW:
# Assumes you have bioclim raster and shape file for range extent
# Produces new rater that is contrained by range extent shape file
##############################################################################

## CLEAR WORKSPACE
rm(list = ls(all.names = TRUE))

#Import libraries
library(rgdal)
library(sf)
library(tmap)
library(tidyverse)
library(rnaturalearth)
library(raster)
library(cowplot)
library(rgeos)

#Transform M. cardinalis distribution raster into a shape file
#sdm <- raster("SDM/Output/ThresholdedEnsemble_Unprojected.grd") #bring in ensemble SDM raster
#card <- rasterToPolygons(sdm, dissolve=TRUE)
#raster::shapefile(card,"SDM/Output/c_range_2.shp",overwrite=TRUE)

# Import M.cardinalis ensamble range extent as sf polygon
c_range <- st_read("SDM/Output/c_range_2.shp")

#Import climate variable, 1981 to 2010 rasters
laea_CMD<-raster("Donor_selection/data/1981_2010/CMD.grd")
laea_MAP<-raster("Donor_selection/data/1981_2010/MAP.grd")
laea_MAT<-raster("Donor_selection/data/1981_2010/MAT.grd")
laea_PAS<-raster("Donor_selection/data/1981_2010/PAS.grd")
laea_EXT<-raster("Donor_selection/data/1981_2010/EXT.grd")

laea_PPT_sm<-raster("Donor_selection/data/1981_2010/PPT_sm.grd")
laea_PPT_wt<-raster("Donor_selection/data/1981_2010/PPT_wt.grd")
laea_Tave_sm<-raster("Donor_selection/data/1981_2010/Tave_sm.grd")
laea_Tave_wt<-raster("Donor_selection/data/1981_2010/Tave_wt.grd")

#re-project all into WGS 1984 (EPSG 4326)
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
wgs_CMD <- projectRaster(laea_CMD, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
wgs_MAP <- projectRaster(laea_MAP, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
wgs_MAT <- projectRaster(laea_MAT, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
wgs_PAS <- projectRaster(laea_PAS, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
wgs_EXT <- projectRaster(laea_EXT, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)

wgs_PPT_sm <- projectRaster(laea_PPT_sm, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
wgs_PPT_wt <- projectRaster(laea_PPT_wt, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
wgs_Tave_sm <- projectRaster(laea_Tave_sm, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
wgs_Tave_wt <- projectRaster(laea_Tave_wt, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
c_range <- st_transform(c_range, crs = 4326) # reproject to WGS 1984 (EPSG 4326)


#CRS match
crs(c_range) 
crs(wgs_CMD)

# Reduce range extent to exact rectangular extent of M. cardinalis species range
CMD.clip <- raster::crop(wgs_CMD, extent(c_range))
MAP.clip <- raster::crop(wgs_MAP, extent(c_range))
MAT.clip <- raster::crop(wgs_MAT, extent(c_range))
PAS.clip <- raster::crop(wgs_PAS, extent(c_range))
EXT.clip <- raster::crop(wgs_EXT, extent(c_range))

PPT_sm.clip <- raster::crop(wgs_PPT_sm, extent(c_range))
PPT_wt.clip <- raster::crop(wgs_PPT_wt, extent(c_range))
Tave_sm.clip <- raster::crop(wgs_Tave_sm, extent(c_range))
Tave_wt.clip <- raster::crop(wgs_Tave_wt, extent(c_range))

# Mask climate variable by exact range extent
CMD.mask <- mask(CMD.clip, c_range)
MAP.mask <- mask(MAP.clip, c_range)
MAT.mask <- mask(MAT.clip, c_range)
PAS.mask <- mask(PAS.clip, c_range)
EXT.mask <- mask(EXT.clip, c_range)

PPT_sm.mask <- mask(PPT_sm.clip, c_range)
PPT_wt.mask <- mask(PPT_wt.clip, c_range)
Tave_sm.mask <- mask(Tave_sm.clip, c_range)
Tave_wt.mask <- mask(Tave_wt.clip, c_range)

#visualize
plot(CMD.mask,main="CMD: Climatic moisture deficit (mm)")
plot(MAP.mask,main="MAP: Precipitation as snow (mm)")
plot(MAT.mask,main="MAT: Mean annual precipitation (mm)")
plot(PAS.mask,main="PAS: recipitation as snow (mm)")
plot(EXT.mask,main="EXT: Extreme temperature over 30 years(°C)")

plot(PPT_sm.mask,main="Tave_wt: Winter mean temperature (°C)")
plot(PPT_wt.mask,main="Tave_sm = Summer mean temperature (°C)")
plot(Tave_sm.mask,main="PPT_wt = Winter precipitation (mm)")
plot(Tave_wt.mask,main="PPT_sm = Summer precipitation (mm)")



#Write out Raster
writeRaster(CMD.mask, file="Donor_selection/data/mask/CMD.mask.grd", overwrite=TRUE)
writeRaster(MAP.mask, file="Donor_selection/data/mask/MAP.mask.grd", overwrite=TRUE)
writeRaster(MAT.mask, file="Donor_selection/data/mask/MAT.mask.grd", overwrite=TRUE)
writeRaster(PAS.mask, file="Donor_selection/data/mask/PAS.mask.grd", overwrite=TRUE)
writeRaster(EXT.mask, file="Donor_selection/data/mask/EXT.mask.grd", overwrite=TRUE)

writeRaster(PPT_sm.mask, file="Donor_selection/data/mask/PPT_sm.mask.grd", overwrite=TRUE)
writeRaster(PPT_wt.mask, file="Donor_selection/data/mask/PPT_wt.mask.grd", overwrite=TRUE)
writeRaster(Tave_sm.mask, file="Donor_selection/data/mask/Tave_sm.mask.grd", overwrite=TRUE)
writeRaster(Tave_wt.mask, file="Donor_selection/data/mask/Tave_wt.mask.grd", overwrite=TRUE)




