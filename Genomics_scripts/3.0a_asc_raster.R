#############################################################################################################
## Downscale asc raster from Climate NA to cover range extent
## Author Daniel Anstett
## 
## 
## Last Modified September 16, 2022
#############################################################################################################
#Import libraries
library(raster)
library(tidyverse)
library(sf)

#############################################################################################################
#Import asc file
#Range extent all of North America
na_asc <- raster("C:/Users/anstett3/Downloads/Climatena_v730/InputFiles/na800.asc")

# Import M.cardinalis ensamble range extent as sf polygon
c_range <- st_read("SDM/Output/c_range_2.shp")

#Crop asc file by M. cardinalis range extent
c_800 <- raster::crop(na_asc, extent(c_range))

#Save c_800 raster as an asc file
writeRaster(c_800, "C:/Users/anstett3/Downloads/Climatena_v730/InputFiles/c_800.asc",
                format="ascii", overwrite = T)