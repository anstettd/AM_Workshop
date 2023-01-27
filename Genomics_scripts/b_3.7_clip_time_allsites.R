##################################################################################
## Extract offset values for each timeseires population
## Author Daniel Anstett
## 
## 
## Last Modified September 19, 2022
###################################################################################

#Library install and import
library(tidyverse) # for data manipulation
library(raster) # for handeling rasters
library(tidyverse)
library(cowplot)
library(sf)
library(rgdal)
library(rnaturalearth)
library(rnaturalearthdata)
#library(devtools)
#devtools::install_github("ropensci/rnaturalearthhires")

##############################################################################

## INPUTS

#Timeseries offset raster
mask_offset_1215 <- raster("Genomics_scripts/Data/offset_1215.tif") #pop data
mask_offset_2011 <- raster("Genomics_scripts/Data/offset_2011.tif") #pop data
mask_offset_2012 <- raster("Genomics_scripts/Data/offset_2012.tif") #pop data
mask_offset_2013 <- raster("Genomics_scripts/Data/offset_2013.tif") #pop data
mask_offset_2014 <- raster("Genomics_scripts/Data/offset_2014.tif") #pop data
mask_offset_2015 <- raster("Genomics_scripts/Data/offset_2015.tif") #pop data
mask_offset_2016 <- raster("Genomics_scripts/Data/offset_2016.tif") #pop data
clim_diff_1215 <- raster("Genomics_scripts/Data/clim_distance.tif") #pop data


#Future cliamte change offset raster
mask_offset_45 <- raster("Genomics_scripts/Data/offset_4.5_peakbf2.tif") #pop data
mask_offset_85 <- raster("Genomics_scripts/Data/offset_8.5_peakbf2.tif") #pop data
mask_offset_45_grain <- raster("Genomics_scripts/Data/offset_4.5_peakbf2_grain.tif") #pop data
mask_offset_85_grain <- raster("Genomics_scripts/Data/offset_8.5_peakbf2_grain.tif") #pop data

#Stack rasters
rasStack_gcc <- stack(mask_offset_45,mask_offset_85)
rasStack_1016 <- stack(mask_offset_1215,mask_offset_2011,mask_offset_2012,mask_offset_2013,
                       mask_offset_2014,mask_offset_2015,mask_offset_2016,clim_diff_1215,
                       mask_offset_45_grain,mask_offset_85_grain)

#Define CRS
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS

#Baseline & Timeseries
demo_pop_raw <- read_csv("Genomics_scripts/Data/offset_pop_lambda_allsites.csv")
demo_pop <- demo_pop_raw %>% dplyr::select(Longitude,Latitude,ID)
coordinates(demo_pop)=cbind(demo_pop$Longitude,demo_pop$Latitude)
proj4string(demo_pop) <- EPSG4326
crs(demo_pop)

rasStack_1016 <- raster::extract(rasStack_1016,demo_pop)
rasValue_gcc <- raster::extract(rasStack_gcc,demo_pop)

#Save offset data in dataframe
offset_demo_pop <- cbind(demo_pop_raw ,rasStack_1016,rasValue_gcc)

colnames(offset_demo_pop)[14] <- "offset_4.5"
colnames(offset_demo_pop)[15] <- "offset_8.5"
colnames(offset_demo_pop)[16] <- "offset_4.5_old"
colnames(offset_demo_pop)[17] <- "offset_8.5_old"

write_csv(offset_demo_pop,"Genomics_scripts/Data/offset_demo_pop.csv")


