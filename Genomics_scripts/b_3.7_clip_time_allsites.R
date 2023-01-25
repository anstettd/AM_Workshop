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
baseline_pop <- read_csv("Genomics_scripts/Data/paper_ID_site_select.csv")
timeseries_pop <- baseline_pop %>% filter(Paper_ID<13) %>% dplyr::select(Long,Lat, Paper_ID,)
#imeseries_pop_sf <- st_as_sf(timeseries_pop,coords=c("Long","Lat"), crs=EPSG4326)
#time_sp <- SpatialPoints(timeseries_pop)
coordinates(timeseries_pop)=cbind(timeseries_pop$Long,timeseries_pop$Lat)
proj4string(timeseries_pop) <- EPSG4326
crs(timeseries_pop)

rasStack_1016 <- raster::extract(rasStack_1016,timeseries_pop)
rasValue_gcc <- raster::extract(rasStack_gcc,timeseries_pop)

#Save offset data in dataframe
timeseries_offset <- baseline_pop %>% filter(Paper_ID<13)
offset_pop <- cbind(timeseries_offset,rasStack_1016,rasValue_gcc)


write_csv(offset_pop,"Genomics_scripts/Data/offset_pop_bf5.csv")


