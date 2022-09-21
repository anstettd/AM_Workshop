##################################################################################
## Gradient forest plots
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
#mask_offset_2011 <- raster("Genomics_scripts/Data/offset_2016.tif") #pop data
#mask_offset_2012 <- raster("Genomics_scripts/Data/offset_2016.tif") #pop data
#mask_offset_2013 <- raster("Genomics_scripts/Data/offset_2016.tif") #pop data
#mask_offset_2014 <- raster("Genomics_scripts/Data/offset_2016.tif") #pop data
#mask_offset_2015 <- raster("Genomics_scripts/Data/offset_2016.tif") #pop data
mask_offset_2016 <- raster("Genomics_scripts/Data/offset_2016.tif") #pop data

#Future cliamte change offset raster
mask_offset_45 <- raster("Genomics_scripts/Data/offset_4.5_peakbf2.tif") #pop data
mask_offset_85 <- raster("Genomics_scripts/Data/offset_8.5_peakbf2.tif") #pop data

#Define CRS
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS

#Mask by M. caridnalis range extent
c_range <- st_read("Shape/c_range50.shp") 
c_range <- st_transform(c_range, crs = 4326) # reproject to WGS 1984 (EPSG 4326)
#mask_offset_2011 <- mask(mask_offset_2011, c_range)
#mask_offset_2012 <- mask(mask_offset_2012, c_range)
#mask_offset_2013 <- mask(mask_offset_2013, c_range)
#mask_offset_2014 <- mask(mask_offset_2014, c_range)
#mask_offset_2015 <- mask(mask_offset_2015, c_range)
mask_offset_2016 <- mask(mask_offset_2016, c_range)
mask_offset_45 <- mask(mask_offset_45, c_range)
mask_offset_85 <- mask(mask_offset_85, c_range)

rasStack_gcc <- stack(mask_offset_45,mask_offset_85)

#Baseline & Timeseries
baseline_pop <- read_csv("Genomics_scripts/Data/paper_ID_site_select.csv")
timeseries_pop <- baseline_pop %>% filter(Paper_ID<13) %>% dplyr::select(Long,Lat, Paper_ID,)
#imeseries_pop_sf <- st_as_sf(timeseries_pop,coords=c("Long","Lat"), crs=EPSG4326)
#time_sp <- SpatialPoints(timeseries_pop)
coordinates(timeseries_pop)=cbind(timeseries_pop$Long,timeseries_pop$Lat)
proj4string(timeseries_pop) <- EPSG4326
crs(timeseries_pop)

rasValue_2016 <- raster::extract(mask_offset_2016,timeseries_pop)
rasValue_gcc <- raster::extract(rasStack_gcc,timeseries_pop)

#Save offset data in dataframe
timeseries_offset <- baseline_pop %>% filter(Paper_ID<13)
offset_pop <- cbind(timeseries_offset,rasValue_2016,rasValue_gcc)

#Plot offset against latitude
#4.5
ggplot(offset_pop, aes(x=Elevation, y=offset_4.5_peakbf2)) + 
  geom_point()+
  geom_smooth(method=lm)+
  scale_y_continuous(name="Genetic Offset RCP4.5")+
  scale_x_continuous(name="Elevation (m)")+
  theme_classic() + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))

ggplot(offset_pop, aes(x=Lat, y=offset_4.5_peakbf2)) + 
  geom_point()+
  geom_smooth(method=lm)+
  scale_y_continuous(name="Genetic Offset RCP4.5")+
  scale_x_continuous(name="Latitude")+
  theme_classic() + theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))

#8.5
ggplot(offset_pop, aes(x=Elevation, y=offset_8.5_peakbf2)) + 
  geom_point()+
  geom_smooth(method=lm)+
  scale_y_continuous(name="Genetic Offset RCP8.5")+
  scale_x_continuous(name="Elevation (m)")+
  theme_classic() + theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))

ggplot(offset_pop, aes(x=Lat, y=offset_8.5_peakbf2)) + 
  geom_point()+
  geom_smooth(method=lm)+
  scale_y_continuous(name="Genetic Offset RCP8.5")+
  scale_x_continuous(name="Latitude")+
  theme_classic() + theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))







