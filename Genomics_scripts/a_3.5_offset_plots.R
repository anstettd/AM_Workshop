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
library(tmap)
library(rnaturalearth)
library(rnaturalearthdata)
#library(devtools)
#devtools::install_github("ropensci/rnaturalearthhires")
library(rnaturalearthhires)
library(rgeos)
library(RColorBrewer)
##############################################################################

## INPUTS

#Timeseries offset raster
#mask_offset_2016 <- raster("Genomics_scripts/Data/offset_2016.tif") #pop data

#Future cliamte change offset raster
mask_offset_45 <- raster("Genomics_scripts/Data/offset_4.5_peakbf2.tif") #pop data
mask_offset_85 <- raster("Genomics_scripts/Data/offset_8.5_peakbf2.tif") #pop data


#Define CRS
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS

#Mask by M. caridnalis range extent
#c_range <- st_read("SDM/Output/c_range_2.shp") 
#c_range <- st_transform(c_range, crs = 4326) # reproject to WGS 1984 (EPSG 4326)
#mask_offset_2016 <- mask(mask_offset_2016, c_range)
#mask_offset_45 <- mask(mask_offset_45, c_range)
#mask_offset_85 <- mask(mask_offset_85, c_range)


# California & Oregon Map Setup
states<-ne_states(country=c("canada","united states of america"),returnclass= "sf")
calo <- states %>%
  filter(name_en=="Oregon" | name_en=="California" | name_en=="Nevada")

#Baseline & Timeseries
baseline_pop <- read_csv("Genomics_scripts/Data/paper_ID_site_select.csv")
timeseries_pop <- baseline_pop %>% filter(Paper_ID<13) %>% dplyr::select(Long,Lat)
timeseries_pop_sf <- st_as_sf(timeseries_pop,coords=c("Long","Lat"), crs=EPSG4326)
baseline_pop <- baseline_pop  %>% dplyr::select(Long,Lat)
baseline_pop_sf <- st_as_sf(baseline_pop,coords=c("Long","Lat"), crs=EPSG4326)
##############################################################################

#2011-2016

off_pallet <- c("#2166AC","#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15")

#Plot offset SSP245 (RCP 4.5)
tmap_mode("plot")
#tmap_mode("view")
offset_2016 <- tm_shape(mask_offset_2016, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = off_pallet)+
  #  tm_raster(palette = rev(brewer.pal(6, "RdBu")))+
  #  tm_raster(palette = "Reds")+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(timeseries_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  #  tm_shape(baseline_pop_sf)+
  #  tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(1.03, 0.73),legend.title.size = 0.001)
offset_2016
#tmap_save(offset_2016, filename = "Offset_graphs/offset2016_peakbf2.pdf",width=4, height=7)


##############################################################################

#Future climate change

off_pallet <- c("#2166AC","#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15")

#Plot offset SSP245 (RCP 4.5)
tmap_mode("plot")
#tmap_mode("view")
offset45 <- tm_shape(mask_offset_45, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = off_pallet)+
#  tm_raster(palette = rev(brewer.pal(6, "RdBu")))+
#  tm_raster(palette = "Reds")+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(timeseries_pop_sf)+
  tm_dots(size=0.5,shape=1,col="black",border.lwd = 2.5)+
#  tm_shape(baseline_pop_sf)+
#  tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.62, 0.48),legend.title.size = 0.001)
offset45
tmap_save(offset45, filename = "Offset_graphs/offset45_peakbf2.pdf",width=4, height=7)

off_pallet2 <- c("#2166AC","#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15","#5c0915")
off_pallet3 <- c("#2166AC","#67A9CF","#D1E5F0","#f7c1c8","#f21836","#A50F15")

#Plot offset SSP585 (RCP 8.5)
tmap_mode("plot")
#tmap_mode("view")
offset85 <- tm_shape(mask_offset_85, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = off_pallet3)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(timeseries_pop_sf)+
  tm_dots(size=0.5,shape=1,col="black",border.lwd = 2.5)+
  #  tm_shape(baseline_pop_sf)+
  #  tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.62, 0.48),legend.title.size = 0.001)
offset85
tmap_save(offset85, filename = "Offset_graphs/offset85_peakbf2.pdf",width=4, height=7)


