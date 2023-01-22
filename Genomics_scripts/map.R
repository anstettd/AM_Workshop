##################################################################################
## Make population mapes of baseline and timeseires
## Author Daniel Anstett
## 
## 
## Last Modified Jan 20, 2023
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
#Define CRS
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS

#Baseline & Timeseries
baseline_pop <- read_csv("Genomics_scripts/Data//Baseline_Timeseries_pops_final2.csv")%>% filter(Paper_ID<55)
timeseries_pop <- baseline_pop %>% filter(Paper_ID<13) %>% dplyr::select(Long,Lat)
timeseries_pop_sf <- st_as_sf(timeseries_pop,coords=c("Long","Lat"), crs=EPSG4326)
baseline_pop <- baseline_pop %>% filter(Paper_ID>12) %>% dplyr::select(Long,Lat)
baseline_pop_sf <- st_as_sf(baseline_pop,coords=c("Long","Lat"), crs=EPSG4326)

# California & Oregon Map Setup
states<-ne_states(country=c("canada","united states of america"),returnclass= "sf")
calo <- states %>%
  filter(name_en=="Oregon" | name_en=="California") #| name_en=="Nevada")

#Plot offset SSP245 (RCP 4.5)
tmap_mode("plot")
#tmap_mode("view")
mim <-
  tm_shape(calo)+
  tm_borders()+
  tm_shape(baseline_pop_sf)+
  tm_dots(size=0.4,shape=21,col="blue")+
  tm_shape(timeseries_pop_sf)+
  tm_dots(size=0.6,shape=21,col="red")+
  tm_layout(legend.position = c(1.03, 0.73),legend.title.size = 0.001,frame.lwd = NA)
mim
tmap_save(mim, filename = "maps/base_time1.png",width=4, height=7)


