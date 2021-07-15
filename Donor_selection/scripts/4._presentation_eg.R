##############################################################################
### SCRIPT PURPOSE: Make additional pots 

# Author: Daniel Anstett 
# last update:  July 15 2021
##############################################################################

## CLEAR WORKSPACE
rm(list = ls(all.names = TRUE))

## LIBRARIES
library(tidyverse) # for data manipulation
library(raster) # for handeling rasters
library(cowplot)
library(sf)
library(tmap)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
#devtools::install_github("ropensci/rnaturalearthhires")
library(rnaturalearthhires)
library(rgeos)
##############################################################################

## INPUTS

#Future Climate Data per Restoration Site
gcc.year <- read_csv("Donor_selection/data/timeseries_lat_2 GCMsY_2041_2071.csv")
gcc.season <- read_csv("Donor_selection/data/timeseries_lat_2 GCMsS_2041_2071.csv")
#Selected wanted variables
gcc.year <- gcc.year %>% dplyr::select(GCM:Elevation,
                                       CMD,MAP,MAT,PAS,EXT)
gcc.season <- gcc.season %>% dplyr::select(PPT_sm,PPT_wt,Tave_sm,Tave_wt,)
gcc.clim <- cbind(gcc.year,gcc.season)

#Set up known points within M. cardinalis into sf object
known <- read_csv("SDM/data_files/presences.csv")
known <- known %>% dplyr::select(Longitude,Latitude)
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
known_sf <- st_as_sf(known,coords=c("Longitude","Latitude"), crs=EPSG4326)
#check data is set up properly
ggplot()+ geom_sf(data = known_sf)+ ggtitle("Known M. cardinnalis Populations")

#Set up focal 55 populations into sf object
pop_var_raw <- read_csv("Genomics_scripts/Data/paper_ID_site_select.csv")
gen_pop <- pop_var_raw %>% dplyr::select(Long,Lat)
gen_pop_sf <- st_as_sf(gen_pop,coords=c("Long","Lat"), crs=EPSG4326)
#check data is set up properly
ggplot()+ geom_sf(data = gen_pop_sf)+ ggtitle("Focal Populations")



# California & Oregon Map Setup
states<-ne_states(country=c("canada","united states of america"),returnclass= "sf")
calo <- states %>%
  filter(name_en=="Oregon" | name_en=="California" | name_en=="Nevada")
#st_crs(calo) in WGS 1984


#Climate Raster
#Annual
MAT.mask <- raster("Donor_selection/data/mask/MAT.mask.grd")
MAP.mask <- raster("Donor_selection/data/mask/MAP.mask.grd")
PAS.mask <- raster("Donor_selection/data/mask/PAS.mask.grd")
EXT.mask <- raster("Donor_selection/data/mask/EXT.mask.grd")
CMD.mask <- raster("Donor_selection/data/mask/CMD.mask.grd")

#Seasonal
PPT_sm.mask <- raster("Donor_selection/data/mask/PPT_sm.mask.grd")
PPT_wt.mask <- raster("Donor_selection/data/mask/PPT_wt.mask.grd")
Tave_sm.mask <- raster("Donor_selection/data/mask/Tave_sm.mask.grd")
Tave_wt.mask <- raster("Donor_selection/data/mask/Tave_wt.mask.grd")

##############################################################################