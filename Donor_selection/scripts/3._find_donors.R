##############################################################################
### SCRIPT PURPOSE: Transform future climate at select sites to bioclim

# Modified from Angert et al. 2018, American Naturalist
# Author: Daniel Anstett & Amy Angert
# last update:  June 24 2021
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

#Known points with M. cardinalis
known <- read_csv("SDM/data_files/presences.csv")
known <- known %>% dplyr::select(Longitude,Latitude)
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
known_sf <- st_as_sf(known,coords=c("Longitude","Latitude"), crs=EPSG4326)
#check data is set up properly
ggplot()+ geom_sf(data = known_sf)+ ggtitle("Known M. cardinnalis Populations")

# California & Oregon Map Setup
states<-ne_states(country=c("canada","united states of america"),returnclass= "sf")
calo <- states %>%
  filter(name_en=="Oregon" | name_en=="California" | name_en=="Nevada")
#st_crs(calo) in WGS 1984


#Climate Raster
#Annual
CMD.mask <- raster("Donor_selection/data/mask/CMD.mask.grd")
MAP.mask <- raster("Donor_selection/data/mask/MAP.mask.grd")
MAT.mask <- raster("Donor_selection/data/mask/MAT.mask.grd")
PAS.mask <- raster("Donor_selection/data/mask/PAS.mask.grd")
EXT.mask <- raster("Donor_selection/data/mask/EXT.mask.grd")
#Seasonal
PPT_sm.mask <- raster("Donor_selection/data/mask/PPT_sm.mask.grd")
PPT_wt.mask <- raster("Donor_selection/data/mask/PPT_wt.mask.grd")
Tave_sm.mask <- raster("Donor_selection/data/mask/Tave_sm.mask.grd")
Tave_wt.mask <- raster("Donor_selection/data/mask/Tave_wt.mask.grd")

##############################################################################


#Site S17, Deep Creek

#CMD

#Some climate mitigation
CMD.S17.ssp245<-CMD.mask
CMD.S17.ssp245[CMD.mask<gcc.clim$CMD[8]-100] <-NA
CMD.S17.ssp245[CMD.mask>gcc.clim$CMD[8]+100] <-NA
#plot(CMD.S17.ssp245, main="CMD at ssp235")

tmap_mode("plot")
tm_shape(CMD.mask, bbox=st_bbox(calo)) + #legal boundires
  tm_raster()+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(known_sf)+
  tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)

tm_shape(CMD.S17.ssp245, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = "#FF00FF",legend.show = FALSE)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(known_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)

#use to swtich between interactive mode ("view") and stationary mode ("plot")
tmap_mode("view")
tmap_mode("plot")









#No climate mitigation
CMD.S17.ssp585<-CMD.mask
plot(CMD.S17.ssp585)
CMD.S17.ssp585[CMD.mask<gcc.clim$CMD[8]-100] <-NA
CMD.S17.ssp585[CMD.mask>gcc.clim$CMD[8]+100] <-NA
plot(CMD.S17.ssp585, main="CMD at ssp585")

#MAP
MAP.S17.ssp245<-MAP.mask
plot(MAP.S17.ssp245)
MAP.S17.ssp245[MAP.mask<gcc.clim$MAP[8]-100] <-NA
MAP.S17.ssp245[MAP.mask>gcc.clim$MAP[8]+100] <-NA
plot(MAP.S17.ssp245, main="MAP at ssp235")

MAP.S17.ssp585<-MAP.mask
plot(MAP.S17.ssp585)
MAP.S17.ssp585[MAP.mask<gcc.clim$MAP[8]-100] <-NA
MAP.S17.ssp585[MAP.mask>gcc.clim$MAP[8]+100] <-NA
plot(MAP.S17.ssp585, main="MAP at ssp585")

#MAT
MAT.S17.ssp245<-MAT.mask
plot(MAT.S17.ssp245)
MAT.S17.ssp245[MAT.mask<gcc.clim$MAT[8]-1] <-NA
MAT.S17.ssp245[MAT.mask>gcc.clim$MAT[8]+1] <-NA
plot(MAT.S17.ssp245, main="MAT at ssp235")

MAT.S17.ssp585<-MAT.mask
plot(MAT.S17.ssp585)
MAT.S17.ssp585[MAT.mask<gcc.clim$MAT[8]-1] <-NA
MAT.S17.ssp585[MAT.mask>gcc.clim$MAT[8]+1] <-NA
plot(MAT.S17.ssp585, main="MAT at ssp585")

#PAS
PAS.S17.ssp245<-PAS.mask
plot(PAS.S17.ssp245)
PAS.S17.ssp245[PAS.mask<gcc.clim$PAS[8]-5] <-NA
PAS.S17.ssp245[PAS.mask>gcc.clim$PAS[8]+5] <-NA
plot(PAS.S17.ssp245, main="PAS at ssp235")

PAS.S17.ssp585<-PAS.mask
plot(PAS.S17.ssp585)
PAS.S17.ssp585[PAS.mask<gcc.clim$PAS[8]-5] <-NA
PAS.S17.ssp585[PAS.mask>gcc.clim$PAS[8]+5] <-NA
plot(PAS.S17.ssp585, main="PAS at ssp585")

#EXT
EXT.S17.ssp245<-EXT.mask
plot(EXT.S17.ssp245)
EXT.S17.ssp245[EXT.mask<gcc.clim$EXT[8]-1] <-NA
EXT.S17.ssp245[EXT.mask>gcc.clim$EXT[8]+1] <-NA
plot(EXT.S17.ssp245, main="EXT at ssp235")

EXT.S17.ssp585<-EXT.mask
plot(EXT.S17.ssp585)
EXT.S17.ssp585[EXT.mask<gcc.clim$EXT[8]-1] <-NA
EXT.S17.ssp585[EXT.mask>gcc.clim$EXT[8]+1] <-NA
plot(EXT.S17.ssp585, main="EXT at ssp585")

#Tave_sm
Tave_sm.S17.ssp245<-Tave_sm.mask
plot(Tave_sm.S17.ssp245)
Tave_sm.S17.ssp245[Tave_sm.mask<gcc.clim$Tave_sm[8]-1] <-NA
Tave_sm.S17.ssp245[Tave_sm.mask>gcc.clim$Tave_sm[8]+1] <-NA
plot(Tave_sm.S17.ssp245, main="Tave_sm at ssp235")

Tave_sm.S17.ssp585<-Tave_sm.mask
plot(Tave_sm.S17.ssp585)
Tave_sm.S17.ssp585[Tave_sm.mask<gcc.clim$Tave_sm[8]-1] <-NA
Tave_sm.S17.ssp585[Tave_sm.mask>gcc.clim$Tave_sm[8]+1] <-NA
plot(Tave_sm.S17.ssp585, main="Tave_sm at ssp585")

#Tave_wt
Tave_wt.S17.ssp245<-Tave_wt.mask
plot(Tave_wt.S17.ssp245)
Tave_wt.S17.ssp245[Tave_wt.mask<gcc.clim$Tave_wt[8]-1] <-NA
Tave_wt.S17.ssp245[Tave_wt.mask>gcc.clim$Tave_wt[8]+1] <-NA
plot(Tave_wt.S17.ssp245, main="Tave_wt at ssp235")

Tave_wt.S17.ssp585<-Tave_wt.mask
plot(Tave_wt.S17.ssp585)
Tave_wt.S17.ssp585[Tave_wt.mask<gcc.clim$Tave_wt[8]-1] <-NA
Tave_wt.S17.ssp585[Tave_wt.mask>gcc.clim$Tave_wt[8]+1] <-NA
plot(Tave_wt.S17.ssp585, main="Tave_wt at ssp585")

#PPT_sm
PPT_sm.S17.ssp245<-PPT_sm.mask
plot(PPT_sm.S17.ssp245)
PPT_sm.S17.ssp245[PPT_sm.mask<gcc.clim$PPT_sm[8]-10] <-NA
PPT_sm.S17.ssp245[PPT_sm.mask>gcc.clim$PPT_sm[8]+10] <-NA
plot(PPT_sm.S17.ssp245, main="PPT_sm at ssp235")

PPT_sm.S17.ssp585<-PPT_sm.mask
plot(PPT_sm.S17.ssp585)
PPT_sm.S17.ssp585[PPT_sm.mask<gcc.clim$PPT_sm[8]-10] <-NA
PPT_sm.S17.ssp585[PPT_sm.mask>gcc.clim$PPT_sm[8]+10] <-NA
plot(PPT_sm.S17.ssp585, main="PPT_sm at ssp585")

#PPT_wt
PPT_wt.S17.ssp245<-PPT_wt.mask
plot(PPT_wt.S17.ssp245)
PPT_wt.S17.ssp245[PPT_wt.mask<gcc.clim$PPT_wt[8]-100] <-NA
PPT_wt.S17.ssp245[PPT_wt.mask>gcc.clim$PPT_wt[8]+100] <-NA
plot(PPT_wt.S17.ssp245, main="PPT_wt at ssp235")

PPT_wt.S17.ssp585<-PPT_wt.mask
plot(PPT_wt.S17.ssp585)
PPT_wt.S17.ssp585[PPT_wt.mask<gcc.clim$PPT_wt[8]-100] <-NA
PPT_wt.S17.ssp585[PPT_wt.mask>gcc.clim$PPT_wt[8]+100] <-NA
plot(PPT_wt.S17.ssp585, main="PPT_wt at ssp585")













