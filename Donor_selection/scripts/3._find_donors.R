##############################################################################
### SCRIPT PURPOSE: Find donor sites using only climate

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
MAT.clip <- raster("Donor_selection/data/clip/MAT.clip.grd")
MAP.clip <- raster("Donor_selection/data/clip/MAP.clip.grd")
PAS.clip <- raster("Donor_selection/data/clip/PAS.clip.grd")
EXT.clip <- raster("Donor_selection/data/clip/EXT.clip.grd")
CMD.clip <- raster("Donor_selection/data/clip/CMD.clip.grd")

#Seasonal
PPT_sm.clip <- raster("Donor_selection/data/clip/PPT_sm.clip.grd")
PPT_wt.clip <- raster("Donor_selection/data/clip/PPT_wt.clip.grd")
Tave_sm.clip <- raster("Donor_selection/data/clip/Tave_sm.clip.grd")
Tave_wt.clip <- raster("Donor_selection/data/clip/Tave_wt.clip.grd")

##############################################################################


#Site S17, Deep Creek

#MAT
MAT.S17.ssp245<-MAT.clip
plot(MAT.S17.ssp245)
MAT.S17.ssp245[MAT.clip<gcc.clim$MAT[8]-1] <-NA
MAT.S17.ssp245[MAT.clip>gcc.clim$MAT[8]+1] <-NA
plot(MAT.S17.ssp245, main="MAT at ssp235")

MAT.S17.ssp585<-MAT.clip
plot(MAT.S17.ssp585)
MAT.S17.ssp585[MAT.clip<gcc.clim$MAT[8]-1] <-NA
MAT.S17.ssp585[MAT.clip>gcc.clim$MAT[8]+1] <-NA
plot(MAT.S17.ssp585, main="MAT at ssp585")

##plot maps

#Climate Layer
tmap_mode("plot")
#tmap_mode("view")
clim_MAT <- tm_shape(MAT.clip, bbox=st_bbox(calo)) + #legal boundires
  tm_raster()+
  tm_shape(calo)+
  tm_borders()+
  #tm_shape(gen_pop_sf)+
  #tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
MAT_clim_layer
tmap_save(clim_MAT, filename = "Graphs/clim_MAT.pdf",width=5, height=6)

#ssp245 matching locations only
tmap_mode("plot")
#tmap_mode("view")
clim_MAT_245 <- tm_shape(MAT.S17.ssp245, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = "#FF00FF",legend.show = FALSE)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
clim_MAT_245
tmap_save(clim_MAT_245, filename = "Graphs/clim_MAT_245.pdf",width=5, height=6)

#ssp585 matching locations only
tmap_mode("plot")
#tmap_mode("view")
clim_MAT_585 <- tm_shape(MAT.S17.ssp585, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = "#FF00FF",legend.show = FALSE)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
clim_MAT_585
tmap_save(clim_MAT_585, filename = "Graphs/clim_MAT_585.pdf",width=5, height=6)

#Additional Exploration
#Use in presentation

#ssp245 not pops
tmap_mode("plot")
#tmap_mode("view")
clim_MAT_245_no_pop <- tm_shape(MAT.S17.ssp245, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = "#FF00FF",legend.show = FALSE)+
  tm_shape(calo)+
  tm_borders()+
  #tm_shape(gen_pop_sf)+
  #tm_dots(size=0.2,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
clim_MAT_245_no_pop
tmap_save(clim_MAT_245_no_pop, filename = "Graphs/clim_MAT_245_no_pop.pdf",width=5, height=6)

#Try 0.5 deg versus 1 deg

MAT.S17.ssp245_0.5 <-MAT.clip
MAT.S17.ssp245_0.5[MAT.clip<gcc.clim$MAT[8]-0.5] <-NA
MAT.S17.ssp245_0.5[MAT.clip>gcc.clim$MAT[8]-0.5] <-NA
plot(MAT.S17.ssp245_0.5, main="MAT at ssp235")

#ssp245 not pops
tmap_mode("plot")
#tmap_mode("view")
clim_MAT_245_0.5 <- tm_shape(MAT.S17.ssp245_0.5, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = "#FF00FF",legend.show = FALSE)+
  tm_shape(calo)+
  tm_borders()+
  #tm_shape(gen_pop_sf)+
  #tm_dots(size=0.2,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
clim_MAT_245_0.5 
tmap_save(clim_MAT_245_0.5 , filename = "Graphs/clim_MAT_245_0.5.pdf",width=5, height=6)

#entire_range using SDM

MAT.mask <- raster("Donor_selection/data/mask/MAT.mask.grd")
#Climate Layer masked by SDM
tmap_mode("plot")
#tmap_mode("view")
smd_ver2 <- tm_shape(MAT.mask, bbox=st_bbox(calo)) + #legal boundires
  tm_raster()+
  tm_shape(calo)+
  tm_borders()+
  #tm_shape(gen_pop_sf)+
  #tm_dots(size=0.1,shape=1)+
  tm_layout(legend.show=FALSE)
smd_ver2
tmap_save(smd_ver2, filename = "Graphs/smd_ver2.pdf",width=5, height=6)



#CMD

#Some climate mitigation
CMD.S17.ssp245<-CMD.clip
CMD.S17.ssp245[CMD.clip<gcc.clim$CMD[8]-100] <-NA
CMD.S17.ssp245[CMD.clip>gcc.clim$CMD[8]+100] <-NA
#plot(CMD.S17.ssp245, main="CMD at ssp235")

#No climate mitigation
CMD.S17.ssp585<-CMD.clip
plot(CMD.S17.ssp585)
CMD.S17.ssp585[CMD.clip<gcc.clim$CMD[8]-100] <-NA
CMD.S17.ssp585[CMD.clip>gcc.clim$CMD[8]+100] <-NA
#plot(CMD.S17.ssp585, main="CMD at ssp585")

#use to swtich between interactive mode ("view") and stationary mode ("plot")
tmap_mode("plot")
tmap_mode("view")

#plot maps
#on sdm
tm_shape(CMD.clip, bbox=st_bbox(calo)) + #legal boundires
  tm_raster()+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)

#some mitigation
tm_shape(CMD.S17.ssp245, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = "#FF00FF",legend.show = FALSE)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)

#no mitigation
tm_shape(CMD.S17.ssp, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = "#FF00FF",legend.show = FALSE)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)



#MAP
MAP.S17.ssp245<-MAP.clip
plot(MAP.S17.ssp245)
MAP.S17.ssp245[MAP.clip<gcc.clim$MAP[8]-100] <-NA
MAP.S17.ssp245[MAP.clip>gcc.clim$MAP[8]+100] <-NA
plot(MAP.S17.ssp245, main="MAP at ssp235")

MAP.S17.ssp585<-MAP.clip
plot(MAP.S17.ssp585)
MAP.S17.ssp585[MAP.clip<gcc.clim$MAP[8]-100] <-NA
MAP.S17.ssp585[MAP.clip>gcc.clim$MAP[8]+100] <-NA
plot(MAP.S17.ssp585, main="MAP at ssp585")

#use to swtich between interactive mode ("view") and stationary mode ("plot")
tmap_mode("plot")
tmap_mode("view")

#plot maps
#on sdm
tm_shape(MAP.clip, bbox=st_bbox(calo)) + #legal boundires
  tm_raster()+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)

#some mitigation
tm_shape(MAP.S17.ssp245, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = "#FF00FF",legend.show = FALSE)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)

#no mitigation
tm_shape(MAP.S17.ssp585, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = "#FF00FF",legend.show = FALSE)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)



#PAS
PAS.S17.ssp245<-PAS.clip
plot(PAS.S17.ssp245)
PAS.S17.ssp245[PAS.clip<gcc.clim$PAS[8]-5] <-NA
PAS.S17.ssp245[PAS.clip>gcc.clim$PAS[8]+5] <-NA
plot(PAS.S17.ssp245, main="PAS at ssp235")

PAS.S17.ssp585<-PAS.clip
plot(PAS.S17.ssp585)
PAS.S17.ssp585[PAS.clip<gcc.clim$PAS[8]-5] <-NA
PAS.S17.ssp585[PAS.clip>gcc.clim$PAS[8]+5] <-NA
plot(PAS.S17.ssp585, main="PAS at ssp585")

#EXT
EXT.S17.ssp245<-EXT.clip
plot(EXT.S17.ssp245)
EXT.S17.ssp245[EXT.clip<gcc.clim$EXT[8]-1] <-NA
EXT.S17.ssp245[EXT.clip>gcc.clim$EXT[8]+1] <-NA
plot(EXT.S17.ssp245, main="EXT at ssp235")

EXT.S17.ssp585<-EXT.clip
plot(EXT.S17.ssp585)
EXT.S17.ssp585[EXT.clip<gcc.clim$EXT[8]-1] <-NA
EXT.S17.ssp585[EXT.clip>gcc.clim$EXT[8]+1] <-NA
plot(EXT.S17.ssp585, main="EXT at ssp585")

#Tave_sm
Tave_sm.S17.ssp245<-Tave_sm.clip
plot(Tave_sm.S17.ssp245)
Tave_sm.S17.ssp245[Tave_sm.clip<gcc.clim$Tave_sm[8]-1] <-NA
Tave_sm.S17.ssp245[Tave_sm.clip>gcc.clim$Tave_sm[8]+1] <-NA
plot(Tave_sm.S17.ssp245, main="Tave_sm at ssp235")

Tave_sm.S17.ssp585<-Tave_sm.clip
plot(Tave_sm.S17.ssp585)
Tave_sm.S17.ssp585[Tave_sm.clip<gcc.clim$Tave_sm[8]-1] <-NA
Tave_sm.S17.ssp585[Tave_sm.clip>gcc.clim$Tave_sm[8]+1] <-NA
plot(Tave_sm.S17.ssp585, main="Tave_sm at ssp585")

#Tave_wt
Tave_wt.S17.ssp245<-Tave_wt.clip
plot(Tave_wt.S17.ssp245)
Tave_wt.S17.ssp245[Tave_wt.clip<gcc.clim$Tave_wt[8]-1] <-NA
Tave_wt.S17.ssp245[Tave_wt.clip>gcc.clim$Tave_wt[8]+1] <-NA
plot(Tave_wt.S17.ssp245, main="Tave_wt at ssp235")

Tave_wt.S17.ssp585<-Tave_wt.clip
plot(Tave_wt.S17.ssp585)
Tave_wt.S17.ssp585[Tave_wt.clip<gcc.clim$Tave_wt[8]-1] <-NA
Tave_wt.S17.ssp585[Tave_wt.clip>gcc.clim$Tave_wt[8]+1] <-NA
plot(Tave_wt.S17.ssp585, main="Tave_wt at ssp585")

#PPT_sm
PPT_sm.S17.ssp245<-PPT_sm.clip
plot(PPT_sm.S17.ssp245)
PPT_sm.S17.ssp245[PPT_sm.clip<gcc.clim$PPT_sm[8]-10] <-NA
PPT_sm.S17.ssp245[PPT_sm.clip>gcc.clim$PPT_sm[8]+10] <-NA
plot(PPT_sm.S17.ssp245, main="PPT_sm at ssp235")

PPT_sm.S17.ssp585<-PPT_sm.clip
plot(PPT_sm.S17.ssp585)
PPT_sm.S17.ssp585[PPT_sm.clip<gcc.clim$PPT_sm[8]-10] <-NA
PPT_sm.S17.ssp585[PPT_sm.clip>gcc.clim$PPT_sm[8]+10] <-NA
plot(PPT_sm.S17.ssp585, main="PPT_sm at ssp585")

#PPT_wt
PPT_wt.S17.ssp245<-PPT_wt.clip
plot(PPT_wt.S17.ssp245)
PPT_wt.S17.ssp245[PPT_wt.clip<gcc.clim$PPT_wt[8]-100] <-NA
PPT_wt.S17.ssp245[PPT_wt.clip>gcc.clim$PPT_wt[8]+100] <-NA
plot(PPT_wt.S17.ssp245, main="PPT_wt at ssp235")

PPT_wt.S17.ssp585<-PPT_wt.clip
plot(PPT_wt.S17.ssp585)
PPT_wt.S17.ssp585[PPT_wt.clip<gcc.clim$PPT_wt[8]-100] <-NA
PPT_wt.S17.ssp585[PPT_wt.clip>gcc.clim$PPT_wt[8]+100] <-NA
plot(PPT_wt.S17.ssp585, main="PPT_wt at ssp585")













