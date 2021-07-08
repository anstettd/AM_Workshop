##################################################################################
## Make % prescence maps of each climate assocaited snp
## Author Daniel Anstett 
## 
##
## Last Modified July 8, 2021
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)
library(cowplot)
library(sf)
library(tmap)
library(rnaturalearth)
library(rnaturalearthdata)
#devtools::install_github("ropensci/rnaturalearthhires")
library(rnaturalearthhires)
library(rgeos)
library(geodist)

###################################################################################
#Import population data
pop_var_raw <- read_csv("Genomics_scripts/Data/paper_ID_site_select.csv")
deg_distance <- pop_var_raw %>% dplyr::select(Long,Lat)
geography <- geodist(deg_distance)
geography_km <- geography/1000

#Import presence/abscence results
binary_1 <- read_csv("Genomics_scripts/Data/freq_binary_1.csv")
binary_2 <- read_csv("Genomics_scripts/Data/freq_binary_2.csv")
binary_3 <- read_csv("Genomics_scripts/Data/freq_binary_3.csv")
binary_4 <- read_csv("Genomics_scripts/Data/freq_binary_4.csv")
binary_5 <- read_csv("Genomics_scripts/Data/freq_binary_5.csv")
binary_6 <- read_csv("Genomics_scripts/Data/freq_binary_6.csv")
binary_7 <- read_csv("Genomics_scripts/Data/freq_binary_7.csv")
binary_8 <- read_csv("Genomics_scripts/Data/freq_binary_8.csv")
binary_9 <- read_csv("Genomics_scripts/Data/freq_binary_9.csv")

# Get Proportion of SNPs Present
psp_1 <- as.data.frame(colMeans(binary_1[5:59],na.rm = TRUE))
psp_2 <- as.data.frame(colMeans(binary_2[5:59],na.rm = TRUE))
psp_3 <- as.data.frame(colMeans(binary_3[5:59],na.rm = TRUE))
psp_4 <- as.data.frame(colMeans(binary_4[5:59],na.rm = TRUE))
psp_5 <- as.data.frame(colMeans(binary_5[5:59],na.rm = TRUE))
psp_6 <- as.data.frame(colMeans(binary_6[5:59],na.rm = TRUE))
psp_7 <- as.data.frame(colMeans(binary_7[5:59],na.rm = TRUE))
psp_8 <- as.data.frame(colMeans(binary_8[5:59],na.rm = TRUE))
psp_9 <- as.data.frame(colMeans(binary_9[5:59],na.rm = TRUE))

#put Proportion of SNPs Present into population dataframe
pop_var <- cbind(pop_var_raw,psp_1,psp_2,psp_3,psp_4,psp_5,psp_6,psp_7,psp_8,psp_9)
colnames(pop_var)[6:14] <- c("MAT","MAP","PAS","EXT","CMD","Tave_wt","Tave_sm","PPT_wt","PPT_sm") 
###################################################################################
#Mapping setup

# California & Oregon Map Setup
states<-ne_states(country=c("canada","united states of america"),returnclass= "sf")
calo <- states %>%
  filter(name_en=="Oregon" | name_en=="California" | name_en=="Nevada")
#st_crs(calo) in WGS 1984


#Population values
MAT_points <- pop_var %>% dplyr::select(Long,Lat,MAT)
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
MAT_sf <- st_as_sf(MAT_points,coords=c("Long","Lat"), crs=EPSG4326)
ggplot()+ geom_sf(data = MAT_sf)+ ggtitle("CMD points") #check data is set up properly

###################################################################################
###################################################################################

# Map Proportion of putatively adaptive SNPs Present

#MAT
tmap_mode("plot")
#tmap_mode("view")
tm_shape(calo)+
  tm_borders()+
  tm_shape(MAT_sf)+
  tm_bubbles(size = 0.15,col="MAT")+ 
  #tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)


###################################################################################

#plot missing variation for P8 Deep Creek
binary1_miss_p8 <- binary_1 %>% filter(P8==0) #select only snps not in P8
psp_1_p8 <- as.data.frame(colMeans(binary1_miss_p8[5:59],na.rm = TRUE))
pop_var_p8 <- cbind(pop_var_raw,psp_1_p8)
colnames(pop_var_p8)[6] <- "MAT" 

#Population values for p8
MAT_points_p8 <- pop_var_p8 %>% dplyr::select(Long,Lat,MAT)
MAT_sf_p8 <- st_as_sf(MAT_points_p8,coords=c("Long","Lat"), crs=EPSG4326)

#MAT
tmap_mode("plot")
tmap_mode("view")
tm_shape(calo)+
  tm_borders()+
  tm_shape(MAT_sf_p8)+
  tm_bubbles(size = 0.15,col="MAT")+ 
  #tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)

###################################################################################
#Distance to putatively adaptive snp









