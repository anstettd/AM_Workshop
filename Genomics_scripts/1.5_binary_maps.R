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
geography_km <- as.data.frame(geography/1000)

#Import abundance results
freq_1 <- read_csv("Genomics_scripts/Data/freq_1.csv")
freq_2 <- read_csv("Genomics_scripts/Data/freq_2.csv")
freq_3 <- read_csv("Genomics_scripts/Data/freq_3.csv")
freq_4 <- read_csv("Genomics_scripts/Data/freq_4.csv")
freq_5 <- read_csv("Genomics_scripts/Data/freq_5.csv")
freq_6 <- read_csv("Genomics_scripts/Data/freq_6.csv")
freq_7 <- read_csv("Genomics_scripts/Data/freq_7.csv")
freq_8 <- read_csv("Genomics_scripts/Data/freq_8.csv")
freq_9 <- read_csv("Genomics_scripts/Data/freq_9.csv")

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

#Define CRS WGS 1984
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS

###################################################################################
###################################################################################
#Map distribution of single snps

# Map of P1, present in most localities
b1_snp2 <- freq_1 %>% filter(SNP==224717) #select abundance data for snp 2
#b1_snp2 <- binary_1 %>% filter(SNP==224717) #select binary data for snp 2, remove comment to see binary data
b1_snp2 <- b1_snp2[5:59] #take only snps values
b1_snp2 <- as.data.frame(t(b1_snp2)) #transpose and put in dataframe
pop_var_snp2 <- cbind(pop_var_raw,b1_snp2) #bind to lat/long population data
colnames(pop_var_snp2)[6] <- "MAT" 

MAT_points_snp2 <- pop_var_snp2 %>% dplyr::select(Long,Lat,MAT) #select relevant data
MAT_sf_snp2 <- st_as_sf(MAT_points_snp2,coords=c("Long","Lat"), crs=EPSG4326)

#Plot MAT associated allele snp2
tmap_mode("plot")
#tmap_mode("view")
snp2 <- tm_shape(calo)+
  tm_borders()+
  tm_shape(MAT_sf_snp2)+
  tm_bubbles(size = 0.15,col="MAT")+ 
  #tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
snp2
tmap_save(snp2, filename = "snp2_freq.pdf",width=5, hight=6)

# Map of P14, missing from most localities
b1_snp14 <- freq_1 %>% filter(SNP==477379) #filter for snp
b1_snp14 <- b1_snp14[5:59] #take only snps values
b1_snp14 <- as.data.frame(t(b1_snp14)) #transpose and put in dataframe
pop_var_snp14 <- cbind(pop_var_raw,b1_snp14) #bind to lat/long population data
colnames(pop_var_snp14)[6] <- "MAT" 

MAT_points_snp14 <- pop_var_snp14 %>% dplyr::select(Long,Lat,MAT) #select relevant data
MAT_sf_snp14 <- st_as_sf(MAT_points_snp14,coords=c("Long","Lat"), crs=EPSG4326)

#Plot MAT associated allele snp14
tmap_mode("plot")
tmap_mode("view")
tm_shape(calo)+
  tm_borders()+
  tm_shape(MAT_sf_snp14)+
  tm_bubbles(size = 0.15,col="MAT")+ 
  #tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)


# Map of P25, missing from most localities
b1_snp25 <- freq_1 %>% filter(SNP==946612) #filter for snp
b1_snp25 <- b1_snp25[5:59] #take only snps values
b1_snp25 <- as.data.frame(t(b1_snp25)) #transpose and put in dataframe
pop_var_snp25 <- cbind(pop_var_raw,b1_snp25) #bind to lat/long population data
colnames(pop_var_snp25)[6] <- "MAT" 

MAT_points_snp25 <- pop_var_snp25 %>% dplyr::select(Long,Lat,MAT) #select relevant data
MAT_sf_snp25 <- st_as_sf(MAT_points_snp25,coords=c("Long","Lat"), crs=EPSG4326)

#Plot MAT associated allele snp25
tmap_mode("plot")
tmap_mode("view")
tm_shape(calo)+
  tm_borders()+
  tm_shape(MAT_sf_snp25)+
  tm_bubbles(size = 0.15,col="MAT")+ 
  tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)


###################################################################################
# Map Proportion of putatively adaptive SNPs Present

#Set up sf objects for snp proportions

#Annual
#MAT
MAT_points <- pop_var %>% dplyr::select(Long,Lat,MAT)
MAT_sf <- st_as_sf(MAT_points,coords=c("Long","Lat"), crs=EPSG4326)
ggplot()+ geom_sf(data = MAT_sf)+ ggtitle("MAT points") #check data is set up properly

#MAP
MAP_points <- pop_var %>% dplyr::select(Long,Lat,MAP)
MAP_sf <- st_as_sf(MAP_points,coords=c("Long","Lat"), crs=EPSG4326)
ggplot()+ geom_sf(data = MAP_sf)+ ggtitle("MAP points") #check data is set up properly

#CMD
CMD_points <- pop_var %>% dplyr::select(Long,Lat,CMD)
CMD_sf <- st_as_sf(CMD_points,coords=c("Long","Lat"), crs=EPSG4326)
ggplot()+ geom_sf(data = CMD_sf)+ ggtitle("CMD points") #check data is set up properly

#Make Maps 
#MAT
tmap_mode("plot")
tmap_mode("view")
tm_shape(calo)+
  tm_borders()+
  tm_shape(MAT_sf)+
  tm_bubbles(size = 0.15,col="MAT")+ 
  #tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)

#MAP
tmap_mode("plot")
tmap_mode("view")
tm_shape(calo)+
  tm_borders()+
  tm_shape(MAP_sf)+
  tm_bubbles(size = 0.15,col="MAP")+ 
  #tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)

#CMD
tmap_mode("plot")
tmap_mode("view")
tm_shape(calo)+
  tm_borders()+
  tm_shape(CMD_sf)+
  tm_bubbles(size = 0.15,col="CMD")+ 
  #tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)




###################################################################################

##plot missing variation for P8 Deep Creek

##Setup P8 MAT
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


##Setup P8 MAP
binary2_miss_p8 <- binary_2 %>% filter(P8==0) #select only snps not in P8
PSP_2_p8 <- as.data.frame(colMeans(binary2_miss_p8[5:59],na.rm = TRUE))
pop_2var_p8 <- cbind(pop_var_raw,PSP_2_p8)
colnames(pop_2var_p8)[6] <- "MAP" 

#Population values for p8
MAP_points_p8 <- pop_2var_p8 %>% dplyr::select(Long,Lat,MAP)
MAP_sf_p8 <- st_as_sf(MAP_points_p8,coords=c("Long","Lat"), crs=EPSG4326)

#MAP
tmap_mode("plot")
tmap_mode("view")
tm_shape(calo)+
  tm_borders()+
  tm_shape(MAP_sf_p8)+
  tm_bubbles(size = 0.15,col="MAP")+ 
  #tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)


##Setup P8 CMD
binary5_miss_p8 <- binary_5 %>% filter(P8==0) #select only snps not in P8
PSP_5_p8 <- as.data.frame(colMeans(binary5_miss_p8[5:59],na.rm = TRUE))
pop_5var_p8 <- cbind(pop_var_raw,PSP_5_p8)
colnames(pop_5var_p8)[6] <- "CMD" 

#Population values for p8
CMD_points_p8 <- pop_5var_p8 %>% dplyr::select(Long,Lat,CMD)
CMD_sf_p8 <- st_as_sf(CMD_points_p8,coords=c("Long","Lat"), crs=EPSG4326)

#CMD
tmap_mode("plot")
tmap_mode("view")
tm_shape(calo)+
  tm_borders()+
  tm_shape(CMD_sf_p8)+
  tm_bubbles(size = 0.15,col="CMD")+ 
  #tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)


###################################################################################
#Distance to punitively adaptive snp

##MAT
dist_p8_MAT <- data.frame() #Generate empty data frame

#select shortest distance to snp punitively adaptive to MAT 
for (i in 1:nrow(binary_1)){
  geography_p8 <- geography_km[,8] #get distance from p8
  b1_t <- as.data.frame(t(binary_1[i,5:59])) #get line i from binary matrix
  geo_b1 <- as.data.frame(cbind(geography_p8,b1_t)) #bind columns together
  colnames(geo_b1) <- c("geo","snp")  #rename
  geo_b1 <-  geo_b1 %>% filter(snp==1) #select for presence of "adaptive" snps
  geo_b1 <-  geo_b1 %>% filter(geo==min(geo)) #select for minimum distance
  dist_p8_MAT[i,1] <- geo_b1[1,1] #write out in distance for i snp
}
colnames(dist_p8_MAT) <- "Distance"

#Make graphic
ggplot(data=dist_p8_MAT, aes(dist_p8_MAT$Distance))+
  geom_histogram(binwidth = 10)+
  theme_classic()


##MAP
dist_p8_MAP <- data.frame() #Generate empty data frame

#select shortest distance to snp punitively adaptive to MAP 
for (i in 1:nrow(binary_2)){
  geography_p8 <- geography_km[,8] #get distance from p8
  b2_t <- as.data.frame(t(binary_2[i,5:59])) #get line i from binary matrix
  geo_b2 <- as.data.frame(cbind(geography_p8,b2_t)) #bind columns together
  colnames(geo_b2) <- c("geo","snp")  #rename
  geo_b2 <-  geo_b2 %>% filter(snp==1) #select for presence of "adaptive" snps
  geo_b2 <-  geo_b2 %>% filter(geo==min(geo)) #select for minimum distance
  dist_p8_MAP[i,1] <- geo_b2[1,1] #write out in distance for i snp
}
colnames(dist_p8_MAP) <- "Distance"

#Make graphic
ggplot(data=dist_p8_MAP, aes(dist_p8_MAP$Distance))+
  geom_histogram(binwidth = 10)+
  theme_classic()



##CMD
dist_p8_CMD <- data.frame() #Generate empty data frame

#select shortest distance to snp punitively adaptive to CMD 
for (i in 1:nrow(binary_5)){
  geography_p8 <- geography_km[,8] #get distance from p8
  b5_t <- as.data.frame(t(binary_5[i,5:59])) #get line i from binary matrix
  geo_b5 <- as.data.frame(cbind(geography_p8,b5_t)) #bind columns together
  colnames(geo_b5) <- c("geo","snp")  #rename
  geo_b5 <-  geo_b5 %>% filter(snp==1) #select for presence of "adaptive" snps
  geo_b5 <-  geo_b5 %>% filter(geo==min(geo)) #select for minimum distance
  dist_p8_CMD[i,1] <- geo_b5[1,1] #write out in distance for i snp
}
colnames(dist_p8_CMD) <- "Distance"

#Make graphic
ggplot(data=dist_p8_CMD, aes(dist_p8_CMD$Distance))+
  geom_histogram(binwidth = 10)+
  theme_classic()

###################################################################################







