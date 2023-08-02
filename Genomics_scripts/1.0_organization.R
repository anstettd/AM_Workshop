##################################################################################
## Merging of data from sequencing centre (NanuQ) with locality data
## Author Daniel Anstett
## 
##
## Last Modified April 28, 2021
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)
library(car)
library(data.table)


#import NanuQ data
NanuQ <- read.csv("Genomics_scripts/Data/NanuQ_all.csv", header=T)
Dylan <- read.csv("Genomics_scripts/Data/Sequenced_Mimulus_NovaSeq.csv", header=T)
Climate_year <- read.csv("Genomics_scripts/Data/Baseline_Timeseries_climate_Normal_1961_1990Y.csv", header=T)
Climate_season <- read.csv("Genomics_scripts/Data/Baseline_Timeseries_climate_Normal_1961_1990S.csv", header=T)

###################################################################################
#Combine NanuQ with Dylan datasets 

victory <- left_join(NanuQ,Dylan,by=c("Name"="ID")) 
#write.csv(victory, "Genomics_scripts/Data/M_caridnalis_genomics_meta_data.csv")
vic <- victory %>% select(Name,Site,Population,Year)
#write.csv(vic, "Genomics_scripts/Data/victory.csv")
vic1 <- victory %>% select(Site,Year)  

###################################################################################
#Produce list of baseline and non-baseline samples
#non-baseline
vic_non_base <- victory %>% filter(Collection == "Timeseries" | Collection == "Yakima")
vic_non_base_names <- vic_non_base$Filename.Prefix
vic_non_base_names<-as.data.frame(vic_non_base_names)
#write_csv(vic_non_base_names, "Genomics_scripts/Data/non_baseline_names.csv")

#Baseline
vic_base <- victory %>% filter(Collection == "Baseline" | Collection == "Both")
vic_base_names <- vic_base$Filename.Prefix
vic_base_names <- gsub("NS","--sample-name NS", vic_base_names) 
vic_base_names <-as.data.frame(vic_base_names)
#write_csv(vic_base_names, "Genomics_scripts/Data/baseline_names.csv")

#Timeseries
vic_time <- victory %>% filter(Collection == "Timeseries" | Collection == "Both")
vic_time_names <- vic_time$Filename.Prefix
vic_time_names <- gsub("NS","--sample-name NS", vic_time_names) 
vic_time_names <-as.data.frame(vic_time_names)
#write_csv(vic_time_names, "Genomics_scripts/Data/timeseries_names.csv")

#Yakima
vic_yakima <- victory %>% filter(Collection == "Yakima")
vic_yakima_names <- vic_yakima$Filename.Prefix
vic_yakima_names <- gsub("NS","--sample-name NS", vic_yakima_names) 
vic_yakima_names <-as.data.frame(vic_yakima_names)
#write_csv(vic_yakima_names, "Genomics_scripts/Data/yakima_names.csv")

#No Yakima
vic_no_yakima <- victory %>% filter(Collection != "Yakima")
vic_no_yakima_names <- vic_no_yakima$Filename.Prefix
vic_no_yakima_names <- gsub("NS","--sample-name NS", vic_no_yakima_names) 
vic_no_yakima_names <-as.data.frame(vic_no_yakima_names)
#write_csv(vic_no_yakima_names, "Genomics_scripts/Data/no_yakima_names.csv")

#Make files for demography populations not in timeseries
#Coast Fork Willamette 
vic_55 <- victory %>% filter(Paper_ID==55)
vic_55<- vic_55$Filename.Prefix
vic_55<-as.data.frame(vic_55)
write_csv(vic_55, "Genomics_scripts/Data/pop_demo_55.csv")

#Rainbow Pool
vic_29 <- victory %>% filter(Paper_ID==29)
vic_29<- vic_29$Filename.Prefix
vic_29<-as.data.frame(vic_29)
write_csv(vic_29, "Genomics_scripts/Data/pop_demo_29.csv")

#Carlon
vic_28 <- victory %>% filter(Paper_ID==28)
vic_28<- vic_28$Filename.Prefix
vic_28<-as.data.frame(vic_28)
write_csv(vic_28, "Genomics_scripts/Data/pop_demo_28.csv")

#Buck Meadows
vic_27 <- victory %>% filter(Paper_ID==27)
vic_27<- vic_27$Filename.Prefix
vic_27<-as.data.frame(vic_27)
write_csv(vic_27, "Genomics_scripts/Data/pop_demo_27.csv")

#Whitewater Cayon
vic_17 <- victory %>% filter(Paper_ID==17)
vic_17<- vic_17$Filename.Prefix
vic_17<-as.data.frame(vic_17)
write_csv(vic_17, "Genomics_scripts/Data/pop_demo_17.csv")

#Kitchen Creek
vic_15 <- victory %>% filter(Paper_ID==15)
vic_15<- vic_15$Filename.Prefix
vic_15<-as.data.frame(vic_15)
write_csv(vic_15, "Genomics_scripts/Data/pop_demo_15.csv")

#Hauser Creek
vic_14 <- victory %>% filter(Paper_ID==14)
vic_14<- vic_14$Filename.Prefix
vic_14<-as.data.frame(vic_14)
write_csv(vic_14, "Genomics_scripts/Data/pop_demo_14.csv")
  

###################################################################################
#Climate Data

#Combine data
Climate_season <- Climate_season %>% select(Tmax_wt:CMI_at)
climate_90 <- cbind(Climate_year,Climate_season)

#Select only 
##### Annual #####
# MAT = Mean annual temperature (°C)
# MAP = Mean annual precipitation (mm)
# PAS = Precipitation as snow (mm) between August in previous year and July in current year
# EXT = Extreme temperature over 30 years
# CMD = Hargreaves climatic moisture deficit (mm)
##### Seasonal #####
# Tave_wt = Winter mean temperature (°C)
# Tave_sm = Summer mean temperature (°C)
# PPT_wt = Winter precipitation (mm)
# PPT_sm = Summer precipitation (mm)

climate <- climate_90 %>% select(Site_Name,Paper_ID,Latitude,Longitude,Elevation,
                                 MAT,MAP,PAS,EXT,CMD,
                                 Tave_wt,Tave_sm,PPT_wt,PPT_sm)
#remove Yakima
climate <- climate %>% filter(Site_Name!="Yakima, WA")
# order by Paper_ID
climate <- climate[order(climate$Paper_ID),]
write_csv(climate, "Donor_selection/Data/climate_pop.csv")

#make site lat file
site_lat <- climate %>% select(Site_Name:Longitude)
#write_csv(site_lat, "Donor_selection/Data/site_lat.csv",col_names=TRUE)

# Format data for BayPass
#Baseline
no_t_env_baypass <- climate %>% select(MAT,MAP,PAS,EXT,CMD,Tave_wt,Tave_sm,PPT_wt,PPT_sm)
no_t_env_baypass.s <- as.data.frame(scale(no_t_env_baypass))
env_baypass <- transpose(no_t_env_baypass.s)
#write_csv(env_baypass, "Genomics_scripts/Data/env_baseline.csv",col_names=FALSE)

###################################################################################
#Produce list of sample ID and Population
#Baseline
vic_base_ID <- vic_base %>% select(Read.Set.Id,Paper_ID)
vic_base_ID <- vic_base_ID[order(vic_base_ID$Paper_ID),] # order by Paper_ID
#write_csv(vic_base_ID, "Genomics_scripts/Data/baseline_pop_id.csv",col_names=FALSE)

#Timeseries
vic_time_ID <- vic_time %>% select(Read.Set.Id,Paper_ID)
vic_time_ID <- vic_time_ID[order(vic_time_ID$Paper_ID),] # order by Paper_ID
#write_csv(vic_time_ID, "Genomics_scripts/Data/time_pop_id.csv",col_names=FALSE)

#Yakima
vic_yakima_ID <- vic_yakima %>% select(Read.Set.Id,Paper_ID)
#write_csv(vic_yakima_ID, "Genomics_scripts/Data/yakima_pop_id.csv",col_names=FALSE)

#Timeseries pop_ID
time_pop <- vic_time %>% select(Read.Set.Id,Paper_ID,Year)
time_pop <- time_pop[order(time_pop$Paper_ID),] # order by Paper_ID
time_pop_f <- time_pop %>% unite(Pop_ID,Paper_ID,Year,sep="_")  ##### This is wrong
write_tsv(time_pop_f, "Genomics_scripts/Data/timeseries_pop_id.txt",col_names=FALSE)


#Timeseries Geo_Region ID
time_Geo_Region <- vic_time %>% select(Read.Set.Id,Geo_Region,Year)
time_Geo_Region <- time_Geo_Region[order(time_Geo_Region$Geo_Region),] # order by Region
time_region_f <- time_Geo_Region %>% unite(Region_ID,Geo_Region,Year,sep="_")  ##### This is wrong
write_tsv(time_region_f, "Genomics_scripts/Data/timeseries_region_id.txt",col_names=FALSE)

###################################################################################
#Make Site Year Combination Table
vic_site <- unique(vic1$Site)
vic_year <- unique(vic1$Year)
vic_year <- vic_year[order(vic_year)]
site_year <- data.frame()
for (i in 1:length(vic_site)){
  for (j in 1:length(vic_year)){
  input<-tally(vic1 %>% filter (Site==as.character(vic_site[i]),Year==vic_year[j]))
  site_year[i,j]<-input
}}

#Replace "0" with NA
site_year[site_year == 0] <- NA

#Lable Row & Col Names
rownames(site_year) <- vic_site
colnames(site_year) <- vic_year
###################################################################################
#Export

write.csv(site_year, "Genomics_scripts/Data/Mimulus_genomics_site_year2.csv")


