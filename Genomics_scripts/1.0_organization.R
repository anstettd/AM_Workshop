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

#import NanuQ data
NanuQ <- read.csv("Genomics_scripts/Data/NanuQ_all.csv", header=T)
Dylan <- read.csv("Genomics_scripts/Data/Final_Mimulus_NovaSeq_150_ID_List_Angert_2020_30X.csv", header=T)

###################################################################################
#Combine NanuQ with Dylan datasets 

victory <- left_join(NanuQ,Dylan,by=c("Name"="ID")) 
vic <- victory %>% select(Name,Site,Population,Year)
vic1 <- victory %>% select(Site,Year)  



###################################################################################
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

write.csv(site_year, "Genomics_scripts/Data/Mimulus_genomics_site_year.csv")


