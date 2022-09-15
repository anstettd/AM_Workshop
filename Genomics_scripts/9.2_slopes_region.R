##################################################################################
## Get SNP slopes over the timeseries
## Author Daniel Anstett
## 
## 
## Last Modified April 11, 2022
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)
library(boot)

#Import timeseries frequencies
freq_MAT <- read_csv("Genomics_scripts/Data/freq_MAT_peakbf2_region.csv")
freq_MAP <- read_csv("Genomics_scripts/Data/freq_MAP_peakbf2_region.csv")
freq_CMD <- read_csv("Genomics_scripts/Data/freq_CMD_peakbf2_region.csv")

###################################################################################
###################################################################################
## Generate SNP slope tables per site

# freq_MAT
freq_MAT_slope<-data.frame()
counter<-1
for (i in 3:dim(freq_MAT)[2]){
  for(j in 1:3){
    chr<-colnames(freq_MAT)[i]
    popSNP <- freq_MAT %>% filter(Site==as.character(unique(freq_MAT$Site)[j])) %>% select(Site,Year,all_of(chr))
    colnames(popSNP)[3]<-"snp_ID"
    
    if(all(is.na(popSNP$snp_ID))==FALSE){
      rSNP <- glm(snp_ID ~ Year, family = binomial, data = popSNP)
      freq_MAT_slope[counter,1]<-unique(popSNP$Site)
      freq_MAT_slope[counter,2]<-chr
      freq_MAT_slope[counter,3]<-rSNP$coefficients[2]
      counter<-counter+1
    } else {
      freq_MAT_slope[counter,1]<-unique(popSNP$Site)
      freq_MAT_slope[counter,2]<-chr
      freq_MAT_slope[counter,3]<-NA
      counter<-counter+1
    }
    
  }
}
colnames(freq_MAT_slope)<-c("Site","snp_ID","Slope")

# freq_MAP
freq_MAP_slope<-data.frame()
counter<-1
for (i in 3:dim(freq_MAP)[2]){
  for(j in 1:3){
    chr<-colnames(freq_MAP)[i]
    popSNP <- freq_MAP %>% filter(Site==as.character(unique(freq_MAP$Site)[j])) %>% select(Site,Year,all_of(chr))
    colnames(popSNP)[3]<-"snp_ID"
    
    if(all(is.na(popSNP$snp_ID))==FALSE){
      rSNP <- glm(snp_ID ~ Year, family = binomial, data = popSNP)
      freq_MAP_slope[counter,1]<-unique(popSNP$Site)
      freq_MAP_slope[counter,2]<-chr
      freq_MAP_slope[counter,3]<-rSNP$coefficients[2]
      counter<-counter+1
    } else {
      freq_MAP_slope[counter,1]<-unique(popSNP$Site)
      freq_MAP_slope[counter,2]<-chr
      freq_MAP_slope[counter,3]<-NA
      counter<-counter+1
    }
    
  }
}
colnames(freq_MAP_slope)<-c("Site","snp_ID","Slope")


# freq_CMD
freq_CMD_slope<-data.frame()
counter<-1
for (i in 3:dim(freq_CMD)[2]){
  for(j in 1:3){
    chr<-colnames(freq_CMD)[i]
    popSNP <- freq_CMD %>% filter(Site==as.character(unique(freq_CMD$Site)[j])) %>% select(Site,Year,all_of(chr))
    colnames(popSNP)[3]<-"snp_ID"
    
    if(all(is.na(popSNP$snp_ID))==FALSE){
      rSNP <- glm(snp_ID ~ Year, family = binomial, data = popSNP)
      freq_CMD_slope[counter,1]<-unique(popSNP$Site)
      freq_CMD_slope[counter,2]<-chr
      freq_CMD_slope[counter,3]<-rSNP$coefficients[2]
      counter<-counter+1
    } else {
      freq_CMD_slope[counter,1]<-unique(popSNP$Site)
      freq_CMD_slope[counter,2]<-chr
      freq_CMD_slope[counter,3]<-NA
      counter<-counter+1
    }
    
  }
}
colnames(freq_CMD_slope)<-c("Site","snp_ID","Slope")

###################################################################################
#Export

write_csv(freq_MAT_slope, "Genomics_scripts/Data/freq_MAT_slope_peakbf2_region.csv")
write_csv(freq_MAP_slope, "Genomics_scripts/Data/freq_MAP_slope_peakbf2_region.csv")
write_csv(freq_CMD_slope, "Genomics_scripts/Data/freq_CMD_slope_peakbf2_region.csv")
