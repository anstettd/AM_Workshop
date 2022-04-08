##################################################################################
## Get SNP slopes over the timeseries
## Author Daniel Anstett
## 
## 
## Last Modified April 8, 2022
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

#Import timeseries frequencies
freq_MAT <- read_csv("Genomics_scripts/Data/freq_MAT.csv")
freq_MAP <- read_csv("Genomics_scripts/Data/freq_MAP.csv")
freq_CMD <- read_csv("Genomics_scripts/Data/freq_CMD.csv")


###################################################################################
###################################################################################
## Generate SNP slope tables per site

# snp1A
snp1A_slope<-data.frame()
counter<-1
for (i in 3:dim(snp1A)[2]){
  for(j in 1:12){
    chr<-colnames(snp1A)[i]
    popSNP <- snp1A %>% filter(Site==j) %>% select(Site,Year,all_of(chr))
    colnames(popSNP)[3]<-"snp_ID"
    
    if(all(is.na(popSNP$snp_ID))==FALSE){
      rSNP <- lm(snp_ID ~ Year, data = popSNP)
      snp1A_slope[counter,1]<-unique(popSNP$Site)
      snp1A_slope[counter,2]<-chr
      snp1A_slope[counter,3]<-abs(rSNP$coefficients[2])
      counter<-counter+1
    } else {
      snp1A_slope[counter,1]<-unique(popSNP$Site)
      snp1A_slope[counter,2]<-chr
      snp1A_slope[counter,3]<-NA
      counter<-counter+1
    }
    
  }
}
colnames(snp1A_slope)<-c("Site","snp_ID","Slope")

# snp2A
snp2A_slope<-data.frame()
counter<-1
for (i in 3:dim(snp2A)[2]){
  for(j in 1:12){
    chr<-colnames(snp2A)[i]
    popSNP <- snp2A %>% filter(Site==j) %>% select(Site,Year,all_of(chr))
    colnames(popSNP)[3]<-"snp_ID"
    
    if(all(is.na(popSNP$snp_ID))==FALSE){
      rSNP <- lm(snp_ID ~ Year, data = popSNP)
      snp2A_slope[counter,1]<-unique(popSNP$Site)
      snp2A_slope[counter,2]<-chr
      snp2A_slope[counter,3]<-abs(rSNP$coefficients[2])
      counter<-counter+1
    } else {
      snp2A_slope[counter,1]<-unique(popSNP$Site)
      snp2A_slope[counter,2]<-chr
      snp2A_slope[counter,3]<-NA
      counter<-counter+1
    }
    
  }
}
colnames(snp2A_slope)<-c("Site","snp_ID","Slope")


# snp5A
snp5A_slope<-data.frame()
counter<-1
for (i in 3:dim(snp5A)[2]){
  for(j in 1:12){
    chr<-colnames(snp5A)[i]
    popSNP <- snp5A %>% filter(Site==j) %>% select(Site,Year,all_of(chr))
    colnames(popSNP)[3]<-"snp_ID"
    
    if(all(is.na(popSNP$snp_ID))==FALSE){
      rSNP <- lm(snp_ID ~ Year, data = popSNP)
      snp5A_slope[counter,1]<-unique(popSNP$Site)
      snp5A_slope[counter,2]<-chr
      snp5A_slope[counter,3]<-abs(rSNP$coefficients[2])
      counter<-counter+1
    } else {
      snp5A_slope[counter,1]<-unique(popSNP$Site)
      snp5A_slope[counter,2]<-chr
      snp5A_slope[counter,3]<-NA
      counter<-counter+1
    }
    
  }
}
colnames(snp5A_slope)<-c("Site","snp_ID","Slope")