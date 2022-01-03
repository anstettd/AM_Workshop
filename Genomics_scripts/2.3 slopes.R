##################################################################################
## Get SNP slopes
## Author Daniel Anstett
## 
## Currently done just for ENV 1,2 and 5 (MAT, MAP and CMD)
## Last Modified Nov 16, 2021
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

snp1A <- read_csv("Genomics_scripts/Data/snp1A.csv")
snp2A <- read_csv("Genomics_scripts/Data/snp2A.csv")
snp5A <- read_csv("Genomics_scripts/Data/snp5A.csv")

snp1A_random <- read_csv("Genomics_scripts/Data/snp1A_random.csv")
snp2A_random <- read_csv("Genomics_scripts/Data/snp2A_random.csv")
snp5A_random <- read_csv("Genomics_scripts/Data/snp5A_random.csv")

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

###################################################################################
## Generate SNP_random slope tables per site

# snp1A_random
snp1A_random_slope<-data.frame()
counter<-1
for (i in 3:dim(snp1A_random)[2]){
  for(j in 1:12){
    chr<-colnames(snp1A_random)[i]
    popSNP <- snp1A_random %>% filter(Site==j) %>% select(Site,Year,all_of(chr))
    colnames(popSNP)[3]<-"snp_ID"
    
    if(all(is.na(popSNP$snp_ID))==FALSE){
      rSNP <- lm(snp_ID ~ Year, data = popSNP)
      snp1A_random_slope[counter,1]<-unique(popSNP$Site)
      snp1A_random_slope[counter,2]<-chr
      snp1A_random_slope[counter,3]<-abs(rSNP$coefficients[2])
      counter<-counter+1
    } else {
      snp1A_random_slope[counter,1]<-unique(popSNP$Site)
      snp1A_random_slope[counter,2]<-chr
      snp1A_random_slope[counter,3]<-NA
      counter<-counter+1
    }
    
  }
}
colnames(snp1A_random_slope)<-c("Site","snp_ID","Slope")

# snp2A_random
snp2A_random_slope<-data.frame()
counter<-1
for (i in 3:dim(snp2A_random)[2]){
  for(j in 1:12){
    chr<-colnames(snp2A_random)[i]
    popSNP <- snp2A_random %>% filter(Site==j) %>% select(Site,Year,all_of(chr))
    colnames(popSNP)[3]<-"snp_ID"
    
    if(all(is.na(popSNP$snp_ID))==FALSE){
      rSNP <- lm(snp_ID ~ Year, data = popSNP)
      snp2A_random_slope[counter,1]<-unique(popSNP$Site)
      snp2A_random_slope[counter,2]<-chr
      snp2A_random_slope[counter,3]<-abs(rSNP$coefficients[2])
      counter<-counter+1
    } else {
      snp2A_random_slope[counter,1]<-unique(popSNP$Site)
      snp2A_random_slope[counter,2]<-chr
      snp2A_random_slope[counter,3]<-NA
      counter<-counter+1
    }
    
  }
}
colnames(snp2A_random_slope)<-c("Site","snp_ID","Slope")

# snp5A_random
snp5A_random_slope<-data.frame()
counter<-1
for (i in 3:dim(snp5A_random)[2]){
  for(j in 1:12){
    chr<-colnames(snp5A_random)[i]
    popSNP <- snp5A_random %>% filter(Site==j) %>% select(Site,Year,all_of(chr))
    colnames(popSNP)[3]<-"snp_ID"
    
    if(all(is.na(popSNP$snp_ID))==FALSE){
      rSNP <- lm(snp_ID ~ Year, data = popSNP)
      snp5A_random_slope[counter,1]<-unique(popSNP$Site)
      snp5A_random_slope[counter,2]<-chr
      snp5A_random_slope[counter,3]<-abs(rSNP$coefficients[2])
      counter<-counter+1
    } else {
      snp5A_random_slope[counter,1]<-unique(popSNP$Site)
      snp5A_random_slope[counter,2]<-chr
      snp5A_random_slope[counter,3]<-NA
      counter<-counter+1
    }
    
  }
}
colnames(snp5A_random_slope)<-c("Site","snp_ID","Slope")

###################################################################################
#Export

write_csv(snp1A_slope, "Genomics_scripts/Data/snp1A_slope.csv")
write_csv(snp2A_slope, "Genomics_scripts/Data/snp2A_slope.csv")
write_csv(snp5A_slope, "Genomics_scripts/Data/snp5A_slope.csv")

write_csv(snp1A_random_slope, "Genomics_scripts/Data/snp1A_rslope.csv")
write_csv(snp2A_random_slope, "Genomics_scripts/Data/snp2A_rslope.csv")
write_csv(snp5A_random_slope, "Genomics_scripts/Data/snp5A_rslope.csv")






