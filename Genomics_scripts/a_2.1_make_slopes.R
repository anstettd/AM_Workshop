##################################################################################
## Get slopes (strength of selection) for all BF
## 
## Author Daniel Anstett
## 
##
## Last Modified September 26, 2022
###################################################################################
#Import libraries
library(tidyverse)

###################################################################################
#Functions

slope_env <- function(df){
  freq_MAT_slope<-data.frame()
  counter<-1
  for (i in 3:dim(df)[2]){
    print(i-2)
    for(j in 1:12){
      chr<-colnames(df)[i]
      popSNP <- df %>% filter(Site==j) %>% select(Site,Year,all_of(chr))
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
  return(freq_MAT_slope)
}

###################################################################################
#Data Import
#freq_MAT_1 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_small.csv")

#Runs much too slowly for 1.8 million SNPS (variables)
slope_MAT_1 <- slope_env(freq_MAT_1)






