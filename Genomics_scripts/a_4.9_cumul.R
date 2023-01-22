##################################################################################
## Get cumulative slope measure
## 
## Author Daniel Anstett
## 
## 
## Last Modified May 17, 2022
###################################################################################
#Import libraries
library(tidyverse)

###################################################################################
#Import files
snp_list_unique <- read.csv("Genomics_scripts/Data/snp_list_unique.csv")
snp_list <- read.csv("Genomics_scripts/Data/snp_list.csv")
snp_list_mat <- snp_list  %>% filter(env=="MAT")
snp_list_map <- snp_list  %>% filter(env=="MAP")
snp_list_cmd <- snp_list  %>% filter(env=="CMD")

offset_pop <- read_csv("Genomics_scripts/Data/offset_pop_lambda.csv")

###################################################################################
#Produce integrative slope

slope.cumul <- data.frame()
cumul_mat <- snp_list_mat %>% group_by(Site) %>% summarise(cumul_slope = sum(Slope))
cumul_map <- snp_list_map %>% group_by(Site) %>% summarise(cumul_slope = sum(Slope))
cumul_cmd <- snp_list_cmd %>% group_by(Site) %>% summarise(cumul_slope = sum(Slope))
cumul_unique <- snp_list_unique %>% group_by(Site) %>% summarise(cumul_slope = sum(Slope))

#Populate blank dataframe 
slope.cumul[1,1] <- cumul_mat$cumul_slope[1]
slope.cumul[2,1] <- cumul_mat$cumul_slope[2]
slope.cumul[3,1] <- cumul_mat$cumul_slope[3]
slope.cumul[4,1] <- cumul_mat$cumul_slope[4]
slope.cumul[5,1] <- cumul_mat$cumul_slope[5]
slope.cumul[6,1] <- cumul_mat$cumul_slope[6]
slope.cumul[7,1] <- 0
slope.cumul[8,1] <- 0
slope.cumul[9,1] <- 0
slope.cumul[10,1] <- 0
slope.cumul[11,1] <- cumul_mat$cumul_slope[7]
slope.cumul[12,1] <- 0

slope.cumul[,2] <- cumul_map$cumul_slope

slope.cumul[1,3] <- cumul_cmd$cumul_slope[1]
slope.cumul[2,3] <- cumul_cmd$cumul_slope[2]
slope.cumul[3,3] <- cumul_cmd$cumul_slope[3]
slope.cumul[4,3] <- cumul_cmd$cumul_slope[4]
slope.cumul[5,3] <- cumul_cmd$cumul_slope[5]
slope.cumul[6,3] <- cumul_cmd$cumul_slope[6]
slope.cumul[7,3] <- 0
slope.cumul[8,3] <- cumul_cmd$cumul_slope[7]
slope.cumul[9,3] <- cumul_cmd$cumul_slope[8]
slope.cumul[10,3] <- cumul_cmd$cumul_slope[9]
slope.cumul[11,3] <- cumul_cmd$cumul_slope[10]
slope.cumul[12,3] <- cumul_cmd$cumul_slope[11]

colnames(slope.cumul) <- c("cumul_MAT","cumul_MAT","cumul_CMD")

write_csv(slope.cumul,"Genomics_scripts/Data/slope.cumul.csv")
write_csv(cumul_unique,"Genomics_scripts/Data/slope.cumul.unique.csv")











