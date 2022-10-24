##################################################################################
## Get positive slope info per population
## Author Daniel Anstett
## 
## 
## Last Modified September 19, 2022
###################################################################################

#Library install and import
library(tidyverse) 

#Bring in dataframe with offset information
offset_pop <- read_csv("Genomics_scripts/Data/offset_pop.csv")

#Import lambda from demography data for all sites
demo_all <- read_csv("Genomics_scripts/Data/lambda_2010_2016.csv")
demo <- demo_all %>% filter(Site =="SweetwaterCreek" | Site == "WestForkMojaveRiver" | 
                              Site =="NorthForkMiddleForkTule" | Site == "RedwoodCreek" |
                              Site == "Wawona" | Site == "OregonCreek" |
                              Site == "LittleJamesonCreek" | Site == "DeepCreek" |
                              Site == "O'NeilCreek" | Site == "DeerCreek" |
                              Site == "RockCreek" | Site == "MillCreek")
demo_1215 <- demo %>% filter(Year == 2012 | Year == 2013 | Year == 2014 | Year == 2015)
site_popID <- read_csv("Genomics_scripts/Data/site_popID.csv")
demo_1215_ID <- left_join(demo_1215,site_popID,by="Site")



