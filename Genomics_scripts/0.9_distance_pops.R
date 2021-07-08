##################################################################################
## Get the geogrpahic distance between populations
## Author Daniel Anstett 
## 
##
## Last Modified July 8, 2021
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)
library(geodist)



###################################################################################
#Import population data
pop_var_raw <- read_csv("Genomics_scripts/Data/paper_ID_site_select.csv")
deg_distance <- pop_var_raw %>% dplyr::select(Long,Lat)
geodist(deg_distance)
