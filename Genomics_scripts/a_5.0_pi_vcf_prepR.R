##################################################################################
## Prepare Timeseries site/year combo names
## Generate text files per site/year with genomics ID only
## Author Daniel Anstett
## 
## 
## Last Modified May 19, 2021
###################################################################################
#Import libraries
library(tidyverse)
###################################################################################


timeseries_pop_id <- read.csv("Genomics_scripts/Data/timeseries_pop_id.csv", header=F)
space_time <- unique(timeseries_pop_id$V2)
  
  for(i in 1:62){
  timeseries_pop_id_select <- timeseries_pop_id %>% filter(V2==space_time[i]) %>% select(V1)
  name_space_time <- paste("Genomics_scripts/Data/space_time/",space_time[i],".csv",sep="")
  
  write_csv(timeseries_pop_id_select , file = name_space_time)
}

