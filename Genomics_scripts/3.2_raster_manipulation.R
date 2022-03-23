##################################################################################
## Set up Rasters for GF
## Author Daniel Anstett
## 
##
## Last Modified March 23, 2021
#################################################################################### Clear environment
rm(list = ls())

# Get this package retrieving function
## This function will automatically load packages that you already have
## and will install packages you don't yet have then load them
ipak <- function(pkg){
  # Function written by Dr. Evan Fricke
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = T)
  sapply(pkg, require, character.only = T)
}

# Define the packages that the script needs
myPackages <- c("tidyverse","raster")

# Load the packages
ipak(myPackages)

###################################################################################
##Import 1981-2010 raster data for West NA
MAT.clip <- raster("Donor_selection/data/clip/MAT.clip.grd")
MAP.clip <- raster("Donor_selection/data/clip/MAP.clip.grd")
PAS.clip <- raster("Donor_selection/data/clip/PAS.clip.grd")
EXT.clip <- raster("Donor_selection/data/clip/EXT.clip.grd")
CMD.clip <- raster("Donor_selection/data/clip/CMD.clip.grd")

#Seasonal
PPT_sm.clip <- raster("Donor_selection/data/clip/PPT_sm.clip.grd")
PPT_wt.clip <- raster("Donor_selection/data/clip/PPT_wt.clip.grd")
Tave_sm.clip <- raster("Donor_selection/data/clip/Tave_sm.clip.grd")
Tave_wt.clip <- raster("Donor_selection/data/clip/Tave_wt.clip.grd")

#Stack Raster
env_wna <- stack(list(MAT=MAT.clip,MAP=MAP.clip,PAS=PAS.clip,EXT=EXT.clip,CMD=CMD.clip,
                      PPT_sm=PPT_sm.clip,PPT_wt=PPT_wt.clip,Tave_sm=Tave_sm.clip,Tave_wt=Tave_wt.clip))
env_wna <- as.data.frame(env_wna, xy=TRUE)
write.csv(env_wna, "Genomics_scripts/Data/env_wna.csv")




