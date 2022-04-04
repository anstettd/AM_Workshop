##################################################################################
## Gradient forest for full SNP dataset
## Author Daniel Anstett
## 
## Modified from Keller & Fitzpatric 2015
## Last Modified Dec 13, 2021
###################################################################################

# Clear environment
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
myPackages <- c("randomForest", "extendedForest", "gradientForest","tidyverse","raster")

# Load the packages
ipak(myPackages)

# Install packages
#install.packages("extendedForest", repos="http://R-Forge.R-project.org")
#install.packages("gradientForest", repos="http://R-Forge.R-project.org")

library(randomForest)
library(gradientForest)
library(tidyverse)
library(dismo) # for biovars function
library(raster) # for raster grids
library(rgdal) # for transforming projections
library(sf)

###################################################################################

##Import Data
snp_clim_bf10 <- read_csv("Genomics_scripts/Data/snp_clim_bay.csv") #pop data
#env_wna <- read_csv("Genomics_scripts/Data/env_wna.csv") #grided WNA climate data

#snp_clim_bf10NA <- snp_clim_bf10 %>%
#  select_if(~ !any(is.na(.)))
####snp_clim_ful <- read_csv("Genomics_scripts/Data/snp_clim_full.csv") # full data


## Generate specific dataframes for GF model
env_site <- snp_clim_bf10 %>% dplyr::select(MAT,MAP,PAS,EXT,CMD,Tave_wt,Tave_sm,PPT_wt,PPT_sm)

#Computationally intensive, so just load resulting filtered file
#test_snp <- snp_clim_bf10[,15:114] %>%
#  select_if(~ !any(is.na(.)))
#write_csv(test_snp,"Genomics_scripts/Data/test_snp.csv")

#Import filtered test snp datatset
test_snp <- read.csv("Genomics_scripts/Data/test_snp.csv")

env_site<-as.data.frame(env_site)
test_snp<-as.data.frame(test_snp)

df_in_1<-data.frame(env_site, test_snp)
pred<-colnames(env_site)
resp<-colnames(test_snp)

#resp<-as.factor(resp)
#resp<-droplevels(resp)


###################################################################################
#Make a climatic 
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
#env_wna <- as.data.frame(env_wna, xy=TRUE)

#colnames(env_wna)[1]<-"X"
#colnames(env_wna)[2]<-"Y"

# The rasters are for all of North America and too large for storage in this repo
# Trim them to study area and save only trimmed files in the repo
# Define extent as 1 degree beyond lat-long extent of points

clim <- read_csv("SDM/data_files/points_Normal_1961_1990MSY.csv")
prj.laea <- "+proj=laea +lon_0=-100 +lat_0=45 +type=crs" ## set laea CRS



ext <- extent(min(clim$Longitude)-1, max(clim$Longitude)+1, min(clim$Latitude)-1, max(clim$Latitude)+1)
bbox = as(ext, "SpatialPolygons") #convert coordinates to a bounding box
prj.wgs = "+proj=longlat + type=crs" #define unprojected coordinate system
proj4string(bbox) <- CRS(prj.wgs) #set projection
bbox.lcc = spTransform(bbox, CRS=CRS(prj.laea)) #re-project to match rasters

#Now we're going to use the bbox info as the y argument to extract info from the raster


env_wna_df <- extract(env_wna,)

extract(env_wna, bbox.lcc,cellnumbers=T)

###################################################################################


gfRef <- gradientForest(cbind(envGF, SNPs_ref), predictor.vars=colnames(envGF),
                        response.vars=colnames(SNPs_ref), ntree=500, 
                        maxLevel=maxLevel, trace=T, corr.threshold=0.50)


# Gradient Forest Model
gf <- gradientForest(df_in_1,
                      predictor.vars = pred, response.vars = resp,
                      ntree = 500, transform = NULL, compact = T,
                      nbin = 201 , corr.threshold = 0.5)


gf <- gradientForest(cbind(env_site, test_snp),predictor.vars = colnames(env_site),
                     response.vars = colnames(test_snp),ntree = 500,
                     maxLevel=maxLevel, trace=T, corr.threshold=0.50)                    

#Importance Plot
plot(gf)


##################################################################################################################

# Mapping spatial genetic variation --------------------------------------------
###### functions to support mapping #####
# builds RGB raster from transformed environment
# snpPreds = dataframe of transformed variables from gf or gdm model
# rast = a raster mask to which RGB values are to be mapped
# cellNums = cell IDs to which RGB values should be assigned
pcaToRaster <- function(snpPreds, rast, mapCells){
  require(raster)
  
  pca <- prcomp(snpPreds, center=TRUE, scale.=FALSE)
  
  ##assigns to colors, edit as needed to maximize color contrast, etc.
  a1 <- pca$x[,1]; a2 <- pca$x[,2]; a3 <- pca$x[,3]
  r <- a1+a2; g <- -a2; b <- a3+a2-a1
  
  ##scales colors
  scalR <- (r-min(r))/(max(r)-min(r))*255
  scalG <- (g-min(g))/(max(g)-min(g))*255
  scalB <- (b-min(b))/(max(b)-min(b))*255
  
  ##assigns color to raster
  rast1 <- rast2 <- rast3 <- rast
  rast1[mapCells] <- scalR
  rast2[mapCells] <- scalG
  rast3[mapCells] <- scalB
  ##stacks color rasters
  outRast <- stack(rast1, rast2, rast3)
  return(outRast)
}

# OK, on to mapping. Script assumes:
# (1) a dataframe named env_trns containing extracted raster data (w/ cell IDs)
# and env. variables used in the models & with columns as follows: cell, bio1, bio2, etc.
#
# (2) a raster mask of the study region to which the RGB data will be written

# transform env using gf models, see ?predict.gradientForest
predRef <- predict(gf, env_wna[,-1]) # remove cell column before transforming

# map continuous variation 
refRGBmap <- pcaToRaster(predRef, mask, env_wna$cell)
plotRGB(refRGBmap)
writeRaster(refRGBmap, "/.../refSNPs_map.tif", format="GTiff", overwrite=TRUE)




################################################################################


# Calculate and map "genetic offset" under climate change ----------------------
# Script assumes:
# (1) a dataframe of transformed env. variables for CURRENT climate 
# (e.g., predGI5 from above).
#
# (2) a dataframe named env_trns_future containing extracted raster data of 
# env. variables for FUTURE a climate scenario, same structure as env_trns

# first transform FUTURE env. variables
projGI5 <- predict(gfGI5, env_trns_future[,-1])

# calculate euclidean distance between current and future genetic spaces  
genOffsetGI5 <- sqrt((projGI5[,1]-predGI5[,1])^2+(projGI5[,2]-predGI5[,2])^2
                     +(projGI5[,3]-predGI5[,3])^2+(projGI5[,4]-predGI5[,4])^2
                     +(projGI5[,5]-predGI5[,5])^2+(projGI5[,6]-predGI5[,6])^2
                     +(projGI5[,7]-predGI5[,7])^2)

# assign values to raster - can be tricky if current/future climate
# rasters are not identical in terms of # cells, extent, etc.
mask[env_trns_future$cell] <- genOffsetGI5
plot(mask)





















##################################################################################################################
#Basic Plots

most_important <- names(importance(gf))[1:9]

##Split Density Plot
#The second plot is the splits density plot (plot.type="S"), which shows binned 
#split importance and location on each gradient (spikes), kernel density of splits 
#(black lines), of observations(red lines) and of splits standardised by 
#observations density (blue lines). Each distribution integrates to predictor 
#importance. These show where important changes in the abundance of
#multiple species are occurring along the gradient; they indicate a composition change rate
plot(gf, plot.type = "S", imp.vars = most_important,leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
     cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))

##Cumulaive Plot
#for each species shows cumulative importance distributions of splits improvement scaled by
#R2 weighted importance, and standardised by density of observations. These show cumulative
#change in abundance of individual species, where changes occur on the gradient, and the species
#changing most on each gradient.
plot(gf, plot.type = "C", imp.vars = most_important,show.overall = F, legend = T, leg.posn = "topleft",
       leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4,cex.axis = 0.6, line.ylab = 0.9, 
       par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0,0.3, 0, 0)))


plot(gf, plot.type = "P", show.names = F, horizontal = F, cex.axis = 1, cex.labels = 0.7, line = 2.5)
















