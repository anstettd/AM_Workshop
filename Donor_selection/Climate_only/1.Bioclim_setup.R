##############################################################################
### SCRIPT PURPOSE: Transform climate rasters to bioclim 
# Modified from Angert et al. 2018, American Naturalist
# Author: Amy Angert & Daniel Anstett
# last update:  5 Jan 2021

## OVERALL WORKFLOW:

# !!! For raster data, this assumes you have downloaded ascii rasters for 48 monthly variables from climateNA (https://adaptwest.databasin.org/pages/adaptwest-climatena)
# Then, read in climateNA output file or rasters with monthly temperature and precipitation records for 1961-1990
# Calculate bioclimatic variables from monthly temperature and precipitation

##############################################################################


############################################################################## 
### LOAD LIBRARIES AND PREPARE INPUTS

## CLEAR WORKSPACE
rm(list = ls(all.names = TRUE))

## LIBRARIES
library(tidyverse) # for data manipulation
library(dismo) # for biovars function
library(raster) # for raster grids
library(rgdal) # for transforming projections

## INPUTS
## (just points files for now; rasters loaded below)
clim <- read_csv("SDM/data_files/points_Normal_1961_1990MSY.csv")

# Ensure no missing values
complete <- clim[complete.cases(clim), ] # should be same dim as clim if there are no missing values

##############################################################################

##############################################################################
### MAKE GRIDS FOR MAPPING BIOCLIM VARIABLES AND MODEL PROJECTIONS

## Read in raw ASCII files for monthly variables
# !!! Use drop-down menu to set working directory to folder with downloaded ASCII files
#setwd("~/Dropbox/AM_Workshop/Climate/NA_6190")
allfiles.list <- list.files(pattern = '.asc') # list of '.asc' files
allfiles <- stack(allfiles.list) # import set of rasters
#plot(allfiles[[2]]) # optional visual check, plotting "PPT02"

# These are in Lambert Conformal Conic projection
prj.lcc <- "+proj=lcc +lon_0=-95 +lat_1=49 +lat_2=77 +type=crs"  #specify projection ('type=crs' is key to using old proj4 syntax)
proj4string(allfiles) <- CRS(prj.lcc) #define projection of rasters
#projection(allfiles) #optional check correct projection applied

# These are for all of North America and too large for storage in this repo
# Trim them to study area and save only trimmed files in the repo
# Define extent as 1 degree beyond lat-long extent of points
ext <- extent(min(clim$Longitude)-1, max(clim$Longitude)+1, min(clim$Latitude)-1, max(clim$Latitude)+1)
bbox = as(ext, "SpatialPolygons") #convert coordinates to a bounding box
prj.wgs = "+proj=longlat + type=crs" #define unprojected coordinate system
proj4string(bbox) <- CRS(prj.wgs) #set projection
bbox.lcc = spTransform(bbox, CRS=CRS(prj.lcc)) #re-project to match rasters

# Crop climate rasters by bounding box
clip <- function(raster, shape) {  #function to crop by bounding box
  a1_crop <- crop(raster, shape)
  step1 <- rasterize(shape, a1_crop)
  a1_crop*step1
}
allfiles.clip <- clip(allfiles, bbox.lcc) #apply the function to the raster stack
names(allfiles.clip) <- names(allfiles) #replace layer names
plot(allfiles.clip[[2]]) #plot clipped layer (PPT02 as an example)
plot(allfiles[[2]]) #compare to unclipped layer

## Calculate bioclim rasters
prec.rast <- allfiles.clip[[1:12]]; names(prec.rast)
tmax.rast <- allfiles.clip[[25:36]]; names(tmax.rast)
tmin.rast <- allfiles.clip[[37:48]]; names(tmin.rast)
b <- biovars(prec.rast, tmin.rast, tmax.rast)
plot(b[[1]])

# Pull out individual variables and transform as needed
bio1 <- (b[["bio1"]])
bio2 <- (b[["bio2"]])
bio3 <- (b[["bio3"]])
bio4 <- (b[["bio4"]])
bio5 <- (b[["bio5"]])
bio6 <- (b[["bio6"]])
bio7 <- (b[["bio7"]])
bio8 <- (b[["bio8"]])
bio9 <- (b[["bio9"]])
bio10 <- (b[["bio10"]])
bio11 <- (b[["bio11"]])
bio12 <- log((b[["bio12"]])+0.5)
bio13 <- log((b[["bio13"]])+0.5)
bio14 <- log((b[["bio14"]])+0.5)
bio15 <- log((b[["bio15"]])+0.5)
bio16 <- log((b[["bio16"]])+0.5)
bio17 <- log((b[["bio17"]])+0.5)
bio18 <- log((b[["bio18"]])+0.5)
bio19 <- log((b[["bio19"]])+0.5)


##############################################################################
# Save individual rasters
# !!! Be sure working directory is set back to project folder
#setwd("~/Dropbox/AM_Workshop/AM_Workshop")

for (i in 1:19) {
  bio = get(paste("bio",i, sep=""))
  writeRaster(bio, file=paste("Donor_selection/data/bio", i, ".grd",sep=""), 
              bylayer = TRUE, datatype = 'INT4S', bandorder = 'grd', 
              overwrite = TRUE)
}
##############################################################################







