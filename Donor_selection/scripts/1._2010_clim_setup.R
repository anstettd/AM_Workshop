##############################################################################
### SCRIPT PURPOSE: Transform climate rasters to bioclim 
# Modified from Angert et al. 2018, American Naturalist
# Author: Amy Angert & Daniel Anstett
# last update:  June 23 2021

## OVERALL WORKFLOW:
# !!! For raster data, this assumes you have downloaded tiff rasters for 36 annual variables 
#from climateNA (https://adaptwest.databasin.org/pages/adaptwest-climatena)
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
#Functions

# A function to crop climate rasters by a bounding box
clip <- function(raster, shape) { 
  a1_crop <- crop(raster, shape)
  step1 <- rasterize(shape, a1_crop)
  return(a1_crop*step1)
}

##############################################################################
### MANIPULATE RASTER FILES TO DECREASE SIZE TO RELEVANT AREA

## Read in TIF files all at once
# !!! Set working directory to folder with downloaded ASCII files
setwd("~/Dropbox/AM_Workshop/Climate/Normal_1981_2010_bioclim/Selected")
allfiles.list <- list.files(pattern = '.tif') # list of '.tif' files
allfiles <- stack(allfiles.list) # import set of rasters
plot(allfiles[[2]]) # optional visual check, plotting "plot second raster"

#CRS manipulation
crs(allfiles[[1]]) # get CRS. In laea projection, Lambert azimuthal equal-area projection
prj.laea <- "+proj=laea +lon_0=-100 +lat_0=45 +type=crs" ## set laea CRS
proj4string(allfiles) <- CRS(prj.laea) #define projection of rasters
projection(allfiles) #optional check correct projection applied

# The rasters are for all of North America and too large for storage in this repo
# Trim them to study area and save only trimmed files in the repo
# Define extent as 1 degree beyond lat-long extent of points
ext <- extent(min(clim$Longitude)-1, max(clim$Longitude)+1, min(clim$Latitude)-1, max(clim$Latitude)+1)
bbox = as(ext, "SpatialPolygons") #convert coordinates to a bounding box
prj.wgs = "+proj=longlat + type=crs" #define unprojected coordinate system
proj4string(bbox) <- CRS(prj.wgs) #set projection
bbox.lcc = spTransform(bbox, CRS=CRS(prj.laea)) #re-project to match rasters

#mask(i.e. clip) files by bounding box
allfiles.clip <- clip(allfiles, bbox.lcc) #apply the function to the raster stack
names(allfiles.clip) <- names(allfiles) #replace layer names
plot(allfiles.clip[[2]]) #plot clipped layer (PPT02 as an example)
plot(allfiles[[2]]) #compare to unclipped layer

##############################################################################
# Save individual rasters
# !!! Be sure working directory is set back to project folder
setwd("~/Dropbox/AM_Workshop/AM_Workshop")

writeRaster(allfiles.clip[[1]],"Donor_selection/data/1981_2010/CMD.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)
writeRaster(allfiles.clip[[2]],"Donor_selection/data/1981_2010/EXT.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)
writeRaster(allfiles.clip[[3]],"Donor_selection/data/1981_2010/MAP.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)
writeRaster(allfiles.clip[[4]],"Donor_selection/data/1981_2010/MAT.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)
writeRaster(allfiles.clip[[5]],"Donor_selection/data/1981_2010/PAS.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)

writeRaster(allfiles.clip[[6]],"Donor_selection/data/1981_2010/PPT_sm.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)
writeRaster(allfiles.clip[[7]],"Donor_selection/data/1981_2010/PPT_wt.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)
writeRaster(allfiles.clip[[8]],"Donor_selection/data/1981_2010/Tave_sm.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)
writeRaster(allfiles.clip[[9]],"Donor_selection/data/1981_2010/Tave_wt.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)
##############################################################################














#Messy work

## get CRS
##test_crs<-raster("~/Dropbox/AM_Workshop/Climate/Normal_1981_2010_bioclim/Selected/Normal_1981_2010_CMD.tif")
##crs(test_crs) # In laea projection, Lambert azimuthal equal-area projection


#re-project rasters into WGS 1984 (EPSG 4326)
##EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
##cmd_raster <- projectRaster(allfiles[[1]], crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)


#what the crs() shows new
##CRS arguments:
#  +proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs 

#what the crs shows old
#CRS arguments:
#  +proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=WGS84 +units=m
#+no_defs

#Inutputed crs
#prj.lcc <- "+proj=lcc +lon_0=-95 +lat_1=49 +lat_2=77 +type=crs" 


