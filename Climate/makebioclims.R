##################################################################################
## Turn Climate WNA into bioclim
## Daniel Anstett
## December 15th
## Fairly Simple using the BIOVARS function in the 'dimso' package
###################################################################################

#Load libraries
library(dismo)
library(raster)


# PART 1: LOAD CLIMATE RASTERS, PROJECT, CHECK (QA/QC) & TRIM. 
# Will load all files based on directory pattern. 
setwd()


allfiles = list.files(); allfiles # will just take the '.asc' files
file <- list.files(pattern = '.asc'); file; length(file) # should be 36.
allfiles = stack(file) # Importing raster, wait ~ 2 - 4 mins.
allfiles; projection(allfiles); dim(allfiles)
#plot(allfiles[[2]]) #OPTIONAL ~ visual check, plotting "PPT02"

# Import layer used to crop/trim rasters (reduce file size). 
# Cropping layer will be 90m DEM cut
# Stored on dropbox > pseudoabsence methods > megadem.cut.tif
dem <- raster('C:/Users/DW/Desktop/bioclim 30s/hydrosheds/megadem.full.tif')

projection(dem); #plot(dem); extent(dem)

# Set trimming extent (North to Portland; South to Tijuana & Coast - Yuma
e <- extent(-125, -114, 32.5, 47) # extent object
dem.crop <- crop(dem, e); # ~ 2min delay; plot(dem.crop)

# Project from lat/long -> projection(allfiles)
projection(dem); projection(allfiles)


e2 <- projectExtent(dem.crop, projection(allfiles))

samp <- allfiles[[2]] # what we want! just choose one at random as sample
dem.proj <- projectRaster(dem.crop, samp, method="bilinear") # dem file will be slightly damaged (do not extract values)
#plot(samp); plot(dem.proj, add=T) # Success!

# Southern Baja to Lower mainland with longitude of entire range 
xmin= -125; xmax = -114; ymin=27;  ymax=49.5
prj.wgs = "+proj=longlat +ellps=WGS84"
prj.aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0"
prj.lcc = "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
e = extent(xmin, xmax, ymin, ymax)
library(rgeos)
bbox = as(e, "SpatialPolygons")
proj4string(bbox) = CRS(prj.wgs)
bbox.lcc = spTransform(bbox, CRS=CRS(prj.lcc))
#Clip function to trim probability rasters by bbox.lcc (extent obj)
clip<-function(raster,shape) {
  a1_crop<-crop(raster,shape)
  step1<-rasterize(shape,a1_crop)
  a1_crop*step1}

# TRIM down files & crop. #untrimmed prior to bioclim
allfiles <- clip(allfiles, bbox.lcc)



