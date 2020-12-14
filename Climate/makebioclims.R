##################################################################################
## Matthew Bayly
## MAKING BIOCLIM GRIDS
## Convert monthly temp & precip grids to bioclimatic variables used for modeling. 
##   
## Fairly Simple using the BIOVARS function in the 'dimso' package
## July 3, 2014
###################################################################################

# Obtain Climate raster:
# Most user friendly download point:
# Climate data for North America, South America, and Europe:
# choose: ClimateWNA & download as 'ASCII' files (do not download CSV files)
# download: 36 Monthly variables (1961-1990 normal period)
# DO NOT DOWNLOAD 36 BIOCLIMATIC variables!!!!!!!!!!!!!!!
# these are not calculated in the same as the ones used for modeling!
#http://www.ualberta.ca/~ahamann/data/climatewna.html
#* Note: Projection time frame = 1961-1990 normal (just for making pretty maps). 

######################################################################
## INDEX:
# PART 1: LOAD CLIMATE RASTERS, PROJECT, CHECK (QA/QC) & TRIM. 
# PART 2: CONVERT MONTHLY RASTERS TO BIOVARS. 
# PART 3: LOG TRANSFORM RASTERS TO MATCH MODELS. 
# PART 4: SAVE & CHECK OUT. 

## load libraries now!
library(dismo)
library(raster)


path.root = "C:/Users/DW/Desktop/temp.aug.31" 
path.dat = paste(path.root, "/datafiles.run", sep="")
path.obj = paste(path.root, "/objects", sep="")
path.gis = paste(path.root, "/path.gis", sep="")
path.fig = paste(path.root, "/figs", sep="")
path.climatewna = "C:/Users/DW/Desktop/bioclim 30s/climwna19601990" 


######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
# PART 1: LOAD CLIMATE RASTERS, PROJECT, CHECK (QA/QC) & TRIM. 
# Clear any junk out of folder or save downloads to new ClimateWNA subfolder 
# Will load all files based on directory pattern. 
setwd(path.climatewna)
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




######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
# PART 2: CONVERT MONTHLY RASTERS TO BIOVARS.

# biovars works with mothly tmin; tmax & PPT. 
setwd(path.gis)

prc <- allfiles[[1:12]]; names(prc)
tmx <- allfiles[[13:24]]; names(tmx)
tmn <- allfiles[[25:36]]; names(tmn)

b <- biovars(prc, tmn, tmx) # time intensive

#Haiku break:
#This autumn
#As reason for growing old
#A cloud and a bird

b; names(b); plot(b[[1]])
setwd(path.gis)
writeRaster(b, filename="bio.grd", bylayer = TRUE, datatype = 'INT4S', bandorder = 'grd', overwrite = TRUE)

# Adjust variables (temp provided when downloaded as oC*10 to save disc same)

bio1 <- (b[["bio1"]])/10
bio2 <- (b[["bio2"]])/10
bio3 <- (b[["bio3"]]) # don't  adjust bio3 (only a % diff)
bio4 <- (b[["bio4"]])/10
bio5 <- (b[["bio5"]])/10
bio6 <- (b[["bio6"]])/10
bio7 <- (b[["bio7"]])/10
bio8 <- (b[["bio8"]])/10
bio9 <- (b[["bio9"]])/10
bio10 <- (b[["bio10"]])/10
bio11 <- (b[["bio11"]])/10
bio12 <- (b[["bio12"]])# precip variables, no adjustment needed. 
bio13 <- (b[["bio13"]])
bio14 <- (b[["bio14"]])
bio15 <- (b[["bio15"]])
bio16 <- (b[["bio16"]])
bio17 <- (b[["bio17"]])
bio18 <- (b[["bio18"]])
bio19 <- (b[["bio19"]])

# save raw rasters
setwd(path.gis)
for (i in 1:19) {
	bio = get(paste("bio",i, sep=""))
	writeRaster(bio, file=paste("bio", i, ".grd",sep=""), bylayer = TRUE, datatype = 'INT4S', bandorder = 'grd', overwrite = TRUE)
	}
	
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
# PART 3: LOG TRANSFORM RASTERS TO MATCH MODELS. 

# predictor set used for modeling
#bio15, lnbio10, lnbio14, lnbio12, bio11, bio4, lnbio3, bio2
bio3 = log(bio3+0.5)
bio10 = log(bio10+0.5)
bio12 = log(bio12+0.5)
bio14 = log(bio14+0.5)


######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
# PART 4: SAVE & CHECK OUT.

# store predictor set in a sub folder to avoid confusion with raw values

dir.create(paste(path.gis, "/pred.set", sep="")) 
path.pred.set <- paste(path.gis, "/pred.set", sep="")
setwd(path.pred.set)

# save rasters used in modelling
#bio15, lnbio10, lnbio14, lnbio12, bio11, bio4, lnbio3, bio2
writeRaster(bio2, file="bio2.grd", bylayer = TRUE, datatype = 'INT4S', bandorder = 'grd', overwrite = TRUE)
writeRaster(bio3, file="bio3.grd", bylayer = TRUE, datatype = 'INT4S', bandorder = 'grd', overwrite = TRUE)
writeRaster(bio4, file="bio4.grd", bylayer = TRUE, datatype = 'INT4S', bandorder = 'grd', overwrite = TRUE)
writeRaster(bio10, file="bio10.grd", bylayer = TRUE, datatype = 'INT4S', bandorder = 'grd', overwrite = TRUE)
writeRaster(bio11, file="bio11.grd", bylayer = TRUE, datatype = 'INT4S', bandorder = 'grd', overwrite = TRUE)
writeRaster(bio12, file="bio12.grd", bylayer = TRUE, datatype = 'INT4S', bandorder = 'grd', overwrite = TRUE)
writeRaster(bio14, file="bio14.grd", bylayer = TRUE, datatype = 'INT4S', bandorder = 'grd', overwrite = TRUE)
writeRaster(bio15, file="bio15.grd", bylayer = TRUE, datatype = 'INT4S', bandorder = 'grd', overwrite = TRUE)


