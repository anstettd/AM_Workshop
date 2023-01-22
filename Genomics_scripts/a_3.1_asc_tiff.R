#############################################################################################################
## Convert asc to TiFF files 
## Author Daniel Anstett
## 
## 
## Last Modified September 16, 2022
#############################################################################################################
#Import libraries
library(raster)
library(tidyverse)
library(sf)
library(rgdal)

#############################################################################################################
#Import asc file
mat_2011 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2011/mat.asc")
map_2011 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2011/map.asc")
cmd_2011 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2011/cmd.asc")

mat_2012 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2012/mat.asc")
map_2012 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2012/map.asc")
cmd_2012 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2012/cmd.asc")

mat_2013 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2013/mat.asc")
map_2013 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2013/map.asc")
cmd_2013 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2013/cmd.asc")

mat_2014 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2014/mat.asc")
map_2014 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2014/map.asc")
cmd_2014 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2014/cmd.asc")

mat_2015 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2015/mat.asc")
map_2015 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2015/map.asc")
cmd_2015 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2015/cmd.asc")

mat_2016 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2016/mat.asc")
map_2016 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2016/map.asc")
cmd_2016 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Year_2016/cmd.asc")

mat_8110 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Normal_1981_2010_800/mat.asc")
map_8110 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Normal_1981_2010_800/map.asc")
cmd_8110 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/Normal_1981_2010_800/cmd.asc")

mat_4170_45 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/13GCMs_ensemble_ssp245_2041-2070Y/mat.asc")
map_4170_45 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/13GCMs_ensemble_ssp245_2041-2070Y/map.asc")
cmd_4170_45 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/13GCMs_ensemble_ssp245_2041-2070Y/cmd.asc")

mat_4170_85 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/13GCMs_ensemble_ssp585_2041-2070Y/mat.asc")
map_4170_85 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/13GCMs_ensemble_ssp585_2041-2070Y/map.asc")
cmd_4170_85 <- raster("C:/Users/anstett3/Documents/Genomics/Large_files/asc/13GCMs_ensemble_ssp585_2041-2070Y/cmd.asc")



#Reproject to WGS 1984 (EPSG4326)
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
crs(mat_2011) <- EPSG4326
crs(map_2011) <- EPSG4326
crs(cmd_2011) <- EPSG4326

crs(mat_2012) <- EPSG4326
crs(map_2012) <- EPSG4326
crs(cmd_2012) <- EPSG4326

crs(mat_2013) <- EPSG4326
crs(map_2013) <- EPSG4326
crs(cmd_2013) <- EPSG4326

crs(mat_2014) <- EPSG4326
crs(map_2014) <- EPSG4326
crs(cmd_2014) <- EPSG4326

crs(mat_2015) <- EPSG4326
crs(map_2015) <- EPSG4326
crs(cmd_2015) <- EPSG4326

crs(mat_2016) <- EPSG4326
crs(map_2016) <- EPSG4326
crs(cmd_2016) <- EPSG4326

crs(mat_8110) <- EPSG4326
crs(map_8110) <- EPSG4326
crs(cmd_8110) <- EPSG4326

crs(mat_4170_45) <- EPSG4326
crs(map_4170_45) <- EPSG4326
crs(cmd_4170_45) <- EPSG4326

crs(mat_4170_85) <- EPSG4326
crs(map_4170_85) <- EPSG4326
crs(cmd_4170_85) <- EPSG4326



# Average 2012 to 2015 rasters
mat_1215 <- (mat_2012 + mat_2013 + mat_2014 + mat_2015)/4 
map_1215 <- (map_2012 + map_2013 + map_2014 + map_2015)/4  
cmd_1215 <- (cmd_2012 + cmd_2013 + cmd_2014 + cmd_2015)/4  




#Export tif file
writeRaster(mat_1215, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_1215/MAT.tif",format="GTiff")
writeRaster(map_1215, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_1215/MAP.tif",format="GTiff")
writeRaster(cmd_1215, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_1215/CMD.tif",format="GTiff")

writeRaster(mat_2011, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2011/MAT.tif",format="GTiff")
writeRaster(map_2011, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2011/MAP.tif",format="GTiff")
writeRaster(cmd_2011, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2011/CMD.tif",format="GTiff")

writeRaster(mat_2012, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2012/MAT.tif",format="GTiff")
writeRaster(map_2012, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2012/MAP.tif",format="GTiff")
writeRaster(cmd_2012, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2012/CMD.tif",format="GTiff")

writeRaster(mat_2013, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2013/MAT.tif",format="GTiff")
writeRaster(map_2013, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2013/MAP.tif",format="GTiff")
writeRaster(cmd_2013, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2013/CMD.tif",format="GTiff")

writeRaster(mat_2014, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2014/MAT.tif",format="GTiff")
writeRaster(map_2014, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2014/MAP.tif",format="GTiff")
writeRaster(cmd_2014, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2014/CMD.tif",format="GTiff")

writeRaster(mat_2015, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2015/MAT.tif",format="GTiff")
writeRaster(map_2015, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2015/MAP.tif",format="GTiff")
writeRaster(cmd_2015, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2015/CMD.tif",format="GTiff")

writeRaster(mat_2016, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2016/MAT.tif",format="GTiff")
writeRaster(map_2016, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2016/MAP.tif",format="GTiff")
writeRaster(cmd_2016, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_2016/CMD.tif",format="GTiff")

writeRaster(mat_4170_45, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_4170_45/MAT.tif",format="GTiff")
writeRaster(map_4170_45, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_4170_45/MAP.tif",format="GTiff")
writeRaster(cmd_4170_45, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_4170_45/CMD.tif",format="GTiff")

writeRaster(mat_4170_85, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_4170_85/MAT.tif",format="GTiff")
writeRaster(map_4170_85, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_4170_85/MAP.tif",format="GTiff")
writeRaster(cmd_4170_85, "C:/Users/anstett3/Documents/Genomics/Large_files/Year_4170_85/CMD.tif",format="GTiff")





#Export tif file
writeRaster(mat_2016, "C:/Users/anstett3/Documents/Genomics/Large_files/Year/mat.grd ",
              format="raster", overwrite = T)

