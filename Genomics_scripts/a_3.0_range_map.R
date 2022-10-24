#############################################################################################################
## Calculate range map
## Author Daniel Anstett
## 
## Modified from Qin Li
## Last Modified September 13, 2022
#############################################################################################################
#Import libraries
library(sp)
library(maptools)
library(rgdal)
library(rgeos)
library(raster)
library(tidyverse)
library(spatstat)

### (1) Buffered Polygon 2 (recommended)

# lat/lon
proj_xy <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

# AEA projection parameters
# parameter for North America/California
# http://spatialreference.org/ref/esri/north-america-albers-equal-area-conic/
aea_NA <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# load shapefile for North America/US boundary (or from other sources)
data(wrld_simpl, package = "maptools")
us_shp <- wrld_simpl[grepl(c("United States"), wrld_simpl$NAME), ]
us_shp_aea <- spTransform(us_shp, CRS(aea_NA))

# load species locality and make projection
#spp_locality = readRDS(file.choose()) # locality data
spp_locality <- read.csv("SDM/data_files/presences.csv")
spp_locality <- spp_locality %>% mutate(sp_name="M_cardinalis")
spp_locality_select <- spp_locality %>% dplyr::select(sp_name,Latitude,Longitude)
spp_shp <- SpatialPoints(spp_locality[,c("Longitude","Latitude")])
#spp_shp <- SpatialPoints(spp_locality[,c("Species","Longitude","Latitude")])
spp_shp <- SpatialPointsDataFrame(spp_shp, 
                    data.frame(ID=1:nrow(spp_locality),Species=spp_locality$sp_name))
proj4string(spp_shp) <- CRS(proj_xy)
spp_shp_aea <- spTransform(spp_shp, CRS(aea_NA))

# MCP
#temp <- mcp(spp_shp_aea, percent=100)
#temp10k <- gBuffer(temp, width = 10000)
#spp_mcp_aea10k <- gIntersection(temp10k,us_shp_aea)

# Buffered Polygon 1: points with buffer 10km
temp10k1 <- gBuffer(spp_shp_aea, width = 10000)
spp_pb_aea10k1 <- gIntersection(temp10k1,us_shp_aea)

# Buffered Polygon 2: points with buffer 50km first and reduce by 40km (to avoid gaps between points)
temp50k <- gBuffer(spp_shp_aea, width = 50000)
spp_pb_aea50k <- gIntersection(temp50k,us_shp_aea)
temp10k2 <- gBuffer(temp50k, width = -20000) # reduce the buffer (-40km)
spp_pb_aea10k2 <- gIntersection(temp10k2,us_shp_aea)
plot(spp_pb_aea10k2)



#Export 50 km buffer
shapefile(x = spp_pb_aea10k2, file = "Shape/c_range50.shp",overwrite=TRUE)
#writeOGR(spp_pb_aea10k2, dsn = '.', layer = 'poly', driver = "ESRI Shapefile")




### (2) Range Polygon 2 by locality density
# ref: 'rangemap' package
# https://cran.r-project.org/web/packages/rangemap/vignettes/rangemap_short_tutorial_I.html
options(rgl.useNULL = TRUE)
library(rangemap)
rgeos::set_RGEOS_CheckValidity(2L)
pts_temp = spp_locality[,c("sp_name","Longitude","Latitude")] # need to be this order

# concave hull
# not quite sure about the parameter setting:
# split_distance = 1500000? buffer_distance = 500000?
sp_range = rangemap_hull(occurrences = pts_temp, hull_type = "concave",
                         buffer_distance = 500000, split = TRUE,
                         cluster_method = "hierarchical",
                         split_distance = 1500000,
                         extent_of_occurrence=F, area_of_occupancy=F)

#rangemap_plot(sp_range) #plot

range_shp = sp_range@species_range
range_shp = aggregate(range_shp, dissolve = TRUE) # perhaps multiple polygons, so merge them
range_shp = spTransform(range_shp, CRS(aea_NA))

#plot(range_shp)
#plot(SpatialPoints(spp_locality[,c("Species","Longitude","Latitude")]), add=T)

er_owin = as.owin(range_shp)

# create a ‘ppp’ (point pattern) object from xy
p_sp = ppp(spp_locality$Longitude, spp_locality$Latitude, window=er_owin)
plot(p_sp)

# compute Kernel Density to show locality density
ds_node_plot = density(p_sp)

library(wesanderson)
# Gradient color
pal <- wes_palette("Zissou1", 100, type = "continuous")

pdf(file = "range_map_with_locality_density.pdf", width=5, height=5)
plot(ds_node_plot, col=pal, add=T)
plot(SpatialPoints(spp_locality[,c("Species","Longitude","Latitude")]), add=T)
dev.off()
