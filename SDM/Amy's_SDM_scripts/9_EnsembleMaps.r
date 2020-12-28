#################################################################################
### SCRIPT PURPOSE: calculate, plot, and save ensemble average of SDMs 
# Modified from Angert et al. 2018, American Naturalist
# Author: Amy Angert
# last update:  28 Dec 2020

##################################################################################### LOAD LIBRARIES AND PREPARE INPUTS

## Clear workspace
rm(list = ls(all.names = TRUE))

## Libraries needed for spatial stuff (others for predict functions called as needed below)
library(raster)
library(maptools)
library(rgdal)
library(rgeos)

## SDM models 
GLM.mod <- get(load("SDM/Output/GLM.mod2.Rda"))
GAM.mod <- get(load("SDM/Output/GAM.mod2.Rda"))
RF.mod <- get(load("SDM/Output/RF.mod2.Rda"))
BRT.mod <- get(load("SDM/Output/BRT.mod3.Rda"))
MAX.mod <- get(load("SDM/Output/MAX.mod.Rda"))

## Projections
prj.wgs = "+proj=longlat + type=crs"
prj.lcc <- "+proj=lcc +lon_0=-95 +lat_1=49 +lat_2=77 +type=crs"

## Bioclim grids
preds.lcc <- stack("SDM/data_files/bio2.grd", # These are in LCC projection
                   "SDM/data_files/bio3.grd", 
                   "SDM/data_files/bio10.grd", 
                   "SDM/data_files/bio11.grd", 
                   "SDM/data_files/bio12.grd", #for some reason these last 3 
                   "SDM/data_files/bio15.grd", #layers lose their names
                   "SDM/data_files/bio17.grd") #upon import
#so rename here; variable names must match models
names(preds.lcc) <- c("bio2", "bio3", "bio10", "bio11", "bio12", "bio15", "bio17")
proj4string(preds.lcc) <- CRS(prj.lcc) #define current lcc projection of rasters
preds = projectRaster(preds.lcc, crs=CRS(prj.wgs)) #transform raster projection to match pres/abs models

## State polygons for pretty maps
# All of USA
sta = readOGR("SDM/data_files/gz_2010_us_040_00_500k/gz_2010_us_040_00_500k.shp")
projection(sta) = CRS(prj.wgs)
# Define extent of study area
clim <- read_csv("SDM/data_files/points_Normal_1961_1990MSY.csv")
ext <- extent(min(clim$Longitude)-1, max(clim$Longitude)+1, min(clim$Latitude)-1, max(clim$Latitude)+1)
bbox = as(ext, "SpatialPolygons") #convert coordinates to a bounding box
# Crop state lines to study area
sta.crop <- crop(sta, bbox)
sta.lcc = spTransform(sta.crop, CRS=CRS(prj.lcc))

## Presence/absence points for pretty maps
presences <- clim %>% filter(ID2==1)
absences <- clim %>% filter(ID2==0)
coordinates(presences) <- ~Longitude+Latitude #convert to spatial data
projection(presences) <- CRS('+proj=longlat') #define projection
presences.lcc <- spTransform(presences, CRS=CRS(prj.lcc)) #transform projection so points layer matches SDM projections
coordinates(absences) <- ~Longitude+Latitude #convert to spatial data
projection(absences) <- CRS('+proj=longlat') #define projection
absences.lcc <- spTransform(absences, CRS=CRS(prj.lcc)) #transform projection so points layer matches SDM projections

################################################################################


################################################################################
### Project SDM model predictions across study area

## GLM
pred.glm <- raster::predict(preds, GLM.mod, type="response")
plot(pred.glm) #this is a really shitty model... does not match cardinalis distribution

## GAM
library(gam) #must load library for correct predict algorithm
pred.gam <- raster::predict(preds, GAM.mod, type="response")
plot(pred.gam) #also a shitty model... does not match cardinalis distribution

## RF
library(randomForest) #must load library for correct predict algorithm
pred.rf <- raster::predict(preds, RF.mod, type="response")
plot(pred.rf) #not sure why this looks binary

## BRT
library(gbm) #must load library for correct predict algorithm
pred.brt <- raster::predict(preds, BRT.mod, type="response")
plot(pred.brt)

## MAX
library(dismo)
pred.max <- raster::predict(preds, MAX.mod, type="response")
plot(pred.max)

################################################################################


################################################################################
### Calculate ensemble average of model predictions

## Stack predictions from each model
pred.all <- stack(pred.glm, pred.gam, pred.rf, pred.brt, pred.max)

## Calculate simple ensemble based on arithmetic average
sim.ensem <- mean(pred.all) #unprojected
plot(sim.ensem)
writeRaster(sim.ensem, file="SDM/Output/UnweightedEnsemble_Unprojected.grd", overwrite=TRUE)
sim.ensem.lcc <- projectRaster(sim.ensem, crs=CRS(prj.lcc)) #convert to Lambers conic projection
plot(sim.ensem.lcc)
writeRaster(sim.ensem.lcc, file="SDM/Output/UnweightedEnsemble_LCCProjection.grd", overwrite=TRUE)

## Calculate weighted ensemble based on AUC score of each model
# (This gives more weight to more accurate models)

# Read in saved AUC scores for each model
glm.accs <- get(load("SDM/Output/GLM.mod2.cvaccs.Rda")) 
glm.auc <- glm.accs$AUC[1]

gam.accs <- get(load("SDM/Output/GAM.mod2.cvaccs.Rda")) 
gam.auc <- gam.accs$AUC[1]

rf.accs <- get(load("SDM/Output/RF.mod2.accs.Rda")) 
rf.auc <- rf.accs$AUC[1]

brt.accs <- get(load("SDM/Output/BRT.mod3.accs.Rda")) 
brt.auc <- brt.accs$AUC[1]

max.accs <- get(load("SDM/Output/MAX.mod.accs.Rda")) 
max.auc <- max.accs$V1[5]

# Compile AUCs into a vector of weights
weights <- c(glm.auc, gam.auc, rf.auc, brt.auc, max.auc)

# Calculate weighted  average
wtd.ensem <- weighted.mean(pred.all, weights)
plot(wtd.ensem) #unprojected
writeRaster(wtd.ensem, file="SDM/Output/WeightedEnsemble_Unprojected.grd", overwrite=TRUE)
wtd.ensem.lcc <- projectRaster(wtd.ensem, crs=CRS(prj.lcc)) #convert to Lambers conic projection
plot(wtd.ensem.lcc)
writeRaster(wtd.ensem.lcc, file="SDM/Output/WeightedEnsemble_LCCProjection.grd", overwrite=TRUE)

################################################################################


################################################################################
### Pretty map

## Set up gridlines & lat/lon labels	
frame.grd = gridlines(sta.crop)
frame.grd.lcc <- spTransform(frame.grd, CRS=CRS(prj.lcc))
gridatt <- gridat(frame.grd, side="EN")
gridat.lcc = spTransform(gridatt, CRS=CRS(prj.lcc))

## Set up color ramp 
library(colorspace)
rbPal <- diverge_hcl(10)

## Save plot
pdf(file="SDM/figures/MAP_EnsemProb.pdf", width=11, height=8.5)
plot(wtd.ensem.lcc, box=FALSE, axes=FALSE, legend=FALSE, col=rbPal) #the ensemble layer
plot(presences.lcc, pch=1, cex=0.8, add=TRUE) #add presence points
plot(absences.lcc, pch=4, col="grey40", cex=0.5, add=TRUE) #add pseudoabsence points
plot(sta.lcc, add=T) #add state lines
plot(frame.grd.lcc, add=TRUE, lty="dashed", col="darkgrey", lwd=1) #add gridlines
text(coordinates(gridat.lcc), labels=parse(text=as.character(gridat.lcc$labels)), pos=gridat.lcc$pos, offset=0.5, col="black", cex=0.7) #add lat-long labels to gridlines
legend("bottomleft", legend="Weighted Ensemble", bty="n", cex=1.5) #add title
plot(wtd.ensem.lcc, legend.only=TRUE, legend.width=1, legend.shrink=0.75, col=rbPal, axis.args=list(at=seq(0, 1, by=0.1), labels=seq(0, 1, by=0.1), cex.axis=0.8)) #add legend for color ramp
dev.off()

################################################################################
