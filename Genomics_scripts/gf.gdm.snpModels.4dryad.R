################################################################################
# Scripts to model & map genetic turnover using gradient forests and generalized
# dissimilarity modeling described in Fitzpatrick & Keller (in press) Ecology Letters
#
# written by MC Fitzpatrick at the Appalachian Lab, Frostburg, MD 2013-2014
# with input on mapping functions from M. Lisk
#
# Code is provided as is, without support 
################################################################################

################################################################################
# BEGIN SCRIPTS FOR GRADIENT FOREST MODELING - GDM SCRIPTS BEGINS @ LINE 148
################################################################################
# load libraries, read in data, etc. -------------------------------------------
library(gradientForest)
library(raster)

# read in data file with minor allele freqs & env/space variables
gfData <- read.csv("poplarSNP.ENV.data.4.gradientForest.csv")
envGF <- gfData[,3:13] # get climate & MEM variables

# build individual SNP datasets
SNPs_ref <- gfData[,grep("REFERENCE",colnames(gfData))] # reference
GI5 <- gfData[,grep("GI5",colnames(gfData))] # GIGANTEA-5 (GI5)
################################################################################


# GRADIENT FOREST MODELING -----------------------------------------------------
maxLevel <- log2(0.368*nrow(envGF)/2) #account for correlations, see ?gradientForest 

# Fit gf models for reference SNPs 
gfRef <- gradientForest(cbind(envGF, SNPs_ref), predictor.vars=colnames(envGF),
                         response.vars=colnames(SNPs_ref), ntree=500, 
                         maxLevel=maxLevel, trace=T, corr.threshold=0.50)

# Fit gf models for GI5 SNPs
gfGI5 <- gradientForest(cbind(envGF, GI5), predictor.vars=colnames(envGF),
                          response.vars=colnames(GI5), ntree=500, 
                          maxLevel=maxLevel, trace=T, corr.threshold=0.50)

# plot output, see ?plot.gradientForest
type = "O"
plot(gfRef, plot.type=type)
plot(gfGI5, plot.type=type)
################################################################################


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

# Function to map difference between spatial genetic predictions
# predMap1 = dataframe of transformed variables from gf or gdm model for first set of SNPs
# predMap2 = dataframe of transformed variables from gf or gdm model for second set of SNPs
# rast = a raster mask to which Procrustes residuals are to be mapped
# mapCells = cell IDs to which Procrustes residuals values should be assigned
RGBdiffMap <- function(predMap1, predMap2, rast, mapCells){
  require(vegan)
  PCA1 <- prcomp(predMap1, center=TRUE, scale.=FALSE)
  PCA2 <- prcomp(predMap2, center=TRUE, scale.=FALSE)
  diffProcrust <- procrustes(PCA1, PCA2, scale=TRUE, symmetrical=FALSE)
  residMap <- residuals(diffProcrust)
  rast[mapCells] <- residMap
  return(list(max(residMap), rast))
}

# OK, on to mapping. Script assumes:
# (1) a dataframe named env_trns containing extracted raster data (w/ cell IDs)
# and env. variables used in the models & with columns as follows: cell, bio1, bio2, etc.
#
# (2) a raster mask of the study region to which the RGB data will be written

# transform env using gf models, see ?predict.gradientForest
predRef <- predict(gfRef, env_trns[,-1]) # remove cell column before transforming
predGI5 <- predict(gfGI5, env_trns[,-1])

# map continuous variation - reference SNPs
refRGBmap <- pcaToRaster(predRef, mask, env_trns$cell)
plotRGB(refRGBmap)
writeRaster(refRGBmap, "/.../refSNPs_map.tif", format="GTiff", overwrite=TRUE)

# map continuous variation - GI5 SNPs
GI5RGBmap <- pcaToRaster(predGI5, mask, env_trns$cell)
plotRGB(refRGBmap)
writeRaster(refRGBmap, "/.../GI5SNPs_map.tif", format="GTiff", overwrite=TRUE)

# Difference between maps (GI5 and reference) 
diffGI5 <- RGBdiffMap(predRef, predGI5, rast=mask, mapCells=env_trns$cell)
plot(diffGI5[[2]])
writeRaster(diffGI5[[2]], "/.../diffRef_GI5.tif", format="GTiff", overwrite=TRUE)
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
################################################################################
# END SCRIPTS FOR GRADIENT FOREST MODELING
################################################################################

#-------------------------------------------------------------------------------

################################################################################
# BEGIN SCRIPTS FOR GENERALIZED DISSIMILARITY MODELING
################################################################################
# load libraries, read in data, etc. -------------------------------------------
library(Gdm01)
library(raster)
# read in data file with Fst & env variables (in GDM site-pair format)
# note Fst values were scaled 0-1 to facilitate model convergence
gdmData <- read.csv("poplarFst.ENV.data.4.GDM.csv")

# build individual SNP datasets
SNPs_ref <- gdmData[,c(1,6:24)] # reference
GI5 <- gdmData[,c(3,6:24)] # GIGANTEA-5 (GI5)
################################################################################


# Fit generalized dissimilarity models -----------------------------------------
GEO=T # use geographic distance as a predictor?
# reference SNPs
gdmRef <- gdm.fit(SNPs_ref, geo=GEO)
gdmRef$explained
gdm.plot(gdmRef, plot.layout=c(3,4))
refSplines <- isplineExtract(gdmRef) # extract spline data for custom plotting

# GI5 SNPs
gdmGI5 <- gdm.fit(GI5, geo=GEO)
gdmGI5$explained
gdm.plot(gdmGI5, plot.layout=c(3,4))
GI5Splines <- isplineExtract(gdmGI5) # extract spline data for custom plotting
################################################################################


# Mapping spatial genetic variation --------------------------------------------

# Follows identical procedure as above for gf, but uses ?gdm.transform instead
# Also note that if geo. dist. is a predictor in the model, x & y rasters
# can be supplied in addition to env. layers.
predRef <- gdm.transform(gdmRef, env_trns[,-1]) # remove cell column
predGI5 <- gdm.transform(gdmGI5, env_trns[,-1])

# map continuous variation - reference SNPs
refRGBmap <- pcaToRaster(predRef, mask, env_trns$cell)
plotRGB(refRGBmap)
writeRaster(refRGBmap, "/.../refSNPs_map.tif", format="GTiff", overwrite=TRUE)

# map continuous variation - GI5 SNPs
GI5RGBmap <- pcaToRaster(predGI5, mask, env_trns$cell)
plotRGB(refRGBmap)
writeRaster(refRGBmap, "/.../GI5SNPs_map.tif", format="GTiff", overwrite=TRUE)

# Difference between maps (GI5 and reference) 
diffGI5 <- RGBdiffMap(predRef, predGI5, rast=mask, mapCells=env_trns$cell)
plot(diffGI5[[2]])
writeRaster(diffGI5[[2]], "/.../diffRef_GI5.tif", format="GTiff", overwrite=TRUE)
################################################################################

# Calculate and map "genetic offset" under climate change ----------------------
# Same approach as above for gf, but uses ?gdm.predict function, which calculates
# distance directly. These values can then be mapped.

genOffsetGI5 <- gdm.predict(gdmGI5, envPred)

# see ?gdm.predict for description of format of envPred.
# Basic structure is (e.g., for a model with geo=T, bio_1, elev):
# Observed  weights x.t1      y.t1    x.t2    y.t2  bio_1.t1  elev.t1   bio_1.t2  elev.t2
#   0.214     1     -152.25   65.42   -152.25 65.42 -5.6      235       -2.1      235   
# where Observed distances are required, but can be dummy values ranging 0-1

# assign values to raster - can be tricky if current/future climate
# rasters are not identical in terms of # cells, extent, etc.
mask[env_trns_future$cell] <- genOffsetGI5
plot(mask)
################################################################################
# END SCRIPTS FOR GENERALIZED DISSIMILARITY MODELING
################################################################################
