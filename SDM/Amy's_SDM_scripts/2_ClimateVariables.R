###################################################################################
## PROGRAM FUNCTIONS: EXRACT CLIMATIC VARIABLES FOR POINT COORDINATES!
###################################################################################
# EXRACT CLIMATIC VARIABLES FOR POINT COORDINATES:
# Extract climate variables from point coordinates from ClimateWNA and convert to bioclim variables
# Climate records are averaged from a set time period ie. 30 years prior 
# to the date of collection. ClimateWNA run initially with the time series function 
# from 1901 - 2012 for all records then collapsed and averaged using the script below.
# CLIMATEWNA run with DEM 90m res obtained from Hydrosheds
# (http://hydrosheds.cr.usgs.gov/index.php)
# A 30m ultra high res DEM is also available from another source, but not necessary here
# June 26th, 2014 Revisions. 


## set pathnames
path.root="/Users/amyangert/Desktop/cardinalis SDM" 
path.dat=paste(path.root,"/data files",sep="") #(normals calculated from date of collection)

## set pathnames - Matthew
path.root = "C:/Users/DW/Desktop/temp.sept.30" 
path.dat = paste(path.root, "/data files", sep="")
path.obj = paste(path.root, "/R objects", sep="")
path.gis = paste(path.dat, "/ecoregions.shp", sep="")


################################################################################	
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
######## START INITIALIZATION

########################### 
## set directories
setwd(path.root)

########################################################
## Source Climate Data according to the year of record 
## ie. 20 year avg from a record collected from 1979 sources 1959 - 1979
## (Matthew Bayly: April 7, 2014)
########################################################
## Use the time series option in ClimateWNA to get annual
## records for the total range of the data (ie. 1947 - 2012)
## in my case the file size was large so it had to be broken
## into two files according to 1/2 of the records 
##
library(plyr)
library(dismo)

#first and second refer to two large files from a size issue, 
#ignore if you have all records in a single dataframe

# OCT 30 2014 ~ ALREADY DONE FOR RECORDS SKIP
# PART 1: Extract elevation values for records
# elevation values for records
#attach elevation values to spreadhseet 
#dem <- raster(file.choose()) # load 90m DEM
#setwd(path.dat)
#all.pres <- read.csv(file='all.pres.csv')
#head(all.pres); summary(all.pres); dim(all.pres)
#coordinates(all.pres) <- ~LON.climate+LAT.climate
#projection(all.pres) <- CRS('+proj=lonlat');
#plot(dem); points(all.pres)# check projection & spread of points
#plot(all.pres) # how many points are left down in Mexico 
#ex <- extract(dem, all.pres) # extract coordinate values from raster stack 
#ex <- data.frame(ex)
#sitevar <- cbind(all.pres, ex) 
# select only records that fall within the bounds of this layer 
# remove mexico records 
#sitevar[!is.na(sitevar$ex),]
#write.csv(sitevar, file = "all.pres.elev.csv") # save file to drive
#########################################################
# Export & run in climateWNA


# PART 2: Load time.series records & temporarily split < 1930 records
# run climateWNA! external to this script.
# since the earliest data from climatewna is 1901 need to average all records 
# from < 1930 to a time period of 1901 - 1930.
# For each record extract time series specific interval 
# first second refer to two XL files, if possible skip next five lines to total.all
#first <- read.csv(file.choose())
#second <- read.csv(file.choose())
#dim(first); dim(second) # do columns match?
#colnames(first[,46:51]); colnames(second[,46:51])
#total.all <- rbind(first, second) # merge data frames
	
# load file from climateWNA
	setwd(path.dat)
	total.all <- read.csv(file=paste(path.dat, "/old/", "to_climateWNA_1901-2012AMT.csv", sep=""))
# have to temporarily add in this annoying catechism otherwise 
# ID's such as herb.4, herb.40 & herb.400 will not be seen as unique by R
	total.all$Merge <- paste(total.all$ID1, "temp", sep=":")
	total.cut <- total.all[which(total.all$ID2 > 1930), ]

dim(total.all) - dim(total.cut) # 1008/112 = 9 records
# we will deal with these records later
# need to come back and revisit records prior to 1930 to average 1900 - 1930 

# PART 3: Use a 'difference' column to get desired time interval
# diff = YEAR OF COLLECTION "-" TIME SERIES YEAR
diff <- (total.cut$ID2 - total.cut$Year) # this will help select the time interval
total.cut <- cbind(total.cut, diff) # 
#tail(total.cut); summary(total.cut$diff)
summary(total.cut$diff)

# For records earlier than 1930, just use the average period from 1902 - 1930
remaining <- total.all[which(total.all$ID2 < 1931), ] # remaining records 
yr30b <- remaining[which(remaining$Year < 1931 & remaining$Year > 1900), ] # select period 1901 - 1930
summary(yr30b$ID2); dim(yr30b) # about seven records


# PART 4: SET THE TIME INTERVAL (ie. 5yr, 15yr, 30yr?)
n <- 30 # how long should the average be
yr30 <-  total.cut[which(total.cut$diff < (n+1) & total.cut$diff > -1), ] # 30 YEAR AVERAGE
summary(yr30$diff); dim(yr30) # should be 0 - interval specified above (30)
#tail(yr30[, c(1:7, 56)], 90) # good visual check

# PART 5: Average down the columns for final datasheet
# CHECK for possible errors prior
# herb.8 not unique from herb.80?
# Average down the columns from the specified time interval for each variable.
library(plyr)
	
	# for all records >1930
	yr30avg <- daply(yr30[, c(1:54)], .(yr30$Merge), colwise(mean)) 
	yr30avg <- data.frame(yr30avg)
	rownamesR <- rownames(yr30avg)
	yr30avg$ID1 <- rownamesR
	yr30avg <- data.frame(lapply(yr30avg, as.numeric), stringsAsFactors=FALSE)
	yr30avg$ID1 <- rownamesR

	# for all records <1930
	yr30avgb <- daply(yr30b[, c(1:54)], .(yr30b$Merge), colwise(mean)) 
	yr30avgb <- data.frame(yr30avgb)
	rownamesR <- rownames(yr30avgb)
	yr30avgb$ID1 <- rownamesR
	yr30avgb <- data.frame(lapply(yr30avgb, as.numeric), stringsAsFactors=FALSE)
	yr30avgb$ID1 <- rownamesR
	
#tail(yr30avg[ ,1:7], 200)
#head(yr30avg[ ,1:7], 200)
dim(yr30avg); dim(yr30avgb)
colnames(yr30avg); colnames(yr30avgb)
final <- rbind(yr30avg, yr30avgb)

# temporarily save output, need to remove comment 
setwd(paste(path.dat, "/old", sep=""))
write.csv(final, file="yr30avg.csv", row.names=FALSE)


#  YOU NOW HAVE YEAR-SPECIFIC CLIMATE AVERAGED Data!
############################################################
############################################################
############################################################
############################################################

#PART 6: Convert monthly variables to Biovars!
# Then convert variables to  the 19 bioclim vars
#?biovars
setwd(paste(path.dat, "/old", sep=""))
climateWNA <- read.csv(file='yr30avg.csv')
# remove na values 8 records prior to 1980
dim(climateWNA)
climateWNA <- climateWNA[complete.cases(climateWNA), ] # there were none
dim(climateWNA)
colnames(climateWNA)

#CHECK dataframe is set up properly!

colnames(climateWNA)
colnames(climateWNA[7]) # should be 'Tmax01'
colnames(climateWNA[18]) # should be 'Tmax12'
tmax <- climateWNA[ , 7:18]
colnames(climateWNA[19]) # should be Tmin01
colnames(climateWNA[30])# should be Tmin12
tmin <- climateWNA[ , 19:30]
colnames(climateWNA[43]) # PPT01
colnames(climateWNA[54]) # PPT12
prec <- climateWNA[ , 43:54]
library(dismo)
tmin <- as.matrix(tmin)
tmax <- as.matrix(tmax)
prec <- as.matrix(prec)
b <- biovars(prec, tmin, tmax)
b <- data.frame(b)
dim(b); colnames(b)

b3 <- cbind(climateWNA, b)
dim(b3); head(b3, 3)

setwd(paste(path.dat, "/old", sep=""))
write.csv(b3, file=paste("bios",n, ".csv", sep=""), row.names=FALSE)
####################################################
####################################################
####################################################

