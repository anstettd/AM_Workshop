################################################################################# 
### SCRIPT PURPOSE: Generate pseudoabsence records 
  # Based on methods used in Angert et al. 2018, Am Nat
  # Author: Amy Angert
  # last update:  11 Dec 2020

## Simplifications and changes made for AM workshop
  # Converted to tidyverse syntax (unless manipulating spatial files)
  # Using relative paths within R project instead of setwd commands
  # Chose 1 thinned presence input file so that everything is not done across i=1:10 replicates (using i=6 from Am Nat paper, chosen at random)

## Still to do
  # Add missing code to draw background points in the first place
  # Build model with full presence dataset instead of dropping occupancy points?
  # Add instructions for downloading climate records
  # Delete merges below where climate variables are pulled in from another pre-existing file

###############################################################################



############################################################################### 
### LOAD LIBRARIES AND INPUTS

## LIBRARIES
library(tidyverse)
library(sp)
library(raster)
library(rgeos)

## INPUTS
# read in presence records
  # these are from herbarium collections, thinned to reduce oversampled areas
    dat <- read_csv("SDM/data_files/occurrences.thinned.csv")
  # rename column
    herb <- dat %>% rename(rowID = Species)
  # convert to spatial points file
    coordinates(herb) <- ~longitude+latitude
    projection(herb) <- CRS('+proj=longlat')
    
  # these include occupancy dataset, all cleaned herbarium records, and 20K pseudoabs
    all = read_csv("SDM/data_files/all.records.aug.31.csv")
  # drop occupancy survey points
    herb.all <- all %>% filter(DATASET=="herb")
  # drop presences, retain only pseudos that were drawn previously 
  # ***THIS NEEDS TO CHANGE***
    pseudo <- all %>% filter(DATASET=="herb" & PRESABS==0)
  # convert to spatial points files
    coordinates(herb.all) = ~Longitude + Latitude
    projection(herb.all) = CRS('+proj=longlat')
    coordinates(pseudo) = ~Longitude + Latitude
    projection(pseudo) = CRS('+proj=longlat')

################################################################################


################################################################################
### DRAW 20K BACKGROUND POINTS
    
## FIGURE OUT WHERE THIS CODE IS...    
################################################################################
    

################################################################################
### SUBSAMPLE FROM BACKGROUND POINTS TO MAKE PSEUDOABSENCES 

## RATIONALE FOR STEP 1: Pseudoabsences should fall in sampled regions yet not in occupied habitat. Here, the zerodist2 function is used to sample background points from within donut-shaped areas surrounding each presence record, so that resulting pseudoabsences are >0.6 km but <80 km from a presence.

# Keep any background points within 80 km of a presence record
  ok = zerodist2(pseudo, herb, 80)
  psu1 = pseudo[unique(ok[,1]),] # not using tidyverse here - doesn't work on spatial data frames

# Drop any background points that are closer then 0.6 km from a presence record
  drop = zerodist2(pseudo, herb, 0.600)
  psu2 = psu1[-drop[,1],]
    
## RATIONALE FOR STEP 2: Different ENM algorithms are optimized for different ratios of presence records relative to pseudoabsences.
    
# 10:1 ratio abs:pres for GLM and GAM algorithms
npsu.glm = 10*length(herb) 
samp.glm = sample(1:nrow(psu2), npsu.glm, replace=FALSE) 	
psu2.glm = psu2[samp.glm,]
assign(paste("pseu","a",sep=""), psu2.glm)

# 4:1 ratio abs:pres for RF and MAX algorithms
npsu.rf = 4*length(herb) 
samp.rf = sample(1:nrow(psu2), npsu.rf, replace=FALSE) 	
psu2.rf = psu2[samp.rf,]
assign(paste("pseu","b",sep=""), psu2.rf)

# 1:1 ratio abs:pres for BRT algorithm
npsu.brt = length(herb) 
samp.brt = sample(1:nrow(psu2), npsu.brt, replace=FALSE) 	
psu2.brt = psu2[samp.brt,]
assign(paste("pseu","c",sep=""), psu2.brt)

################################################################################

## stopped here 12/11
################################################################################
## MERGE PRESENCES & PSEUDOABSENCES

herb.clim = base::merge(herb.all, herb, by.x="MASTER.ID", by.y="rowID", all.x=FALSE, all.y=TRUE) #associate with climatic variables
	herb.clim = herb.clim[,c(1:75)] #drop cols not present in pseu
	#names(herb.clim)[1] = "ID1"
	#names(herb.clim)[2] = "ID2" 
	herb.clim$PRESABS = 1
	pseu = get(paste("pseu",i,sep="")) #get matching pseu set
	pseu = as.data.frame(pseu) #remove spatial class so lat and lon merge
#Daniel: Change Column Names From Longitude and Latitude to MASTER.ID and DATASET 	
	names(pseu)[1] = "MASTER.ID"
	names(pseu)[2] = "DATASET"
	pseu[,c(3,4,5,6,7,1,2,8:75)]
	pseu$PRESABS = 0
	pseu$MASTER.ID <- pseu$ID 
	order_cols <- colnames(herb.clim)
	pseu <- pseu[,c(order_cols)]
	presabs = rbind(herb.clim, pseu) #bind pres & psedos into one frame
	presabs$bio3 = log(presabs$bio3+0.5) #ln-transform
	presabs$bio10 = log(presabs$bio10+0.5)
	presabs$bio12 = log(presabs$bio12+0.5)
	presabs$bio14 = log(presabs$bio14+0.5)
	assign(paste("dat",i,sep=""), presabs) #save
	write.csv(presabs, file=paste("SDM/Output/dat",i,".csv", sep=""))

#check numbers
table(dat1$PRESABS)
-table(dat5$PRESABS)
table(dat9$PRESABS)
#summary(dat9)


# for the 1:4 ratio
#setwd(path.dat)
for (i in 1:10) {
	herb = get(paste("herb",i,sep="")) #get thinned herb set
	herb.clim = base::merge(herb.all, herb, by.x="MASTER.ID", by.y="rowID", all.x=FALSE, all.y=TRUE) #associate with climatic variables
	herb.clim = herb.clim[,c(1:75)] #drop cols not present in pseu
	#names(herb.clim)[1] = "ID1"
	#names(herb.clim)[2] = "ID2" 
	herb.clim$PRESABS = 1
	pseu = get(paste("pseu",i,"b",sep="")) #get matching pseu set
	pseu = as.data.frame(pseu) #remove spatial class so lat and lon merge
	names(pseu)[1] = "MASTER.ID"
	names(pseu)[2] = "DATASET"
	pseu[,c(3,4,5,6,7,1,2,8:75)]	
	pseu$PRESABS = 0
	pseu$MASTER.ID <- pseu$ID 
	order_cols <- colnames(herb.clim)
	pseu <- pseu[,c(order_cols)]
	presabs = rbind(herb.clim, pseu) #bind pres & psedos into one frame
	presabs$bio3 = log(presabs$bio3+0.5) #ln-transform
	presabs$bio10 = log(presabs$bio10+0.5)
	presabs$bio12 = log(presabs$bio12+0.5)
	presabs$bio14 = log(presabs$bio14+0.5)
	assign(paste("dat",i,"b",sep=""), presabs) #save
	write.csv(presabs, file=paste("SDM/Output/dat",i,"b.csv", sep=""))
	}
#check numbers
table(dat1b$PRESABS)
-table(dat5b$PRESABS)
table(dat9b$PRESABS)
#summary(dat9b)



# LAST for the 1:1 ratio
#setwd(path.dat)
for (i in 1:10) {
	herb = get(paste("herb",i,sep="")) #get thinned herb set
	herb.clim = base::merge(herb.all, herb, by.x="MASTER.ID", by.y="rowID", all.x=FALSE, all.y=TRUE) #associate with climatic variables
	herb.clim = herb.clim[,c(1:75)] #drop cols not present in pseu
	#names(herb.clim)[1] = "ID1"
	#names(herb.clim)[2] = "ID2" 
	herb.clim$PRESABS = 1
	pseu = get(paste("pseu",i,"c",sep="")) #get matching pseu set
	pseu = as.data.frame(pseu) #remove spatial class so lat and lon merge
	names(pseu)[1] = "MASTER.ID"
	names(pseu)[2] = "DATASET"
	pseu$PRESABS = 0
	pseu$MASTER.ID <- pseu$ID 
	order_cols <- colnames(herb.clim)
	pseu <- pseu[,c(order_cols)]
	presabs = rbind(herb.clim, pseu) #bind pres & psedos into one frame
	presabs$bio3 = log(presabs$bio3+0.5) #ln-transform
	presabs$bio10 = log(presabs$bio10+0.5)
	presabs$bio12 = log(presabs$bio12+0.5)
	presabs$bio14 = log(presabs$bio14+0.5)
	assign(paste("dat",i,"c",sep=""), presabs) #save
	write.csv(presabs, file=paste("SDM/Output/dat",i,"c.csv", sep=""))
	}
#check numbers
table(dat1c$PRESABS)
-table(dat5c$PRESABS)
table(dat9c$PRESABS)
-table(dat9c$PRESABS)


#summary(dat9c)






######## END MERGE PRES & PSEUDOABS
################################################################################


######################
#######Check for errors or other problems
################################################################################
#OPTIONAL 
#OPTIONAL 
#OPTIONAL 
	for (i in 1:10) {
	dat <- read.csv(paste("dat",i,".csv", sep=""))
	coordinates(dat) <- ~Longitude+Latitude
	projection(dat) <- CRS('+proj=longlat')
	assign(paste("dat",i, sep=""), dat)
	}
	
	for (i in 1:10) {
	dat <- read.csv(paste("dat",i,"b.csv", sep=""))
	coordinates(dat) <- ~Longitude+Latitude
	projection(dat) <- CRS('+proj=longlat')
	assign(paste("dat",i,"b", sep=""), dat)
	}
	
	for (i in 1:10) {
	dat <- read.csv(paste("dat",i,"c.csv", sep=""))
	coordinates(dat) <- ~Longitude+Latitude
	projection(dat) <- CRS('+proj=longlat')
	assign(paste("dat",i,"c", sep=""), dat)
	}
	
	
	## plot points to visualize
	#initial set of all records
	herb.pres = herb.all[herb.all$PRESABS==1,] #
	plot(herb.pres, pch=19, cex=0.5, col='green')
	allx <- mean(herb.pres$Longitude); ally <- mean(herb.pres$Latitude)	
	points(allx, ally, pch = 3, cex=3, col='green', lwd=4)
	
	par(mfrow=c(2,5))
	for (i in 1:10) {
	dat = get(paste("dat",i,"c", sep=""))
	pres = dat[dat$PRESABS==1,]
	abs = dat[dat$PRESABS==0,]
	plot(abs, pch=19, cex=0.5, col="black")
	points(pres, pch=21, cex=0.5, col="red")
	px <- mean(pres$Longitude); ax <- mean(abs$Longitude)
	py <- mean(pres$Latitude); ay <- mean(abs$Latitude)
	points(px, py, pch = 3, cex=3, col='red', lwd=4); points(ax, ay, pch = 3, cex=3, col='black', lwd=4)
	points(allx, ally, pch = 3, cex=3, col='green', lwd=4)
	}
	par(mfrow=c(1,1))
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	