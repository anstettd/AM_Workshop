###################################################################################
## Amy Angert
## PROGRAM FUNCTIONS: 
	# Thin pseudoabsence records to match thinned replicates of herbarium dataset
## last update:  29 May 2014

## Amy notes 11 Dec 2020: simplifications and changes needed for AM workshop
  # Delete all setwd paths to use relative paths within R project
  # Add missing code to draw pseudoabsences in the first place
  # Randomly select 1 thinned presence file so that everything is not done across i=1:10 replicates
  # No need to drop occupancy records; model can be built with full presence dataset
  # Add instructions for downloading climate records and delete merges below where climate variables are pulled in from another file

###################################################################################



################################################################################
######## START INITIALIZATION

## set pathnames
#path.root = "/Users/amylauren/Desktop/cardinalis SDM May 2014" 
#path.dat = paste(path.root, "/data files", sep="")
#path.obj = paste(path.root, "/R objects", sep="")
#path.gis = paste(path.dat, "/ecoregions.shp", sep="")

## set pathnames - Matthew
#path.root = "C:/Users/DW/Desktop/temp.sept.30" 
#path.dat = paste(path.root, "/data files", sep="")
#path.obj = paste(path.root, "/R objects", sep="")
#path.gis = paste(path.dat, "/ecoregions.shp", sep="")


## Import localities
#setwd(path.dat)
library(sp)
library(raster)
library(rgeos)

for (i in 1:10) {
	dat <- read.csv(paste("SDM/data_files/occurrences.thinned",i,".csv", sep=""))
	coordinates(dat) <- ~longitude+latitude
	projection(dat) <- CRS('+proj=longlat')
	dat$rowID <- dat$Species
	assign(paste("herb",i, sep=""), dat)
	}
plot(herb1, pch=19, cex=0.5, col="red")
points(herb2, pch=19, cex=0.5, col="blue")
points(herb3, pch=19, cex=0.5, col="green")
plot(herb8, pch=19, cex=0.5, col="red")
points(herb9, pch=19, cex=0.5, col="green")
points(herb10, pch=19, cex=0.5, col="blue")
	
all = read.csv("SDM/data_files/all.records.aug.31.csv") #includes occupancy dataset, cleaned herbarium records, and 20K pseudoabs
herb.all = all[all$DATASET=="herb",] #drop occupancy points
pseudo = herb.all[herb.all$PRESABS==0,] #pull out pseudos
dim(pseudo); table(pseudo$PRESABS)
coordinates(pseudo) = ~Longitude + Latitude
projection(pseudo) = CRS('+proj=longlat')

coordinates(herb.all) = ~Longitude + Latitude
projection(herb.all) = CRS('+proj=longlat')


######## END INITIALIZATION
################################################################################




################################################################################
######## START SUBSAMPLE PSEUDOS
library(sp)

## 10:1 ratio abs:pres for GLM and GAM
for (i in 1:10) {
	herb = get(paste("herb",i,sep=""))
	npsu = 10*length(herb) 
	ok = zerodist2(pseudo, herb, 80)
	drop = zerodist2(pseudo, herb, 0.600)
	psu1 = pseudo[unique(ok[,1]),]
	psu2 = psu1[-drop[,1],]
	samp = sample(1:nrow(psu2), npsu, replace=FALSE) 	
	psu2 = psu2[samp,]
	assign(paste("pseu",i,sep=""), psu2)
	}

## 4:1 ratio abs:pres for RF and MAX
for (i in 1:10) {
	herb = get(paste("herb",i,sep=""))
	npsu = 4*length(herb) 
	ok = zerodist2(pseudo, herb, 80)
	drop = zerodist2(pseudo, herb, 0.600)
	psu1 = pseudo[unique(ok[,1]),]
	psu2 = psu1[-drop[,1],]
	samp = sample(1:nrow(psu2), npsu, replace=FALSE) 	
	psu2 = psu2[samp,]
	assign(paste("pseu",i,"b",sep=""), psu2)
	}

## 1:1 ratio abs:pres for BRT
for (i in 1:10) {
	herb = get(paste("herb",i,sep=""))
	npsu = length(herb) 
	ok = zerodist2(pseudo, herb, 80)
	drop = zerodist2(pseudo, herb, 0.600)
	psu1 = pseudo[unique(ok[,1]),]
	psu2 = psu1[-drop[,1],]
	samp = sample(1:nrow(psu2), npsu, replace=FALSE) 	
	psu2 = psu2[samp,]
	assign(paste("pseu",i,"c",sep=""), psu2)
	}

######## END SUBSAMPLE PSEUDOS
################################################################################

################################################################################
######## START MERGE PRES & PSEUDOABS
#setwd(path.dat)
for (i in 1:10) {
	herb = get(paste("herb",i,sep="")) #get thinned herb set
#Daniel: merge is a function in base R, sp, and raster, I assume you wanted merge in base R
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
	}
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
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	