###################################################################################
## Amy Angert
## Modified from Thomas C. Edwards, Jr., "Species Distribution Modelling Using R"
## PROGRAM FUNCTIONS: 
##   Build set of MAXENT models with presence and user-input pseudoabsences 
##   Perform 5-fold cross validation and assess accuracy
##   MUST have downloaded maxent.jar from http://www.cs.princeton.edu/~schapire/maxent/
##     AND placed into C:/Users/YourUserNameHere/Documents/R/win-library/2.14/dismo/java/
##   MUST have x64 java from http://www.oracle.com/technetwork/java/javase/downloads/jdk7-downloads-1880260.html
##     AND modified Path in W7 system variables to add path to x64 java
##   See:  MakingMAXENTWorkInR.pdf for further details
## last update:  22 April 2016
###################################################################################

## load libraries now if desired; loaded below when needed
#library(dismo)
#library(raster)
#library(PresenceAbsence)
#library(maptools)

## set pathnames
path.root="/Users/amylauren/Google Drive/Occupancy AmNat Amy" 
path.dat=paste(path.root, "/data files", sep="")
path.obj=paste(path.root, "/R objects", sep="")
path.fig=paste(path.root, "/figures", sep="")
path.cod=paste(path.root, "/R code", sep="")

## set pathnames - Matthew
path.root = "C:/Users/DW/Desktop/temp.sept.30" 
path.dat = paste(path.root, "/data files", sep="")
path.dat.fix = paste(path.root, "/data files", sep="") # older files relocated to another directory
path.obj = paste(path.root, "/R objects", sep="")
path.eco = paste(path.obj, "/ecoregions.shp", sep="")
path.bio = paste(path.obj, "/wc0.5", sep="")
path.cod=paste(path.root, "/R code", sep="")
path.fig=paste(path.root, "/figures", sep="")

################################################################################
######## START INITIALIZATION
## variables chosen based on dev and collinearity: bio15, lnbio10, lnbio14, lnbio12, bio11, bio4, lnbio3, bio2
## replicate thinned datasets with pseudo:pres at 10:1
## some variables are already ln-transformed
## see file 'RCode_ThinPseudoAbs'
## now just call in saved .csv files

setwd(path.dat) 

#for (i in 1:10) {
#	dat = read.csv(paste("dat",i,".csv", sep=""))
#	dat = dat[,c(4,59:61,67:69,71:72)] #column indices changed 9/4/14
#	assign(paste("dat",i, sep=""), dat)
#	}
#rm(dat)

#for (i in 1:10) {
#	dat = read.csv(paste("dat",i,"b.csv", sep=""))
#	dat = dat[,c(4,59:61,67:69,71:72)] #column indices changed 9/4/14
#	assign(paste("dat",i, sep=""), dat)
#	}
#rm(dat)

##use these (1:1 ratio)
for (i in 1:10) {
	dat = read.csv(paste("dat",i,"c.csv", sep=""))
	dat = dat[,c(4,59:61,67:69,71:72)] #column indices changed 9/4/14
	assign(paste("dat",i, sep=""), dat)
	}
rm(dat)

## examine training data
dim(dat1); dim(dat10)
table(dat1$PRES); table(dat10$PRES)     
head(dat1); head(dat10)                     
str(dat1); str(dat10)                      

## testing data
all = read.csv("all.records.aug.31.csv") #includes occupancy dataset, cleaned herbarium records, and 20K pseudoabs drawn to match envir space of true absences
ext = all[all$DATASET=="occ",] #pull out occupancy dataset
ext$bio3 = log(ext$bio3+0.5) #make needed ln-transforms of predictors
ext$bio10 = log(ext$bio10+0.5)
ext$bio12 = log(ext$bio12+0.5)
ext$bio14 = log(ext$bio14+0.5)


######## END INITIALIZATION
################################################################################




################################################################################
######## {Do not re-run} START MAXENT MODEL: TRUE PRESENCES & BUFFERED-POINTS PSEUDOABSENCES
# Temporarily give MaxEnt more memory; options(java.parameters = "-Xmx1g" )
library(dismo)                         

# From MaxEnt sensitivity analysis (Amy -  Feb 2015): 
	# Opt for generality, maximizing AUC/AIC differences 
	# ENMeval run and from sensitivty tests, decided to adjust settings from default
		# FEATURES = LQHP (linear, quadratic, hinge, product)
		# REGULARIZATION MULTIPLIER = 2.5
		
setwd(path.obj)
for (i in 1:10) {
	## call up replicate training data
	dat = get(paste("dat",i, sep=""))
	# maxent settings adjusted 
	mod.MAX=maxent(dat[,c(2:9)],dat$PRESABS, args=c("betamultiplier=2.5", "threshold=FALSE")) 
	
	assign(paste("MAX.mod1.",i, sep=""), mod.MAX)
	save(mod.MAX, file=paste("MAX.mod1.",i,".pseudo11.Rda", sep=""))
	}

#plot(mod.MAX) # var importance plot
#response(mod.MAX) # prediction vs. var plot

######## END MAXENT MODEL: TRUE PRESENCES & BUFFERED-POINTS PSEUDOABSENCES
################################################################################




################################################################################
######## {START HERE} LOAD FINAL MODELS 

library(dismo)                         
setwd(path.obj)

for (i in 1:10) {
	mod = get(load(paste("MAX.mod1.",i,".pseudo11.Rda", sep="")))
	assign(paste("MAX.mod1.",i, sep=""), mod)
	}

######## LOAD FINAL MODELS
################################################################################




################################################################################
######## START RESUBSITUTION ACCURACY, MODEL = MAXENT

library(PresenceAbsence) 

setwd(path.cod)
source("accuracy.R")

accs=c()
for (i in 1:10) {
	## call up replicate training data and predictions
	dat = get(paste("dat",i, sep=""))
	mod = get(paste("MAX.mod1.",i, sep=""))
	modl = "mod.MAX" # add var to keep track of model
	temp = accuracy.max(dat, mod, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "MAX.mod1"
	accs = rbind(accs, temp)	
}	
setwd(path.obj)
save(accs, file="MAX.mod1.accs.pseudo11.Rda")
	
######## END INTERNAL ACCURACY & CLASSIFICATION 
################################################################################





################################################################################
######## START RESUBSTITUTION RELIABILITY CALCULATIONS, MODEL=MAX

setwd(path.cod)
source("calibration.R")

x=seq(0,1,0.05)
y=seq(0,1,0.05)
cal.MAX.training = as.data.frame(matrix(NA, 10,4))
names(cal.MAX.training) = c("int", "slope", "p_int", "p_slope")
setwd(path.fig)
pdf(file="MAX_CalPlots_Training_pseudo11.pdf", width=11, height=8.5)
par(mfrow=c(3,4))
setwd(path.obj)
for (i in 1:10) {
	## pull in replicate data and model predictions	
	dat = get(paste("dat",i, sep=""))
	mod = get(paste("MAX.mod1.",i, sep=""))
	preds = predict(mod, dat)
	cal = calib.mod(dat$PRESABS, preds)
	cal.mod = glm(dat$PRESABS ~ log((preds)/(1-preds)), family=binomial)
	plot(predict(cal.mod, type="response") ~ preds, xlim=c(0,1), ylim=c(0,1), xlab="predicted", ylab="observed")
	lines(y~x, lty="dashed")
	assign(paste("cal.",i, sep=""), cal)
	cal.MAX.training[i,1] = cal$calib.coeffs[1]
	cal.MAX.training[i,2] = cal$calib.coeffs[2]
	cal.MAX.training[i,3] = cal$'testa0|b1'
	cal.MAX.training[i,4] = cal$'testb1|a'
	}
dev.off()

setwd(path.obj)
save(cal.MAX.training, file="MAX.mod1.cal.training.Rda")

######## END RESUBSTITUTION RELIABILITY CALCULATIONS, MODEL=MAX
################################################################################





################################################################################
######## START CROSS-VALIDATION METRICS FROM MAXENT

library(PresenceAbsence) 

setwd(path.cod)
source("accuracy.R")

cv.accs=c()
modl="mod.MAX" # var placeholder
x.fold=5
n.col=9
for (i in 1:10) {
	## call up replicate training data and predictions
	dat = get(paste("dat",i, sep=""))
	mod = get(paste("MAX.mod1.",i, sep=""))
	temp = cv.accuracy.max(dat, mod, modl, x.fold, n.col)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "MAX.mod1"
	cv.accs=rbind(cv.accs, temp)
	}
setwd(path.obj)
save(cv.accs, file="MAX.mod1.cvaccs.pseudo11.Rda")

######## END CROSS-VALIDATION METRICS FROM MAXENT
################################################################################





################################################################################
######## START EXTERNAL ACCURACY METRICS FROM MAXENT
## predict to independent occupancy dataset

library(PresenceAbsence)   

setwd(path.cod)
source("accuracy.R")

ext.accs=c()
modl = "mod.MAX" # var placeholder
for (i in 1:10) {
	## pull in replicate model
	mod = get(paste("MAX.mod1.",i, sep=""))
	temp = ext.accuracy.max(mod, ext, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "MAX.mod1"
	ext.accs=rbind(ext.accs, temp)
	}
setwd(path.obj)
save(ext.accs, file="MAX.mod1.extaccs.pseudo11.Rda")

######## END EXTERNAL ACCURACY METRICS FROM MAXENT
################################################################################




################################################################################
######## START EXTERNAL RELIABILITY METRICS FROM MAXENT

setwd(path.cod)
source("calibration.R")

setwd(path.fig)
pdf(file="MAX_CalPlots_Testing_pseudo11.pdf", width=11, height=8.5)
par(mfrow=c(3,4))
x=seq(0,1,0.05)
y=seq(0,1,0.05)
cal.MAX.testing = as.data.frame(matrix(NA, 10,4))
names(cal.MAX.testing) = c("int", "slope", "p_int", "p_slope")
setwd(path.obj)
for (i in 1:10) {
	## pull in replicate data and model predictions	
	mod = get(paste("MAX.mod1.",i, sep=""))
	preds = predict(mod, ext)
	preds<-preds+0.00001 #jitter preds of 0 or 1
	preds[preds>=1]<-0.99999
	cal = calib.mod(ext$PRESABS, preds)
	cal.mod = glm(ext$PRESABS ~ log((preds)/(1-preds)), family=binomial)
	plot(predict(cal.mod, type="response") ~ preds, xlim=c(0,1), ylim=c(0,1), xlab="predicted", ylab="observed")
	lines(y~x, lty="dashed")
	assign(paste("cal.",i, sep=""), cal)
	cal.MAX.testing[i,1] = cal$calib.coeffs[1]
	cal.MAX.testing[i,2] = cal$calib.coeffs[2]
	cal.MAX.testing[i,3] = cal$'testa0|b1'
	cal.MAX.testing[i,4] = cal$'testb1|a'
	}
dev.off()

setwd(path.obj)
save(cal.MAX.testing, file="MAX.mod1.cal.testing.Rda")

######## END EXTERNAL RELIABILITY METRICS FROM MAXENT
################################################################################




################################################################################
######## START SOME EXTRA STUFF

## compare model variable importance
setwd(path.obj)
library(dismo)
par(ask=T)
par(mfrow=c(1,1))
for (i in 1:10) {
	mod = get(load(paste("MAX.mod1.",i,".Rda", sep="")))   
	plot(mod) 
	}

## compare variable responses
setwd(path.obj)
library(dismo)
par(ask=T)
par(mfrow=c(2,4))
for (i in 1:10) {
	mod = get(load(paste("MAX.mod1.",i,".Rda", sep="")))   
	response(mod)
	}
	
######## END SOME EXTRA STUFF
################################################################################
