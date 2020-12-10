
###################################################################################
## Amy Angert
## Modified from Thomas C. Edwards, Jr., "Species Distribution Modelling Using R"
## PROGRAM FUNCTIONS: 
##   For each of 10 replicate training data sets (thinned herbarium presences + matching pseudoabsences)
##		Build full and reduced GAM models
##   	Perform resubstitution, cross-validation, and external validation to assess accuracy
##		Calculate Miller's calibration statistics to assess reliability
## last update:  22 April 2016 (Amy triple checking everything to make final figs and tables)
###################################################################################

## load libraries now if desired; loaded below when needed
#library(PresenceAbsence)
#library(DAAG)
#library(gam)
#library(car)

## set pathnames
#path.root="/Users/amylauren/Google Drive/Occupancy AmNat Amy" 
#path.dat=paste(path.root, "/data files", sep="")
#path.obj=paste(path.root, "/R objects", sep="")
#path.fig=paste(path.root, "/figures", sep="")
#path.cod=paste(path.root, "/R code", sep="")

## set pathnames - Matthew
#path.root = "C:/Users/DW/Desktop/temp.sept.30" 
#path.dat = paste(path.root, "/data files", sep="")
#path.dat.fix = paste(path.root, "/data files", sep="") # older files relocated to another directory
#path.obj = paste(path.root, "/R objects", sep="")
#path.eco = paste(path.obj, "/ecoregions.shp", sep="")
#path.bio = paste(path.obj, "/wc0.5", sep="")
#path.cod=paste(path.root, "/R code", sep="")
#path.fig=paste(path.root, "/figures", sep="")


################################################################################
######## START INITIALIZATION
## import data
## slim down variables to conform to structure required by mod.form function	
## variables chosen based on dev and collinearity: bio15, lnbio10, lnbio14, lnbio12, bio11, bio4, lnbio3, bio2
## variables requiring ln-transform are already transformed (see RCode_ThinPseudos)

#setwd(path.dat) 

#for (i in 1:10) {
#	dat = read.csv(paste("dat",i,".csv", sep=""))
#	dat = dat[,c(4,59:61,67:69,71:72)] #column indices changed 9/4/14
#	dat = na.omit(dat) #for some reason many presences are now missing values for predictor variables
#	assign(paste("dat",i, sep=""), dat)
#	}
#rm(dat)

#for (i in 1:10) {
#	dat = read.csv(paste("dat",i,"b.csv", sep=""))
#	dat = dat[,c(4,59:61,67:69,71:72)] #column indices changed 9/4/14
#	dat = na.omit(dat) #for some reason many presences are now missing values for predictor variables
#	assign(paste("dat",i, sep=""), dat)
#	}
#rm(dat)

## use these (1:1 ratio)
for (i in 1:10) {
	dat = read.csv(paste("SDM/Output/dat",i,"c.csv", sep=""))
	dat = dat[,c(4,59:61,67:69,71:72)] #column indices changed 9/4/14
	assign(paste("dat",i, sep=""), dat)
	}
rm(dat)

## examine training data
dim(dat1); dim(dat10)
table(dat1$PRESABS); table(dat10$PRESABS)     
head(dat1); head(dat10)                     
str(dat1); str(dat10)                      

## testing data
all = read.csv("SDM/data_files/all.records.aug.31.csv") #includes occupancy dataset, cleaned herbarium records, and 20K pseudoabs drawn to match envir space of true absences
ext = all[all$DATASET=="occ",] #pull out occupancy dataset
ext$bio3 = log(ext$bio3+0.5) #make needed ln-transforms of predictors
ext$bio10 = log(ext$bio10+0.5)
ext$bio12 = log(ext$bio12+0.5)
ext$bio14 = log(ext$bio14+0.5)

######## END INITIALIZATION
################################################################################




################################################################################
######## {Do not rerun} START COMPARATIVE GAM MODELLING 

#setwd(path.cod)
source("SDM/R_Code/modforms.R")

#setwd(path.obj)
library(gam)

for (i in 1:10) {
	## call up replicate training data
	dat = get(paste("dat",i, sep=""))
	
	## GAM model 1: all smoothers using defaults df (smoothers=4)
	mod1.GAM=gam(mod.form.4(dat,1,2), family=binomial, data=dat)
	#mod1.fit=100*(1-mod1.GAM$deviance/mod1.GAM$null.deviance) 
	assign(paste("GAM.mod1.",i, sep=""), mod1.GAM)
	#assign(paste("GAM.mod1.",i,".fit", sep=""), mod1.fit)

	## GAM model 2: all smoothers using specified df (smoothers=3)
	mod2.GAM=gam(mod.form.3(dat,1,2), family=binomial, data=dat)
	#mod2.fit=100*(1-mod2.GAM$deviance/mod2.GAM$null.deviance) 
	assign(paste("GAM.mod2.",i, sep=""), mod2.GAM)
	#assign(paste("GAM.mod2.",i,".fit", sep=""), mod2.fit)

	## GAM model 3: all smoothers using specified df (smoothers=2)
	mod3.GAM=gam(mod.form.2(dat,1,2), family=binomial, data=dat)
	#mod3.fit=100*(1-mod3.GAM$deviance/mod3.GAM$null.deviance)  
	assign(paste("GAM.mod3.",i, sep=""), mod3.GAM)
	#assign(paste("GAM.mod3.",i,".fit", sep=""), mod3.fit)

	## GAM model 4: stepwise w/diff scopes
	mod4.GAM=step.gam(mod2.GAM,scope=list(
		"bio2"	=~1+ bio2	+s(bio2,2) +s(bio2,3)	+s(bio2,4),
		"bio3"	=~1+ bio3 	+s(bio3,2) 	+s(bio3,3) 	+s(bio3,4),
		"bio4"	=~1+ bio4 	+s(bio4,2) 	+s(bio4,3) 	+s(bio4,4),
		"bio10"	=~1+ bio10 	+s(bio10,2) +s(bio10,3) +s(bio10,4),
		"bio11"	=~1+ bio11 	+s(bio11,2) +s(bio11,3) +s(bio11,4),
		"bio12"	=~1+ bio12 	+s(bio12,2) +s(bio12,3)	+s(bio12,4),
		"bio14"	=~1+ bio14 	+s(bio14,2) +s(bio14,3)	+s(bio14,4),
		"bio15"	=~1+ bio15 	+s(bio15,2) +s(bio15,3)	+s(bio15,4)),
		trace=F)
	#mod4.fit=100*(1-mod4.GAM$deviance/mod4.GAM$null.deviance)  
	assign(paste("GAM.mod4.",i, sep=""), mod4.GAM)
	#assign(paste("GAM.mod4.",i,".fit", sep=""), mod4.fit)
	save(mod4.GAM, file=paste("GAM.mod4.",i,".pseudo11.Rda", sep=""))

	#mods.fit=cbind(mod1.fit,mod2.fit,mod3.fit,mod4.fit)  # all model fits
	#assign(paste("GAM.mods.",i,".fit", sep=""), mods.fit)
	}
	
######## END MODEL COMPARATIVE GAM MODELLING
################################################################################



################################################################################
######## {Start here} START LOAD SAVED MODEL OBJECTS AND PREDICTIONS

#setwd(path.obj)
library(gam)

for (i in 1:10) {
	mod = get(load(paste("GAM.mod4.",i,".pseudo11.Rda", sep="")))
	assign(paste("GAM.mod4.",i, sep=""), mod)
	pred = predict(mod, type="response")
	assign(paste("GAM.mod4.",i,".pred", sep=""), pred)
	}
	
######## END LOAD SAVED MODEL OBJECTS AND PREDICTIONS
################################################################################



################################################################################
######## START RESUBSTITUTION ACCURACY COMPARISONS, MODEL = STEP GAM

library(PresenceAbsence)   
#setwd(path.cod)
source("accuracy.R")

accs = c()
for (i in 1:10) {
	dat = get(paste("dat",i, sep=""))
	mod = get(paste("GAM.mod4.",i, sep=""))
	pred=predict(mod, type="response") # predict by model
	modl="mod4.GAM" # add var to keep track of model
	temp = accuracy(dat, pred, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "GAM.mod4"
	accs = rbind(accs, temp)
	}
#setwd(path.obj)
save(accs, file="GAM.mod4.accs.pseudo11.Rda")

######## END RESUBSTITUTION ACCURACY CALCULATIONS, MODEL= STEP GAM
################################################################################


	
################################################################################
######## START RESUBSTITUTION RELIABILITY CALCULATIONS, MODEL= STEP GAM

#setwd(path.cod)
source("calibration.R")

#setwd(path.fig)
pdf(file="GAM_CalPlots_Training.pseudo11.pdf", width=11, height=8.5)
par(mfrow=c(3,4))
x=seq(0,1,0.05)
y=seq(0,1,0.05)
cal.GAM.training = as.data.frame(matrix(NA, 10,4))
names(cal.GAM.training) = c("int", "slope", "p_int", "p_slope")
setwd(path.obj)
for (i in 1:10) {
	## pull in replicate data and model predictions	
	dat = get(paste("dat",i, sep=""))
	mod = get(paste("GAM.mod4.",i, sep=""))
	preds = predict(mod, dat, type="response")
	cal = calib.mod(dat$PRESABS, preds)
	cal.mod = glm(dat$PRESABS ~ log((preds)/(1-preds)), family=binomial)
	plot(predict(cal.mod, type="response") ~ preds, xlim=c(0,1), ylim=c(0,1), xlab="predicted", ylab="observed")
	lines(y~x, lty="dashed")
	assign(paste("cal.",i, sep=""), cal)
	cal.GAM.training[i,1] = cal$calib.coeffs[1]
	cal.GAM.training[i,2] = cal$calib.coeffs[2]
	cal.GAM.training[i,3] = cal$'testa0|b1'
	cal.GAM.training[i,4] = cal$'testb1|a'
	}
dev.off()

setwd(path.obj)
save(cal.GAM.training, file="GAM.mod4.cal.training.Rda")


######## END RESUBSTITUTION RELIABILITY CALCULATIONS, MODEL= STEP GAM
################################################################################




################################################################################
######## START CROSS-VALIDATION ACCURACY CALCULATIONS, MODEL=STEP GAM

## (will vary slightly if re-run due to random assignment to folds)

library(PresenceAbsence)	
library(DAAG)          
setwd(path.cod)
source("accuracy.R")

cv.accs = c()
for (i in 1:10) {
	## pull in replicate data and model
	mod = get(paste("GAM.mod4.",i, sep=""))
	dat = get(paste("dat",i, sep=""))
	modl="mod4.GAM" # assign model to varname
	temp = cv.accuracy(mod, dat, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "GAM.mod4"
	cv.accs = rbind(cv.accs, temp)
	}
setwd(path.obj)
save(cv.accs, file="GAM.mod4.cvaccs.pseudo11.Rda")

######## END CROSS-VALIDATION ACCURACY CALCULATIONS, MODEL= STEP GAM
################################################################################




################################################################################
######## START EXTERNAL VALIDATION CALCULATIONS, MODEL= STEP GAM
## predict to independent occupancy dataset

library(PresenceAbsence)
setwd(path.cod)
source("accuracy.R")

ext.accs = c()
for (i in 1:10) {
	## pull in replicate model
	mod = get(paste("GAM.mod4.",i, sep=""))
	modl="mod4.GAM"  # assign model to varname
	temp = ext.accuracy(mod, ext, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "GAM.mod4"
	ext.accs = rbind(ext.accs, temp)
  }
setwd(path.obj)
save(ext.accs, file="GAM.mod4.extaccs.pseudo11.Rda")

######## END EXTERNAL VALIDATION CALCULATIONS, MODEL=STEP GAM
################################################################################




################################################################################
######## START EXTERNAL RELIABILITY  CALCULATIONS, MODEL=STEP GAM
setwd(path.cod)
source("calibration.R")

setwd(path.fig)
pdf(file="GAM_CalPlots_Testing.pseudo11.pdf", width=11, height=8.5)
par(mfrow=c(3,4))
x=seq(0,1,0.05)
y=seq(0,1,0.05)
cal.GAM.testing = as.data.frame(matrix(NA, 10,4))
names(cal.GAM.testing) = c("int", "slope", "p_int", "p_slope")
setwd(path.obj)
for (i in 1:10) {
	## pull in replicate data and model predictions	
	mod = get(paste("GAM.mod4.",i, sep=""))
	preds = predict(mod, ext, type="response")
	cal = calib.mod(ext$PRESABS, preds)
	cal.mod = glm(ext$PRESABS ~ log((preds)/(1-preds)), family=binomial)
	plot(predict(cal.mod, type="response") ~ preds, xlim=c(0,1), ylim=c(0,1), xlab="predicted", ylab="observed")
	lines(y~x, lty="dashed")
	assign(paste("cal.",i, sep=""), cal)
	cal.GAM.testing[i,1] = cal$calib.coeffs[1]
	cal.GAM.testing[i,2] = cal$calib.coeffs[2]
	cal.GAM.testing[i,3] = cal$'testa0|b1'
	cal.GAM.testing[i,4] = cal$'testb1|a'
	}
dev.off()

setwd(path.obj)
save(cal.GAM.testing, file="GAM.mod4.cal.testing.Rda")

######## END EXTERNAL RELIABILITY CALCULATIONS, MODEL= STEP GAM
################################################################################


