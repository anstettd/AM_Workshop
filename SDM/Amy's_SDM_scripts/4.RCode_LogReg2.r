###################################################################################
## Amy Angert
## Modified from Thomas C. Edwards, Jr., "Species Distribution Modelling Using R"
## PROGRAM FUNCTIONS: 
##   For each of 10 replicate training data sets (thinned herbarium presences + matching pseudoabsences)
##		Build full and reduced model logistic GLMs
##   	Perform resubstitution, cross-validation, and external validation to assess accuracy
##		Calculate Miller's calibration statistics to assess reliability
## last update:  22 April 2016 (Amy going through all code to triple-check for errors and make final figures and tables)
###################################################################################

## load libraries now if desired; loaded below when needed
#library(PresenceAbsence)
#library(DAAG)

## set pathnames - Amy Angert
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
#	dat = na.omit(dat) 
#	assign(paste("dat",i, sep=""), dat)
#	}
#rm(dat)

#for (i in 1:10) {
#	dat = read.csv(paste("dat",i,"b.csv", sep=""))
#	dat = dat[,c(4,59:61,67:69,71:72)] #column indices changed 9/4/14
#	dat = na.omit(dat) 
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
######## {DO NOT RERUN} START BUILD FULL AND REDUCED MODELS 

#setwd(path.cod)
source("SDM/R_code/modforms.R")

#setwd(path.obj)

for (i in 1:10) {
	## build initial model all variables: mod1.i.LR
	dat = get(paste("dat",i, sep=""))
	mod1.LR=glm(mod.form.quad(dat,1,2),family=binomial,data=dat)
	assign(paste("LR.mod1.",i, sep=""), mod1.LR)
	#mod1.fit=100*(1-mod1.LR$deviance/mod1.LR$null.deviance)
	mod1.pred=predict(mod1.LR,type="response") # model prediction
	assign(paste("LR.mod1.",i,".pred", sep=""), mod1.pred)
	#summary(mod1.LR)        # full model summary stats
	save(mod1.LR, file=paste("SDM/Output/LR.mod1.",i,".pseudo11.Rda", sep="")) # save model object

	## build parsimonious model w/ backwards variable reduction: mod2.i.LR
	mod2.LR=step(mod1.LR,trace=F) # backwards stepwise variable reduction
	assign(paste("LR.mod2.",i, sep=""), mod2.LR)
	#mod2.fit=100*(1-mod2.LR$deviance/mod2.LR$null.deviance) # model fit
	mod2.pred=predict(mod2.LR,type="response") # model prediction
	assign(paste("LR.mod2.",i,".pred", sep=""), mod2.pred)
	#summary(mod2.LR)  # reduced model summary	
	save(mod2.LR, file=paste("SDM/Output/LR.mod2.",i,".pseudo11.Rda", sep="")) # save model object
	}

######## END BUILD FULL AND REDUCED MODELS
################################################################################




################################################################################
######## {START HERE} LOAD FINAL REDUCED MODELS AND THEIR PREDICTIONS 

#setwd(path.obj)

for (i in 1:10) {
	mod = get(load(paste("SDM/Output/LR.mod2.",i,".pseudo11.Rda", sep="")))
	assign(paste("LR.mod2.",i, sep=""), mod)
	pred=predict(mod,type="response") # model prediction
	assign(paste("LR.mod2.",i,".pred", sep=""), pred)
	}

######## LOAD FINAL REDUCED MODELS AND THEIR PREDICTIONS
################################################################################




################################################################################
######## START RESUBSTITUTION ACCURACY CALCULATIONS, MODEL=LOGISTIC GLM
## requires pkg PresenceAbsence
## build testing dataframes using mod2.i predictions

library(PresenceAbsence)   
#setwd(path.cod)
source("SDM/R_code/accuracy.R")

accs = c()
for (i in 1:10) {
	dat = get(paste("dat",i, sep=""))
	pred = get(paste("LR.mod2.",i,".pred", sep=""))	
	modl="mod2.LR"
	temp = accuracy(dat, pred, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "LR.mod2"
	accs = rbind(accs, temp)
	}
#setwd(path.obj)
save(accs, file="SDM/Output/LR.mod2.accs.pseudo11.Rda")

######## END RESUBSTITUTION ACCURACY CALCULATIONS, MODEL=LOGISTIC GLM
################################################################################




################################################################################
######## START RESUBSTITUTION RELIABILITY CALCULATIONS, MODEL=LOGISTIC GLM
#setwd(path.cod)
source("SDM/R_code/calibration.R")

#setwd(path.fig)
pdf(file="SDM/Output/GLM_CalPlots_Training_pseudo11.pdf", width=11, height=8.5)
par(mfrow=c(3,4))
x=seq(0,1,0.05)
y=seq(0,1,0.05)
cal.LR.training = as.data.frame(matrix(NA, 10,4))
names(cal.LR.training) = c("int", "slope", "p_int", "p_slope")
for (i in 1:10) {
	## pull in replicate data and model predictions	
	dat = get(paste("dat",i, sep=""))
	mod = get(paste("LR.mod2.",i, sep=""))
	preds = predict(mod, dat, type="response")
	cal = calib.mod(dat$PRESABS, preds)
	cal.mod = glm(dat$PRESABS ~ log((preds)/(1-preds)), family=binomial)
	plot(predict(cal.mod, type="response") ~ preds, xlim=c(0,1), ylim=c(0,1), xlab="predicted", ylab="observed")
	lines(y~x, lty="dashed")
	assign(paste("cal.",i, sep=""), cal)
	cal.LR.training[i,1] = cal$calib.coeffs[1]
	cal.LR.training[i,2] = cal$calib.coeffs[2]
	cal.LR.training[i,3] = cal$'testa0|b1'
	cal.LR.training[i,4] = cal$'testb1|a'
	}
dev.off()

#setwd(path.obj)
save(cal.LR.training, file="SDM/Output/LR.mod2.cal.training.Rda")

######## END RESUBSTITUTION RELIABILITY CALCULATIONS, MODEL=LOGISTIC GLM
################################################################################




################################################################################
######## START CROSS-VALIDATION ACCURACY CALCULATIONS, MODEL=LOGISTIC GLM
## perform  k-fold cross-validation; requires pkg DAAG
## results will vary slight if re-run due to random division into folds

library(PresenceAbsence)  
library(DAAG)               
#setwd(path.cod)
source("SDM/R_code/accuracy.R")

cv.accs = c()
for (i in 1:10) {
	## pull in replicate data and model
	mod = get(paste("LR.mod2.",i, sep=""))
	dat = get(paste("dat",i, sep=""))
	modl="mod2.LR" # assign model to varname
	temp = cv.accuracy(mod, dat, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "LR.mod2"
	cv.accs = rbind(cv.accs, temp)
	}
#setwd(path.obj)
save(cv.accs, file="SDM/Output/LR.mod2.cvaccs.pseudo11.Rda")

######## END CROSS-VALIDATION ACCURACY CALCULATIONS, MODEL=LOGISTIC GLM
################################################################################




################################################################################
######## START EXTERNAL ACCURACY CALCULATIONS, MODEL=LOGISTIC GLM
## predict to independent occupancy dataset

library(PresenceAbsence)
#setwd(path.cod)
source("SDM/R_code/accuracy.R")

ext.accs = c()
for (i in 1:10) {
	## pull in replicate model
	mod = get(paste("LR.mod2.",i, sep=""))
	modl="mod2.LR" # add var to keep track of model
	temp = ext.accuracy(mod, ext, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "LR.mod2"
	ext.accs = rbind(ext.accs, temp)
	}
#setwd(path.obj)
save(ext.accs, file="SDM/Output/LR.mod2.extaccs.pseudo11.Rda")

######## END EXTERNAL VALIDATION CALCULATIONS, MODEL=LOGISTIC GLM
################################################################################




################################################################################
######## START EXTERNAL RELIABILITY  CALCULATIONS, MODEL=LOGISTIC GLM
#setwd(path.cod)
source("SDM/R_code/calibration.R")

setwd(path.fig)
pdf(file="SDM/Output/GLM_CalPlots_Testing_pseudo11.pdf", width=11, height=8.5)
par(mfrow=c(3,4))
x=seq(0,1,0.05)
y=seq(0,1,0.05)
cal.LR.testing = as.data.frame(matrix(NA, 10,4))
names(cal.LR.testing) = c("int", "slope", "p_int", "p_slope")
#setwd(path.obj)
for (i in 1:10) {
	## pull in replicate data and model predictions	
	mod = get(paste("LR.mod2.",i, sep=""))
	preds = predict(mod, newdata=ext, type="response")
	cal = calib.mod(ext$PRESABS, preds)
	cal.mod = glm(ext$PRESABS ~ log((preds)/(1-preds)), family=binomial)
	plot(predict(cal.mod, type="response") ~ preds, xlim=c(0,1), ylim=c(0,1), xlab="predicted", ylab="observed")
	lines(y~x, lty="dashed")
	assign(paste("cal.",i, sep=""), cal)
	cal.LR.testing[i,1] = cal$calib.coeffs[1]
	cal.LR.testing[i,2] = cal$calib.coeffs[2]
	cal.LR.testing[i,3] = cal$'testa0|b1'
	cal.LR.testing[i,4] = cal$'testb1|a'
	}
dev.off()

#setwd(path.obj)
save(cal.LR.testing, file="SDM/Output/LR.mod2.cal.testing.Rda")

######## END EXTERNAL RELIABILITY CALCULATIONS, MODEL=LOGISTIC GLM
################################################################################




################################################################################
######## START SPATIAL VARIATION IN RESIDUALS
#setwd(path.dat)
for (i in 1:10) {
	dat = read.csv(paste("SDM/Output/dat",i,"c.csv", sep=""))
	dat = dat[,c(4,7,59:61,67:69,71:72)] #column indices changed 9/4/14
	assign(paste("dat",i, sep=""), dat)
	}
rm(dat)

for (i in 1:10) {
	dat = get(paste("dat",i, sep=""))
	pred = get(paste("LR.mod2.",i,".pred", sep=""))
	dat = cbind(dat, pred)
	resid = dat$PRESABS - dat$pred
	dat = cbind(dat, resid)
	resid.pres = lm(dat$resid[dat$PRESABS==1]~dat$Latitude[dat$PRESABS==1])
	resid.abs = lm(dat$resid[dat$PRESABS==0]~dat$Latitude[dat$PRESABS==0])
	par(mfrow=c(1,2))
	plot(dat$resid[dat$PRESABS==1] ~ dat$Latitude[dat$PRESABS==1])
	plot(dat$resid[dat$PRESABS==0] ~ dat$Latitude[dat$PRESABS==0])
	}
## ---> these show that residual variance decreases to the north
		# presences and absences in the south are in suitable and unsuitable habitat (small and large residuals)
		# presences and absences in the north are in uniformly unsuitable habitat (large residuals for presences, small residuals for absences)
		
######## END SPATIAL VARIATION IN RESIDUALS
################################################################################

