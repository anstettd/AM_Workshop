###################################################################################
### SCRIPT PURPOSE: SDM based on random forest models of presence/absence ~ bioclimatic predictors
# Modified from Angert et al. 2018, American Naturalist
# Author: Amy Angert
# last update:  14 Dec 2020

## OVERALL WORKFLOW:
#	Build RF model
#	Perform resubstitution and cross-validation to assess accuracy

###################################################################################
### LOAD LIBRARIES AND PREPARE INPUTS

## LIBRARIES
library(randomForest) # for RF models
library(PresenceAbsence) # for accuracy stats
library(DAAG) # for cross-validation sampling

## INPUTS
# read in file
dat <- read_csv('SDM/data_files/sdm_input.csv')

# slim dataframe to conform to structure required by mod.form function
# (see script modforms.R)
dat.input <- dat %>% 
  select(presabs, bio15, bio10, bio14, bio12, bio11, bio4, bio3, bio2)

################################################################################


################################################################################
### RANDOM FORESTS MODEL

source("SDM/R_code/modforms.R")

mod1.RF <- randomForest(mod.form.lin(dat.input, 1, 2), importance=T, keep.forest=T, data=dat.input)           
save(mod1.RF, file="RF.mod1.Rda")	

################################################################################


################################################################################
######## {If necessary} LOAD SAVED MODEL AND CALCULATE PREDICTIONS

mod1.RF = get(load("RF.mod1.Rda"))
mod1.pred = predict(mod1.RF, type="prob")[,2]

################################################################################


################################################################################
######## START RESUBSTITUTION ACCURACY CALCULATIONS, MODEL=RF

library(PresenceAbsence)   # PresenceAbsence for accuracy metrics
#setwd(path.cod)
source("SDM/R_code/accuracy.R")

accs = c()
for (i in 1:10) {
	dat = get(paste("dat",i, sep=""))
	pred = get(paste("RF.mod1.",i,".pred", sep=""))
	modl = "mod1.RF"                     # add var to keep track of model
	temp = accuracy(dat, pred, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "RF.mod1"
	accs = rbind(accs, temp)
	}
#setwd(path.obj)
save(accs, file="SDM/Output/RF.mod1.accs.pseudo11.Rda")

######## END RESUBSTITUTION ACCURACY CALCULATIONS, MODEL=RF
################################################################################



################################################################################
######## START RESUBSTITUTION RELIABILITY CALCULATIONS, MODEL= RF

#setwd(path.cod)
source("SDM/R_code/calibration.R")

#setwd(path.fig)
pdf(file="SDM/Output/RF_CalPlots_Training.pseudo11.pdf", width=11, height=8.5)
par(mfrow=c(3,4))
x=seq(0,1,0.05)
y=seq(0,1,0.05)
cal.RF.training = as.data.frame(matrix(NA, 10,4))
names(cal.RF.training) = c("int", "slope", "p_int", "p_slope")
#setwd(path.obj)
for (i in 1:10) {
	## pull in replicate data and model predictions	
	dat = get(paste("dat",i, sep=""))
	mod = get(paste("RF.mod1.",i, sep=""))
	preds = predict(mod, type="prob")[,2]
	cal = calib.mod(dat$PRESABS, preds)
	preds<-preds+0.00001 #first jitter preds of 0 or 1
	preds[preds>=1]<-0.99999
	cal.mod = glm(dat$PRESABS ~ log((preds)/(1-preds)), family=binomial)
	plot(predict(cal.mod, type="response") ~ preds, xlim=c(0,1), ylim=c(0,1), xlab="predicted", ylab="observed")
	lines(y~x, lty="dashed")
	assign(paste("cal.",i, sep=""), cal)
	cal.RF.training[i,1] = cal$calib.coeffs[1]
	cal.RF.training[i,2] = cal$calib.coeffs[2]
	cal.RF.training[i,3] = cal$'testa0|b1'
	cal.RF.training[i,4] = cal$'testb1|a'
	}
dev.off()

#setwd(path.obj)
save(cal.RF.training, file="SDM/Output/RF.mod1.cal.training.Rda")

######## END RESUBSTITUTION RELIABILITY CALCULATIONS, MODEL= RF
################################################################################




################################################################################
######## START EXTERNAL VALIDATION, MODEL=RF
## predict to independent occupancy dataset

library(PresenceAbsence)  
#setwd(path.cod)
source("SDM/R_code/accuracy.R")

ext.accs = c()
for (i in 1:10) {
	mod = get(paste("RF.mod1.",i, sep=""))
	modl="mod1.RF"                     # add var to keep track of model
	temp = ext.accuracy.rf(mod, ext, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "RF.mod1"
	ext.accs = rbind(ext.accs, temp)
	}
#setwd(path.obj)
save(ext.accs, file="SDM/Output/RF.mod1.extaccs.pseudo11.Rda")

######## END EXTERNAL VALIDATION CALCULATIONS, MODEL=RF
################################################################################




################################################################################
######## START EXTERNAL RELIABILITY  CALCULATIONS, MODEL=RF
#setwd(path.cod)
source("SDM/R_code/calibration.R")

#setwd(path.fig)
pdf(file="SDM/Output/RF_CalPlots_Testing.pseudo11.pdf", width=11, height=8.5)
par(mfrow=c(3,4))
x=seq(0,1,0.05)
y=seq(0,1,0.05)
cal.RF.testing = as.data.frame(matrix(NA, 10,4))
names(cal.RF.testing) = c("int", "slope", "p_int", "p_slope")
setwd(path.obj)
for (i in 1:10) {
	## pull in replicate data and model predictions	
	mod = get(paste("RF.mod1.",i, sep=""))
	preds = predict(mod, newdata=ext, type="prob")[,2]
	cal = calib.mod(ext$PRESABS, preds)
	preds<-preds+0.00001 #first jitter preds of 0 or 1
	preds[preds>=1]<-0.99999
	cal.mod = glm(ext$PRESABS ~ log((preds)/(1-preds)), family=binomial)
	plot(predict(cal.mod, type="response") ~ preds, xlim=c(0,1), ylim=c(0,1), xlab="predicted", ylab="observed")
	lines(y~x, lty="dashed")
	assign(paste("cal.",i, sep=""), cal)
	cal.RF.testing[i,1] = cal$calib.coeffs[1]
	cal.RF.testing[i,2] = cal$calib.coeffs[2]
	cal.RF.testing[i,3] = cal$'testa0|b1'
	cal.RF.testing[i,4] = cal$'testb1|a'
	}
dev.off()

#setwd(path.obj)
save(cal.RF.testing, file="SDM/Output/RF.mod1.cal.testing.Rda")

######## END EXTERNAL RELIABILITY, MODEL=RF
################################################################################




################################################################################
######## START EXTRA STUFF

## variable importance plots
#setwd(path.obj)
library(randomForest)
par(ask=T)
for (i in 1:10) {
	mod = get(load(paste("RF.mod1.",i,".Rda", sep="")))
	varImpPlot(mod,main="Variable Importance Plots")
	}

## partial plots
#setwd(path.obj)
library(randomForest)
par(mfrow=c(2,4))
for (i in 1:10) {
	mod = get(load(paste("RF.mod1.",i,".Rda", sep="")))
	dat = get(paste("dat",i, sep=""))
	partialPlot(mod,dat, bio2, which.class="1", main="Partial Plot: Bio 2",ylab="Logit(Presence)") 
	partialPlot(mod,dat, bio3, which.class="1", main="Partial Plot: Bio 3",ylab="Logit(Presence)") 
	partialPlot(mod,dat, bio4, which.class="1", main="Partial Plot: Bio 4",ylab="Logit(Presence)") 
	partialPlot(mod,dat, bio10, which.class="1", main="Partial Plot: Bio 10",ylab="Logit(Presence)") 
	partialPlot(mod,dat, bio11, which.class="1", main="Partial Plot: Bio 11",ylab="Logit(Presence)") 
	partialPlot(mod,dat, bio12, which.class="1", main="Partial Plot: Bio 12",ylab="Logit(Presence)") 
	partialPlot(mod,dat, bio14, which.class="1", main="Partial Plot: Bio 14",ylab="Logit(Presence)") 
	partialPlot(mod,dat, bio15, which.class="1", main="Partial Plot: Bio 15",ylab="Logit(Presence)") 
	}
	
######## END EXTRA STUFF
################################################################################
