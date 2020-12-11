###################################################################################
## Amy Angert
## Modified from Thomas C. Edwards, Jr., "Species Distribution Modelling Using R"
## PROGRAM FUNCTIONS: 
##   For each of 10 replicate training data sets (thinned herbarium presences + matching pseudoabsences)
##   	Build a random forest model
## 	 	Attempt to reduce overfitting by 
##	   		Using 1:1 ratio of pseudoabs:presence
##	   		Using different random sample of pseudos for replicate model runs
##   	Perform cross and external validation
## last update:  22 April 2016 (Amy triple checking everything for final tables and figs)
###################################################################################

## load libraries now if desired; loaded below when needed
#library(randomForest)
#library(PresenceAbsence)

## set pathnames
#path.root="/Users/amylauren/Google Drive/Occupancy AmNat Amy" 
#path.dat=paste(path.root, "/data files", sep="")
#path.dat.fix = paste(path.dat, "/Fixed datafiles - Matthew Sept 4", sep="")
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

## variables chosen based on dev and collinearity: bio15, lnbio10, lnbio14, lnbio12, bio11, bio4, lnbio3, bio2
## replicate thinned datasets with pseudo:pres at 4:1
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

## use these (1:1 ratio)
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
######## {Do not re-run} START RF MODEL

library(randomForest) 
#setwd(path.cod)
source("SDM/R_code/modforms.R")

#setwd(path.obj)
for (i in 1:10) {
	dat = get(paste("dat", i, sep=""))
	mod1.RF = randomForest(mod.form.lin(dat,1,2), importance=T, keep.forest=T, data=dat)           
	assign(paste("RF.mod1.",i, sep=""), mod1.RF)
	save(mod1.RF, file=paste("RF.mod1.",i,".pseudo11.Rda", sep=""))	
	mod1.pred = predict(mod1.RF, type="prob")[,2] # predict from model
	assign(paste("RF.mod1.",i,".pred", sep=""), mod1.pred)
	}

######## END RF MODEL
################################################################################




################################################################################
######## {Start here} LOAD SAVED MODEL AND PREDICTIONS

library(randomForest) 
#setwd(path.obj)

for (i in 1:10) {
	dat = get(paste("dat", i, sep=""))
	mod1.RF = get(load(paste("SRF.mod1.",i,".pseudo11.Rda", sep="")))
	assign(paste("RF.mod1.",i, sep=""), mod1.RF)
	mod1.pred = predict(mod1.RF, type="prob")[,2]
	assign(paste("RF.mod1.",i,".pred", sep=""), mod1.pred)
	}


######## END RF MODEL
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
