###################################################################################
## Amy Angert
## Modified from Thomas C. Edwards, Jr., "Species Distribution Modelling Using R"
## PROGRAM FUNCTIONS: 
##   For 10 replicate training datasets of thinned herbarium records and matching pseudoabsences:
##   Build basic BRT model
## 	 Attempt to reduce overfitting by using 1:1 ratio of pseudoabs:presence
##   Examine output and adjust LR / TC as needed
##   Estimate accuracies of final model
## last update:  22 April 2016 (Amy triple checking everything for final figures and tables)
###################################################################################

## load libraries now if desired; loaded below when needed
#library(gbm)
#library(dismo)
#library(PresenceAbsence)

## set pathnames
#path.root="/Users/amylauren/Google Drive/Occupancy AmNat Amy"
#path.dat = paste(path.root, "/data files", sep="") 
#path.obj = paste(path.root, "/R objects", sep="")
#path.fig = paste(path.root, "/figures", sep="")
#path.cod = paste(path.root, "/R code", sep="")

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
## replicate thinned datasets with pseudo:pres at 1:1
## some variables are already ln-transformed
## see file 'RCode_ThinPseudoAbs'
## now just call in saved .csv files

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

##use these (1:1 ratio)
for (i in 1:10) {
	dat = read.csv(paste("SDM/Output/dat",i,"c.csv", sep=""))
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
all = read.csv("SDM/data_files//all.records.aug.31.csv") #includes occupancy dataset, cleaned herbarium records, and 20K pseudoabs drawn to match envir space of true absences
ext = all[all$DATASET=="occ",] #pull out occupancy dataset
ext$bio3 = log(ext$bio3+0.5) #make needed ln-transforms of predictors
ext$bio10 = log(ext$bio10+0.5)
ext$bio12 = log(ext$bio12+0.5)
ext$bio14 = log(ext$bio14+0.5)

######## END INITIALIZATION
################################################################################





################################################################################
######## {Do not re-run) START BRT MODELS
library(gbm)    # load gbm package for BRT
library(dismo)  # for BRT calls per Elith et al (2008) JAnimalEcol 77:802-813

#setwd(path.obj)
for (i in 1:10) {
	## call up replicate training data
	dat = get(paste("dat", i, sep=""))
	resp=paste("as.factor(",colnames(dat[1]),")",sep="")  # assign response to column number
	n.col=ncol(dat)                                       # number of columns
	pred=2:n.col     				                       # assign predictors to column numbers
	
	## basic BRT model
	mod1.BRT=gbm.step(data=dat, gbm.x=pred, gbm.y=1, family="bernoulli", 
		tree.complexity=2, learning.rate=0.01, bag.fraction=0.75, n.folds=5, 
		n.trees=50, plot.main=TRUE, keep.fold.fit=TRUE, step.size=15) 
 	assign(paste("BRT.mod1.",i, sep=""), mod1.BRT)
 	save(mod1.BRT, file=paste("SDM/Output/BRT.mod1.",i,".pseudo11.Rda", sep=""))
 		
	## basic BRT model - LR adjusted up 
	mod2.BRT=gbm.step(data=dat, gbm.x=pred, gbm.y=1, family="bernoulli",
		tree.complexity=2, learning.rate=0.001, bag.fraction=0.75, n.folds=5,
	    plot.main=TRUE, keep.fold.fit=TRUE, step.size=15)
 	assign(paste("BRT.mod2.",i, sep=""), mod2.BRT)
 	save(mod2.BRT, file=paste("SDM/Output/BRT.mod2.",i,".pseudo11.Rda", sep=""))

	## basic BRT model - LR adjusted down
	mod3.BRT=gbm.step(data=dat, gbm.x=pred, gbm.y=1, family="bernoulli",
	    tree.complexity=2, learning.rate=0.1, bag.fraction=0.75, n.folds=5,
	    plot.main=TRUE, keep.fold.fit=TRUE, step.size=15)
 	assign(paste("BRT.mod3.",i, sep=""), mod3.BRT)
 	save(mod3.BRT, file=paste("SDM/Output/BRT.mod3.",i,".pseudo11.Rda", sep=""))

	## basic BRT model - now bump TC up
	mod4.BRT=gbm.step(data=dat, gbm.x=pred, gbm.y=1, family="bernoulli",
	    tree.complexity=3, learning.rate=0.1, bag.fraction=0.75, n.folds=5,
	    plot.main=TRUE, keep.fold.fit=TRUE, step.size=5)
 	assign(paste("BRT.mod4.",i, sep=""), mod4.BRT)
 	save(mod4.BRT, file=paste("SDM/Output/BRT.mod4.",i,".pseudo11.Rda", sep=""))

	## basic BRT model - keep TC up and bump LR back up
	mod5.BRT=gbm.step(data=dat, gbm.x=pred, gbm.y=1, family="bernoulli",
	    tree.complexity=3, learning.rate=0.01, bag.fraction=0.75, n.folds=5,
	    plot.main=TRUE, keep.fold.fit=TRUE, step.size=5)
 	assign(paste("BRT.mod5.",i, sep=""), mod5.BRT)
 	save(mod5.BRT, file=paste("SDM/Output/BRT.mod5.",i,".pseudo11.Rda", sep=""))
	
	## back to min LR, decrease TC further
	mod6.BRT=gbm.step(data=dat, gbm.x=pred, gbm.y=1, family="bernoulli",
	   tree.complexity=1, learning.rate=0.1, bag.fraction=0.75, n.folds=5,
	   plot.main=TRUE, keep.fold.fit=TRUE, step.size=5)
 	assign(paste("BRT.mod6.",i, sep=""), mod6.BRT)
 	save(mod6.BRT, file=paste("SDM/Output/BRT.mod6.",i,".pseudo11.Rda", sep=""))
	}

## Notes about interpreting deviance graphs
# want to reach plateau in deviance (red line)
# want deviance to be low (i.e., low residual variance)
# want fewest trees that get you to plateau (green line) (i.e., parsimony)
# in practice, play with LR and TC to "optimize" (really, tinker) with these outcomes
# might face trade-off between large model with low deviance and small model with higher deviance


# Any failed models? - YES 
		for (j in 3:5) {
			for (i in 1:10) {
		mod = get(paste("BRT.mod",j,".",i, sep=""))
		print(paste("BRT.mod",j,".",i, sep=""))
		print(mod$cv.statistics[[1]])
		}}

# re-run failed models (BRT.mod4 reps - )
	# remove step size argument for failed models
	failed <- c(5,6,7)
	for(i in 1:length(failed)){
	dat = get(paste("dat", failed[i], sep="")); resp=paste("as.factor(",colnames(dat[1]),")",sep=""); n.col=ncol(dat); pred=2:n.col
	mod4.BRT=gbm.step(data=dat, gbm.x=pred, gbm.y=1, family="bernoulli", tree.complexity=3, learning.rate=0.1, bag.fraction=0.75, n.folds=5, plot.main=TRUE, keep.fold.fit=TRUE)
	assign(paste("BRT.mod4.",failed[i], sep=""), mod4.BRT)
 	save(mod4.BRT, file=paste("SDM/Output/BRT.mod4.",failed[i],".pseudo11.Rda", sep=""))}
	
	
## examine BRT output
dev = matrix(ncol=6, nrow=11)
ntrees = matrix(ncol=6, nrow=11)
for (j in 3:5) {
	for (i in 1:10) {
		mod = get(paste("BRT.mod",j,".",i, sep=""))
				if((is.numeric(mod$cv.statistics[[1]]))==T){
				dev[i,j] = mod$cv.statistics[[1]]
				ntrees[i,j] = mod$n.trees
				} else {
				dev[i,j] = "NA" # failed model 
				}
}}

names(dev) = c("mod1", "mod2","mod3", "mod4", "mod5", "mod6")
for (j in 3:5) {
	dev[11,j] = mean(as.numeric(as.character(dev[1:10,j])), na.rm=TRUE)}
dev = as.data.frame(dev)
	dev

ntrees = as.data.frame(ntrees)
names(ntrees) = c("mod1", "mod2","mod3", "mod4", "mod5", "mod6")
	for (j in 3:5) {ntrees[11,j] = mean(as.numeric(as.character(ntrees[1:10,j])), na.rm=TRUE)}
	ntrees

### Selecting best models for pseudo10:1 
	## train rep 1
		# lowest deviance: mod5 < mod4 < mod6 #
		# fewest trees: mod4 < mod3 < mod6 #
		# best: mod4, mod6
	## train rep 2
		# lowest deviance: mod3 < mod5 < mod4
		# fewest trees: mod4 < mod3 < mod5
		# best: mod3, mod4
	## train rep 3
		# lowest deviance: mod6 < mod4 < mod3
		# fewest trees: mod4 < mod3 < mod5
		# best: mod4, mod3
	## train rep 4
		# lowest deviance: mod3 < mod1 < mod5
		# fewest trees: mod4 < mod3 < mod6
		# best: mod3, mod4
	## train rep 5
		# lowest deviance: mod4 < mod3 < mod6
		# fewest trees: mod4 < mod3 < mod6
		# best: mod4, mod3
	## train rep 6
		# lowest deviance: mod5 < mod4 < mod3
		# fewest trees: mod4 < mod3 < mod6
		# best: mod4, mod3
	## train rep 7
		# lowest deviance: mod3 < mod5 < mod4
		# fewest trees: mod4 < mod3 < mod6
		# best: mod3, mod4
	## train rep 8
		# lowest deviance: mod4 < mod5 < mod1
		# fewest trees: mod4 < mod3 < mod6
		# best: mod4, mod5
	## train rep 9
		# lowest deviance: mod3 < mod5 < mod4
		# fewest trees: mod4 < mod3 < mod6
		# best: mod3, mod4
	## train rep 10
		# lowest deviance: mod6 < mod4 < mod3
		# fewest trees: mod4 < mod3 < mod5
		# best: mod4, mod3
	## overall 
		# mod1 1st x 0 + 2nd x 0
		# mod2 1st x 0 + 2nd x 0
		# mod3 1st x 4 + 2nd x 4 --> second best choice
		# mod4 1st x 6 + 2nd x 4  --> best choice
		# mod5 1st x 0 + 2nd x 1  
		# mod6 1st x 0 + 2nd x 1

### Selecting best models for pseudo4:1 
	## train rep 1
		# lowest deviance: mod5 < mod4 < mod3
		# fewest trees: mod4 < mod3 < mod6
		# best: mod4, mod3
	## train rep 2
		# lowest deviance: mod5 < mod4 < mod1
		# fewest trees: mod4 < mod3 < mod6
		# best: mod4, mod5
	## train rep 3
		# lowest deviance: mod4 < mod5 < mod6
		# fewest trees: mod4 < mod3 < mod6
		# best: mod4, mod6
	## train rep 4
		# lowest deviance: mod4 < mod1 < mod3
		# fewest trees: mod4 < mod3 < mod6
		# best: mod4, mod3
	## train rep 5
		# lowest deviance: mod5 < mod6 < mod3
		# fewest trees: mod4 < mod3 < mod6
		# best: mod3/mod6
	## train rep 6
		# lowest deviance: mod5 < mod4 < mod1
		# fewest trees: mod4 < mod3 < mod6
		# best: mod4, mod5
	## train rep 7
		# lowest deviance: mod1 < mod3 < mod6
		# fewest trees: mod4 < mod3 < mod6
		# best: mod3, mod6
	## train rep 8
		# lowest deviance: mod5 < mod4 < mod3
		# fewest trees: mod4 < mod3 < mod6
		# best: mod4, mod3
	## train rep 9
		# lowest deviance: mod3 < mod5 < mod4
		# fewest trees: mod4 < mod3 < mod6
		# best: mod3, mod4
	## train rep 10
		# lowest deviance: mod5 < mod4 < mod6
		# fewest trees: mod4 < mod3 < mod5
		# best: mod4, mod5
	## overall 
		# mod1 1st x 0 + 2nd x 0
		# mod2 1st x 0 + 2nd x 0
		# mod3 1st x 1.5 + 2nd x 3.5 --> second best choice
		# mod4 1st x 7 + 2nd x 1  --> best choice
		# mod5 1st x 0 + 2nd x 3  
		# mod6 1st x 0.5 + 2nd x 2.5

### Selecting best models for pseudo1:1 
	## train rep 1
		# lowest deviance: mod4 < mod3 < mod5
		# fewest trees: mod5 < mod3 < mod4
		# best: mod5, mod5
	## train rep 2
		# lowest deviance: mod4 < mod5 < mod3
		# fewest trees: mod3 < mod4 < mod5
		# best: mod4, mod3
	## train rep 3
		# lowest deviance: mod1 < mod4 < mod3
		# fewest trees: mod4 < mod3 < mod6
		# best: mod4, mod3
	## train rep 4
		# lowest deviance: mod5 < mod1 < mod3 < mod4
		# fewest trees: mod4 < mod3 < mod6 < mod5
		# best: mod3, mod4/mod5
	## train rep 5
		# lowest deviance: mod4 < mod5 < mod1 < mod2
		# fewest trees: mod3 < mod4 < mod6 < mod5
		# best: mod4, mod5
	## train rep 6
		# lowest deviance: mod5 < mod6 < mod1
		# fewest trees: mod3 < mod6 < mod5
		# best: mod6, mod5
	## train rep 7
		# lowest deviance: mod4 < mod5 < mod1 < mod3
		# fewest trees: mod4 < mod3 < mod6 < mod5
		# best: mod4, mod5/mod3
	## train rep 8
		# lowest deviance: mod5 < mod3 < mod1
		# fewest trees: mod3 < mod6 < mod5
		# best: mod3, mod5
	## train rep 9
		# lowest deviance: mod4 < mod1 < mod5 < mod6
		# fewest trees: mod4 < mod3 < mod6 < mod5
		# best: mod4, mod5/mod6
	## train rep 10
		# lowest deviance: mod5 < mod4 < mod1 < mod3
		# fewest trees: mod4 < mod3 < mod6 < mod5
		# best: mod4, mod5
	## overall 
		# mod1 1st x 0 + 2nd x 0
		# mod2 1st x 0 + 2nd x 0
		# mod3 1st x 2 + 2nd x 2.5 --> decent choice (prioritizes low ntrees)
		# mod4 1st x 7 + 2nd x 0.5  --> best overall choice, but it crashes for some training sets
		# mod5 1st x 0 + 2nd x 7.5 --> best runner-up (prioritizes low dev)
		# mod6 1st x 1 + 2nd x 0.5
		
######## END BRT MODELS
################################################################################




################################################################################
######## {START HERE} LOAD FINAL MODELS 

#setwd(path.obj)

for (i in 1:10) {
	mod = get(load(paste("SDM/Output/BRT.mod4.",i,".pseudo11.Rda", sep="")))
	assign(paste("BRT.mod4.",i, sep=""), mod)
	}

######## LOAD FINAL REDUCED MODELS AND THEIR PREDICTIONS
################################################################################




################################################################################
######## START INTERNAL ACCURACY CALCULATIONS, MODEL=BRT

library(PresenceAbsence)   
#setwd(path.cod)
source("SDM/R_code/accuracy.R")

## be sure to adjust model #
accs = c()
for (i in 1:10) {
	dat = get(paste("dat", i, sep=""))
	mod = get(paste("BRT.mod4.",i, sep=""))
	#pred = mod$fold.fit
		pred = mod$fitted
		modl="BRT.mod4"              
		temp = accuracy.brt(dat, pred, modl)
		temp$rep = i
		temp$thresh = c("SensSpec", "Kappa")
		temp$model = "BRT.mod4"
		accs = rbind(accs, temp)	
}
#setwd(path.obj)
save(accs, file="SDM/Output/BRT.mod4.resubaccs.pseudo11.Rda")

######## END INTERNAL ACCURACY CALCULATIONS, MODEL=RF
################################################################################




################################################################################
######## START RESUBSTITUTION RELIABILITY CALCULATIONS, MODEL= BRT

#setwd(path.cod)
source("SDM/R_code/calibration.R")

#setwd(path.fig)
pdf(file="SDM/Output/BRT_CalPlots_Training.mod4.pseudo11.pdf", width=11, height=8.5)
par(mfrow=c(3,4))
x=seq(0,1,0.05)
y=seq(0,1,0.05)
cal.BRT.training = as.data.frame(matrix(NA, 10,4))
names(cal.BRT.training) = c("int", "slope", "p_int", "p_slope")
#setwd(path.obj)
for (i in 1:10) {
	## pull in replicate data and model predictions	
	dat = get(paste("dat",i, sep=""))
	mod = get(paste("BRT.mod4.",i, sep=""))
	preds = mod$fitted
	cal = calib.mod(dat$PRESABS, preds)
	cal.mod = glm(dat$PRESABS ~ log((preds)/(1-preds)), family=binomial)
	plot(predict(cal.mod, type="response") ~ preds, xlim=c(0,1), ylim=c(0,1), xlab="predicted", ylab="observed")
	lines(y~x, lty="dashed")
	assign(paste("cal.",i, sep=""), cal)
	cal.BRT.training[i,1] = cal$calib.coeffs[1]
	cal.BRT.training[i,2] = cal$calib.coeffs[2]
	cal.BRT.training[i,3] = cal$'testa0|b1'
	cal.BRT.training[i,4] = cal$'testb1|a'
	}
dev.off()

#setwd(path.obj)
save(cal.BRT.training, file="SDM/Output/BRT.mod4.cal.training.Rda")

######## END RESUBSTITUTION RELIABILITY CALCULATIONS, MODEL= BRT
################################################################################




################################################################################
######## START EXTERNAL ACCURACY CALCULATIONS, MODEL=BRT
## predict to independent occupancy dataset

library(PresenceAbsence)   
#setwd(path.cod)
source("SDM/R_code/accuracy.R")

## be sure to adjust model #
extaccs = c()
for (i in 1:10) {
	mod = get(paste("BRT.mod4.",i, sep=""))
	modl = "BRT.mod4"
	temp = ext.accuracy.brt(mod, ext, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "BRT.mod4"
	extaccs = rbind(extaccs, temp)	
	}
#setwd(path.obj)
save(extaccs, file="SDM/Output/BRT.mod4.extaccs.pseudo11.Rda")

######## END EXTERNAL ACCURACY CALCULATIONS, MODEL=BRT
################################################################################




################################################################################
######## START EXTERNAL RELIABILITY CALCULATIONS, MODEL=BRT

#setwd(path.cod)
source("SDM/R_code/calibration.R")

#setwd(path.fig)
pdf(file="SDM/Output/BRT_CalPlots_Testing.mod4.pseudo11.pdf", width=11, height=8.5)
par(mfrow=c(3,4))
x=seq(0,1,0.05)
y=seq(0,1,0.05)
cal.BRT.testing = as.data.frame(matrix(NA, 10,4))
names(cal.BRT.testing) = c("int", "slope", "p_int", "p_slope")
#setwd(path.obj)
for (i in 1:10) {
	## pull in replicate data and model predictions	
	mod = get(paste("BRT.mod4.",i, sep=""))
	preds = predict(mod, newdata=ext, type="response", n.trees=50)	
	cal = calib.mod(ext$PRESABS, preds)
	cal.mod = glm(ext$PRESABS ~ log((preds)/(1-preds)), family=binomial)
	plot(predict(cal.mod, type="response") ~ preds, xlim=c(0,1), ylim=c(0,1), xlab="predicted", ylab="observed")
	lines(y~x, lty="dashed")
	assign(paste("cal.",i, sep=""), cal)
	cal.BRT.testing[i,1] = cal$calib.coeffs[1]
	cal.BRT.testing[i,2] = cal$calib.coeffs[2]
	cal.BRT.testing[i,3] = cal$'testa0|b1'
	cal.BRT.testing[i,4] = cal$'testb1|a'
	}
dev.off()

#setwd(path.obj)
save(cal.BRT.testing, file="SDM/Output/BRT.mod4.cal.testing.Rda")

######## END EXTERNAL RELIBIALITY CALCULATIONS, MODEL=BRT
################################################################################




################################################################################
######## START SOME EXTRA STUFF

## compare relative variable importance
BRT.mod3.1$contributions

## examine response:predictor plots
par(ask=T)
for (i in 1:10) {
	mod = get(paste("BRT.mod3.",i, sep=""))
	gbm.plot(mod, n.plots=7) # response:predictor plots
	}
for (i in 1:10) {
	mod = get(paste("BRT.mod5.",i, sep=""))
	gbm.plot(mod, n.plots=7) # response:predictor plots
	}
	
## search for & examine interactions
mod.int=gbm.interactions(mod) # examine pairwise interactions
mod.int$rank.list # matrix of 5 top interactions

## plot top pairwise interactions
#par(mfrow=c(1,2))
gbm.perspec(mod, mod.int$rank.list[1,1], mod.int$rank.list[1,3], theta=30)

######## END SOME EXTRA STUFF
################################################################################

