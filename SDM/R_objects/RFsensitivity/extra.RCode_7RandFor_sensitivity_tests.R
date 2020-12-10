###################################################################################
## RANDOM FOREST SENSITIVTY TESTS
## load libraries now if desired; loaded below when needed
library(randomForest)
library(PresenceAbsence)

## set pathnames
path.root="/Users/amyangert/Dropbox/SDM with Matt/Amy work August 2014" 
path.dat=paste(path.root, "/data files", sep="")
path.dat.fix = paste(path.dat, "/Fixed datafiles - Matthew Sept 4", sep="")
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
## replicate thinned datasets with pseudo:pres at 1:1
## some variables are already ln-transformed
## see file 'RCode_ThinPseudoAbs'
## now just call in saved .csv files

setwd(path.dat.fix)


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









#============================================================================================#
# RANDOM FOREST POST-HOC SENSITIVITY TESTS
# March 16, 2015 (Matthew)
#============================================================================================#

library(randomForest) 
setwd(path.cod)
source("modforms.R")
setwd(path.obj)
dir.create(paste(path.obj, "/RFsensitivity", sep=""))				# create a new folder to hold models from sensitivity tests
path.RFsens <- paste(path.obj, "/RFsensitivity", sep="")
setwd(path.RFsens)

# SENSITIVITY TESTS
# MODEL COMPARISONS:
	# (X3) Portion of data split for cross validation 
			# argument: 
			# 90% train & 10% test - 10 fold cv
			# 80% train & 20% test - 5 fold cv
			# 66% train & 33% test - 3 fold cv
	# (X4) Number of trees (X4)
			# testing alternatives 50, 500, 1000, 1500
	# (X2) mtry - Number of variables randomly sampled as candidates at each split
			# default or with RFtune
	# 3X4X3 = 36 alternative models
	# * 10 replicate datasets = 360 models 
	
# MODEL EVALUATION:
	# AUC score from cross validation 
	# external AUC score 
	# Calibration slope estimate 
		# & p-value (is different from 1?)
	# Calibration intercept estimate
		# & p-value (is different from 0?)
	
	
# Test alternative models with nested for loops	

trees <- c("50", "500","1000", "1500") # trees 
cv_folds <- c(10, 5, 3) 
my_mtry <- c(2,3,4) # default would be around 2.8 & differences marginal from exploratory work +/- 0.5

for(l in 1:length(my_mtry)){
	this_mtry <- my_mtry[l]
	for(m in 1:length(cv_folds)){
			this_cv_fold = cv_folds[m]
		for(j in 1:length(trees)){
				this_tree = trees[j]; this_tree <- as.numeric(this_tree)	
				for (i in 1:10) {	
					setwd(path.obj)
					dat = get(paste("dat", i, sep=""))						# get data pseudo replicate
					# run model
					mod1.RF = randomForest(mod.form.lin(dat,1,2), mtry=this_mtry, 
						importance=T, keep.forest=T, data=dat, ntree=this_tree, cv.fold=this_cv_fold)           
				# assign and save model objects			
					assign(paste("RF_cv",this_cv_fold,"_ntree",this_tree,"_mtry", this_mtry, "_rep", i, sep=""), mod1.RF)
					setwd(path.RFsens)
				# save models to file
					save(mod1.RF, file=paste("RF_cv",this_cv_fold,"_ntree",this_tree,"_mtry", this_mtry, "_rep", i, ".rda", sep=""))	
				# predict models to data
					mod1.pred = predict(mod1.RF, type="prob")[,2] 					# predict from model
					assign(paste("RF_cv",this_cv_fold,"_ntree",this_tree,"_mtry", this_mtry, "_rep", i,".pred", sep=""), mod1.pred)
					
					}
					}
					}
					}
######## END RF SENSITIVITY MODELS
###############################################################################
######## START RESUBSTITUTION ACCURACY CALCULATIONS, MODEL=RF
		library(PresenceAbsence)   # PresenceAbsence for accuracy metrics
		setwd(path.cod)
		source("accuracy.R")
		accs = c()
		setwd(path.RFsens)
		models <- dir()
			for(l in 1:length(my_mtry)){
				this_mtry <- my_mtry[l]
					for(m in 1:length(cv_folds)){
						this_cv_fold = cv_folds[m]
							for(j in 1:length(trees)){
								this_tree = trees[j]; this_tree <- as.numeric(this_tree)	
									for (i in 1:10) {
										dat = get(paste("dat",i, sep=""))
										pred = get(paste("RF_cv",this_cv_fold,"_ntree",this_tree,"_mtry", this_mtry, "_rep", i,".pred", sep=""))
										# add var to keep track of model
										modl = paste("RF_cv",this_cv_fold,"_ntree",this_tree,"_mtry", this_mtry, "_rep", i, sep="")                    
										temp = accuracy(dat, pred, modl)	
										temp$rep = i
										temp$thresh = c("SensSpec")
										temp$model = modl
											temp$mtry <- my_mtry[l]
											temp$cvfold <- cv_folds[m]
											temp$ntree <- trees[j]
										accs = rbind(accs, temp)
										}
									setwd(path.RFsens)
									#save(accs, file=paste("accs_cv",this_cv_fold,"_ntree",this_tree,"_mtry", this_mtry, "_rep", i,".rda",sep=""))
									}}}
######## END
#setwd(path.RFsens); write.csv(accs, file="accs.csv")
################################################################################
######## START EXTERNAL VALIDATION, MODEL=RF
## predict to independent occupancy dataset
library(PresenceAbsence)  
setwd(path.cod)
source("accuracy.R")
ext.accs = c()
			for(l in 1:length(my_mtry)){
				this_mtry <- my_mtry[l]
					for(m in 1:length(cv_folds)){
						this_cv_fold = cv_folds[m]
							for(j in 1:length(trees)){
								this_tree = trees[j]; this_tree <- as.numeric(this_tree)	
									for (i in 1:10) {
											mod = get(paste("RF_cv",this_cv_fold,"_ntree",this_tree,"_mtry", this_mtry, "_rep", i, sep=""))
											modl = paste("RF_cv",this_cv_fold,"_ntree",this_tree,"_mtry", this_mtry, "_rep", i, sep="")                    
											temp = ext.accuracy.rf(mod, ext, modl)
											temp$rep = i
											temp$thresh = c("SensSpec")
											temp$model = "RF.mod1"
											temp$mtry <- my_mtry[l]
											temp$cvfold <- cv_folds[m]
											temp$ntree <- trees[j]
											temp$Unique = paste("RF_cv",this_cv_fold,"_ntree",this_tree,"_mtry", this_mtry, sep="")                    
											ext.accs = rbind(ext.accs, temp)
											setwd(path.RFsens)
											#save(ext.accs, file=paste("ex_accs_cv",this_cv_fold,"_ntree",this_tree,"_mtry", this_mtry, "_rep", i,".rda",sep=""))
}}}}

full <- cbind(accs, ext.accs)
accuracy <- cbind(accs, ext.accs)
colnames(accuracy)[7] <- "AUC_cv"
colnames(accuracy)[20] <- "AUC_ext"

library(plyr)
AUC_cv <- ddply(accuracy, ~Unique, summarize, mean=mean(AUC_cv))
AUC_ext <- ddply(accuracy, ~Unique, summarize, mean=mean(AUC_ext))
accuracy <- cbind(AUC_cv, AUC_ext[,2])
colnames(accuracy) <- c("MODEL", "AUC_cv", "AUC_ext")

setwd(path.RFsens); write.csv(accuracy, file="aucs_models.csv")
######## END EXTERNAL VALIDATION CALCULATIONS, MODEL=RF
################################################################################
################################################################################
######## START EXTERNAL RELIABILITY  CALCULATIONS, MODEL=RF
setwd(path.cod)
source("calibration.R")
x=seq(0,1,0.05)
y=seq(0,1,0.05)
frammy <- c(); frammy <- data.frame(frammy)
setwd(path.RFsens)
			for(l in 1:length(my_mtry)){
				this_mtry <- my_mtry[l]
					for(m in 1:length(cv_folds)){
						this_cv_fold = cv_folds[m]
							for(j in 1:length(trees)){
								this_tree = trees[j]; this_tree <- as.numeric(this_tree)	
									for (i in 1:10) {
										mod = get(paste("RF_cv",this_cv_fold,"_ntree",this_tree,"_mtry", this_mtry, "_rep", i, sep=""))
										preds = predict(mod, newdata=ext, type="prob")[,2]
										modl = paste("RF_cv",this_cv_fold,"_ntree",this_tree,"_mtry", this_mtry, "_rep", i, sep="")                    
										cal = calib.mod(ext$PRESABS, preds)
										preds<-preds+0.00001 #first jitter preds of 0 or 1
										preds[preds>=1]<-0.99999
										cal.mod = glm(ext$PRESABS ~ log((preds)/(1-preds)), family=binomial)
										assign(paste("cal.",i, sep=""), cal)
										tempr <- c(cal$calib.coeffs[1], cal$calib.coeffs[2], cal$'testa0|b1', cal$'testb1|a', modl)
										tempr <- data.frame(tempr); tempr <- t(tempr)
										colnames(tempr) <- c("int", "slope", "p_int", "p_slope", "model")
										frammy <- rbind(frammy, tempr); frammy <- rbind(frammy, tempr) # twice to match above

	}}}}
	
final <- cbind(full, frammy)
setwd(path.RFsens); write.csv(final, file="RF_tune.csv")


### double check mtry

ext_vars <- ext[, c("bio2", "bio3", "bio4", "bio10", "bio11", "bio12", "bio14", "bio15")]
tuneRF(ext_vars, as.factor(ext$PRESABS), mtryStart=1.5, ntreeTry=500, stepFactor=0.5, improve=0.001, trace=TRUE, plot=TRUE, doBest=FALSE)
tuneRF(dat5[,2:9], as.factor(dat5$PRESABS), mtryStart=1.5, ntreeTry=500, stepFactor=0.5, improve=0.001, trace=TRUE, plot=TRUE, doBest=FALSE)



######## END POST HOC-RF SENSITIVITY TEST
################################################################################


