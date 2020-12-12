###################################################################################
## Amy Angert
## Modified from Thomas C. Edwards, Jr., "Species Distribution Modelling Using R"
## PROGRAM FUNCTIONS: 
##   Generate descriptive statistics for selected predictor variables
##   Examine correlations among selected predictor varaibles
##   Estimate variable importance using a simple GLM
## last update:  1 May 2014
###################################################################################


## load libraries
#modeest    # for determining data mode

## set pathnames
#path.root="/Users/amyangert/Desktop/cardinalis SDM" 
#path.dat=paste(path.root,"/data files",sep="") #(normals calculated from date of collection)

## set pathnames - Matthew
#path.root = "C:/Users/DW/Desktop/temp.sept.30" 
#path.dat = paste(path.root, "/data files", sep="")
#path.obj = paste(path.root, "/R objects", sep="")
#path.gis = paste(path.dat, "/ecoregions.shp", sep="")


################################################################################
######## START INITIALIZATION
## make training data dataframe (old, before Matt merged files)
#setwd(path.dat)
#ab = read.csv("pseudo.absence.herb.30yr.csv") #load pseudo-absences 
#pr1 = read.csv("30yr.avg.bioclim.herb.csv") #load cleaned herbarium records + bioclim
#pr2 = read.csv("30yr.avg.monthly.herb.csv") #load cleaned herbarium records + bioclim
#dim(ab); dim(pr1); dim(pr2)
#pr = merge(pr1, pr2, by="ID") #merge montly and bioclim variables
#ab.cull = ab[,c(1,3:4,6:72)] # drop unnecessary columns
#pr.cull = pr[,c(1,3:23,29:76)]
#ab.cull$PRES = 0 #add column for pres/abs
#pr.cull$PRES = 1
#names(ab.cull); names(pr.cull) #rename to match
#names(pr.cull)[1] = "ID1"
#names(pr.cull)[2] = "Latitude"
#names(pr.cull)[3] = "Longitude"
#dat = rbind(ab.cull, pr.cull) #bind pres and abs data into one frame

## read in training data dataframe
all = read.csv("SDM/data_files/all.records.aug.31.csv") #includes occupancy dataset, cleaned herbarium records, and 20K pseudoabs drawn to match envir space of true absences
dat = all[all$DATASET=="herb",] #cull out occupancy dataset
npres = sum(dat$PRESABS); npres #subsample for 10:1 ratio of pseudoabs:pres
pres = dat[dat$PRESABS==1,]
pseud = dat[dat$PRESABS==0,]
pseud = pseud[sample(nrow(pseud), 10*npres, replace=FALSE),]
dat = rbind(pres,pseud)

######## END INITIALIZATION 	
################################################################################



################################################################################
######## START DISTRIBUTION OF PREDICTORS

par(ask=T)
par(mfrow=c(1,2))
# Daniel: This for loop requires hitting enter ~20 times manually the code below doesn't actually generate the output file
for (i in 45:dim(dat)[2]) {
	hist(dat[,i], main=names(dat)[i])
	hist(log(dat[,i]+0.5), main=paste("log",names(dat)[i]))
	}
#substantially more normal when log-transformed: 
	#bio3 (subtle), bio10 (subtle) bio12, bio13 (questionable), bio14, bio16, bio17, bio18, bio19 (questionable)
	# PPT01 & PPT02 (both questionable), PPT05-PPT12
	# (all monthly temperature variables fine)
dat$lnbio3 = log(dat$bio3+0.5)
dat$lnbio10 = log(dat$bio10+0.5)
dat$lnbio12 = log(dat$bio12+0.5)
dat$lnbio13 = log(dat$bio13+0.5)
dat$lnbio14 = log(dat$bio14+0.5)
dat$lnbio16 = log(dat$bio16+0.5)
dat$lnbio17 = log(dat$bio17+0.5)
dat$lnbio18 = log(dat$bio18+0.5)
dat$lnbio19 = log(dat$bio19+0.5)
dat$lnPPT05 = log(dat$PPT05+0.5)
dat$lnPPT06 = log(dat$PPT06+0.5)
dat$lnPPT07 = log(dat$PPT07+0.5)
dat$lnPPT08 = log(dat$PPT08+0.5)
dat$lnPPT09 = log(dat$PPT08+0.5)
dat$lnPPT10 = log(dat$PPT10+0.5)
dat$lnPPT11 = log(dat$PPT11+0.5)
dat$lnPPT12 = log(dat$PPT12+0.5)

#replace raw with ln-transforms (in preliminary investigations, ln-var almost always have higher adj dev and lower AIC than untransformed var)
dat.ln = dat[,c(1:48,85:92,57:58,76,60:65,77,67,78:80,71,81:84)]
write.csv(dat.ln, "SDM/Output/herb.dat.lnpreds.csv")
######## END EXAMINE DISTRIBUTION OF PRED
################################################################################


################################################################################
######## START DESCRIPTIVE STATS AMONG PREDICTORS
library(modeest) # package modeest for mode
t1=aggregate(dat[,c(9:92)],list(dat$PRES),FUN=mean,na.rm=T); t1 # mean of bioclim by P/A
t2=aggregate(dat[,c(9:92)],list(dat$PRES),FUN=median,na.rm=T); t2 # median of bioclim by P/A
t3=aggregate(dat[,c(9:92)],list(dat$PRES),FUN=mfv,na.rm=T); t3 # mode of bioclim by P/A
# t3 trips ~100 warnings. Not sure if this is a problem.
#argument 'na.rm' is soft-deprecated, please start using 'na_rm' instead

######## END DESCRIPTIVE STATS AMONG PRED PREDICTORS
################################################################################



################################################################################
######## START GRAPHICAL DESCRIPTIVE STATS AMONG PREDICTORS
## basic boxplots by spp presence:absence
par(mfrow=c(1,3))
boxplot(dat$bio1~dat$PRES,xlab="Presence:Absence",ylab="MAT")
boxplot(dat$bio4~dat$PRES,xlab="Presence:Absence",ylab="Temperature Seasonality")
boxplot(dat$bio5~dat$PRES,xlab="Presence:Absence",ylab="Max Temp Warmest Month")
boxplot(dat$bio6~dat$PRES,xlab="Presence:Absence",ylab="Min Temp Coldest Month")
boxplot(dat$bio8~dat$PRES,xlab="Presence:Absence",ylab="Mean Temp Wettest Q")
boxplot(dat$bio9~dat$PRES,xlab="Presence:Absence",ylab="Mean Temp Driest Q")
boxplot(dat$bio12~dat$PRES,xlab="Presence:Absence",ylab="Ann Precip")
boxplot(dat$lnbio18~dat$PRES,xlab="Presence:Absence",ylab="log Precip of Warmest Q")
boxplot(dat$bio19~dat$PRES,xlab="Presence:Absence",ylab="Precip of Coldest Q")
boxplot(dat$Tave03~dat$PRES,xlab="Presence:Absence",ylab="Mar Ave T")
boxplot(dat$Tave03~dat$PRES,xlab="Presence:Absence",ylab="Mar Ave T")
boxplot(dat$Tave03~dat$PRES,xlab="Presence:Absence",ylab="Mar Ave T")
boxplot(dat$Tave06~dat$PRES,xlab="Presence:Absence",ylab="Jun Ave T")
boxplot(dat$Tave09~dat$PRES,xlab="Presence:Absence",ylab="Sep Ave T")
boxplot(dat$Tave12~dat$PRES,xlab="Presence:Absence",ylab="Dec Ave T")
boxplot(dat$PPT03~dat$PRES,xlab="Presence:Absence",ylab="Mar Ppt")
boxplot(dat$lnPPT06~dat$PRES,xlab="Presence:Absence",ylab="log Jun Ppt")
boxplot(dat$lnPPT09~dat$PRES,xlab="Presence:Absence",ylab="log Sep Ppt")
boxplot(dat$lnPPT12~dat$PRES,xlab="Presence:Absence",ylab="log Dec Ppt")

######## END GRAPHICAL DESCRIPTIVE STATS AMONG PRED PREDICTORS
################################################################################



################################################################################
######## START CORRELATIONS AMONG PREDICTORS
## numeric correlations
cut.point=.7  # set cutpoint for correlation
c1=cor(dat.ln[,c(8:20)],use="pairwise.complete.obs",method="pearson");c1 # est. correlation among monthly tmax
c1=cor(dat.ln[,c(21:32)],use="pairwise.complete.obs",method="pearson");c1 # est. correlation among monthly tmin
c1=cor(dat.ln[,c(33:44)],use="pairwise.complete.obs",method="pearson");c1 # est. correlation among monthly tave
c1=cor(dat.ln[,c(45:56)],use="pairwise.complete.obs",method="pearson");c1 # est. correlation among monthly ppt
c1=cor(dat.ln[,c(57:75)],use="pairwise.complete.obs",method="pearson");c1 # est. correlation among bioclim
c1=cor(dat.ln[,c(9:56)],use="pairwise.complete.obs",method="pearson");c1 # est. correlation among all monthlies
c2=subset(c1> cut.point | c1< -cut.point);c2 # matrix of cor>cutpoint (FALSE indicates cor<cutpoint)


## START MODIFIED panel.cor CORRELATION FUNCTION		   
##   determine correlations among predictor variables: modified from 
##   http://addictedtor.free.fr/graphiques/graphcode.php?graph=137
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
  { usr <- par("usr"); on.exit(par(usr)) 
    par(usr = c(0, 1, 0, 1)) 
    r <- abs(cor(x, y, use = "pairwise.complete.obs")) 
    txt <- format(c(r, 0.123456789), digits=digits)[1] 
    txt <- paste(prefix, txt, sep="") 
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 

    test <- cor.test(x,y) 
    # borrowed from printCoefmat
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " ")) 

    text(0.5, 0.5, txt, cex = cex * r)
    text(.8, .8, Signif, cex=cex, col=2) 
  }
## END MODIFIED panel.cor CORRELATION FUNCTION	

## plot correlations using modified panel.cor function
pairs(dat.ln[,c(57:75)],lower.panel=panel.smooth, upper.panel=panel.cor, 
      main="Bioclim Variables")
pairs(dat.ln[,c(33:44)],lower.panel=panel.smooth, upper.panel=panel.cor, 
      main="Monthly Ave Temperature Variables")
pairs(dat.ln[,c(9:20)],lower.panel=panel.smooth, upper.panel=panel.cor,
      main="Monthly Max Temperature Variables")
#Daniel: Error: unexpected ')' in "      main="Monthly Ave Temperature Variables")"
pairs(dat.ln[,c(21:32)],lower.panel=panel.smooth, upper.panel=panel.cor, 
      main="Monthly Min Temperature Variables")
pairs(dat.ln[,c(45:56)],lower.panel=panel.smooth, upper.panel=panel.cor, 
      main="Monthly Precip Variables")
#Daniel: Error: unexpected ')' in "      main="Monthly Min Temperature Variables")"

######## END CORRELATIONS AMONG PRED PREDICTORS
################################################################################



################################################################################
######## START VARIABLE IMPORTANCE AMONG PREDICTORS
## START function variable importance
## Cheerfully borrowed and modified from:  Niklaus E. Zimmermann
## This function evaluates the predictive power of individual
##   predictor variables as available in sequence within a data frame

#### START function variable importance
varimp.glm=function(tr.spp,tr.var,pres,pf,pl) {
  tmp.mat=matrix(ncol=2,nrow=(pl-pf+1))
  for (i in pf:pl) {
    ## option: linear+quadratic; linear only
    tmp=glm(tr.spp[,pres]~tr.var[,i]+I((tr.var[,i])^2),na.action=na.omit,family=binomial)
    #tmp=glm(tr.spp[,pres]~tr.var[,i],na.action=na.omit,family=binomial) # linear only glm
    tmp.mat[(i-pf+1),1]=tmp$aic
    tmp.mat[(i-pf+1),2]=(1-(tmp$deviance/tmp$null.deviance))
    }
    return(tmp.mat)
}
#Daniel: Warning message: In doTryCatch(return(expr), name, parentenv, handler) :reached elapsed time limit


#### END function variable importance

## estimate VIP values => AIC & Adj deviance
tr.vip=dat[,c(3,9:92)] # keep P/A & ALL predictors 
tr.vip=dat.ln[,c(3,9:75)] # keep P/A & all predictors 
tr.vip=dat.ln[,c(3,57:75)] # keep only P/A & bioclim predictors
tr.vip=dat.ln[,c(3,9:56)] #keep only P/A and monthly predictors
tr.vip=dat.ln[,c(3,33:44)] # keep only P/A & ave temp predictors
tr.vip=dat.ln[,c(3,9:20)] # keep only P/A & max temp predictors
tr.vip=dat.ln[,c(3,21:32)] # keep only P/A & min temp predictors
tr.vip=dat.ln[,c(3,45:56)] # keep only P/A & precip predictors
pres=1                     # column for presence:absence
v.start=2                  # column start predictor variables
v.stop=ncol(tr.vip)        # last column predictor variables
v.num=v.stop-1             # number predictor variables
dev.fit=varimp.glm(tr.vip,tr.vip,pres,v.start,v.stop) # call VIP function
x.labs=as.data.frame(names(tr.vip[2:v.stop]))
dev.fit =cbind(as.data.frame(dev.fit), x.labs)                 # output matrix; col=1 AIC, col=2 Adj deviance
names(dev.fit) = c("AIC", "AdjDev", "Var")
dev.fit = dev.fit[order(dev.fit$AIC),]
dev.fit

## built basic barplot if desired
d.max=ceiling(signif(max(dev.fit[,2]),2)*10)/10 # convoluted; max of y-axis
ylim.r=range(0,d.max)                           # range y-axis
x.labs=names(tr.vip[2:v.stop])                  # x-axis labels
barplot(dev.fit[,2],col="darkgreen",ylim=ylim.r,
  main="pred VIPs",ylab="adj.D2",names=x.labs)
abline(h=0); abline(mean(dev.fit[,2]),0,lt=3) # ref lines; dash=mean adj.dev

######## END VARIABLE IMPORTANCE AMONG PREDICTORS
################################################################################


################################################################################
######## START VARIABLE SELECTION BASED ON IMPORTANCE AND COLLINEARITY

### Models with bioclim only
# Choose best: bio15
# bio15 has no |r|>0.7 with other bioclim
# Choose 2nd best: lnbio10
# If bio10 is included, then bio1, bio5, bio8, and bio9 should be excluded (|r|>0.7)
# Choose 3rd best: lnbio14
# If lnbio14 is included, then bio17 and bio18 should be excluded
# Skip bio1, lnbio17, bio5
# Choose 7th best: lnbio12
# If lnbio12 is included, then lnbio13, bio16, and bio19 should be excluded
# Choose next best: bio11
# If bio11 is included, then bio1, bio6, and bio8 should be excluded
# Skip bio9, bio8, lnbio18, lnbio17, lnbio13, lnbio19
# Choose next best: bio4
# If bio4 is included, then bio7 should be excluded
# Choose next best of remaining: lnbio3
# lnbio3 has no |r|>0.7 with other bioclim
# Choose last remaining: bio2
# Final include list:
	#bio15, bio10, bio14, bio12, bio11, bio4, bio3, bio2

### Models with monthlies only (guided only by AIC/dev rankings)
# Choose best: Tave06
# Tave06 knocks out Tmax1-12, Tmin3-10, Tave2-12, PPT05
# Choose next best: lnPPT06
# lnPPT06 knocks out lnPPT4,10-12
# Choose next best: Tmin10
# Tmin10 knocks out Tave01, Tmin1-2,11-12
# Choose next best: lnPPT07
# lnPPT07 knocks out lnPPT08-09
# Choose next best: PPT01
# PPT01 knocks out PPT2-3
# Final include list
	# Tave06, lnPPT06, Tmin10, lnPPT07, PPT01

### Models with monthlies only (guided by categories then rankings)
# Goal is to have at least 1 Tmax, Tmin, and PPT
# Rankings of maxes: 11, 6, 12, 10, 1, 5, 7, 4, 9, 8, 2, 3
# Rankings of mins: 8, 7, 9, 6, 10, 5, 4, 11, 12, 1, 3, 2
# Rankings of ppt: 6, 5, 11, 10, 12, 7, 8, 9, 1, 4, 3, 2
# Choose best max: Tmax11
# Tmax11 not compatible with any Tmins (all |r| > 0.7)
# Discard Tmax11 and choose Tmax6
# Tmax6 compatible with Tmin8 (barely) or Tmin9 (definitely)
# Both temps compatible with PPT6
# Possible to include any other Tmax? No. Tmax01 is ok with Tmax06, but not Tmin08. No others compatible with Tmax06.
# Possible to include any other Tmin? Yes. Tmin01 is ok with Tmin08, Tmax06, and lnPPT06.
# Possible to include any other Ppt? Yes. Many Ppt are compatible with PPT06. Next best are all summer months also. Move to best winter/fall month, PPT01.
# Final include list:
	# Tmax06, Tmin08, lnPPT06, Tmin01, PPT01

######## END VARIABLE SELECTION BASED ON IMPORTANCE AND COLLINEARITY
################################################################################


################################################################################
####### Save files with different predictor sets
dat.bio = dat.ln[,c(1:8,58:60,66:68,70:71)]
dat.mo1 = dat.ln[,c(1:8,38,30,45,50,51)]
dat.mo2 = dat.ln[,c(1:8,14,21,28,45,50)]
save(dat.bio, dat.mo1, dat.mo2, file="SDM/Output/herb.dat.preds.Rda")








