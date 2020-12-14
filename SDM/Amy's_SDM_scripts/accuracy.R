accuracy = function(dat, pred, modl) {
	## make dataframe with 
		## col1 = model label (required by optimal.thresholds)
		## col2 = presence/absence
		## col3 = model prediction	
	datx=cbind(modl,dat[1],pred) # build dataframe w/model predictions

	## determine best threshold
	mod.cut=optimal.thresholds(datx,opt.methods=c("MaxSens+Spec","MaxKappa")) 

	## generate confusion matrix
	mod.cfmat.sensspec=table(datx[[2]],factor(as.numeric(datx$pred>=mod.cut$pred[1])))
	mod.cfmat.kappa=table(datx[[2]],factor(as.numeric(datx$pred>=mod.cut$pred[2])))

	## calculate model accuracies with standard deviation=F
	mod.acc=presence.absence.accuracy(datx,threshold=mod.cut$pred,st.dev=F)
	tss=mod.acc$sensitivity+mod.acc$specificity-1 # code TSS metric
	mod.acc=cbind(mod.acc[1:7],tss) # bind all metrics
}

cv.accuracy = function(mod, dat, modl) {
	## make cross-validated predictions
	cv.mod=CVbinary(mod, nfolds=5, print.details=F) # crossval predict
	cv.pred=cv.mod$cvhat   # assign new name to estimates
	
	## bind into dataframe
	cv.dat=cbind(modl,dat[1],cv.pred) # build obs and prediction dataframe

	## determine best threshold
	cv.cut=optimal.thresholds(cv.dat,opt.methods=c("MaxSens+Spec", "MaxKappa")) 
	
	## generate confusion matrices
	cv.cfmat.sensspec=table(cv.dat[[2]],factor(as.numeric(cv.dat$cv.pred>cv.cut$cv.pred[1]))) 
	cv.cfmat.kappa=table(cv.dat[[2]],factor(as.numeric(cv.dat$cv.pred>cv.cut$cv.pred[2]))) 

	## calculate accuracies
	cv.acc=presence.absence.accuracy(cv.dat,threshold=cv.cut$cv.pred,st.dev=F) 
	cv.tss=cv.acc$sensitivity+cv.acc$specificity-1 # code TSS metric
	cv.acc=cbind(cv.acc[1:7],cv.tss) # bind all metrics
}

ext.accuracy = function(mod, ext, modl) {	
	## predict to new data
	mod.epred=predict(mod, newdata=ext, type="response")
	
	## build testing dataframes using model predictions
	datx=cbind(modl, ext[3], mod.epred) # build dataframe w/mod predictions

	## determine best threshold
	mod.ecut=optimal.thresholds(datx,opt.methods=c("MaxSens+Spec","MaxKappa")) 

	## generate confusion matrices
	ext.cfmat.sensspec=table(datx[[2]],factor(as.numeric(datx$mod.epred>=mod.ecut$mod.epred[1])))
	ext.cfmat.kappa=table(datx[[2]],factor(as.numeric(datx$mod.epred>=mod.ecut$mod.epred[2])))
		
	## calculate model accuracies with standard deviation=F
	ext.acc=presence.absence.accuracy(datx,threshold=mod.ecut$mod.epred,st.dev=F)
	tss=ext.acc$sensitivity+ext.acc$specificity-1 # code TSS metric
	ext.acc=cbind(ext.acc[1:7],tss) # bind all metrics
}

accuracy.brt = function(dat, pred, modl) {
	## make dataframe with 
		## col1 = model label (required by optimal.thresholds)
		## col2 = presence/absence
		## col3 = model prediction	
	datx=cbind(modl,dat[1],pred) # build dataframe w/model predictions
	names(datx)[3] = "cvpred"
	datx$cvpred = exp(datx$cvpred)/(1+exp(datx$cvpred)) #convert from logit
	
	## determine best threshold
	mod.cut=optimal.thresholds(datx,opt.methods=c("MaxSens+Spec","MaxKappa")) 

	## generate confusion matrix
	mod.cfmat.sensspec=table(datx[[2]],factor(as.numeric(datx$cvpred>=mod.cut$cvpred[1])))
	mod.cfmat.kappa=table(datx[[2]],factor(as.numeric(datx$cvpred>=mod.cut$cvpred[2])))

	## calculate model accuracies with standard deviation=F
	mod.acc=presence.absence.accuracy(datx,threshold=mod.cut$cvpred,st.dev=F)
	tss=mod.acc$sensitivity+mod.acc$specificity-1 # code TSS metric
	mod.acc=cbind(mod.acc[1:7],tss) # bind all metrics
}

ext.accuracy.rf = function(mod, ext, modl) {	
	## predict to new data
	mod.epred=predict(mod, newdata=ext, type="prob")[,2]
	
	## build testing dataframes using model predictions
	datx=cbind(modl, ext[3], mod.epred) # build dataframe w/mod predictions

	## determine best threshold
	mod.ecut=optimal.thresholds(datx,opt.methods=c("MaxSens+Spec","MaxKappa")) 

	## generate confusion matrices
	ext.cfmat.sensspec=table(datx[[2]],factor(as.numeric(datx$mod.epred>=mod.ecut$mod.epred[1])))
	ext.cfmat.kappa=table(datx[[2]],factor(as.numeric(datx$mod.epred>=mod.ecut$mod.epred[2])))
		
	## calculate model accuracies with standard deviation=F
	ext.acc=presence.absence.accuracy(datx,threshold=mod.ecut$mod.epred,st.dev=F)
	tss=ext.acc$sensitivity+ext.acc$specificity-1 # code TSS metric
	ext.acc=cbind(ext.acc[1:7],tss) # bind all metrics
}

ext.accuracy.brt = function(mod, ext, modl) {	
	## predict to new data
	mod.epred=predict(mod, newdata=ext, type="response", n.trees=50)
	
	## build testing dataframes using model predictions
	datx=cbind(modl, ext[3], mod.epred) # build dataframe w/mod predictions

	## determine best threshold
	mod.ecut=optimal.thresholds(datx,opt.methods=c("MaxSens+Spec","MaxKappa")) 

	## generate confusion matrices
	ext.cfmat.sensspec=table(datx[[2]],factor(as.numeric(datx$mod.epred>=mod.ecut$mod.epred[1])))
	ext.cfmat.kappa=table(datx[[2]],factor(as.numeric(datx$mod.epred>=mod.ecut$mod.epred[2])))
		
	## calculate model accuracies with standard deviation=F
	ext.acc=presence.absence.accuracy(datx,threshold=mod.ecut$mod.epred,st.dev=F)
	tss=ext.acc$sensitivity+ext.acc$specificity-1 # code TSS metric
	ext.acc=cbind(ext.acc[1:7],tss) # bind all metrics
}

accuracy.max = function(dat, mod, modl) {
	## get model predictions
	mod.pred = predict(mod,dat)
	
	## make dataframe with 
		## col1 = model label (required by optimal.thresholds)
		## col2 = presence/absence
		## col3 = model prediction	
	datx=cbind(modl,dat[1],mod.pred) # build dataframe w/model predictions

	## evaluate model
	mod.val=evaluate(p=dat[dat$PRESABS==1,2:9], a=dat[dat$PRESABS==0,2:9], mod.MAX) 

	## determine thresholds
	mod.cut=threshold(mod.val) # view maxent thresholds
	mod.cut = as.data.frame(c(mod.cut$spec_sens, mod.cut$kappa))
	mod.cut$method = c("sensspec", "kappa")
	names(mod.cut)[1] = "pred"

	## generate confusion matrices
	mod.cfmat.specsens=table(datx[[2]], factor(as.numeric(datx$mod.pred>=mod.cut$pred[1])))
	mod.cfmat.kappa=table(datx[[2]], factor(as.numeric(datx$mod.pred>=mod.cut$pred[2])))

	## calculate model accuracies with standard deviation=F
	mod.acc=presence.absence.accuracy(datx, threshold=mod.cut$pred, st.dev=F) 
	tss=mod.acc$sensitivity + mod.acc$specificity - 1 # code TSS metric
	mod.acc=cbind(mod.acc[1:7],tss) # bind all metrics	
}

cv.accuracy.max = function(dat, mod, modl, x.fold, n.col) {
	#assign to cross folds and predict
	dat.xf = sample(rep(c(1:5), length=nrow(dat))) # vector of random xfolds
	mod.predXF = rep(0, length=nrow(dat)) # empty vector of 0
	for (j in 1:x.fold) {
		tr = dat[dat.xf!=j, ] # training not eq. i
		te = dat[dat.xf==j, ] # test eq. i
		mx = maxent(tr[2:n.col], tr[1]) # maxent model on training
		mod.predXF[dat.xf==j] = predict(mx, te) # predict to test
		}
	dat <- cbind(dat, mod.predXF)

	## evaluate model, determine optimal threshold
	mod.val=evaluate(p=dat[dat$PRESABS==1,2:9], a=dat[dat$PRESABS==0,2:9], mod) 
	mod.cutXF = threshold(mod.val)
	mod.cutXF = as.data.frame(c(mod.cutXF$spec_sens, mod.cutXF$kappa))
	mod.cutXF$method = c("sensspec", "kappa")
	names(mod.cutXF)[1] = "cvpred"

	## build testing dataframe using model predictions
	dat2XF=cbind(modl, dat[1], mod.predXF) 

	## calculate model accuracies with standard deviation=F
	mod.accXF = presence.absence.accuracy(dat2XF, threshold=mod.cutXF$cvpred, st.dev=F)
	tss = mod.accXF$sensitivity + mod.accXF$specificity - 1 # code TSS metric
	mod.accXF = cbind(mod.accXF[1:7], tss) # bind all metrics
}

ext.accuracy.max = function(mod, ext, modl) {
	## predict to new data
	mod.epred = predict(mod, ext) 

	## build testing dataframes using model predictions
	datx = cbind(modl, ext[3], mod.epred) # bind obs and predictions

	## evaluate model and find thresholds
	mod.eval = evaluate(p=ext[ext$PRESABS==1,c(58:60,66:68,70:71)], a=ext[ext$PRESABS==0,c(58:60,66:68, 70:71)], mod) 
	mod.ecut = threshold(mod.eval)
	mod.ecut = as.data.frame(c(mod.ecut$spec_sens, mod.ecut$kappa))
	mod.ecut$method = c("sensspec", "kappa")
	names(mod.ecut)[1] = "epred"

	## generate confusion matrix
	mod.ecfmat.sensspec = table(datx[[2]], factor(as.numeric(datx$mod.epred >= mod.ecut$epred[1])))
	mod.ecfmat.kappa = table(datx[[2]], factor(as.numeric(datx$mod.epred >= mod.ecut$epred[2])))

	## calculate model accuracies with standard deviation=F
	mod.eacc = presence.absence.accuracy(datx, threshold=mod.ecut$epred, st.dev=F) 
	tss = mod.eacc$sensitivity + mod.eacc$specificity - 1 # code TSS metric
	mod.eacc=cbind(mod.eacc[1:7], tss) # bind all metrics	
}