## Calibration functions to test for spatial variation in occupancy
## From Wintle et al. 2005 (see wintle tutorial file for other attributions and citations)

"calib.mod" <- function(obs,preds){
	pred<-preds+0.00001 #first jitter preds of 0 or 1
	pred[pred>=1]<-0.99999
	mod<-glm(obs ~ log((pred)/(1-(pred))), family=binomial)
	lp<-log((pred)/(1-(pred))) #ie forcing a to 0 and b to 1
	a0b1<-glm(obs ~ offset(lp) - 1, family = binomial)
	miller1<- 1 - pchisq(a0b1$deviance - mod$deviance, 2)
	ab1<-glm(obs ~ offset(lp), family = binomial) #this allows an intercept
	miller2<- 1 - pchisq(a0b1$deviance - ab1$deviance, 1)
	miller3<- 1 - pchisq(ab1$deviance - mod$deviance, 1)
	calib.res <- list(mod,mod$coef,miller1,miller2,miller3) # gather results
	names(calib.res) <- c("model","calib.coeffs","testa0b1","testa0|b1","testb1|a")
	return(calib.res)
	}


"calib" <- function(obs, preds) {
		pred <- preds + 1e-005
		pred[pred >= 1] <- 0.99999
		mod <- glm(obs ~ log((pred)/(1 - (pred))), family = binomial)
		return(summary(mod))
		return(mod$coef[2])
		}


"graph.calib" <- function(obs, preds, x) {
		cut.prob.1 <- cut(preds, x)
		prob.bin.1 <- as.vector(tapply(preds,cut.prob.1,mean))
		prob.aves.1 <- as.vector(tapply(obs,cut.prob.1,sum)/table(cut.prob.1))
		se.1 <- sqrt(prob.aves.1*(1-prob.aves.1)/table(cut.prob.1))
		up.1 <- prob.aves.1 + as.vector(se.1)
		down.1 <- prob.aves.1 - as.vector(se.1)
		plot(obs ~ preds)
		points(prob.aves.1 ~ prob.bin.1, cex=2, col="red")
		for (i in 1:x){
			lines(c(prob.bin.1[i],prob.bin.1[i]),c(up.1[i],down.1[i]), col="red")
			}
		}

