library(tidyverse)
library(gradientForest)

require(gradientForest)
load("C:/Users/Julia/Desktop/GZ.sps.mat.Rdata")

load("C:/Users/Julia/Desktop/GZ.phys.site.Rdata")

nSites <- dim(Sp_mat)[1]
nSpecs <- dim(Sp_mat)[2]
lev <- floor(log2(nSites * 0.368/2))
lev

df_in<-cbind(Phys_site, Sp_mat)


gf <- gradientForest(data = cbind(Phys_site, Sp_mat), predictor.vars = colnames(Phys_site),
                     response.vars = colnames(Sp_mat), ntree = 500, transform = NULL,
                     maxLevel = lev, corr.threshold = 0.5, compact = T, nbin = 201)


data(CoMLsimulation)
preds <- colnames(Xsimulation)
specs <- colnames(Ysimulation)

df_in_2<-data.frame(Ysimulation,Xsimulation)

f1 <- gradientForest(data.frame(Ysimulation,Xsimulation), preds, specs, ntree=10)
f1