##################################################################################
## Get positive slope info per population
## Author Daniel Anstett
## 
## 
## Last Modified Jan 25, 2023
###################################################################################

#Library install and import
library(tidyverse) 

#Import lambda from demography data for all sites
demo_all <- read_csv("Genomics_scripts/Data/lambda_2010_2016.csv")
demo_select <- demo_all %>% filter(Site !="BearCreek" & Site != "CherokeeCreek" & 
                              Site !="ChinoCreek" & Site != "CorralCreek" &
                              Site != "FiddleCreek" & Site != "HauserCreek" &
                              Site != "MossCreekMercedDrainage" & Site != "NorthForkTuolumne" &
                              Site != "ParadiseCreek" & Site != "SweetwaterCreek" &
                              Site != "WillowCreek" & Site != "KitchenCreek" &
                              Site != "WhitewaterCanyon") 
demo_select <- demo_select %>% filter(!Year%in%c(2010,2011))
demo_select <- na.omit(demo_select)


demo_order <- demo_select[order(demo_select$Latitude),]
demo_ID <- data.frame(unique(demo_order$Site))
demo_ID$ID<-1:nrow(demo_ID)
names(demo_ID) <- c("Site","ID")
demo_select_ID <- left_join(demo_select,demo_ID,by="Site")
demo_select_raw <- demo_select %>% select(-Year,-lambda,-SiteYear)
demo_ID_join <- left_join(demo_ID,unique(demo_select_raw),by="Site")

###################################################################################
#Calc lambda slopes over time

lambda_pop <- data.frame()
for (i in 1:nrow(demo_ID)){
  demo_pop <- demo_select_ID %>% filter(ID == i)
  lambda_pop[i,1] <- summary(lm(lambda~Year,data=demo_pop))$coefficients[2]
}

offset_pop_lambda <- cbind(demo_ID_join,lambda_pop)

colnames(offset_pop_lambda)[8] <- "lambda_slope"

#write_csv(offset_pop_lambda,"Genomics_scripts/Data/offset_pop_lambda_allsites.csv")





