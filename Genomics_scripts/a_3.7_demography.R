##################################################################################
## Get positive slope info per population
## Author Daniel Anstett
## 
## 
## Last Modified September 19, 2022
###################################################################################

#Library install and import
library(tidyverse) 

#Bring in dataframe with offset information
offset_pop <- read_csv("Genomics_scripts/Data/offset_pop_bf5.csv")
offset_pop <- offset_pop %>% mutate(Region = ifelse(Lat >= 40, "North", 
                                                        ifelse((Lat >35) & (Lat <40), "Centre","South")))

#Import lambda from demography data for all sites
demo_all <- read_csv("Genomics_scripts/Data/lambda_2010_2016.csv")
demo <- demo_all %>% filter(Site =="SweetwaterRiver" | Site == "WestForkMojaveRiver" | 
                              Site =="NorthForkMiddleForkTule" | Site == "RedwoodCreek" |
                              Site == "Wawona" | Site == "OregonCreek" |
                              Site == "LittleJamesonCreek" | Site == "DeepCreek" |
                              Site == "O'NeilCreek" | Site == "DeerCreek" |
                              Site == "RockCreek" | Site == "MillCreek")
site_popID <- read_csv("Genomics_scripts/Data/site_popID.csv")
demo_ID <- left_join(demo,site_popID,by="Site")
demo_ID <- demo_ID %>% filter(!Year%in%c(2010,2011))


###################################################################################
#Calc lambda slopes over time

lambda_pop <- data.frame()
for (i in 1:12){
  demo_pop <- demo_ID %>% filter(Paper_ID == i)
  lambda_pop[i,1] <- summary(lm(lambda~Year,data=demo_pop))$coefficients[2]
}

offset_pop_lambda <- cbind(offset_pop,lambda_pop)
colnames(offset_pop_lambda)[14] <- "offset_4.5"
colnames(offset_pop_lambda)[15] <- "offset_8.5"
colnames(offset_pop_lambda)[16] <- "offset_4.5_old"
colnames(offset_pop_lambda)[17] <- "offset_8.5_old"
colnames(offset_pop_lambda)[19] <- "lambda_slope"

write_csv(offset_pop_lambda,"Genomics_scripts/Data/offset_pop_lambda.csv")





