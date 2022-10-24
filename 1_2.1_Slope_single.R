##################################################################################
## Get slopes (strength of selection) for all BF
## 
## Author Jose Lazaro
## 
## Last Modified October 6, 2022
###################################################################################

library(reshape2)
library(dplyr)
library(tidyverse)

###################################################################################
#Data Import
freq_MAT_1 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_1.csv")

freq_MAT_1.melted = reshape2::melt(freq_MAT_1, id.vars = c("Site", "Year"))
colnames(freq_MAT_1.melted)[3] <- "snp_ID"
colnames(freq_MAT_1.melted)[4] <- "freq"

freq_MAT_1.melted_mod <- na.omit(freq_MAT_1.melted)

freq_MAT_1.slope <- group_by(freq_MAT_1.melted_mod, Site, snp_ID) %>%
  arrange(Site, Year) %>%
  summarize(Slope = glm(freq ~ Year, family = binomial)$coefficients[2])

#write_csv(freq_MAT_1.slope, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_1.csv")





