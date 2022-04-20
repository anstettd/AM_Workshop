##################################################################################
## Make SNP slopes histograms
## Author Daniel Anstett
## 
## Currently done just for ENV 1,2 and 5 (MAT, MAP and CMD)
## Last Modified Nov 16, 2021
###################################################################################
#Import libraries
library(tidyverse)

###################################################################################

snp1A_slope <- read_csv("Genomics_scripts/Data/freq_MAT_slope.csv")
snp2A_slope <- read_csv("Genomics_scripts/Data/freq_MAP_slope.csv")
snp5A_slope <- read_csv("Genomics_scripts/Data/freq_CMD_slope.csv")

# Singleton exampel
snp1A_slope_1 <- snp1A_slope %>% filter(Site==1)
wilcox.test(snp1A_slope_1$Slope, mu = 0, alternative = "greater")

#Loop
wil_p <- data.frame()

for (i in 1:12){
  snp1A_slope_1 <- snp1A_slope %>% filter(Site==i)
  snp2A_slope_1 <- snp2A_slope %>% filter(Site==i)
  snp5A_slope_1 <- snp5A_slope %>% filter(Site==i)
 will_test1 <- wilcox.test(snp1A_slope_1$Slope, mu = 0, alternative = "greater")
 will_test2 <- wilcox.test(snp2A_slope_1$Slope, mu = 0, alternative = "greater")
 will_test5 <- wilcox.test(snp5A_slope_1$Slope, mu = 0, alternative = "greater")
 wil_p[i,1] <-will_test1$p.value
 wil_p[i,2] <-will_test2$p.value
 wil_p[i,3] <-will_test5$p.value
}

colnames(wil_p) <- c("MAT","MAP","CMD")

write_csv(wil_p, "Genomics_scripts/Data/wilcox_slope_p-values.csv")


