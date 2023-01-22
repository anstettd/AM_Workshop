##################################################################################
## Kurtosis
## 
## Author Daniel Anstett
## 
## 
## Last Modified Dec 13, 2022
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)
library(moments)
library(reshape2)
library(qvalue)

###################################################################################
#function
emP <- function(obs,example){
  vec_l<-length(example[example<=obs])
  per_pval2 <- vec_l/length(example)[1] 
  return(per_pval2)
}






###################################################################################
#Import observed slopes
mat_slope <- read_csv("Genomics_scripts/Data/freq_MAT_slope_peakbf5.csv") %>% filter(Slope >-2.5, Slope<2.5)
map_slope <- read_csv("Genomics_scripts/Data/freq_MAP_slope_peakbf5.csv") %>% filter(Slope >-2.5, Slope<2.5)
cmd_slope <- read_csv("Genomics_scripts/Data/freq_cmd_slope_peakbf5.csv") %>% filter(Slope >-2.5, Slope<2.5)

setwd("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files")

#Import large files with 1000 X random slopes
mat_bars <- read_csv("rand_slope_mat_multi_peakbf5.csv") %>% filter(Slope >-2.5, Slope<2.5)
map_bars <- read_csv("rand_slope_map_multi_peakbf5.csv") %>% filter(Slope >-2.5, Slope<2.5)
cmd_bars <- read_csv("rand_slope_cmd_multi_peakbf5.csv") %>% filter(Slope >-2.5, Slope<2.5)

setwd("/Users/daniel_anstett/Dropbox/AM_Workshop/AM_Workshop")

#Kurtosis neutral set
kurt_mat <- data.frame()
kurt_map <- data.frame()
kurt_cmd <- data.frame()
for(j in 1:12){
  mat_bars_pop <- mat_bars %>% filter(Site==j)
  map_bars_pop <- map_bars %>% filter(Site==j)
  cmd_bars_pop <- cmd_bars %>% filter(Site==j)
  for (i in 1:1000){
    rand_snp_mat <- mat_bars_pop %>% filter(Seed_ID==i)
    rand_snp_map <- map_bars_pop %>% filter(Seed_ID==i)
    rand_snp_cmd <- cmd_bars_pop %>% filter(Seed_ID==i)
    #kurt[i,j] <- mean(rand_snp$Slope)
    kurt_mat[i,j] <- kurtosis(rand_snp_mat$Slope)
    kurt_map[i,j] <- kurtosis(rand_snp_map$Slope)
    kurt_cmd[i,j] <- kurtosis(rand_snp_cmd$Slope)
  }
}
colnames(kurt_mat) <- c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12")
colnames(kurt_map) <- c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12")
colnames(kurt_cmd) <- c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12")
#Not working
#kurt.melt <- reshape2::melt(kurt, 
#            id.vars = c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12"))


#Kurtosis observed
kurt_obs <- data.frame()
for (k in 1:12){
  obs_slopes_mat <- mat_slope %>% filter(Site==k)
  obs_slopes_map <- map_slope %>% filter(Site==k)
  obs_slopes_cmd <- cmd_slope %>% filter(Site==k)
  kurt_obs[k,1] <- kurtosis(obs_slopes_mat$Slope)
  kurt_obs[k,2] <- kurtosis(obs_slopes_map$Slope)
  kurt_obs[k,3] <- kurtosis(obs_slopes_cmd$Slope)
}

colnames(kurt_obs)  <- c("MAT","MAP","CMD")


############
#Empirical p-values
p_data <- data.frame()
for (i in 1:12){
  p_data[i,1] <- emP(kurt_obs[i,3],kurt_cmd[,i])
  p_data[i,2] <- emP(kurt_obs[i,2],kurt_map[,i])
  p_data[i,3] <- emP(kurt_obs[i,1],kurt_mat[,i])
}

colnames(p_data) <- c("CMD","MAP","MAT")

write_csv(p_data, "Genomics_scripts/Data/kurtosis_pval.csv")



############
#Histograms

mat_hist <- ggplot(kurt_cmd,aes(x=P3))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  geom_vline(xintercept=kurt_obs$CMD[3])+
  scale_y_continuous(name="Density")+
  scale_x_continuous(name="Kurtosis")+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist


mat_hist <- ggplot(kurt_cmd,aes(x=P4))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  geom_vline(xintercept=kurt_obs$CMD[4])+
  scale_y_continuous(name="Density")+
  scale_x_continuous(name="Kurtosis")+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist







