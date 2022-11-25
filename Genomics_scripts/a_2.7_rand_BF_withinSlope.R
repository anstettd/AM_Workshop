##################################################################################
## Select random SNPs within certain slope range
## Get BF mean of 1000 X and compare to observed data
## 
## Author Daniel Anstett
## 
## Last Modified October 12, 2022
###################################################################################
library(tidyverse)

#Import All Slopes
slope_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/slope_mat.csv")
#slope_map <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/slope_map.csv")
#slope_cmd <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/slope_cmd.csv")
slope_mat <- slope_mat %>% mutate(log_p=-log(approx_p,base=10))
slope_mat2 <- slope_mat %>% filter(Slope>0) %>% filter(Slope<2.5)

#Import slopes for SNP double filtered SNP set
mat_all <- read_csv("Genomics_scripts/Data/freq_MAT_slope_peakbf5.csv")
#map_all <- read_csv("Genomics_scripts/Data/freq_MAP_slope_peakbf5.csv")
#cmd_all <- read_csv("Genomics_scripts/Data/freq_CMD_slope_peakbf5.csv")

#Left join to get BF and window info into double filtered SNP set
mat_all <- mat_all %>% select(-Slope)
mat_all_join <- left_join(mat_all,slope_mat,by=c("Site"="Site","snp_ID"="chr_snp"))

############################################################################################
#Calculate BF Mean for all slope between 0 and 2.5

mat_mean <- data.frame() #Setup dataframe
mat_mean_obs <- data.frame()

#Test for 1 pop

for(j in 1:12){
  mat_all_p1 <- mat_all_join %>% filter(Site==j)
  mat_mean_obs[j,1] <- mean(mat_all_p1$BF,na.rm=TRUE)
  for(i in 1:1000){
    set.seed(i)
    index <- sample(1:nrow(slope_mat), dim(mat_all_p1)[1], replace=FALSE)
    rand_sample <- slope_mat[index,]
    mat_mean[i,j] <- mean(rand_sample$BF)
  }
}

#full range
#slope_mat2_p1 <- slope_mat %>% filter(Site==2)
bf_hist_p1 <- ggplot(mat_mean, aes(x=V1)) +
  geom_histogram()+
  geom_vline(xintercept = mat_mean_obs[1,1])+
  scale_y_continuous(name="Log Bayes Factor")+
  scale_x_continuous(name="SNP Count",breaks = seq(-10, 20, by = 10))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
bf_hist_p1
ggsave("Graphs_Oct_22/bf_hist_p1.pdf",plot=bf_hist_p1, width=8, height = 6, units = "in")







###########################################################################################
#Single Case
set.seed(1)
index <- sample(1:nrow(slope_mat), dim(mat_all_p1)[1], replace=FALSE)
rand_sample <- slope_mat[index,]
mean(rand_sample$BF)






