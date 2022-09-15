##################################################################################
## Regression plots of BF>20 SNPs for random and observed
## Author Daniel Anstett
## 
## 
## Last Modified May 19, 2021
###################################################################################
#Import libraries
library(tidyverse)

#Import timeseries frequencies
freq_mat <- read_csv("Genomics_scripts/Data/freq_MAT_peakbf2_region.csv")
freq_map <- read_csv("Genomics_scripts/Data/freq_MAP_peakbf2_region.csv")
freq_cmd <- read_csv("Genomics_scripts/Data/freq_CMD_peakbf2_region.csv")

#Gather data frames (go from wide to long)
freq_mat <- freq_mat %>% gather(SNP_ID,SNP_Freq,3:dim(freq_mat)[2])
freq_map <- freq_map %>% gather(SNP_ID,SNP_Freq,3:dim(freq_map)[2])
freq_cmd <- freq_cmd %>% gather(SNP_ID,SNP_Freq,3:dim(freq_cmd)[2])

#Import one stratified random permutation 
rand_mat <- read_csv("Genomics_scripts/Data/rand_gathered_MAT_peakbf2_region.csv")
rand_map <- read_csv("Genomics_scripts/Data/rand_gathered_MAP_peakbf2_region.csv")
rand_cmd <- read_csv("Genomics_scripts/Data/rand_gathered_CMD_peakbf2_region.csv")

#Add Type
freq_mat$Type <- "Climate Associated"
freq_map$Type <- "Climate Associated"
freq_cmd$Type <- "Climate Associated"

rand_mat$Type <- "Random"
rand_map$Type <- "Random"
rand_cmd$Type <- "Random"

#Merge Climate Associated and Random
joint_mat <- rbind(freq_mat,rand_mat)
joint_map <- rbind(freq_map,rand_map)
joint_cmd <- rbind(freq_cmd,rand_cmd)

#Add env info
joint_mat$env <- "mat"
joint_map$env <- "map"
joint_cmd$env <- "cmd"

#Join
joint_env <- rbind(joint_mat,joint_map,joint_cmd)

#Make categorical variable for frequency bin
joint_env <- joint_env %>% 
  group_by(Site) %>% 
  mutate(Base_Year = min(Year)) 

freq_bin_calc <- joint_env %>% 
  group_by(Site, Type, env) %>% 
  filter(Year==Base_Year) %>% 
  mutate(SNP_Freq_Bin = cut(SNP_Freq,5, labels=c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0")))

joint_env <- left_join(joint_env, freq_bin_calc) %>% 
  group_by(Site, Type, env, SNP_ID) %>% 
  mutate(SNP_Freq_Bin2 = first(SNP_Freq_Bin))

###################################################################################
#Setup plots
joint_env_p1 <- joint_env %>% filter(Site=="North")
joint_env_p2 <- joint_env %>% filter(Site=="Centre")
joint_env_p3 <- joint_env %>% filter(Site=="South")


###################################################################################
#For Paper
env_labs<-c("cmd"="A CMD",
            "map"="B MAP",
            "mat"="C MAT")
#North
ggplot(data=joint_env_p1,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F, alpha=.25,cex=0.6) + theme_classic() +
  labs(y="SNP Frequency",x="Year (North)") + scale_color_manual(values=c("blue","firebrick2")) + theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14, face = "bold", angle = 45,hjust = 1, vjust = 1), 
    axis.title = element_text(size =16, face = "bold"), 
    axis.text.y = element_text(size = 14, face = "bold")) + 
  facet_wrap(.~env,labeller = labeller(env=env_labs)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), 
        strip.text.x=element_text(size=14,face="bold",hjust=0,vjust=-1.2))
ggsave("Graphs_CI_peak/2_1spagetti_peakbf2_north.pdf",width=12, height = 6, units = "in")

#Center
ggplot(data=joint_env_p2,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F, alpha=.25,cex=0.6) + theme_classic() +
  labs(y="SNP Frequency",x="Year (Centre)") + scale_color_manual(values=c("blue","firebrick2")) + theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14, face = "bold", angle = 45,hjust = 1, vjust = 1), 
    axis.title = element_text(size =16, face = "bold"), 
    axis.text.y = element_text(size = 14, face = "bold")) + 
  facet_wrap(.~env,labeller = labeller(env=env_labs)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), 
        strip.text.x=element_text(size=14,face="bold",hjust=0,vjust=-1.2))
ggsave("Graphs_CI_peak/2_2spagetti_peakbf2_centre.pdf",width=12, height = 6, units = "in")

#South
ggplot(data=joint_env_p3,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F, alpha=.25,cex=0.6) + theme_classic() +
  labs(y="SNP Frequency",x="Year (South)") + scale_color_manual(values=c("blue","firebrick2")) + theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14, face = "bold", angle = 45,hjust = 1, vjust = 1), 
    axis.title = element_text(size =16, face = "bold"), 
    axis.text.y = element_text(size = 14, face = "bold")) + 
  facet_wrap(.~env,labeller = labeller(env=env_labs)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), 
        strip.text.x=element_text(size=14,face="bold",hjust=0,vjust=-1.2))
ggsave("Graphs_CI_peak/2_3spagetti_peakbf2_south.pdf",width=12, height = 6, units = "in")




