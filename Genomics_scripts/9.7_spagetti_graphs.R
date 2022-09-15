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
freq_mat <- read_csv("Genomics_scripts/Data/freq_MAT_region.csv")
freq_map <- read_csv("Genomics_scripts/Data/freq_MAP_region.csv")
freq_cmd <- read_csv("Genomics_scripts/Data/freq_CMD_region.csv")

#Gather data frames (go from wide to long)
freq_mat <- freq_mat %>% gather(SNP_ID,SNP_Freq,3:275)
freq_map <- freq_map %>% gather(SNP_ID,SNP_Freq,3:571)
freq_cmd <- freq_cmd %>% gather(SNP_ID,SNP_Freq,3:304)

#Import one stratified random permutation 
rand_mat <- read_csv("Genomics_scripts/Data/rand_gathered_MAT_region.csv")
rand_map <- read_csv("Genomics_scripts/Data/rand_gathered_MAP_region.csv")
rand_cmd <- read_csv("Genomics_scripts/Data/rand_gathered_CMD_region.csv")

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
  labs(y="SNP Frequency",x="Year") + scale_color_manual(values=c("blue","firebrick2")) + theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14, face = "bold", angle = 45,hjust = 1, vjust = 1), 
    axis.title = element_text(size =16, face = "bold"), 
    axis.text.y = element_text(size = 14, face = "bold")) + 
  facet_wrap(.~env,labeller = labeller(env=env_labs)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), 
        strip.text.x=element_text(size=14,face="bold",hjust=0,vjust=-1.2))
ggsave("Graphs_CI/Frequency/env_both_freqchange_paper_North.pdf",width=12, height = 6, units = "in")

#Center
ggplot(data=joint_env_p2,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F, alpha=.25,cex=0.6) + theme_classic() +
  labs(y="SNP Frequency",x="Year") + scale_color_manual(values=c("blue","firebrick2")) + theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14, face = "bold", angle = 45,hjust = 1, vjust = 1), 
    axis.title = element_text(size =16, face = "bold"), 
    axis.text.y = element_text(size = 14, face = "bold")) + 
  facet_wrap(.~env,labeller = labeller(env=env_labs)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), 
        strip.text.x=element_text(size=14,face="bold",hjust=0,vjust=-1.2))
ggsave("Graphs_CI/Frequency/env_both_freqchange_paper_Center.pdf",width=12, height = 6, units = "in")

#p3
ggplot(data=joint_env_p3,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F, alpha=.25,cex=0.6) + theme_classic() +
  labs(y="SNP Frequency",x="Year") + scale_color_manual(values=c("blue","firebrick2")) + theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14, face = "bold", angle = 45,hjust = 1, vjust = 1), 
    axis.title = element_text(size =16, face = "bold"), 
    axis.text.y = element_text(size = 14, face = "bold")) + 
  facet_wrap(.~env,labeller = labeller(env=env_labs)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), 
        strip.text.x=element_text(size=14,face="bold",hjust=0,vjust=-1.2))
ggsave("Graphs_CI/Frequency/env_both_freqchange_paper_South.pdf",width=12, height = 6, units = "in")


###################################################################################
##Separate graphs
##mat

#Site 1
freq_map_p1 <- freq_map %>% filter(Site==4)
rand_map_p1 <- rand_map %>% filter(Site==4)
joint_map_p1 <- joint_map %>% filter(Site==4)


#Random
ggplot(data=rand_map_p1,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F, alpha=.4,cex=0.6) + theme_classic() +
  labs(y="SNP Frequency Random",x="Year") + scale_color_manual(values=c("firebrick2")) + 
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
#ggsave("Graphs_CI/map_random_freqchange_p4.pdf",width=10, height = 6, units = "in")

#Climape Associated
ggplot(data=freq_map_p1,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F, alpha=.4,cex=0.6) + theme_classic() +
  labs(y="SNP Frequency Climape",x="Year") + scale_color_manual(values=c("blue")) + 
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
#ggsave("Graphs_CI/map_obs_freqchange_p4.pdf",width=10, height = 6, units = "in")

#Both
ggplot(data=joint_map_p1,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F, alpha=.2,cex=0.6) + theme_classic() +
  labs(y="SNP Frequency",x="Year") + scale_color_manual(values=c("blue","firebrick2")) + theme(
    axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
    axis.title = element_text(size = 20, face = "bold"), 
    axis.text.y = element_text(size = 16, face = "bold"))
#ggsave("Graphs_CI/map_both_freqchange_p4.pdf",width=10, height = 6, units = "in")



###################################################################################
#For Poster
#p4
ggplot(data=joint_env_p1,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F, alpha=.25,cex=0.6) + theme_classic() +
  #facet_grid(SNP_Freq_Bin) +
  labs(y="SNP Frequency",x="Year") + scale_color_manual(values=c("blue","firebrick2")) + theme(
    legend.position = "none",
    axis.text.x = element_text(size = 20, face = "bold", angle = 45,hjust = 1, vjust = 1), 
    axis.title = element_text(size =0, face = "bold"), 
    axis.text.y = element_text(size = 20, face = "bold")) + 
  facet_wrap(~env)
##ggsave("Graphs_CI/Frequency/Poster/env_both_freqchange_p4.pdf",width=12, height = 6, units = "in")

# Inspecting differences by allele frequency starting points
ggplot(data=joint_env_p1[joint_env_p1$env=="mat",], aes(Year,SNP_Freq,group=SNP_ID, color=Type)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F, alpha=.25,cex=0.6) + theme_classic() +
  labs(y="SNP Frequency",x="Year") + scale_color_manual(values=c("blue","firebrick2")) + theme(
    legend.position = "none",
    axis.text.x = element_text(size = 20, face = "bold", angle = 45,hjust = 1, vjust = 1), 
    axis.title = element_text(size =0, face = "bold"), 
    axis.text.y = element_text(size = 20, face = "bold")) + 
  facet_wrap(~SNP_Freq_Bin2)






