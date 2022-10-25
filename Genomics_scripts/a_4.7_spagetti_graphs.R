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
freq_mat <- read_csv("Genomics_scripts/Data/freq_MAT_peakbf5.csv")
freq_map <- read_csv("Genomics_scripts/Data/freq_MAP_peakbf5.csv")
freq_cmd <- read_csv("Genomics_scripts/Data/freq_CMD_peakbf5.csv")

#Gather data frames
freq_mat <- freq_mat %>% gather(SNP_ID,SNP_Freq,3:dim(freq_mat)[2])
freq_map <- freq_map %>% gather(SNP_ID,SNP_Freq,3:dim(freq_map)[2])
freq_cmd <- freq_cmd %>% gather(SNP_ID,SNP_Freq,3:dim(freq_cmd)[2])

#Import gathered rand (seed=1) data frame from 6.3
rand_mat <- read_csv("Genomics_scripts/Data/rand_gathered_mat_peakbf5.csv")
rand_map <- read_csv("Genomics_scripts/Data/rand_gathered_map_peakbf5.csv")
rand_cmd <- read_csv("Genomics_scripts/Data/rand_gathered_cmd_peakbf5.csv")


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
joint_env_p1 <- joint_env %>% filter(Site==1)
joint_env_p2 <- joint_env %>% filter(Site==2)
joint_env_p3 <- joint_env %>% filter(Site==3)
joint_env_p4 <- joint_env %>% filter(Site==4)
joint_env_p5 <- joint_env %>% filter(Site==5)
joint_env_p6 <- joint_env %>% filter(Site==6)
joint_env_p7 <- joint_env %>% filter(Site==7)
joint_env_p8 <- joint_env %>% filter(Site==8)
joint_env_p9 <- joint_env %>% filter(Site==9)
joint_env_p10 <- joint_env %>% filter(Site==10)
joint_env_p11 <- joint_env %>% filter(Site==11)
joint_env_p12 <- joint_env %>% filter(Site==12)


###################################################################################
#For Paper
env_labs<-c("cmd"="A CMD",
            "map"="B MAP",
            "mat"="C MAT")
#p1
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
ggsave("Graphs_Oct_22/2_01env_both_freqchange_paper_p1.pdf",width=12, height = 6, units = "in")


#p12
ggplot(data=joint_env_p12,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
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
ggsave("Graphs_Oct_22/2_02env_both_freqchange_paper_p12.pdf",width=12, height = 6, units = "in")



#p2
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
ggsave("Graphs_Oct_22/2_03env_both_freqchange_paper_p2.pdf",width=12, height = 6, units = "in")

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
ggsave("Graphs_Oct_22/2_04env_both_freqchange_paper_p3.pdf",width=12, height = 6, units = "in")

#p4
ggplot(data=joint_env_p4,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
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
ggsave("Graphs_Oct_22/2_05env_both_freqchange_paper_p4.pdf",width=12, height = 6, units = "in")

#p5
ggplot(data=joint_env_p5,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
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
ggsave("Graphs_Oct_22/2_06env_both_freqchange_paper_p5.pdf",width=12, height = 6, units = "in")

#p6
ggplot(data=joint_env_p6,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
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
ggsave("Graphs_Oct_22/2_07env_both_freqchange_paper_p6.pdf",width=12, height = 6, units = "in")

#p7
ggplot(data=joint_env_p7,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
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
ggsave("Graphs_Oct_22/2_08env_both_freqchange_paper_p7.pdf",width=12, height = 6, units = "in")

#p8
ggplot(data=joint_env_p8,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
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
ggsave("Graphs_Oct_22/2_09env_both_freqchange_paper_p8.pdf",width=12, height = 6, units = "in")

#p9
ggplot(data=joint_env_p9,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
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
ggsave("Graphs_Oct_22/2_10env_both_freqchange_paper_p9.pdf",width=12, height = 6, units = "in")

#p10
ggplot(data=joint_env_p10,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
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
ggsave("Graphs_Oct_22/2_11env_both_freqchange_paper_p10.pdf",width=12, height = 6, units = "in")

#p11
ggplot(data=joint_env_p11,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
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
ggsave("Graphs_Oct_22/2_12env_both_freqchange_paper_p11.pdf",width=12, height = 6, units = "in")












