##################################################################################
## Make Single SNP Plots
## 
## Author Daniel Anstett
## 
## Last Modified October 12, 2022
###################################################################################
library(reshape2)
library(dplyr)
library(tidyverse)
###################################################################################
#Import Slopes
slope_CMD <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/slope_CMD.csv")
slope_CMD_1 <- slope_CMD %>% filter(Site==2)
#rm(slope_CMD)

slope_CMDn50 <- head(slope_CMD_1 %>% filter(Slope < -25),1)
slope_CMDn25 <- head(slope_CMD_1 %>% filter(Slope < -25) %>% filter(Slope > -30),1)
slope_CMDn5 <- head(slope_CMD_1 %>% filter(Slope < -5) %>% filter(Slope > -7),1)
slope_CMDn2 <- head(slope_CMD_1 %>% filter(Slope < -2) %>% filter(Slope > -5),1)

slope_CMD50 <- head(slope_CMD_1 %>% filter(Slope > 25),1)
slope_CMD25 <- head(slope_CMD_1 %>% filter(Slope > 25) %>% filter(Slope < 30),1)
slope_CMD5 <- head(slope_CMD_1 %>% filter(Slope > 5) %>% filter(Slope < 7),1)
slope_CMD2 <- head(slope_CMD_1 %>% filter(Slope > 2) %>% filter(Slope < 5),1)

#rm(slope_CMD_1)

#Get SNP ID for first of each
slope_eg <- data.frame()
slope_eg <- as.data.frame(rbind(slope_CMDn50,slope_CMDn25,slope_CMDn5,slope_CMDn2,
                  slope_CMD50,slope_CMD25,slope_CMD5,slope_CMD2)$chr_snp)
slope_eg[,2] <- c("1_n50","2_n25","3_n5","4_n2","5_p50","6_p25","7_p5","8_p2")

colnames(slope_eg) <- c("chr_snp","position")

#Use Unix to select columns from .csv files with all SNPS 
#See filter_in_unix.txt

#Import FAT table for selected columns
freq_CMD <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_filter_test.csv")
freq_CMD.melted = reshape2::melt(freq_CMD, id.vars = c("Site", "Year"))
colnames(freq_CMD.melted) <- c("Site","Year","chr_snp","Frequency")
freq_CMD.melted_pos <- left_join(freq_CMD.melted, slope_eg, by="chr_snp")
freq_CMD_2 <- freq_CMD.melted_pos %>% filter(Site==2)

#p2
ggplot(data=freq_CMD_2,aes(Year,Frequency)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F, alpha=.25,cex=0.6) + theme_classic() +
  geom_point()+
  labs(y="SNP Frequency",x="Year") + 
  theme(legend.position = "none",
    axis.text.x = element_text(size = 14, face = "bold", angle = 45,hjust = 1, vjust = 1), 
    axis.title = element_text(size =16, face = "bold"), 
    axis.text.y = element_text(size = 14, face = "bold")) + 
  facet_wrap(.~position)
             #,labeller = labeller(env=env_labs)) +
#  theme(legend.title = element_blank(),
#        legend.text = element_text(size=12,face="bold"),
#        strip.background = element_blank(), 
#        strip.text.x=element_text(size=14,face="bold",hjust=0,vjust=-1.2))
ggsave("Graphs_CI/eg_freq_change_CMD.pdf",width=10, height = 6, units = "in")
















