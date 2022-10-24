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
slope_MAP <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/slope_map.csv")
slope_MAP_1 <- slope_MAP %>% filter(Site==2)
#rm(slope_MAP)

slope_MAPn50 <- head(slope_MAP_1 %>% filter(Slope < -25),1)
slope_MAPn25 <- head(slope_MAP_1 %>% filter(Slope < -25) %>% filter(Slope > -30),1)
slope_MAPn5 <- head(slope_MAP_1 %>% filter(Slope < -5) %>% filter(Slope > -7),1)
slope_MAPn2 <- head(slope_MAP_1 %>% filter(Slope < -2) %>% filter(Slope > -5),1)

slope_MAP50 <- head(slope_MAP_1 %>% filter(Slope > 25),1)
slope_MAP25 <- head(slope_MAP_1 %>% filter(Slope > 25) %>% filter(Slope < 30),1)
slope_MAP5 <- head(slope_MAP_1 %>% filter(Slope > 5) %>% filter(Slope < 7),1)
slope_MAP2 <- head(slope_MAP_1 %>% filter(Slope > 2) %>% filter(Slope < 5),1)

#rm(slope_MAP_1)

#Get SNP ID for first of each
slope_eg <- data.frame()
slope_eg <- as.data.frame(rbind(slope_MAPn50,slope_MAPn25,slope_MAPn5,slope_MAPn2,
                  slope_MAP50,slope_MAP25,slope_MAP5,slope_MAP2)$chr_snp)
slope_eg[,2] <- c("1_n50","2_n25","3_n5","4_n2","5_p50","6_p25","7_p5","8_p2")

colnames(slope_eg) <- c("chr_snp","position")

#Use Unix to select columns from .csv files with all SNPS 
#See filter_in_unix.txt

#Import FAT table for selected columns
freq_MAP <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_filter_test.csv")
freq_MAP.melted = reshape2::melt(freq_MAP, id.vars = c("Site", "Year"))
colnames(freq_MAP.melted) <- c("Site","Year","chr_snp","Frequency")
freq_MAP.melted_pos <- left_join(freq_MAP.melted, slope_eg, by="chr_snp")
freq_MAP_2 <- freq_MAP.melted_pos %>% filter(Site==2)

#p2
ggplot(data=freq_MAP_2,aes(Year,Frequency)) + 
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
ggsave("Graphs_CI/eg_freq_change_MAP.pdf",width=10, height = 6, units = "in")
















