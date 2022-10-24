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
slope_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/slope_mat.csv")
slope_mat_1 <- slope_mat %>% filter(Site==2)
rm(slope_mat)

slope_matn50 <- head(slope_mat_1 %>% filter(Slope < -25),1)
slope_matn25 <- head(slope_mat_1 %>% filter(Slope < -25) %>% filter(Slope > -30),1)
slope_matn5 <- head(slope_mat_1 %>% filter(Slope < -5) %>% filter(Slope > -7),1)
slope_matn2 <- head(slope_mat_1 %>% filter(Slope < -2) %>% filter(Slope > -5),1)

slope_mat50 <- head(slope_mat_1 %>% filter(Slope > 25),1)
slope_mat25 <- head(slope_mat_1 %>% filter(Slope > 25) %>% filter(Slope < 30),1)
slope_mat5 <- head(slope_mat_1 %>% filter(Slope > 5) %>% filter(Slope < 7),1)
slope_mat2 <- head(slope_mat_1 %>% filter(Slope > 2) %>% filter(Slope < 5),1)

rm(slope_mat_1)

#Get SNP ID for first of each
slope_eg <- data.frame()
slope_eg <- as.data.frame(rbind(slope_matn50,slope_matn25,slope_matn5,slope_matn2,
                  slope_mat50,slope_mat25,slope_mat5,slope_mat2)$chr_snp)
slope_eg[,2] <- c("1_n50","2_n25","3_n5","4_n2","5_p50","6_p25","7_p5","8_p2")

colnames(slope_eg) <- c("chr_snp","position")

#Use Unix to select columns from .csv files with all SNPS 
#See filter_in_unix.txt

#sed -n "1 s/CE10_chr1_231269.*//p" freq_MAT_slopes.csv | sed 's/[^,*]//g' | wc -c
#sed -n "1 s/CE10_chr2_41689666.*//p" freq_MAT_slopes.csv | sed 's/[^,*]//g' | wc -c	
#sed -n "1 s/CE10_chr1_305074.*//p" freq_MAT_slopes.csv | sed 's/[^,*]//g' | wc -c
#sed -n "1 s/CE10_chr1_177602.*//p" freq_MAT_slopes.csv | sed 's/[^,*]//g' | wc -c
#sed -n "1 s/CE10_chr1_232342.*//p" freq_MAT_slopes.csv | sed 's/[^,*]//g' | wc -c
#sed -n "1 s/CE10_chr1_2389716.*//p" freq_MAT_slopes.csv | sed 's/[^,*]//g' | wc -c
#sed -n "1 s/CE10_chr1_235316.*//p" freq_MAT_slopes.csv | sed 's/[^,*]//g' | wc -c
#sed -n "1 s/CE10_chr1_232854.*//p" freq_MAT_slopes.csv | sed 's/[^,*]//g' | wc -c
#cut -d "," -f 1,2,847,1244077,1918,327,873,21414,969,887 freq_MAT_slopes.csv > freq_MAT_filter_test.csv

#Import FAT table for selected columns
freq_MAT <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_filter_test.csv")
freq_MAT.melted = reshape2::melt(freq_MAT, id.vars = c("Site", "Year"))
colnames(freq_MAT.melted) <- c("Site","Year","chr_snp","Frequency")
freq_MAT.melted_pos <- left_join(freq_MAT.melted, slope_eg, by="chr_snp")
freq_MAT_2 <- freq_MAT.melted_pos %>% filter(Site==2)

#p2
ggplot(data=freq_MAT_2,aes(Year,Frequency)) + 
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
ggsave("Graphs_CI/eg_freq_change_mat.pdf",width=10, height = 6, units = "in")
















