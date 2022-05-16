##################################################################################
## Make SNP slopes histograms
## Author Daniel Anstett
## 
## Currently done just for ENV 1,2 and 5 (MAT, MAP and CMD)
## Last Modified Nov 16, 2021
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

snp1A_slope <- read_csv("Genomics_scripts/Data/freq_MAT_slope.csv")
snp2A_slope <- read_csv("Genomics_scripts/Data/freq_MAP_slope.csv")
snp5A_slope <- read_csv("Genomics_scripts/Data/freq_CMD_slope.csv")

setwd("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files")

snp1A_bars <- read_csv("rand_slope_mat_multi.csv")
snp2A_bars <- read_csv("rand_slope_map_multi.csv")
snp5A_bars <- read_csv("rand_slope_cmd_multi.csv")

setwd("/Users/daniel_anstett/Dropbox/AM_Workshop/AM_Workshop")

###################################################################################
## Test case
###################################################################################


#Pop 1
snp1A6 <- snp1A_slope %>% filter(Site==1) %>% mutate(climate="(A) MAT")

#graph
snp1A_hist <- ggplot(snp1A6,aes(x=Slope))+
  geom_histogram(position="identity",binwidth = 0.2,color="black",fill="white")+
#  scale_fill_manual(values = c("Observed"="deepskyblue","Random"="hotpink"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (S02 Sweetwater)",limits = c(-2.5,2.5))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist
#ggsave("Graphs_overlap/Single_trim/trim_pop1_overlap.pdf",width=10, height = 5, units = "in")



# ALL pops

#snp1A 
snp1A_hist <- ggplot(snp1A_slope,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",fill="white")+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (MAT)",limits=c(-2.5,2.5))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~Site)
snp1A_hist
#ggsave("Graphs_overlap/trim_MAT_overlap.pdf",width=10, height = 7.5, units = "in")








