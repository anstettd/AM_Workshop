##################################################################################
## Make SNP slopes histograms
## Author Daniel Anstett
## 
## Currently done just for ENV 1,2 and 5 (MAT, MAP and CMD)
## Last Modified April 11, 2022
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

#Import frequency Table
freq_MAT_slope <- read_csv("Genomics_scripts/Data/freq_MAT_slope.csv")
freq_MAP_slope <- read_csv("Genomics_scripts/Data/freq_MAP_slope.csv")
freq_CMD_slope <- read_csv("Genomics_scripts/Data/freq_CMD_slope.csv")

#Filter per site
freq_MAT_1 <- freq_MAT_slope %>% filter (Site==1)



###################################################################################
#MAT Slope histogram
MAT_hist <- ggplot(freq_MAT_slope,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Slope")+
  theme_classic()
MAT_hist <- MAT_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
MAT_hist <- MAT_hist + facet_wrap(.~Site)
MAT_hist
ggsave("Slope_graphs/MAT_all.pdf",width=11, height = 7, units = "in")


#MAP Slope histogram
MAP_hist <- ggplot(freq_MAP_slope,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Slope")+
  theme_classic()
MAP_hist <- MAP_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
MAP_hist <- MAP_hist + facet_wrap(.~Site)
MAP_hist
ggsave("Slope_graphs/MAP_all.pdf",width=11, height = 7, units = "in")


#CMD Slope histogram
CMD_hist <- ggplot(freq_CMD_slope,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Slope")+
  theme_classic()
CMD_hist <- CMD_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
CMD_hist <- CMD_hist + facet_wrap(.~Site)
CMD_hist
ggsave("Slope_graphs/CMD_all.pdf",width=11, height = 7, units = "in")






