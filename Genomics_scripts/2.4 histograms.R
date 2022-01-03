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

snp1A_slope <- read_csv("Genomics_scripts/Data/snp1A_slope.csv")
snp2A_slope <- read_csv("Genomics_scripts/Data/snp2A_slope.csv")
snp5A_slope <- read_csv("Genomics_scripts/Data/snp5A_slope.csv")

snp1A_rslope <- read_csv("Genomics_scripts/Data/snp1A_rslope.csv")
snp2A_rslope <- read_csv("Genomics_scripts/Data/snp2A_rslope.csv")
snp5A_rslope <- read_csv("Genomics_scripts/Data/snp5A_rslope.csv")

###################################################################################
#Add Dummy Variable for data type
snp1A_slope <- snp1A_slope %>% mutate(data_type="Observed")
snp2A_slope <- snp2A_slope %>% mutate(data_type="Observed")
snp5A_slope <- snp5A_slope %>% mutate(data_type="Observed")

snp1A_rslope <- snp1A_rslope %>% mutate(data_type="Random")
snp2A_rslope <- snp2A_rslope %>% mutate(data_type="Random")
snp5A_rslope <- snp5A_rslope %>% mutate(data_type="Random")

#Make combined datasets
snp1A_all <- rbind(snp1A_slope,snp1A_rslope)
snp2A_all <- rbind(snp2A_slope,snp2A_rslope)
snp5A_all <- rbind(snp5A_slope,snp5A_rslope)

###################################################################################
## Make Figure 4
snp1A6 <- snp1A_all %>% filter(Site==6) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==6) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==6) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#Fig 4A
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue","Random"="hotpink"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Slope Through Time",limits = c(-0.01,0.2))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
ggsave("Graphs/snp6_all.pdf",width=10, height = 5, units = "in")


## Make Figure 5
snp5A1 <- snp5A_all %>% filter(Site==1) 
snp5A5 <- snp5A_all %>% filter(Site==5)
snp5A8 <- snp5A_all %>% filter(Site==8)
snp5A_158 <- rbind(snp5A1,snp5A5,snp5A8)

#Fig 5
snp5A_158_hist <- ggplot(snp5A_158,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue","Random"="hotpink"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Slope Through Time",limits = c(-0.01,0.2))+
  theme_classic()
snp5A_158_hist <- snp5A_158_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp5A_158_hist <- snp5A_158_hist + facet_wrap(.~Site)
snp5A_158_hist
ggsave("Graphs/snp5A_158.pdf",width=10, height = 5, units = "in")




###################################################################################
##Generate overlapped histograms

#snp1A 
snp1A_hist <- ggplot(snp1A_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="blue","Random"="red"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Slope Through Time")+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~Site)
snp1A_hist
ggsave("Graphs/snp1A_MAT_overlap.pdf",width=10, height = 7.5, units = "in")

#snp2A 
snp2A_hist <- ggplot(snp2A_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="blue","Random"="red"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Slope Through Time")+
  theme_classic()
snp2A_hist <- snp2A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp2A_hist <- snp2A_hist + facet_wrap(.~Site)
snp2A_hist
ggsave("Graphs/snp2A_MAP_overlap.pdf",width=10, height = 7.5, units = "in")


#snp5A 
snp5A_hist <- ggplot(snp5A_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="blue","Random"="red"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Slope Through Time",limits=c(-0.01,0.5))+
  theme_classic()
snp5A_hist <- snp5A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp5A_hist <- snp5A_hist + facet_wrap(.~Site)
snp5A_hist
ggsave("Graphs/snp5A_CMD_overlap.pdf",width=10, height = 7.5, units = "in")


###################################################################################
##Generate non overlaping histograms

#Generate snp1A histogram
snp1A_hist <- ggplot(snp1A_slope, aes(X=Slope))+
  geom_histogram(aes(Slope),bins=40,color="black", fill="skyblue")+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Slope Through Time")+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~Site)
snp1A_hist
ggsave("Graphs/snp1A_MAT.pdf",width=10, height = 7.5, units = "in")


#Generate snp2A histogram
snp2A_hist <- ggplot(snp2A_slope, aes(X=Slope))+
  geom_histogram(aes(Slope),bins=40,color="black", fill="skyblue")+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Slope Through Time")+
  theme_classic()
snp2A_hist <- snp2A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp2A_hist <- snp2A_hist + facet_wrap(.~Site)
snp2A_hist
ggsave("Graphs/snp2A_MAP.pdf",width=10, height = 7.5, units = "in")

#Generate snp5A histogram
snp5A_hist <- ggplot(snp5A_slope, aes(X=Slope))+
  geom_histogram(aes(Slope),bins=40,color="black", fill="skyblue")+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Slope Through Time")+
  theme_classic()
snp5A_hist <- snp5A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp5A_hist <- snp5A_hist + facet_wrap(.~Site)
snp5A_hist
ggsave("Graphs/snp5A_CMD.pdf",width=10, height = 7.5, units = "in")



###################################################################################
##Generate Random histograms

#Generate snp1A histogram
snp1A_rhist <- ggplot(snp1A_rslope, aes(X=Slope))+
  geom_histogram(aes(Slope),bins=40,color="black", fill="skyblue")+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Slope Through Time")+
  theme_classic()
snp1A_rhist <- snp1A_rhist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_rhist <- snp1A_rhist + facet_wrap(.~Site)
snp1A_rhist
ggsave("Graphs/snp1A_MAT_random.pdf",width=10, height = 7.5, units = "in")


#Generate snp2A histogram
snp2A_rhist <- ggplot(snp2A_rslope, aes(X=Slope))+
  geom_histogram(aes(Slope),bins=40,color="black", fill="skyblue")+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Slope Through Time")+
  theme_classic()
snp2A_rhist <- snp2A_rhist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp2A_rhist <- snp2A_rhist + facet_wrap(.~Site)
snp2A_rhist
ggsave("Graphs/snp2A_MAP_random.pdf",width=10, height = 7.5, units = "in")

#Generate snp5A histogram
snp5A_rhist <- ggplot(snp5A_rslope, aes(X=Slope))+
  geom_histogram(aes(Slope),bins=40,color="black", fill="skyblue")+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Slope Through Time")+
  theme_classic()
snp5A_rhist <- snp5A_rhist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp5A_rhist <- snp5A_rhist + facet_wrap(.~Site)
snp5A_rhist
ggsave("Graphs/snp5A_CMD_random.pdf",width=10, height = 7.5, units = "in")











