##################################################################################
## Make SNP slopes histograms for regions for just observed data (no random)
## Author Daniel Anstett
## 
## Currently done just for ENV 1,2 and 5 (MAT, MAP and CMD)
## Last Modified July 28, 2028
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

snp1A_slope <- read_csv("Genomics_scripts/Data/freq_MAT_slope_region.csv")
snp2A_slope <- read_csv("Genomics_scripts/Data/freq_MAP_slope_region.csv")
snp5A_slope <- read_csv("Genomics_scripts/Data/freq_CMD_slope_region.csv")

#snp1A_rslope <- read_csv("Genomics_scripts/Data/rand_slope_mat.csv")
#snp2A_rslope <- read_csv("Genomics_scripts/Data/rand_slope_map.csv")
#snp5A_rslope <- read_csv("Genomics_scripts/Data/rand_slope_cmd.csv")

###################################################################################
#Add Dummy Variable for data type
snp1A_all <- snp1A_slope %>% mutate(data_type="Observed")
snp2A_all <- snp2A_slope %>% mutate(data_type="Observed")
snp5A_all <- snp5A_slope %>% mutate(data_type="Observed")

#snp1A_rslope <- snp1A_rslope %>% mutate(data_type="Random")
#snp2A_rslope <- snp2A_rslope %>% mutate(data_type="Random")
#snp5A_rslope <- snp5A_rslope %>% mutate(data_type="Random")

###################################################################################
## All population full graphs
###################################################################################

#snp1A 
snp1A_hist <- ggplot(snp1A_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (MAT)")+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~Site)
snp1A_hist
#ggsave("Graphs_overlap/full_MAT_overlap.pdf",width=10, height = 7.5, units = "in")

#snp2A 
snp2A_hist <- ggplot(snp2A_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.35)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (MAP)")+
  theme_classic()
snp2A_hist <- snp2A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp2A_hist <- snp2A_hist + facet_wrap(.~Site)
snp2A_hist
#ggsave("Graphs_overlap/full_MAP_overlap.pdf",width=10, height = 7.5, units = "in")


#snp5A 
snp5A_hist <- ggplot(snp5A_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (CMD)")+
  #,limits=c(-0.01,0.5))+
  theme_classic()
snp5A_hist <- snp5A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp5A_hist <- snp5A_hist + facet_wrap(.~Site)
snp5A_hist
#ggsave("Graphs_overlap/full_CMD_overlap.pdf",width=10, height = 7.5, units = "in")


###################################################################################
## All population trim 2 graphs
###################################################################################

#snp1A 
snp1A_hist <- ggplot(snp1A_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="#440154FF"))+
  #scale_fill_manual(values = c("Observed"="deepskyblue"))+
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

#snp2A 
snp2A_hist <- ggplot(snp2A_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="#440154FF"))+
  #scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (MAP)",limits=c(-2.5,2.5))+
  theme_classic()
snp2A_hist <- snp2A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp2A_hist <- snp2A_hist + facet_wrap(.~Site)
snp2A_hist
#ggsave("Graphs_overlap/trim_MAP_overlap.pdf",width=10, height = 7.5, units = "in")


#snp5A 
snp5A_hist <- ggplot(snp5A_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="#440154FF"))+
  #scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (CMD)",limits=c(-1,1))+
  theme_classic()
snp5A_hist <- snp5A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp5A_hist <- snp5A_hist + facet_wrap(.~Site)
snp5A_hist
#ggsave("Graphs_overlap/trim_CMD_overlap.pdf",width=10, height = 7.5, units = "in")




###################################################################################
## Per population trim 2 graphs
###################################################################################
#Pop 1
snp1A6 <- snp1A_all %>% filter(Site=="North") %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site=="North") %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site=="North") %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (North)",limits = c(-2.5,2.5))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single_trim/trim_pop1_overlap.pdf",width=10, height = 5, units = "in")

#Pop 2
snp1A6 <- snp1A_all %>% filter(Site=="Centre") %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site=="Centre") %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site=="Centre") %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Center)",limits = c(-2.5,2.5))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single_trim/trim_pop2_overlap.pdf",width=10, height = 5, units = "in")


#Pop 3
snp1A6 <- snp1A_all %>% filter(Site=="South") %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site=="South") %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site=="South") %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (South)",limits = c(-2.5,2.5))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single_trim/trim_pop3_overlap.pdf",width=10, height = 5, units = "in")




###################################################################################
## Per population trim 1 graphs
###################################################################################
#Pop 1
snp1A6 <- snp1A_all %>% filter(Site==1) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==1) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==1) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (S02 Sweetwater)",limits = c(-1,1))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single_trim1/trim1_pop1_overlap.pdf",width=10, height = 5, units = "in")

#Pop 2
snp1A6 <- snp1A_all %>% filter(Site==2) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==2) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==2) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (S07 WF Mojave)",limits = c(-1,1))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single_trim1/trim1_pop2_overlap.pdf",width=10, height = 5, units = "in")


#Pop 3
snp1A6 <- snp1A_all %>% filter(Site==3) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==3) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==3) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (S10 NFMF Tule)",limits = c(-1,1))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single_trim1/trim1_pop3_overlap.pdf",width=10, height = 5, units = "in")


#Pop 4
snp1A6 <- snp1A_all %>% filter(Site==4) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==4) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==4) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (S08 Readwoods)",limits = c(-1,1))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single_trim1/trim1_pop4_overlap.pdf",width=10, height = 5, units = "in")


#Pop 5
snp1A6 <- snp1A_all %>% filter(Site==5) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==5) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==5) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (S32 Wawona)",limits = c(-1,1))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single_trim1/trim1_pop5_overlap.pdf",width=10, height = 5, units = "in")

#Pop 6
snp1A6 <- snp1A_all %>% filter(Site==6) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==6) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==6) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.7)+
  scale_fill_manual(values = c("Observed"="#440154FF"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (S29 Oregon Creek)",limits = c(-1,1))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single_trim1/trim1_pop6_overlap.pdf",width=10, height = 5, units = "in")

#Pop 7
snp1A6 <- snp1A_all %>% filter(Site==7) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==7) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==7) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Little Jamison S18)",limits = c(-1,1))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single_trim1/trim1_pop7_overlap.pdf",width=10, height = 5, units = "in")

#Pop 8
snp1A6 <- snp1A_all %>% filter(Site==8) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==8) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==8) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (S17 Deep Creek)",limits = c(-1,1))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single_trim1/trim1_pop8_overlap.pdf",width=10, height = 5, units = "in")

#Pop 9
snp1A6 <- snp1A_all %>% filter(Site==9) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==9) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==9) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (S16 O'Neil Creek)",limits = c(-1,1))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single_trim1/trim1_pop9_overlap.pdf",width=10, height = 5, units = "in")


#Pop 10
snp1A6 <- snp1A_all %>% filter(Site==10) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==10) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==10) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (S36 Deer Creek)",limits = c(-1,1))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single_trim1/trim1_pop10_overlap.pdf",width=10, height = 5, units = "in")


#Pop 11
snp1A6 <- snp1A_all %>% filter(Site==11) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==11) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==11) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (S15 Rock Creek)",limits = c(-1,1))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single_trim1/trim1_pop11_overlap.pdf",width=10, height = 5, units = "in")

#Pop 12
snp1A6 <- snp1A_all %>% filter(Site==12) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==12) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==12) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (S11 Mill Creek)",limits = c(-1,1))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single_trim1/trim1_pop12_overlap.pdf",width=10, height = 5, units = "in")




###################################################################################
## Per population full graphs
###################################################################################

#Pop 1
snp1A6 <- snp1A_all %>% filter(Site==1) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==1) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==1) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
                     #,limits = c(-0.01,0.2))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single/pop1_overlap.pdf",width=10, height = 5, units = "in")

#Pop 2
snp1A6 <- snp1A_all %>% filter(Site==2) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==2) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==2) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  #,limits = c(-0.01,0.2))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single/pop2_overlap.pdf",width=10, height = 5, units = "in")


#Pop 3
snp1A6 <- snp1A_all %>% filter(Site==3) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==3) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==3) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  #,limits = c(-0.01,0.2))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single/pop3_overlap.pdf",width=10, height = 5, units = "in")


#Pop 4
snp1A6 <- snp1A_all %>% filter(Site==4) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==4) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==4) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  #,limits = c(-0.01,0.2))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single/pop4_overlap.pdf",width=10, height = 5, units = "in")


#Pop 5
snp1A6 <- snp1A_all %>% filter(Site==5) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==5) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==5) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  #,limits = c(-0.01,0.2))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single/pop5_overlap.pdf",width=10, height = 5, units = "in")

#Pop 6
snp1A6 <- snp1A_all %>% filter(Site==6) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==6) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==6) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  #,limits = c(-0.01,0.2))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single/pop6_overlap.pdf",width=10, height = 5, units = "in")

#Pop 7
snp1A6 <- snp1A_all %>% filter(Site==7) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==7) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==7) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  #,limits = c(-0.01,0.2))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single/pop7_overlap.pdf",width=10, height = 5, units = "in")

#Pop 8
snp1A6 <- snp1A_all %>% filter(Site==8) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==8) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==8) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  #,limits = c(-0.01,0.2))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single/pop8_overlap.pdf",width=10, height = 5, units = "in")

#Pop 9
snp1A6 <- snp1A_all %>% filter(Site==9) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==9) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==9) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  #,limits = c(-0.01,0.2))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single/pop9_overlap.pdf",width=10, height = 5, units = "in")


#Pop 10
snp1A6 <- snp1A_all %>% filter(Site==10) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==10) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==10) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  #,limits = c(-0.01,0.2))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single/pop10_overlap.pdf",width=10, height = 5, units = "in")


#Pop 11
snp1A6 <- snp1A_all %>% filter(Site==11) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==11) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==11) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (S15 Rock Creek)")+
  #,limits = c(-0.01,0.2))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single/pop11_overlap.pdf",width=10, height = 5, units = "in")

#Pop 12
snp1A6 <- snp1A_all %>% filter(Site==12) %>% mutate(climate="(A) MAT")
snp2A6 <- snp2A_all %>% filter(Site==12) %>% mutate(climate="(B) MAP")
snp5A6 <- snp5A_all %>% filter(Site==12) %>% mutate(climate="(C) CMD")
snpA6_all <- rbind(snp1A6,snp2A6,snp5A6)

#graph
snp1A_hist <- ggplot(snpA6_all,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  #,limits = c(-0.01,0.2))+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~climate)
snp1A_hist
#ggsave("Graphs_overlap/Single/pop12_overlap.pdf",width=10, height = 5, units = "in")





###################################################################################
#Not as useful
###################################################################################
## Make Figure 5
snp5A1 <- snp5A_all %>% filter(Site==1) 
snp5A5 <- snp5A_all %>% filter(Site==5)
snp5A8 <- snp5A_all %>% filter(Site==8)
snp5A_158 <- rbind(snp5A1,snp5A5,snp5A8)

#Fig 5
snp5A_158_hist <- ggplot(snp5A_158,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_fill_manual(values = c("Observed"="deepskyblue"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
#,limits = c(-0.01,0.2))+
  theme_classic()
snp5A_158_hist <- snp5A_158_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp5A_158_hist <- snp5A_158_hist + facet_wrap(.~Site)
snp5A_158_hist
##ggsave("Graphs/snp5A_158.pdf",width=10, height = 5, units = "in")
 

###################################################################################
##Generate non overlaping histograms

#Generate snp1A histogram
snp1A_hist <- ggplot(snp1A_slope, aes(X=Slope))+
  geom_histogram(aes(Slope),bins=40,color="black", fill="skyblue")+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  theme_classic()
snp1A_hist <- snp1A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_hist <- snp1A_hist + facet_wrap(.~Site)
snp1A_hist
##ggsave("Graphs/snp1A_MAT.pdf",width=10, height = 7.5, units = "in")


#Generate snp2A histogram
snp2A_hist <- ggplot(snp2A_slope, aes(X=Slope))+
  geom_histogram(aes(Slope),bins=40,color="black", fill="skyblue")+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  theme_classic()
snp2A_hist <- snp2A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp2A_hist <- snp2A_hist + facet_wrap(.~Site)
snp2A_hist
##ggsave("Graphs/snp2A_MAP.pdf",width=10, height = 7.5, units = "in")

#Generate snp5A histogram
snp5A_hist <- ggplot(snp5A_slope, aes(X=Slope))+
  geom_histogram(aes(Slope),bins=40,color="black", fill="skyblue")+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  theme_classic()
snp5A_hist <- snp5A_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp5A_hist <- snp5A_hist + facet_wrap(.~Site)
snp5A_hist
##ggsave("Graphs/snp5A_CMD.pdf",width=10, height = 7.5, units = "in")



###################################################################################
##Generate Random histograms

#Generate snp1A histogram
snp1A_rhist <- ggplot(snp1A_rslope, aes(X=Slope))+
  geom_histogram(aes(Slope),bins=40,color="black", fill="skyblue")+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  theme_classic()
snp1A_rhist <- snp1A_rhist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp1A_rhist <- snp1A_rhist + facet_wrap(.~Site)
snp1A_rhist
##ggsave("Graphs/snp1A_MAT_random.pdf",width=10, height = 7.5, units = "in")


#Generate snp2A histogram
snp2A_rhist <- ggplot(snp2A_rslope, aes(X=Slope))+
  geom_histogram(aes(Slope),bins=40,color="black", fill="skyblue")+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  theme_classic()
snp2A_rhist <- snp2A_rhist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp2A_rhist <- snp2A_rhist + facet_wrap(.~Site)
snp2A_rhist
##ggsave("Graphs/snp2A_MAP_random.pdf",width=10, height = 7.5, units = "in")

#Generate snp5A histogram
snp5A_rhist <- ggplot(snp5A_rslope, aes(X=Slope))+
  geom_histogram(aes(Slope),bins=40,color="black", fill="skyblue")+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  theme_classic()
snp5A_rhist <- snp5A_rhist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
snp5A_rhist <- snp5A_rhist + facet_wrap(.~Site)
snp5A_rhist
##ggsave("Graphs/snp5A_CMD_random.pdf",width=10, height = 7.5, units = "in")




rand = as.data.frame(rnorm(1000, mean=0, sd=0.1))
rand$data_type="Random"
colnames(rand)=c("Slope", "data_type")
obs = as.data.frame(rnorm(1000, mean=0.2, sd=0.15))
obs$data_type="Observed"
colnames(obs)=c("Slope", "data_type")
dummy.dat <- rbind(rand, obs)

demo_hist <- ggplot(dummy.dat,aes(x=Slope,fill=data_type))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.7)+
  scale_fill_manual(values = c("Observed"="#440154FF"))+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  theme_classic()
demo_hist <- demo_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
demo_hist




