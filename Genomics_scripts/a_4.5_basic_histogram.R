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

mat_all <- read_csv("Genomics_scripts/Data/freq_MAT_slope_peakbf5.csv")
map_all <- read_csv("Genomics_scripts/Data/freq_MAP_slope_peakbf5.csv")
cmd_all <- read_csv("Genomics_scripts/Data/freq_CMD_slope_peakbf5.csv")

mat_all[,4] <- "MAT"
map_all[,4] <- "MAP"
cmd_all[,4] <- "CMD"
colnames(mat_all)[4] <- "env"
colnames(map_all)[4] <- "env"
colnames(cmd_all)[4] <- "env"
env_all <- rbind(mat_all,map_all,cmd_all)
env_1 <- env_all %>% filter(Site==1)
env_2 <- env_all %>% filter(Site==2)
env_3 <- env_all %>% filter(Site==3)
env_4 <- env_all %>% filter(Site==4)
env_5 <- env_all %>% filter(Site==5)
env_6 <- env_all %>% filter(Site==6)
env_7 <- env_all %>% filter(Site==7)
env_8 <- env_all %>% filter(Site==8)
env_9 <- env_all %>% filter(Site==9)
env_10 <- env_all %>% filter(Site==10)
env_11 <- env_all %>% filter(Site==11)
env_12 <- env_all %>% filter(Site==12)

###################################################################################
## All population full graphs
###################################################################################

#mat 
mat_hist <- ggplot(mat_all,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (MAT)")+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~Site)
mat_hist
ggsave("Graphs_Oct_22/full_mat_overlap.pdf",width=10, height = 7.5, units = "in")

#map 
map_hist <- ggplot(map_all,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (map)")+
  theme_classic()
map_hist <- map_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
map_hist <- map_hist + facet_wrap(.~Site)
map_hist
ggsave("Graphs_Oct_22/full_map_overlap.pdf",width=10, height = 7.5, units = "in")

#cmd
cmd_hist <- ggplot(cmd_all,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (cmd)")+
  theme_classic()
cmd_hist <- cmd_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
cmd_hist <- cmd_hist + facet_wrap(.~Site)
cmd_hist
ggsave("Graphs_Oct_22/full_cmd_overlap.pdf",width=10, height = 7.5, units = "in")



###################################################################################
## Full graphs per population
###################################################################################

#p1 
mat_hist <- ggplot(env_1,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 1)")+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_50_p1.pdf",width=10, height = 7.5, units = "in")

#p2 
mat_hist <- ggplot(env_2,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 2)")+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_50_p2.pdf",width=10, height = 7.5, units = "in")

#p3 
mat_hist <- ggplot(env_3,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 3)")+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_50_p3.pdf",width=10, height = 7.5, units = "in")

#p4 
mat_hist <- ggplot(env_4,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 4)")+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_50_p4.pdf",width=10, height = 7.5, units = "in")

#p5 
mat_hist <- ggplot(env_5,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 5)")+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_50_p5.pdf",width=10, height = 7.5, units = "in")

#p6 
mat_hist <- ggplot(env_6,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 6)")+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_50_p6.pdf",width=10, height = 7.5, units = "in")

#p7 
mat_hist <- ggplot(env_7,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 7)")+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_50_p7.pdf",width=10, height = 7.5, units = "in")

#p8 
mat_hist <- ggplot(env_8,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 8)")+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_50_p8.pdf",width=10, height = 7.5, units = "in")

#p9 
mat_hist <- ggplot(env_9,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 9)")+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_50_p9.pdf",width=10, height = 7.5, units = "in")

#p10 
mat_hist <- ggplot(env_10,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 10)")+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_50_p10.pdf",width=10, height = 7.5, units = "in")

#p11 
mat_hist <- ggplot(env_11,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 11)")+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_50_p11.pdf",width=10, height = 7.5, units = "in")

#p12 
mat_hist <- ggplot(env_12,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 12)")+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_50_p12.pdf",width=10, height = 7.5, units = "in")


###################################################################################
## Full graphs per population, -5 to 5
###################################################################################

#p1 
mat_hist <- ggplot(env_1,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 1)",limits=c(-5,5))+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_5_p1.pdf",width=10, height = 7.5, units = "in")

#p2 
mat_hist <- ggplot(env_2,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 2)",limits=c(-5,5))+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_5_p2.pdf",width=10, height = 7.5, units = "in")

#p3 
mat_hist <- ggplot(env_3,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 3)",limits=c(-5,5))+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_5_p3.pdf",width=10, height = 7.5, units = "in")

#p4 
mat_hist <- ggplot(env_4,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 4)",limits=c(-5,5))+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_5_p4.pdf",width=10, height = 7.5, units = "in")

#p5 
mat_hist <- ggplot(env_5,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 5)",limits=c(-5,5))+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_5_p5.pdf",width=10, height = 7.5, units = "in")

#p6 
mat_hist <- ggplot(env_6,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 6)",limits=c(-5,5))+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_5_p6.pdf",width=10, height = 7.5, units = "in")

#p7 
mat_hist <- ggplot(env_7,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 7)",limits=c(-5,5))+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_5_p7.pdf",width=10, height = 7.5, units = "in")

#p8 
mat_hist <- ggplot(env_8,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 8)",limits=c(-5,5))+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_5_p8.pdf",width=10, height = 7.5, units = "in")

#p9 
mat_hist <- ggplot(env_9,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 9)",limits=c(-5,5))+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_5_p9.pdf",width=10, height = 7.5, units = "in")

#p10 
mat_hist <- ggplot(env_10,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 10)",limits=c(-5,5))+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_5_p10.pdf",width=10, height = 7.5, units = "in")

#p11 
mat_hist <- ggplot(env_11,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 11)",limits=c(-5,5))+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_5_p11.pdf",width=10, height = 7.5, units = "in")

#p12 
mat_hist <- ggplot(env_12,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 12)",limits=c(-5,5))+
  theme_classic()
mat_hist <- mat_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
mat_hist <- mat_hist + facet_wrap(.~env)
mat_hist
ggsave("Graphs_Oct_22/obs_5_p12.pdf",width=10, height = 7.5, units = "in")





