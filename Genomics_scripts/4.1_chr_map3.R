##################################################################################
## Climate associated SNP distribution across chromosomes
## Author Daniel Anstett
## 
## 
## Last Modified Marc 22, 2021
###################################################################################

##### Annual #####
# env 1 is MAT = Mean annual temperature (°C)
# env 2 is MAP = Mean annual precipitation (mm)
# env 5 is CMD = Hargreaves climatic moisture deficit (mm)


###################################################################################
#Import libraries
library(tidyverse)

#Import BF>10 SNP data
env1 <- read_csv("Genomics_scripts/Data/env1_adapt.csv")
env2 <- read_csv("Genomics_scripts/Data/env2_adapt.csv")
env5 <- read_csv("Genomics_scripts/Data/env5_adapt.csv")
env5$Env <- 3


#Add identifier
env_all <- rbind(env1,env2,env5)
env_20 <- env_all %>% filter(BF>20)

env_chr14 <- env_20 %>% filter(Chromosome =="CE10_chr1" | Chromosome =="CE10_chr2" |
                                 Chromosome =="CE10_chr3" | Chromosome =="CE10_chr4")
env_chr58 <- env_20 %>% filter(Chromosome =="CE10_chr5" | Chromosome =="CE10_chr6" |
                                 Chromosome =="CE10_chr7" | Chromosome =="CE10_chr8")

#Chromosome 1-4 BF>20
chr1_graph_14 <-ggplot(env_chr14, aes(x=SNP, y=Env)) +
  geom_point(size=2, shape=4,color="black") +
  scale_y_continuous(breaks = seq(1:3),labels=c("1" = "MAT", "2" = "MAP","3" = "CMD"))+
  scale_x_continuous(limits = c(0, 5.4e+07), breaks=c(0,1e+07,2e+07,3e+07,4e+07,5e+07))+
  xlab("Chromosome Position")+ theme_classic() + theme(axis.title.y=element_blank())+
  facet_wrap(~ Chromosome, ncol = 1)+
  #theme(strip.background = element_blank())
  theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5,))+
  theme(panel.spacing = unit(0.4, "lines"))
chr1_graph_14
ggsave("Graphs_overlap/chr14.pdf",width=10, height = 5, units = "in")

#Chromosome 5-8 BF>20
chr1_graph_58 <-ggplot(env_chr58, aes(x=SNP, y=Env)) +
  geom_point(size=2, shape=4,color="black") +
  scale_y_continuous(breaks = seq(1:3),labels=c("1" = "MAT", "2" = "MAP","3" = "CMD"))+
  scale_x_continuous(limits = c(0, 6.2e+07), breaks=c(0,1e+07,2e+07,3e+07,4e+07,5e+07,6e+07))+
  xlab("Chromosome Position")+ theme_classic() + theme(axis.title.y=element_blank())+
  facet_wrap(~ Chromosome, ncol = 1)+
#theme(strip.background = element_blank())
  theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))+
  theme(panel.spacing = unit(0.4, "lines"))
chr1_graph_58
ggsave("Graphs_overlap/chr58.pdf",width=10, height = 5, units = "in")



