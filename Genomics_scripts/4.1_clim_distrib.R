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
# env 3 is PAS = Precipitation as snow (mm) between August in previous year and July in current year
# env 4 is EXT = Extreme temperature over 30 years
# env 5 is CMD = Hargreaves climatic moisture deficit (mm)
##### Seasonal #####
# env 6 is Tave_wt = Winter mean temperature (°C)
# env 7 is Tave_sm = Summer mean temperature (°C)
# env 8 is PPT_wt = Winter precipitation (mm)
# env 9 is PPT_sm = Summer precipitation (mm)

###################################################################################
#Import libraries
library(tidyverse)

#Import BF>10 SNP data
env1 <- read_csv("Genomics_scripts/Data/env1_adapt.csv")
env2 <- read_csv("Genomics_scripts/Data/env2_adapt.csv")
env3 <- read_csv("Genomics_scripts/Data/env3_adapt.csv")
env4 <- read_csv("Genomics_scripts/Data/env4_adapt.csv")
env5 <- read_csv("Genomics_scripts/Data/env5_adapt.csv")
env6 <- read_csv("Genomics_scripts/Data/env6_adapt.csv")
env7 <- read_csv("Genomics_scripts/Data/env7_adapt.csv")
env8 <- read_csv("Genomics_scripts/Data/env8_adapt.csv")
env9 <- read_csv("Genomics_scripts/Data/env9_adapt.csv")

#Add identifier
env_all <- rbind(env1,env2,env3,env4,env5,env6,env7,env8,env9)
env_chr1 <- env_all %>% filter(Chromosome =="CE10_chr1")
env_chr1_20 <- env_chr1 %>% filter(BF>20)

#Chromosome 1 BF>10
chr1_graph <-ggplot(env_chr1, aes(x=SNP, y=Env)) +
  geom_point(size=1, shape=4,color="black") +
  scale_y_continuous(breaks = seq(1:9),labels=c("1" = "MAT", "2" = "MAP","3" = "PAS","4" = "EXT","5" = "CMD",
                            "6" = "T_winter","7" = "T_summer","8" = "Precip_Winter","9" = "Precip_summer"))+
  xlab("Chromosome 1 Position")+ theme_classic() + theme(axis.title.y=element_blank())
chr1_graph

#Chromosome 1 BF>20
chr1_graph_20 <-ggplot(env_chr1_20, aes(x=SNP, y=Env)) +
  geom_point(size=1, shape=4,color="black") +
  scale_y_continuous(breaks = seq(1:9),labels=c("1" = "MAT", "2" = "MAP","3" = "PAS","4" = "EXT","5" = "CMD",
                                                "6" = "T_winter","7" = "T_summer","8" = "Precip_Winter","9" = "Precip_summer"))+
  xlab("Chromosome 1 Position")+ theme_classic() + theme(axis.title.y=element_blank())
chr1_graph_20





