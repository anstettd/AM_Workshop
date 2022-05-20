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

#Import BF>20 SNP data
snp1_filter <- read_csv("Genomics_scripts/Data/snp1_filter.csv")
snp2_filter <- read_csv("Genomics_scripts/Data/snp2_filter.csv")
snp5_filter <- read_csv("Genomics_scripts/Data/snp5_filter.csv")

snp1_random <- read_csv("Genomics_scripts/Data/snp1_random.csv")
snp2_random <- read_csv("Genomics_scripts/Data/snp2_random.csv")
snp5_random <- read_csv("Genomics_scripts/Data/snp5_random.csv")

#Add env variable identifier 
chr_obs_mat <- snp1_filter %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)
chr_obs_map <- snp2_filter %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)
chr_obs_cmd <- snp5_filter %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)

chr_rand_mat <- snp1_random %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)
chr_rand_map <- snp2_random %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)
chr_rand_cmd <- snp5_random %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)

#Add env identfier
chr_obs_mat$env <- "MAT"
chr_obs_map$env <- "MAP"
chr_obs_cmd$env <- "CMD"

chr_rand_mat$env <- "MAT"
chr_rand_map$env <- "MAP"
chr_rand_cmd$env <- "CMD"

#Add data type
chr_obs_mat$type <- "Climate Associated"
chr_obs_map$type <- "Climate Associated"
chr_obs_cmd$type <- "Climate Associated"

chr_rand_mat$type <- "Random"
chr_rand_map$type <- "Random"
chr_rand_cmd$type <- "Random"

#Join Climate Associated and Random
chr_mat <- rbind(chr_obs_mat,chr_rand_mat)
chr_map <- rbind(chr_obs_map,chr_rand_map)
chr_cmd <- rbind(chr_obs_cmd,chr_rand_cmd)


#############################################################################################
#Make plots
#MAT
chr1_graph_14 <-ggplot(chr_mat, aes(x=Position, y=type)) +
  geom_point(size=2, shape=4,color="black") +
  #  scale_y_continuous(breaks = seq(1:3),labels=c("1" = "MAT", "2" = "MAP","3" = "CMD"))+
  #  scale_x_continuous(limits = c(0, 5.4e+07), breaks=c(0,1e+07,2e+07,3e+07,4e+07,5e+07))+
  xlab("Chromosome Position")+ theme_classic() + theme(axis.title.y=element_blank())+
  facet_wrap(~ chr, ncol = 1)+
  #theme(strip.background = element_blank())
  theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5,))+
  theme(panel.spacing = unit(0.4, "lines"))
chr1_graph_14
ggsave("Graphs_CI/chr_mat.pdf",width=10, height = 6, units = "in")

#MAP
chr1_graph_14 <-ggplot(chr_map, aes(x=Position, y=type)) +
  geom_point(size=2, shape=4,color="black") +
  #  scale_y_continuous(breaks = seq(1:3),labels=c("1" = "MAT", "2" = "MAP","3" = "CMD"))+
  #  scale_x_continuous(limits = c(0, 5.4e+07), breaks=c(0,1e+07,2e+07,3e+07,4e+07,5e+07))+
  xlab("Chromosome Position")+ theme_classic() + theme(axis.title.y=element_blank())+
  facet_wrap(~ chr, ncol = 1)+
  #theme(strip.background = element_blank())
  theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5,))+
  theme(panel.spacing = unit(0.4, "lines"))
chr1_graph_14
ggsave("Graphs_CI/chr_map.pdf",width=10, height = 6, units = "in")

#CMD
chr1_graph_14 <-ggplot(chr_cmd, aes(x=Position, y=type)) +
  geom_point(size=2, shape=4,color="black") +
  #  scale_y_continuous(breaks = seq(1:3),labels=c("1" = "MAT", "2" = "MAP","3" = "CMD"))+
  #  scale_x_continuous(limits = c(0, 5.4e+07), breaks=c(0,1e+07,2e+07,3e+07,4e+07,5e+07))+
  xlab("Chromosome Position")+ theme_classic() + theme(axis.title.y=element_blank())+
  facet_wrap(~ chr, ncol = 1)+
  #theme(strip.background = element_blank())
  theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5,))+
  theme(panel.spacing = unit(0.4, "lines"))
chr1_graph_14
ggsave("Graphs_CI/chr_cmd.pdf",width=10, height = 6, units = "in")






