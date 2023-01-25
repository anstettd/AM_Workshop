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
snp1_filter <- read_csv("Genomics_scripts/Data/win_bf_mat30_5.csv")
snp2_filter <- read_csv("Genomics_scripts/Data/win_bf_map30_5.csv")
snp5_filter <- read_csv("Genomics_scripts/Data/win_bf_cmd30_5.csv")

#Add env variable identifier 
chr_obs_mat <- snp1_filter %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)
chr_obs_map <- snp2_filter %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)
chr_obs_cmd <- snp5_filter %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)

#Add env identfier
chr_obs_mat$env <- "MAT"
chr_obs_map$env <- "MAP"
chr_obs_cmd$env <- "CMD"


#Add data type
chr_obs_mat$type <- "Climate Associated"
chr_obs_map$type <- "Climate Associated"
chr_obs_cmd$type <- "Climate Associated"

#Join env variables
chr_all <- rbind(chr_obs_mat,chr_obs_map,chr_obs_cmd)
chr_all$Position <- as.integer(chr_all$Position)


#############################################################################################
#Plot all 
chr1_graph_14 <-ggplot(chr_all, aes(x=Position, y=env)) +
  geom_point(size=2, shape=4,color="black") +
  scale_x_continuous(name="Position",  breaks=c(0,1e+07,2e+07,3e+07,4e+07,5e+07,6e+07)) +
  theme_classic() + 
  #theme(axis.title.y=element_blank())+
  facet_wrap(~ chr, ncol = 1)+
  theme(strip.text.x = element_text(size = 11,face="bold"),
        axis.text.x = element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
        axis.title.y=element_blank())
chr1_graph_14
ggsave("Graphs_CI/chr_map_double_filter.pdf",width=8, height = 8, units = "in")




