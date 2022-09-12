#############################################################################################################
## Get SNPs with highest Bayes factor from peak WZA windows
## Author Daniel Anstett
## 
## 
## 
## Last Modified Sept 2, 2022
#############################################################################################################
#Import libraries
library(tidyverse)

#Import files
snps_peak_mat <- read_csv("Genomics_scripts/Data/snps_peak_mat.csv")
snps_peak_map <- read_csv("Genomics_scripts/Data/snps_peak_map.csv")
snps_peak_cmd <- read_csv("Genomics_scripts/Data/snps_peak_cmd.csv")


#Import snp env associations
env1 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_1_trim.tsv",header=F, sep=" ")
env2 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_2_trim.tsv",header=F, sep=" ")
env5 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_5_trim.tsv",header=F, sep=" ")

#Name Columns
colnames(env1) <- c("Chromosome","SNP","Env","BF")
colnames(env2) <- c("Chromosome","SNP","Env","BF")
colnames(env5) <- c("Chromosome","SNP","Env","BF")

env1_united <- env1 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env2_united <- env2 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env5_united <- env5 %>% unite(chr_snp,"Chromosome","SNP",sep="_")

#Filter Bayes factor by peak windows
peak_bf_mat <- env1_united %>% filter(chr_snp %in% as.character(snps_peak_mat$chr_snp))
peak_bf_map <- env2_united %>% filter(chr_snp %in% as.character(snps_peak_map$chr_snp))
peak_bf_cmd <- env5_united %>% filter(chr_snp %in% as.character(snps_peak_cmd$chr_snp))

peak_bf10_mat <- peak_bf_mat %>% filter(BF>=10)
peak_bf10_map <- peak_bf_map %>% filter(BF>=10)
peak_bf10_cmd <- peak_bf_cmd %>% filter(BF>=10)

write_csv(peak_bf10_mat,"Genomics_scripts/Data/win_bf_mat.csv")
write_csv(peak_bf10_map,"Genomics_scripts/Data/win_bf_map.csv")
write_csv(peak_bf10_cmd,"Genomics_scripts/Data/win_bf_cmd.csv")







