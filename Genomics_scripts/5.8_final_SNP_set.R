#############################################################################################################
## Generate Final SNP Set
## BF > 30
## Peak Boferroni window, BF >10
## 
## 
## Last Modified Sept 28, 2022
#############################################################################################################
#Import libraries
library(tidyverse)

#Import Large Loci Win
loci_win <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/loci_win.csv")                                                                                 
loci_win <- loci_win %>% unite(col="chr_snp", c("chr","snp"), sep="_")

#Import files
snps_peak_mat <- read_csv("Genomics_scripts/Data/snps_peak_mat.csv")
snps_peak_map <- read_csv("Genomics_scripts/Data/snps_peak_map.csv")
snps_peak_cmd <- read_csv("Genomics_scripts/Data/snps_peak_cmd.csv")

#Import snp env associations (Baseline)
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

#############################################################################################################
#Filter baseline by BF=>30

##Filter Bayes factor by peak windows
env1_united_bf30 <- env1_united %>% filter(BF>30)
env2_united_bf30 <- env2_united %>% filter(BF>30)
env5_united_bf30 <- env5_united %>% filter(BF>30)


#Filter Bayes factor by peak windows
snp_mat_peakbf <- env1_united %>% filter(BF<=30) %>% filter(chr_snp %in% as.character(snps_peak_mat$chr_snp))
snp_map_peakbf <- env2_united %>% filter(BF<=30) %>% filter(chr_snp %in% as.character(snps_peak_map$chr_snp))
snp_cmd_peakbf <- env5_united %>% filter(BF<=30) %>% filter(chr_snp %in% as.character(snps_peak_cmd$chr_snp))

#Filter peak windows by BF >10
snp_mat_peakbf10 <- snp_mat_peakbf %>% filter(BF>5)
snp_map_peakbf10 <- snp_map_peakbf %>% filter(BF>5)
snp_cmd_peakbf10 <- snp_cmd_peakbf %>% filter(BF>5)

#Merge peakbf10 with bf30
snp_mat_peakbf30 <- rbind(snp_mat_peakbf10,env1_united_bf30)
snp_map_peakbf30 <- rbind(snp_map_peakbf10,env2_united_bf30)
snp_cmd_peakbf30 <- rbind(snp_cmd_peakbf10,env5_united_bf30)

#Get window IDs for all selected SNPs
snp_mat_peakbf30_win <- left_join(snp_mat_peakbf30,loci_win, by="chr_snp")
snp_map_peakbf30_win <- left_join(snp_map_peakbf30,loci_win, by="chr_snp")
snp_cmd_peakbf30_win <- left_join(snp_cmd_peakbf30,loci_win, by="chr_snp")

colnames(snp_mat_peakbf30_win)[5] <- "win"
colnames(snp_map_peakbf30_win)[5] <- "win"
colnames(snp_cmd_peakbf30_win)[5] <- "win"


write_csv(snp_mat_peakbf30_win,"Genomics_scripts/Data/win_bf_mat30_5.csv")
write_csv(snp_map_peakbf30_win,"Genomics_scripts/Data/win_bf_map30_5.csv")
write_csv(snp_cmd_peakbf30_win,"Genomics_scripts/Data/win_bf_cmd30_5.csv")


