##################################################################################
## Filter out BF>0 SNPs
## Author Daniel Anstett & Julia Anstett
## 
##
## Last Modified July 6, 2021
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)


#Import snp env associations
env1 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_1_trim.tsv",header=F, sep=" ")
env2 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_2_trim.tsv",header=F, sep=" ")
env5 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_5_trim.tsv",header=F, sep=" ")


#Name Columns
colnames(env1) <- c("Chromosome","SNP","Env","BF")
colnames(env2) <- c("Chromosome","SNP","Env","BF")
colnames(env5) <- c("Chromosome","SNP","Env","BF")

#Filter Bayes factor > |10|
env1_filter <- env1 %>% filter(BF >= 0)
env2_filter <- env2 %>% filter(BF >= 0)
env5_filter <- env5 %>% filter(BF >= 0)


#Make Chromosome SNP variable for easier left_joining later
env1_united <- env1_filter %>% unite(chr_snp,Chromosome,SNP)
env1_united <- env1_united %>% select(chr_snp)
env1_filter_united <- cbind(env1_united,env1_filter)

env2_united <- env2_filter %>% unite(chr_snp,Chromosome,SNP)
env2_united <- env2_united %>% select(chr_snp)
env2_filter_united <- cbind(env2_united,env2_filter)

env5_united <- env5_filter %>% unite(chr_snp,Chromosome,SNP)
env5_united <- env5_united %>% select(chr_snp)
env5_filter_united <- cbind(env5_united,env5_filter)


#Export "adaptive" snps
write_csv(env1_filter_united, "Genomics_scripts/Data/env1_BF0.csv",col_names=TRUE)
write_csv(env2_filter_united, "Genomics_scripts/Data/env2_BF0.csv",col_names=TRUE)
write_csv(env5_filter_united, "Genomics_scripts/Data/env5_BF0.csv",col_names=TRUE)



