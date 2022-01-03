##################################################################################
## Filter Timeseries for BF >20
## Author Daniel Anstett
## 
## Currently done just for ENV 1,2 and 5 (MAT, MAP and CMD)
## Last Modified Nov 15, 2021
###################################################################################


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

#Filter to BF>20 data
env120 <- env1 %>% filter(BF>20)
env220 <- env2 %>% filter(BF>20)
env320 <- env3 %>% filter(BF>20)
env420 <- env4 %>% filter(BF>20)
env520 <- env5 %>% filter(BF>20)
env620 <- env6 %>% filter(BF>20)
env720 <- env7 %>% filter(BF>20)
env820 <- env8 %>% filter(BF>20)
env920 <- env9 %>% filter(BF>20)

#Import full snp table for timeseries
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
colnames(loci) <- c("Chromosome","SNP")
loci_united <- loci %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp <-cbind(loci_united,snp) #add snp lables to rows

#Select Timeseries SNPs with BF>20
snp1_filter <-loci_snp %>% filter (chr_snp %in% as.character(env120$chr_snp)) #MAT
snp2_filter <-loci_snp %>% filter (chr_snp %in% as.character(env320$chr_snp)) #MAP
snp5_filter <-loci_snp %>% filter (chr_snp %in% as.character(env520$chr_snp)) #CMD

#Randomly sample the same number of SNPs for each environmental variable
set.seed(1)
rand1 <- loci_snp[sample(nrow(loci_snp), dim(snp1_filter)[1]), ]
snp1_random <-loci_snp %>% filter (chr_snp %in% as.character(rand1$chr_snp)) #MAT

set.seed(2)
rand2 <- loci_snp[sample(nrow(loci_snp), dim(snp2_filter)[1]), ]
snp2_random <-loci_snp %>% filter (chr_snp %in% as.character(rand2$chr_snp)) #MAT

set.seed(5)
rand5 <- loci_snp[sample(nrow(loci_snp), dim(snp5_filter)[1]), ]
snp5_random <-loci_snp %>% filter (chr_snp %in% as.character(rand5$chr_snp)) #MAT



#Export Files
write_csv(snp1_filter, "Genomics_scripts/Data/snp1_filter.csv")
write_csv(snp2_filter, "Genomics_scripts/Data/snp2_filter.csv")
write_csv(snp5_filter, "Genomics_scripts/Data/snp5_filter.csv")

write_csv(snp1_random, "Genomics_scripts/Data/snp1_random.csv")
write_csv(snp2_random, "Genomics_scripts/Data/snp2_random.csv")
write_csv(snp5_random, "Genomics_scripts/Data/snp5_random.csv")









