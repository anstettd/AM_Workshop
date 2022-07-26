##################################################################################
## Filter timeseries SNPs so that they only include the top SNP candiates from baseline (BF>20) 
## Author Daniel Anstett
## 
## Currently done just for ENV 1,2 and 5 (MAT, MAP and CMD)
## Last Modified Nov 15, 2021
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

##### Annual #####
# env 1 is MAT = Mean annual temperature (°C)
# env 2 is MAP = Mean annual precipitation (mm)
# env 5 is CMD = Hargreaves climatic moisture deficit (mm)


#Import BF>10 Baseline SNP data 
env1 <- read_csv("Genomics_scripts/Data/env1_adapt.csv")
env2 <- read_csv("Genomics_scripts/Data/env2_adapt.csv")
env5 <- read_csv("Genomics_scripts/Data/env5_adapt.csv")


#Filter to BF>20 data
env120 <- env1 %>% filter(BF>20)
env220 <- env2 %>% filter(BF>20)
env520 <- env5 %>% filter(BF>20)

#Import full snp table for timeseries
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
colnames(loci) <- c("Chromosome","SNP")
loci_united <- loci %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp <-cbind(loci_united,snp) #add snp lables to rows

#Select Timeseries SNPs with baseline SNPs at BF>20
snp1_filter <-loci_snp %>% filter (chr_snp %in% as.character(env120$chr_snp)) #MAT
snp2_filter <-loci_snp %>% filter (chr_snp %in% as.character(env220$chr_snp)) #MAP
snp5_filter <-loci_snp %>% filter (chr_snp %in% as.character(env520$chr_snp)) #CMD


write_csv(snp1_filter, "Genomics_scripts/Data/snp1_filter.csv")
write_csv(snp2_filter, "Genomics_scripts/Data/snp2_filter.csv")
write_csv(snp5_filter, "Genomics_scripts/Data/snp5_filter.csv")
