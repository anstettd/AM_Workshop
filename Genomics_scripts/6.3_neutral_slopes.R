##################################################################################
## Generate stratified distribution of random non-climate associated slopes
## Author Daniel Anstett
## 
## 
## Last Modified April 15, 2022
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

#Import timeseries frequencies
freq_MAT <- read_csv("Genomics_scripts/Data/freq_MAT.csv")
freq_MAP <- read_csv("Genomics_scripts/Data/freq_MAP.csv")
freq_CMD <- read_csv("Genomics_scripts/Data/freq_CMD.csv")

#Import full snp table for timeseries
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
colnames(loci) <- c("Chromosome","SNP")
loci_united <- loci %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp <-cbind(loci_united,snp) #add snp lables to rows

#Filter full snp table to remove climate associated SNPs
snp_swiss <-loci_snp %>% filter (!chr_snp %in% as.character(freq_MAT$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(freq_MAP$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(freq_CMD$chr_snp))

#What ever comes next




