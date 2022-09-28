##################################################################################
## Prep work for generate SNP proportion matrix per region
## Filter timeseries by SNPS from peak window SNPS (WZA) with BF (BayPass) > 2
## Author Daniel Anstett
## 
## Currently done just for ENV 1,2 and 5 (MAT, MAP and CMD)
## Last Modified September 15, 2022
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

#Import Peak window, BF>2 SNP IDs
win_bf_mat  <- read_csv("Genomics_scripts/Data/win_bf_mat.csv")
win_bf_map  <- read_csv("Genomics_scripts/Data/win_bf_map.csv")
win_bf_cmd  <- read_csv("Genomics_scripts/Data/win_bf_cmd.csv")

#Import full snp table for timeseries
pop_order_time<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp_time<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci_time<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
colnames(loci_time) <- c("Chromosome","SNP")
loci_united_time <- loci_time %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp_time <-cbind(loci_united_time,snp_time) #add snp lables to rows


###################################################################################
#Import Baseline data
#import full snp table
pop_order_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", 
                           header=F, sep="\t")
snp_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", 
                     header=F, sep=" ")
loci_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", 
                      header=F, sep="\t")
colnames(loci_base) <- c("Chromosome","SNP")
loci_united_base <- loci_base %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp_base <-cbind(loci_united_base,snp_base) #add snp lables to rows


#Filter timeseries by peak window, BF>2 SNP IDs
peak_bf2_time_mat <-loci_snp_time %>% filter (chr_snp %in% as.character(win_bf_mat$chr_snp))
peak_bf2_time_map <-loci_snp_time %>% filter (chr_snp %in% as.character(win_bf_map$chr_snp))
peak_bf2_time_cmd <-loci_snp_time %>% filter (chr_snp %in% as.character(win_bf_cmd$chr_snp))


#Filter baseline by peak window, BF>2 SNP IDs
peak_bf2_base_mat <-loci_snp_base %>% filter (chr_snp %in% as.character(win_bf_mat$chr_snp))
peak_bf2_base_map <-loci_snp_base %>% filter (chr_snp %in% as.character(win_bf_map$chr_snp))
peak_bf2_base_cmd <-loci_snp_base %>% filter (chr_snp %in% as.character(win_bf_cmd$chr_snp))


write_csv(peak_bf2_time_mat,"Genomics_scripts/Data/peak_bf2_time_mat.csv")
write_csv(peak_bf2_time_map,"Genomics_scripts/Data/peak_bf2_time_map.csv")
write_csv(peak_bf2_time_cmd,"Genomics_scripts/Data/peak_bf2_time_cmd.csv")

write_csv(peak_bf2_base_mat,"Genomics_scripts/Data/peak_bf2_base_mat.csv")
write_csv(peak_bf2_base_map,"Genomics_scripts/Data/peak_bf2_base_map.csv")
write_csv(peak_bf2_base_cmd,"Genomics_scripts/Data/peak_bf2_base_cmd.csv")







