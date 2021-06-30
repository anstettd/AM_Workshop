##################################################################################
## Calculating Mean snp and indel coverage
## Author Daniel Anstett
## 
##
## Last Modified May 6, 2021
###################################################################################
#Import libraries
library(tidyverse)
###################################################################################

snp_1 <- read.table("Genomics_scripts/Data/cov_snp2_chr1.tsv", sep = "\t", fill = TRUE , header = TRUE)
snp_2 <- read.table("Genomics_scripts/Data/cov_snp2_chr2.tsv", sep = "\t", fill = TRUE , header = TRUE)
snp_3 <- read.table("Genomics_scripts/Data/cov_snp2_chr3.tsv", sep = "\t", fill = TRUE , header = TRUE)
snp_4 <- read.table("Genomics_scripts/Data/cov_snp2_chr4.tsv", sep = "\t", fill = TRUE , header = TRUE)
snp_5 <- read.table("Genomics_scripts/Data/cov_snp2_chr5.tsv", sep = "\t", fill = TRUE , header = TRUE)
snp_6 <- read.table("Genomics_scripts/Data/cov_snp2_chr6.tsv", sep = "\t", fill = TRUE , header = TRUE)
snp_7 <- read.table("Genomics_scripts/Data/cov_snp2_chr7.tsv", sep = "\t", fill = TRUE , header = TRUE)
snp_8 <- read.table("Genomics_scripts/Data/cov_snp2_chr8.tsv", sep = "\t", fill = TRUE , header = TRUE)

indel_1 <- read.table("Genomics_scripts/Data/cov_indel2_chr1.tsv", sep = "\t", fill = TRUE , header = TRUE)
indel_2 <- read.table("Genomics_scripts/Data/cov_indel2_chr2.tsv", sep = "\t", fill = TRUE , header = TRUE)
indel_3 <- read.table("Genomics_scripts/Data/cov_indel2_chr3.tsv", sep = "\t", fill = TRUE , header = TRUE)
indel_4 <- read.table("Genomics_scripts/Data/cov_indel2_chr4.tsv", sep = "\t", fill = TRUE , header = TRUE)
indel_5 <- read.table("Genomics_scripts/Data/cov_indel2_chr5.tsv", sep = "\t", fill = TRUE , header = TRUE)
indel_6 <- read.table("Genomics_scripts/Data/cov_indel2_chr6.tsv", sep = "\t", fill = TRUE , header = TRUE)
indel_7 <- read.table("Genomics_scripts/Data/cov_indel2_chr7.tsv", sep = "\t", fill = TRUE , header = TRUE)
indel_8 <- read.table("Genomics_scripts/Data/cov_indel2_chr8.tsv", sep = "\t", fill = TRUE , header = TRUE)

snp_cov_summary <- data.frame()
indel_cov_summary <- data.frame()

#rbind
snp<-rbind(snp_1,snp_2,snp_3,snp_4,snp_5,snp_6,snp_7,snp_8)
indel<-rbind(indel_1,indel_2,indel_3,indel_4,indel_5,indel_6,indel_7,indel_8)

#calculating mean and stdev of DP pre chromsome
snp_cov_mean <- snp %>% group_by(CHROM) %>% summarise(avg=mean(DP)) %>% ungroup()
snp_cov_stdev <- snp %>% group_by(CHROM) %>% summarise(stdev=sd(DP)) %>% ungroup()
snp_cov_summary <-left_join(snp_cov_mean,snp_cov_stdev,by="CHROM")
snp_cov_summary <- snp_cov_summary  %>% mutate(d1=avg+stdev,d2=avg-stdev) 
snp_cov_summary <- snp_cov_summary %>% select(CHROM,d1,d2)

indel_cov_mean <- indel %>% group_by(CHROM) %>% summarise(avg=mean(DP)) %>% ungroup()
indel_cov_stdev <- indel %>% group_by(CHROM) %>% summarise(stdev=sd(DP)) %>% ungroup()
indel_cov_summary <-left_join(indel_cov_mean,indel_cov_stdev,by="CHROM")
indel_cov_summary <- indel_cov_summary  %>% mutate(d1=avg+stdev,d2=avg-stdev) 
indel_cov_summary <- indel_cov_summary %>% select(CHROM,d1,d2)

#Export
write_csv(snp_cov_summary,"Genomics_scripts/Data/snp_cov2_summary.csv")
write_csv(indel_cov_summary,"Genomics_scripts/Data/indel_cov2_summary.csv")







