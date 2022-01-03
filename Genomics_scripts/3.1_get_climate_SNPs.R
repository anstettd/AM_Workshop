##################################################################################
## Gradient forest data prep
## Author Daniel Anstett
## Filter climate SNPs for baseline preBaypass SNP table 
## Also merge data from 9 climate SNPs into one table
## Generate snpA table for both Bay and full data sets and merge with climate data
## 
## Last Modified Dec 13, 2021
###################################################################################
##Libraries
library(tidyverse)

##Import Data

climate <- read_csv("Donor_selection/Data/climate_pop.csv")

#Import BF>10 SNP data for environmental variable
env1 <- read_csv("Genomics_scripts/Data/env1_adapt.csv")
env2 <- read_csv("Genomics_scripts/Data/env2_adapt.csv")
env3 <- read_csv("Genomics_scripts/Data/env3_adapt.csv")
env4 <- read_csv("Genomics_scripts/Data/env4_adapt.csv")
env5 <- read_csv("Genomics_scripts/Data/env5_adapt.csv")
env6 <- read_csv("Genomics_scripts/Data/env6_adapt.csv")
env7 <- read_csv("Genomics_scripts/Data/env7_adapt.csv")
env8 <- read_csv("Genomics_scripts/Data/env8_adapt.csv")
env9 <- read_csv("Genomics_scripts/Data/env9_adapt.csv")

#Merge into 1 dataframe that includes every unique allele
env_merge <- rbind(env1,env2,env3,env4,env5,env6,env7,env8,env9)
env_all <- env_merge[,1]
env_snp <- unique(env_all)

#Import full snp table for baseline
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
colnames(loci) <- c("Chromosome","SNP")
loci_united <- loci %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp <-cbind(loci_united,snp) #add snp lables to rows

#Select BayPass SNPs from merged dataframe
loci_env <- loci_snp %>% filter (chr_snp %in% as.character(env_snp$chr_snp))

#Generate snpA table for BayPass SNPs
snp_bay<-data.frame()
counter<-1
for (i in seq (2,dim(loci_env)[2]-1,2)){
  for(j in 1:dim(loci_env)[1]){
    tmp_total<-as.numeric(loci_env[j,i]) + as.numeric(loci_env[j,i+1])
    snp_bay[j,counter]<-as.numeric(loci_env[j,i])/tmp_total
  }
  counter<-counter+1
}
rownames(snp_bay)<- loci_env$chr_snp
snp_bay_T <- as.data.frame(t(snp_bay))
snp_clim_bay <- cbind(climate,snp_bay_T)
write_csv(snp_clim_bay, "Genomics_scripts/Data/snp_clim_bay.csv")

#Remove NA's
snp_clim_bf10NA <- snp_clim_bay %>%
  select_if(~ !any(is.na(.)))
write_csv(snp_clim_bay, "Genomics_scripts/Data/snp_clim_bf10NA")


#Generate snpA table for full baseline snp datatset #Does not run locally
snp_full<-data.frame()
counter<-1
for (i in seq (2,dim(loci_snp)[2]-1,2)){
  for(j in 1:dim(loci_snp)[1]){
    tmp_total<-as.numeric(loci_snp[j,i]) + as.numeric(loci_snp[j,i+1])
    snp_full[j,counter]<-as.numeric(loci_snp[j,i])/tmp_total
  }
  counter<-counter+1
}
rownames(snp_full)<- loci_snp$chr_snp
snp_full_T <- as.data.frame(t(snp_full))
snp_clim_full <- cbind(climate,snp_full_T)
#write_csv(snp_clim_full, "Genomics_scripts/Data/snp_clim_full.csv")


