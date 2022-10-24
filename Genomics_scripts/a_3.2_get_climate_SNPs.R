##################################################################################
## Gradient forest data prep
## Author Daniel Anstett
## Filter climate SNPs for baseline preBaypass SNP table 
## Also merge data from 3 climate SNPs into one table
## Generate snpA table for both Bay and full data sets and merge with climate data
## 
## Last Modified September 19, 2022
###################################################################################
##Libraries
library(tidyverse)

##Import Data
climate <- read_csv("Donor_selection/Data/climate_pop.csv")

#Import peak bf>2 SNP data for each environmental variable
env1_peakbf2 <- read_csv("Genomics_scripts/Data/win_bf_mat30_5.csv")
env2_peakbf2 <- read_csv("Genomics_scripts/Data/win_bf_map30_5.csv")
env5_peakbf2 <- read_csv("Genomics_scripts/Data/win_bf_cmd30_5.csv")

#Merge into 1 dataframe that includes every unique allele
env_merge <- rbind(env1_peakbf2,env2_peakbf2,env5_peakbf2)
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
snp_clim_bayNA <- cbind(climate,snp_bay_T)
write_csv(snp_clim_bayNA, "Genomics_scripts/Data/snp_clim_peakbf5NA.csv")

#Remove NA's
snp_clim_bay_noNA <- snp_clim_bayNA %>% select_if(~ !any(is.na(.)))
write_csv(snp_clim_bay_noNA, "Genomics_scripts/Data/snp_clim_peakbf5_noNA.csv")


#Generate snpA table for full baseline snp datatset #Does not run locally
#snp_full<-data.frame()
#counter<-1
#for (i in seq (2,dim(loci_snp)[2]-1,2)){
#  for(j in 1:dim(loci_snp)[1]){
#    tmp_total<-as.numeric(loci_snp[j,i]) + as.numeric(loci_snp[j,i+1])
#    snp_full[j,counter]<-as.numeric(loci_snp[j,i])/tmp_total
#  }
#  counter<-counter+1
#}
#rownames(snp_full)<- loci_snp$chr_snp
#snp_full_T <- as.data.frame(t(snp_full))
#snp_clim_full <- cbind(climate,snp_full_T)
#write_csv(snp_clim_full, "Genomics_scripts/Data/snp_clim_full.csv")


