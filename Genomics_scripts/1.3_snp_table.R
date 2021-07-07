##################################################################################
## Set up snp table
## Author Daniel Anstett & Julia Anstett
## 
##
## Last Modified July 6, 2021
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

#import full snp table
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", 
                      header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", 
                      header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", 
                 header=F, sep="\t")
colnames(loci) <- c("Chromosome","SNP")
loci_united <- loci %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp <-cbind(loci_united,snp) #add snp lables to rows

#Import snp env associations
env1 <- read_csv("Genomics_scripts/Data/env1_adapt.csv")
env2 <- read_csv("Genomics_scripts/Data/env2_adapt.csv")
env3 <- read_csv("Genomics_scripts/Data/env3_adapt.csv")
env4 <- read_csv("Genomics_scripts/Data/env4_adapt.csv")
env5 <- read_csv("Genomics_scripts/Data/env5_adapt.csv")
env6 <- read_csv("Genomics_scripts/Data/env6_adapt.csv")
env7 <- read_csv("Genomics_scripts/Data/env7_adapt.csv")
env8 <- read_csv("Genomics_scripts/Data/env8_adapt.csv")
env9 <- read_csv("Genomics_scripts/Data/env9_adapt.csv")

#left_join snp table to each environemntal varaible to get loci data for each population
env1_loci <- left_join(env1,loci_snp, by="chr_snp")
env2_loci <- left_join(env2,loci_snp, by="chr_snp")
env3_loci <- left_join(env3,loci_snp, by="chr_snp")
env4_loci <- left_join(env4,loci_snp, by="chr_snp")
env5_loci <- left_join(env5,loci_snp, by="chr_snp")
env6_loci <- left_join(env6,loci_snp, by="chr_snp")
env7_loci <- left_join(env7,loci_snp, by="chr_snp")
env8_loci <- left_join(env8,loci_snp, by="chr_snp")
env9_loci <- left_join(env9,loci_snp, by="chr_snp")

env_all <- rbind(env1_loci,env2_loci,env3_loci,env4_loci,env5_loci,
                 env6_loci,env7_loci,env8_loci,env9_loci)
length(unique(env_all$chr_snp))

#Write out snp table
write_csv(env1_loci, "Genomics_scripts/Data/env1_loci.csv",col_names=TRUE)
write_csv(env2_loci, "Genomics_scripts/Data/env2_loci.csv",col_names=TRUE)
write_csv(env3_loci, "Genomics_scripts/Data/env3_loci.csv",col_names=TRUE)
write_csv(env4_loci, "Genomics_scripts/Data/env4_loci.csv",col_names=TRUE)
write_csv(env5_loci, "Genomics_scripts/Data/env5_loci.csv",col_names=TRUE)
write_csv(env6_loci, "Genomics_scripts/Data/env6_loci.csv",col_names=TRUE)
write_csv(env7_loci, "Genomics_scripts/Data/env7_loci.csv",col_names=TRUE)
write_csv(env8_loci, "Genomics_scripts/Data/env8_loci.csv",col_names=TRUE)
write_csv(env9_loci, "Genomics_scripts/Data/env9_loci.csv",col_names=TRUE)













