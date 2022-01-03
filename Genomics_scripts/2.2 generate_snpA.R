##################################################################################
## Generate SNP A proportion matrix
## Author Daniel Anstett
## 
## Currently done just for ENV 1,2 and 5 (MAT, MAP and CMD)
## Last Modified Nov 15, 2021
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

#Import BF>10 SNP data
snp1_filter <- read_csv("Genomics_scripts/Data/snp1_filter.csv")
snp2_filter <- read_csv("Genomics_scripts/Data/snp2_filter.csv")
snp5_filter <- read_csv("Genomics_scripts/Data/snp5_filter.csv")

snp1_random <- read_csv("Genomics_scripts/Data/snp1_random.csv")
snp2_random <- read_csv("Genomics_scripts/Data/snp2_random.csv")
snp5_random <- read_csv("Genomics_scripts/Data/snp5_random.csv")

pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")

###################################################################################
##Generate SNP A proportion matrix

#snp1_filter for SNP A
snp1A<-data.frame()
counter<-1
for (i in seq (2,dim(snp1_filter)[2]-1,2)){
  for(j in 1:dim(snp1_filter)[1]){
    tmp_total<-as.numeric(snp1_filter[j,i]) + as.numeric(snp1_filter[j,i+1])
    snp1A[j,counter]<-as.numeric(snp1_filter[j,i])/tmp_total
  }
  counter<-counter+1
}
colnames(snp1A)<- pop_order[,1] #name each pop/time combination
rownames(snp1A)<- snp1_filter$chr_snp

#snp2_filter for SNP A
snp2A<-data.frame()
counter<-1
for (i in seq (2,dim(snp2_filter)[2]-1,2)){
  for(j in 1:dim(snp2_filter)[1]){
    tmp_total<-as.numeric(snp2_filter[j,i]) + as.numeric(snp2_filter[j,i+1])
    snp2A[j,counter]<-as.numeric(snp2_filter[j,i])/tmp_total
  }
  counter<-counter+1
}
colnames(snp2A)<- pop_order[,1] #name each pop/time combination
rownames(snp2A)<- snp2_filter$chr_snp

#snp5_filter for SNP A
snp5A<-data.frame()
counter<-1
for (i in seq (2,dim(snp5_filter)[2]-1,2)){
  for(j in 1:dim(snp5_filter)[1]){
    tmp_total<-as.numeric(snp5_filter[j,i]) + as.numeric(snp5_filter[j,i+1])
    snp5A[j,counter]<-as.numeric(snp5_filter[j,i])/tmp_total
  }
  counter<-counter+1
}
colnames(snp5A)<- pop_order[,1] #name each pop/time combination
rownames(snp5A)<- snp5_filter$chr_snp


###################################################################################
##Generate SNP A random proportion matrix

#snp1_random for SNP A
snp1A_random<-data.frame()
counter<-1
for (i in seq (2,dim(snp1_random)[2]-1,2)){
  for(j in 1:dim(snp1_random)[1]){
    tmp_total<-as.numeric(snp1_random[j,i]) + as.numeric(snp1_random[j,i+1])
    snp1A_random[j,counter]<-as.numeric(snp1_random[j,i])/tmp_total
  }
  counter<-counter+1
}
colnames(snp1A_random)<- pop_order[,1] #name each pop/time combination
rownames(snp1A_random)<- snp1_random$chr_snp

#snp2_random for SNP A
snp2A_random<-data.frame()
counter<-1
for (i in seq (2,dim(snp2_random)[2]-1,2)){
  for(j in 1:dim(snp2_random)[1]){
    tmp_total<-as.numeric(snp2_random[j,i]) + as.numeric(snp2_random[j,i+1])
    snp2A_random[j,counter]<-as.numeric(snp2_random[j,i])/tmp_total
  }
  counter<-counter+1
}
colnames(snp2A_random)<- pop_order[,1] #name each pop/time combination
rownames(snp2A_random)<- snp2_random$chr_snp

#snp5_random for SNP A
snp5A_random<-data.frame()
counter<-1
for (i in seq (2,dim(snp5_random)[2]-1,2)){
  for(j in 1:dim(snp5_random)[1]){
    tmp_total<-as.numeric(snp5_random[j,i]) + as.numeric(snp5_random[j,i+1])
    snp5A_random[j,counter]<-as.numeric(snp5_random[j,i])/tmp_total
  }
  counter<-counter+1
}
colnames(snp5A_random)<- pop_order[,1] #name each pop/time combination
rownames(snp5A_random)<- snp5_random$chr_snp

#Transpose and separate site_year
snp1A_T <- as.data.frame(t(snp1A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
snp2A_T <- as.data.frame(t(snp2A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
snp5A_T <- as.data.frame(t(snp5A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))

snp1A_random_T <- as.data.frame(t(snp1A_random)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
snp2A_random_T <- as.data.frame(t(snp2A_random)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
snp5A_random_T <- as.data.frame(t(snp5A_random)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))

#Export Files
write_csv(snp1A_T, "Genomics_scripts/Data/snp1A.csv")
write_csv(snp2A_T, "Genomics_scripts/Data/snp2A.csv")
write_csv(snp5A_T, "Genomics_scripts/Data/snp5A.csv")

write_csv(snp1A_random_T, "Genomics_scripts/Data/snp1A_random.csv")
write_csv(snp2A_random_T, "Genomics_scripts/Data/snp2A_random.csv")
write_csv(snp5A_random_T, "Genomics_scripts/Data/snp5A_random.csv")







###################################################################################
#Test Case
snp_props<-data.frame()
snp_props_a<-data.frame()
snp_props_b<-data.frame()
counter<-1
for (i in seq (2,dim(snp1_filter)[2]-1,2)){
  for(j in 1:dim(snp1_filter)[1]){
    tmp_total<-as.numeric(snp1_filter[j,i]) + as.numeric(snp1_filter[j,i+1])
    #print(tmp_total)
    snp_props[j,i-1]<-as.numeric(snp1_filter[j,i])/tmp_total
    snp_props[j,i]<-as.numeric(snp1_filter[j,i+1])/tmp_total
    
    snp_props_a[j,counter]<-as.numeric(snp1_filter[j,i])/tmp_total
    snp_props_b[j,counter]<-as.numeric(snp1_filter[j,i+1])/tmp_total
  }
  counter<-counter+1
}
   