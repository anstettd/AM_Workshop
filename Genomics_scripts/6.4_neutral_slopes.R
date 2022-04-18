##################################################################################
## Generate stratified distribution of random non-climate associated slopes
## Author Daniel Anstett
## 
## 
## Last Modified April 18, 2022
###################################################################################
#Import libraries
library(tidyverse)

#Import SNP A & B tables
rand_mat <- read_csv("Genomics_scripts/Data/rand_mat.csv")
rand_map <- read_csv("Genomics_scripts/Data/rand_map.csv")
rand_cmd <- read_csv("Genomics_scripts/Data/rand_cmd.csv")

# Make SNP A table
rand_mat_A <- data.frame()
counter<-1
for (i in seq (2,dim(snp1_random)[2]-1,2)){
  for(j in 1:dim(snp1_random)[1]){
    tmp_total<-as.numeric(snp1_random[j,i]) + as.numeric(snp1_random[j,i+1])
    rand_mat_A[j,counter]<-as.numeric(snp1_random[j,i])/tmp_total
  }
  counter<-counter+1
}
colnames(snp1A_random)<- pop_order[,1] #name each pop/time combination
rownames(snp1A_random)<- snp1_random$chr_snp
