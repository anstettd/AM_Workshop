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

#Import population order
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
pop_order_2 <- data.frame()
pop_order_2 [1,1] <- "chr_shp"
pop_order_2 [2,1] <- "Site"
pop_order_2 <- rbind(pop_order_2,pop_order)



###################################################################################
# Make SNP A table
rand_mat_A <- data.frame()
counter<-1
for (i in seq (3,dim(rand_mat)[2]-1,2)){
  for(j in 1:dim(rand_mat)[1]){
    tmp_total<-as.numeric(rand_mat[j,i]) + as.numeric(rand_mat[j,i+1])
    rand_mat_A[j,counter]<-as.numeric(rand_mat[j,i])/tmp_total
  }
  counter<-counter+1
}
colnames(rand_mat_A)<- pop_order[,1] #name each pop/time combination
#rownames(rand_mat_A)<- rand_mat$chr_snp

rand_map_A <- data.frame()
counter<-1
for (i in seq (3,dim(rand_map)[2]-1,2)){
  for(j in 1:dim(rand_map)[1]){
    tmp_total<-as.numeric(rand_map[j,i]) + as.numeric(rand_map[j,i+1])
    rand_map_A[j,counter]<-as.numeric(rand_map[j,i])/tmp_total
  }
  counter<-counter+1
}
colnames(rand_map_A)<- pop_order[,1] #name each pop/time combination
#rownames(rand_map_A)<- rand_map$chr_snp

rand_cmd_A <- data.frame()
counter<-1
for (i in seq (3,dim(rand_cmd)[2]-1,2)){
  for(j in 1:dim(rand_cmd)[1]){
    tmp_total<-as.numeric(rand_cmd[j,i]) + as.numeric(rand_cmd[j,i+1])
    rand_cmd_A[j,counter]<-as.numeric(rand_cmd[j,i])/tmp_total
  }
  counter<-counter+1
}
colnames(rand_cmd_A)<- pop_order[,1] #name each pop/time combination
#rownames(rand_cmd_A)<- rand_cmd$chr_snp





