##################################################################################
## Tabulate proportion of allele associated with environmental varaible 
## Author Daniel Anstett & Julia Anstett
## 
##
## Last Modified July 7, 2021
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

#Import snp env associations
env1_loci <- read_csv("Genomics_scripts/Data/env1_loci.csv")
env2_loci <- read_csv("Genomics_scripts/Data/env2_loci.csv")
env3_loci <- read_csv("Genomics_scripts/Data/env3_loci.csv")
env4_loci <- read_csv("Genomics_scripts/Data/env4_loci.csv")
env5_loci <- read_csv("Genomics_scripts/Data/env5_loci.csv")
env6_loci <- read_csv("Genomics_scripts/Data/env6_loci.csv")
env7_loci <- read_csv("Genomics_scripts/Data/env7_loci.csv")
env8_loci <- read_csv("Genomics_scripts/Data/env8_loci.csv")
env9_loci <- read_csv("Genomics_scripts/Data/env9_loci.csv")

#Import Climate
climate <- read_csv("Donor_selection/Data/climate_pop.csv")

###################################################################################

#Generate frequency matrix for prop A & prop B

#Mat
snp_prop_A<-env1_loci %>% select (chr_snp, Chromosome, SNP, Env)
counter=5
pop_num=1
for (i in seq(1,110,2)) {
  
  snpA<-paste("V",i, sep="")
  snpB<-paste("V",i+1, sep="")
  P<-paste("P", pop_num, sep="")
  tmp<-env1_loci %>% select(snpA,snpB)
  
  colnames(tmp)<-c("A", "B")
  
   #tmp <-tmp %>% mutate(Prop= A/(A+B))
   snp_prop_A[,counter]<-tmp$A/(tmp$A + tmp$B)
   colnames (snp_prop_A)[counter]<-P
   
   counter<-counter+1
   pop_num<-pop_num+1
}
  
#snp_prop_A[is.na(snp_prop_A)] <- 0

#B Mat
snp_prop_B<-env1_loci %>% select (chr_snp, Chromosome, SNP, Env)
counter=5
pop_num=1

for (i in seq(1,110,2)) {
  
  snpA<-paste("V",i, sep="")
  snpB<-paste("V",i+1, sep="")
  P<-paste("P", pop_num, sep="")
  tmp<-env1_loci %>% select(snpA,snpB)
  
  colnames(tmp)<-c("A", "B")
  
  #tmp <-tmp %>% mutate(Prop= A/(A+B))
  snp_prop_B[,counter]<-tmp$B/(tmp$A + tmp$B)
  colnames (snp_prop_B)[counter]<-P
  
  counter<-counter+1
  pop_num<-pop_num+1
}

#snp_prop_B[is.na(snp_prop_B)] <- 0
###################################################################################

#set up frequency association table
freq.temp <- data.frame()

for (i in 1:dim(snp_prop_A)[1]) {
  
env_pop <- climate %>% select(Paper_ID,MAT) #Select climate data
snp_prop_A_tmp<-snp_prop_A %>% filter(SNP==snp_prop_A$SNP[i]) %>% 
  select(-chr_snp, -Chromosome, -SNP, -Env) #filter prop A data 
snp_prop_B_tmp<-snp_prop_B %>% filter(SNP==snp_prop_B$SNP[i]) %>%
  select(-chr_snp, -Chromosome, -SNP, -Env) #filter prop B data 

env_pop <- cbind(env_pop,t(snp_prop_A_tmp), t(snp_prop_B_tmp))
colnames(env_pop)[3] <- "prop_A"
colnames(env_pop)[4]<-"prop_B"

lm.temp_A <- lm(env_pop$prop_A~env_pop$MAT) # save lm of cliamte predicting prop A
lm.temp_B <- lm(env_pop$prop_B~env_pop$MAT) # save lm of cliamte predicting prop B


if (isTRUE(lm.temp_A$coefficents[2]>0) && isTRUE(lm.temp_B$coefficents[2]>0)){
  
  print(i)
  print(lm.temp_A) 
  print(lm.temp_B) 

  }else if(lm.temp_A$coefficients[2]>0){
  #print("A")
  tmp_in<-snp_prop_A %>% filter(SNP==snp_prop_A$SNP[i])
  freq.temp<-rbind(freq.temp,tmp_in)
  
  } else if (lm.temp_B$coefficients[2]>0){
  # print("B")
     tmp_in<-snp_prop_B %>% filter(SNP==snp_prop_B$SNP[i])
     freq.temp<-rbind(freq.temp,tmp_in)
     
  }else{
    print(i)
    print(lm.temp_A) 
    print(lm.temp_B) 
  }
}

###################################################################################


# Test for concordance betweeen signs of regression for frequency A and frequency B matrixes

env_pop <- climate %>% select(Paper_ID,MAT)

snp_prop_A_tmp<-snp_prop_A %>% filter(SNP==snp_prop_A$SNP[3661]) %>%
  select(-chr_snp, -Chromosome, -SNP, -Env)

snp_prop_B_tmp<-snp_prop_B %>% filter(SNP==snp_prop_B$SNP[3661]) %>%
  select(-chr_snp, -Chromosome, -SNP, -Env)

env_pop <- cbind(env_pop,t(snp_prop_A_tmp), t(snp_prop_B_tmp))
colnames(env_pop)[3] <- "prop_A"
colnames(env_pop)[4]<-"prop_B"

lm(env_pop$prop_A~env_pop$MAT)
lm(env_pop$prop_B~env_pop$MAT)

plot(x=env_pop$MAT,y=env_pop$prop_A)
plot(x=env_pop$MAT,y=env_pop$prop_B)







for (i in 1:dim(snp_prop_A)[1])
  snp_prop_A_tmp<-snp_prop_A %>% filter(SNP==snp_prop_A$SNP[i]) %>%
    select(-chr_snp, -Chromosome, -SNP, -Env)
  

env_pop


  



