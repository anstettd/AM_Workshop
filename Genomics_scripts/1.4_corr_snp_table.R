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
###################################################################################
##setup functions

#Generate frequency matrix for prop A 
prop_A <- function(snp_table) {
snp_prop_A<- snp_table %>% select (chr_snp, Chromosome, SNP, Env)
counter=5
pop_num=1
for (i in seq(1,110,2)) {
  
  snpA<-paste("V",i, sep="") #sets up string for snpA for pop num
  snpB<-paste("V",i+1, sep="") #sets up string for snpA for pop num
  P<-paste("P", pop_num, sep="") #sets up paper_ID (population ID)
  tmp<-snp_table %>% select(snpA,snpB) #graps all snps per paper_ID
  
  colnames(tmp)<-c("A", "B") #renames column headers
  
   snp_prop_A[,counter]<-tmp$A/(tmp$A + tmp$B) #calc proportion for A, assign to output
   colnames (snp_prop_A)[counter]<-P 
  
   counter<-counter+1
   pop_num<-pop_num+1
}
#snp_prop_A[is.na(snp_prop_A)] <- 0
return(snp_prop_A)
}

#Generate frequency matrix for prop B
prop_B <- function(snp_table){
#B Mat
snp_prop_B<-snp_table %>% select (chr_snp, Chromosome, SNP, Env)
counter=5
pop_num=1

for (i in seq(1,110,2)) {
  
  snpA<-paste("V",i, sep="")
  snpB<-paste("V",i+1, sep="")
  P<-paste("P", pop_num, sep="")
  tmp<-snp_table %>% select(snpA,snpB)
  
  colnames(tmp)<-c("A", "B")
  
  snp_prop_B[,counter]<-tmp$B/(tmp$A + tmp$B)
  colnames (snp_prop_B)[counter]<-P
  
  counter<-counter+1
  pop_num<-pop_num+1
}
#snp_prop_B[is.na(snp_prop_B)] <- 0
return(snp_prop_B)
}


###################################################################################
##set up frequency association table
#Positive assocation with climate
FAT_p <- function(snp_table,climate_table,env_in){
  snp_prop_A_in<-prop_A(snp_table) # call function
  snp_prop_B_in<-prop_B(snp_table) # call function

freq.temp <- data.frame()
for (i in 1:dim(snp_prop_A_in)[1]) { #for each row
  
env_pop <- climate_table %>% select(Paper_ID,as.character(env_in)) #Select climate data
snp_prop_A_tmp<-snp_prop_A_in %>% filter(SNP==snp_prop_A_in$SNP[i]) %>% 
  select(-chr_snp, -Chromosome, -SNP, -Env) #filter prop A data 
snp_prop_B_tmp<-snp_prop_B_in %>% filter(SNP==snp_prop_B_in$SNP[i]) %>%
  select(-chr_snp, -Chromosome, -SNP, -Env) #filter prop B data 

env_pop <- cbind(env_pop,t(snp_prop_A_tmp), t(snp_prop_B_tmp))
colnames(env_pop)[3] <-"prop_A"
colnames(env_pop)[4]<-"prop_B"

lm.temp_A <- lm(env_pop$prop_A~env_pop[,2]) # save lm of cliamte predicting prop A
lm.temp_B <- lm(env_pop$prop_B~env_pop[,2]) # save lm of cliamte predicting prop B

#decides if A or B is positively associated 
if (isTRUE(lm.temp_A$coefficents[2]>0) && isTRUE(lm.temp_B$coefficents[2]>0)){
  
  print(i)
  print(lm.temp_A) 
  print(lm.temp_B) 
  
  } else if(lm.temp_A$coefficients[2]>0){ 
    #print("A")
    tmp_in<-snp_prop_A_in %>% filter(SNP==snp_prop_A_in$SNP[i])
    freq.temp<-rbind(freq.temp,tmp_in)
  
  } else if (lm.temp_B$coefficients[2]>0){
    # print("B")
     tmp_in<-snp_prop_B_in %>% filter(SNP==snp_prop_B_in$SNP[i])
     freq.temp<-rbind(freq.temp,tmp_in)
     
  }else{
    #print(i)
    #print(lm.temp_A) 
    #print(lm.temp_B) 
  }
}
return(freq.temp)
}


#Negative assocation with climate
FAT_n <- function(snp_table,climate_table,env_in){
  snp_prop_A_in<-prop_A(snp_table) # call function
  snp_prop_B_in<-prop_B(snp_table) # call function
  
  freq.temp <- data.frame()
  for (i in 1:dim(snp_prop_A_in)[1]) { #for each row
    
    env_pop <- climate_table %>% select(Paper_ID,as.character(env_in)) #Select climate data
    snp_prop_A_tmp<-snp_prop_A_in %>% filter(SNP==snp_prop_A_in$SNP[i]) %>% 
      select(-chr_snp, -Chromosome, -SNP, -Env) #filter prop A data 
    snp_prop_B_tmp<-snp_prop_B_in %>% filter(SNP==snp_prop_B_in$SNP[i]) %>%
      select(-chr_snp, -Chromosome, -SNP, -Env) #filter prop B data 
    
    env_pop <- cbind(env_pop,t(snp_prop_A_tmp), t(snp_prop_B_tmp))
    colnames(env_pop)[3] <-"prop_A"
    colnames(env_pop)[4]<-"prop_B"
    
    lm.temp_A <- lm(env_pop$prop_A~env_pop[,2]) # save lm of cliamte predicting prop A
    lm.temp_B <- lm(env_pop$prop_B~env_pop[,2]) # save lm of cliamte predicting prop B
    
    #decides if A or B is positively associated 
    if (isTRUE(lm.temp_A$coefficents[2]<0) && isTRUE(lm.temp_B$coefficents[2]<0)){
      
      print(i)
      print(lm.temp_A) 
      print(lm.temp_B) 
      
    } else if(lm.temp_A$coefficients[2]<0){ 
      #print("A")
      tmp_in<-snp_prop_A_in %>% filter(SNP==snp_prop_A_in$SNP[i])
      freq.temp<-rbind(freq.temp,tmp_in)
      
    } else if (lm.temp_B$coefficients[2]<0){
      # print("B")
      tmp_in<-snp_prop_B_in %>% filter(SNP==snp_prop_B_in$SNP[i])
      freq.temp<-rbind(freq.temp,tmp_in)
      
    }else{
      #print(i)
      #print(lm.temp_A) 
      #print(lm.temp_B) 
    }
  }
  return(freq.temp)
}

###################################################################################
###################################################################################

#Climate variables considered
##### Annual #####
# MAT = Mean annual temperature (°C)
# MAP = Mean annual precipitation (mm)
# PAS = Precipitation as snow (mm) between August in previous year and July in current year
# EXT = Extreme temperature over 30 years
# CMD = Hargreaves climatic moisture deficit (mm)
##### Seasonal #####
# Tave_wt = Winter mean temperature (°C)
# Tave_sm = Summer mean temperature (°C)
# PPT_wt = Winter precipitation (mm)
# PPT_sm = Summer precipitation (mm)

#Get Frequency Association Table (FAT) for each environmental variable
freq_1 <- FAT_p(env1_loci,climate,"MAT")
freq_2 <- FAT_n(env2_loci,climate,"MAP")
freq_3 <- FAT_n(env3_loci,climate,"PAS")
freq_4 <- FAT_p(env4_loci,climate,"EXT")
freq_5 <- FAT_p(env5_loci,climate,"CMD")

freq_6 <- FAT_p(env6_loci,climate,"Tave_wt")
freq_7 <- FAT_p(env7_loci,climate,"Tave_sm")
freq_8 <- FAT_n(env8_loci,climate,"PPT_wt")
freq_9 <- FAT_n(env9_loci,climate,"PPT_sm")

#Set up binary table for presence/absence of "adaptive" snp

freq_lable <- freq_1 %>% select(chr_snp:Env)
part1 <- freq_1 %>% select(P1:P55)
part1[part1 > 0] <- 1
freq1_p_a <- cbind(freq_lable,part1)

freq_lable <- freq_2 %>% select(chr_snp:Env)
part2 <- freq_2 %>% select(P1:P55)
part2[part2 > 0] <- 1
freq2_p_a <- cbind(freq_lable,part2)

freq_lable <- freq_3 %>% select(chr_snp:Env)
part3 <- freq_3 %>% select(P1:P55)
part3[part3 > 0] <- 1
freq3_p_a <- cbind(freq_lable,part3)

freq_lable <- freq_4 %>% select(chr_snp:Env)
part4 <- freq_4 %>% select(P1:P55)
part4[part4 > 0] <- 1
freq4_p_a <- cbind(freq_lable,part4)

freq_lable <- freq_5 %>% select(chr_snp:Env)
part5 <- freq_5 %>% select(P1:P55)
part5[part5 > 0] <- 1
freq5_p_a <- cbind(freq_lable,part5)

freq_lable <- freq_6 %>% select(chr_snp:Env)
part6 <- freq_6 %>% select(P1:P55)
part6[part6 > 0] <- 1
freq6_p_a <- cbind(freq_lable,part6)

freq_lable <- freq_7 %>% select(chr_snp:Env)
part7 <- freq_7 %>% select(P1:P55)
part7[part7 > 0] <- 1
freq7_p_a <- cbind(freq_lable,part7)

freq_lable <- freq_8 %>% select(chr_snp:Env)
part8 <- freq_8 %>% select(P1:P55)
part8[part8 > 0] <- 1
freq8_p_a <- cbind(freq_lable,part8)

freq_lable <- freq_9 %>% select(chr_snp:Env)
part9 <- freq_9 %>% select(P1:P55)
part9[part9 > 0] <- 1
freq9_p_a <- cbind(freq_lable,part9)

#Export

write_csv(freq_1, "Genomics_scripts/Data/freq_1.csv")
write_csv(freq_2, "Genomics_scripts/Data/freq_2.csv")
write_csv(freq_3, "Genomics_scripts/Data/freq_3.csv")
write_csv(freq_4, "Genomics_scripts/Data/freq_4.csv")
write_csv(freq_5, "Genomics_scripts/Data/freq_5.csv")
write_csv(freq_6, "Genomics_scripts/Data/freq_6.csv")
write_csv(freq_7, "Genomics_scripts/Data/freq_7.csv")
write_csv(freq_8, "Genomics_scripts/Data/freq_8.csv")
write_csv(freq_9, "Genomics_scripts/Data/freq_9.csv")


write_csv(freq1_p_a, "Genomics_scripts/Data/freq_binary_1.csv")
write_csv(freq2_p_a, "Genomics_scripts/Data/freq_binary_2.csv")
write_csv(freq3_p_a, "Genomics_scripts/Data/freq_binary_3.csv")
write_csv(freq4_p_a, "Genomics_scripts/Data/freq_binary_4.csv")
write_csv(freq5_p_a, "Genomics_scripts/Data/freq_binary_5.csv")
write_csv(freq6_p_a, "Genomics_scripts/Data/freq_binary_6.csv")
write_csv(freq7_p_a, "Genomics_scripts/Data/freq_binary_7.csv")
write_csv(freq8_p_a, "Genomics_scripts/Data/freq_binary_8.csv")
write_csv(freq9_p_a, "Genomics_scripts/Data/freq_binary_9.csv")
  



