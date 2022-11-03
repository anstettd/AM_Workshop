##################################################################################
## Generate SNP proportion matrix per population
## For peak window SNPS (WZA) with BF (BayPass) > 10 and all BF>30
## Done for all alleles positively associated with climate change
## Select if A or B is positively associated with climate change
## Author Daniel Anstett
## 
## Currently done just for ENV 1,2 and 5 (MAT, MAP and CMD)
## Last Modified September 15, 2022
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

#Import peak bf>5 Timeseries data
snp1_time <- read_csv("Genomics_scripts/Data/peak_bf5_time_mat.csv")
snp2_time <- read_csv("Genomics_scripts/Data/peak_bf5_time_map.csv")
snp5_time <- read_csv("Genomics_scripts/Data/peak_bf5_time_cmd.csv")

###################################################################################
#Import peak snp Baseline data bf>5, in order to get basetime frequencies
#basetime = Spacetime combinations (population at moment in time) shared in common between the baseline and timeseries data
env1_base <- read_csv("Genomics_scripts/Data/peak_bf5_base_mat.csv")
env2_base <- read_csv("Genomics_scripts/Data/peak_bf5_base_map.csv")
env5_base <- read_csv("Genomics_scripts/Data/peak_bf5_base_cmd.csv")

#Ensure baseline has same SNPs as timeseries and remove any not in the timeseries 
env1_base <- env1_base %>% filter (chr_snp %in% as.character(snp1_time$chr_snp))
env2_base <- env2_base %>% filter (chr_snp %in% as.character(snp2_time$chr_snp))
env5_base <- env5_base %>% filter (chr_snp %in% as.character(snp5_time$chr_snp))


#Import Baseline Climate 
climate <- read_csv("Donor_selection/Data/climate_pop.csv")
climate <- climate %>% select(Site_Name:MAP,CMD) #Select relevant climate variables

#Import pop names (site_year names)
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")

###################################################################################
###################################################################################
##setup functions

#Generate frequency matrix for prop A 
prop_A <- function(snp_table) {
  snp_prop_A<- snp_table %>% select (chr_snp)
  counter=2
  pop_num=1
  for (i in seq(1,dim(snp_table)[2]-1,2)) {
    
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
  snp_prop_B<-snp_table %>% select (chr_snp)
  counter=2
  pop_num=1
  
  for (i in seq(1,dim(snp_table)[2]-1,2)) {
    
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
FAT_p <- function(snp_base,snp_time,climate_table,env_in){
  snp_prop_A_in<-prop_A(snp_base) # call function
  snp_prop_B_in<-prop_B(snp_base) # call function
  
  snp_prop_A_in_time<-prop_A(snp_time)
  snp_prop_B_in_time<-prop_B(snp_time)
  
  freq.temp <- data.frame()
  for (i in 1:dim(snp_prop_A_in)[1]) { #for each row
    
    env_pop <- climate_table %>% select(Paper_ID,as.character(env_in)) #Select climate data
    snp_prop_A_tmp<-snp_prop_A_in %>% filter(chr_snp==snp_prop_A_in$chr_snp[i]) %>% 
      select(-chr_snp) #filter prop A data 
    snp_prop_B_tmp<-snp_prop_B_in %>% filter(chr_snp==snp_prop_B_in$chr_snp[i]) %>%
      select(-chr_snp) #filter prop B data 
    
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
      tmp_in<-snp_prop_A_in_time %>% filter(chr_snp==snp_prop_A_in$chr_snp[i])
      freq.temp<-rbind(freq.temp,tmp_in)
      
    } else if (lm.temp_B$coefficients[2]>0){
      # print("B")
      tmp_in<-snp_prop_B_in_time %>% filter(chr_snp==snp_prop_B_in$chr_snp[i])
      freq.temp<-rbind(freq.temp,tmp_in)
      
    }
  }
  return(freq.temp)
}


#Negative assocation with climate
FAT_n <- function(snp_base,snp_time,climate_table,env_in){
  snp_prop_A_in<-prop_A(snp_base) # call function
  snp_prop_B_in<-prop_B(snp_base) # call function
  
  snp_prop_A_in_time<-prop_A(snp_time)
  snp_prop_B_in_time<-prop_B(snp_time)
  
  freq.temp <- data.frame()
  for (i in 1:dim(snp_prop_A_in)[1]) { #for each row
    
    env_pop <- climate_table %>% select(Paper_ID,as.character(env_in)) #Select climate data
    snp_prop_A_tmp<-snp_prop_A_in %>% filter(chr_snp==snp_prop_A_in$chr_snp[i]) %>% 
      select(-chr_snp) #filter prop A data 
    snp_prop_B_tmp<-snp_prop_B_in %>% filter(chr_snp==snp_prop_B_in$chr_snp[i]) %>%
      select(-chr_snp) #filter prop B data 
    
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
      tmp_in<-snp_prop_A_in_time %>% filter( chr_snp==snp_prop_A_in$chr_snp[i])
      freq.temp<-rbind(freq.temp,tmp_in)
      
    } else if (lm.temp_B$coefficients[2]<0){
      # print("B")
      tmp_in<-snp_prop_B_in_time %>% filter(chr_snp==snp_prop_B_in$chr_snp[i])
      freq.temp<-rbind(freq.temp,tmp_in)
    }
  }
  return(freq.temp)
}

########################################################################################################
########################################################################################################
## Make table with climate change associated SNPs for timeseries using baseline climate
freq_MAT_1 <- FAT_p(env1_base,snp1_time,climate,"MAT")
freq_MAP_1 <- FAT_n(env2_base,snp2_time,climate,"MAP")
freq_CMD_1 <- FAT_p(env5_base,snp5_time,climate,"CMD")

#Make into dataframe from tibble
freq_MAT<- as.data.frame(freq_MAT_1)
freq_MAP<- as.data.frame(freq_MAP_1)
freq_CMD<- as.data.frame(freq_CMD_1)

#Make first column a row name
rownames(freq_MAT)<- as.vector(freq_MAT$chr_snp)
rownames(freq_MAP)<- as.vector(freq_MAP$chr_snp)
rownames(freq_CMD)<- as.vector(freq_CMD$chr_snp)

#Remove 1st column
freq_MAT <- freq_MAT %>% select(-chr_snp)
freq_MAP <- freq_MAP %>% select(-chr_snp)
freq_CMD <- freq_CMD %>% select(-chr_snp)

# Add in row names to SNP frequency tables 
colnames(freq_MAT)<- pop_order[,1] #name each pop/time combination
colnames(freq_MAP)<- pop_order[,1] #name each pop/time combination
colnames(freq_CMD)<- pop_order[,1] #name each pop/time combination

#Transpose and split up site_year
freq_MAT_T <- as.data.frame(t(freq_MAT)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_MAP_T <- as.data.frame(t(freq_MAP)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_CMD_T <- as.data.frame(t(freq_CMD)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))

# Write out climate change correlated timeseries frequency table 
write_csv(freq_MAT_T, "Genomics_scripts/Data/freq_MAT_peakbf5.csv")
write_csv(freq_MAP_T, "Genomics_scripts/Data/freq_MAP_peakbf5.csv")
write_csv(freq_CMD_T, "Genomics_scripts/Data/freq_CMD_peakbf5.csv")




  

