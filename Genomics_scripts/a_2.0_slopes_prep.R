##################################################################################
## Get slopes (strength of selection) for all BF
## 
## Author Daniel Anstett
## 
##
## Last Modified September 26, 2022
###################################################################################


###################################################################################

#Import libraries
library(tidyverse)

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
    print(i)
    #decides if A or B is positively associated 
    if (is.na(lm.temp_A$coefficients[2]) || is.na(lm.temp_B$coefficients[2])) {
      print("skip")
    } else if (isTRUE(lm.temp_A$coefficents[2]>0) && isTRUE(lm.temp_B$coefficents[2]>0)){
      
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
  
  #Make into dataframe from tibble
  freq <- as.data.frame(freq.temp)
  
  #Make first column a row name
  rownames(freq)<- as.vector(freq$chr_snp)
  
  #Remove 1st column
  freq <- freq %>% select(-chr_snp)
  
  # Add in row names to SNP frequency tables 
  colnames(freq)<- pop_order[,1] #name each pop/time combination
  
  #Transpose and split up site_year
  freq_T <- as.data.frame(t(freq)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
  
  return(freq_T)
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
    print(i)
    #decides if A or B is positively associated 
    if (is.na(lm.temp_A$coefficients[2]) || is.na(lm.temp_B$coefficients[2])) {
      print("skip")
    } else if (isTRUE(lm.temp_A$coefficents[2]>0) && isTRUE(lm.temp_B$coefficents[2]>0)){
      
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
  
  #Make into dataframe from tibble
  freq <- as.data.frame(freq.temp)
  
  #Make first column a row name
  rownames(freq)<- as.vector(freq$chr_snp)
  
  #Remove 1st column
  freq <- freq %>% select(-chr_snp)
  
  # Add in row names to SNP frequency tables 
  colnames(freq)<- pop_order[,1] #name each pop/time combination
  
  #Transpose and split up site_year
  freq_T <- as.data.frame(t(freq)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))

  return(freq_T)
}

########################################################################################################
########################################################################################################

#Import full snp table for timeseries
pop_order_time<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp_time<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci_time<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
colnames(loci_time) <- c("Chromosome","SNP")
loci_united_time <- loci_time %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp_time <-as.data.frame(cbind(loci_united_time,snp_time)) #add snp lables to rows
#loci_snp_time <- loci_snp_time[24575,]
loci_snp_time_lable <- as.data.frame(loci_snp_time$chr_snp)
colnames(loci_snp_time_lable) <- "chr_snp"

#Import Baseline data
#import full snp table
pop_order_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", 
                      header=F, sep="\t")
snp_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", 
                header=F, sep=" ")
loci_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", 
                 header=F, sep="\t")
colnames(loci_base) <- c("Chromosome","SNP")
loci_united_base <- loci_base %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp_base <- as.data.frame(cbind(loci_united_base,snp_base)) #add snp lables to rows
loci_snp_base_lable <- as.data.frame(loci_snp_base$chr_snp)
colnames(loci_snp_base_lable) <- "chr_snp"

#Ensure baseline has same SNPs as timeseries and remove any not in the timeseries 
loci_snp_base_join <- left_join(loci_snp_time_lable,loci_snp_base,by="chr_snp")
loci_snp_base_filter <- na.omit(loci_snp_base_join)

#Ensure timeseries has same SNPs as baseline and remove any not in the baseline 
loci_snp_time_join <- left_join(loci_snp_base_lable,loci_snp_time,by="chr_snp")
loci_snp_time_filter <- na.omit(loci_snp_time_join)

#Import Baseline Climate 
climate <- read_csv("Donor_selection/Data/climate_pop.csv")
climate <- climate %>% select(Site_Name:MAP,CMD) #Select relevant climate variables

#Import pop names (site_year names)
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")


########################################################################################################
#Breakup into seqments
#loci_snp_time_1 <- loci_snp_time_filter[1:100000,]
#loci_snp_base_1 <- loci_snp_base_filter[1:100000,]
#loci_snp_time_2 <- loci_snp_time_filter[100001:200000,]
#loci_snp_base_2 <- loci_snp_base_filter[100001:200000,]
#loci_snp_time_3 <- loci_snp_time_filter[200001:300000,]
#loci_snp_base_3 <- loci_snp_base_filter[200001:300000,]
#loci_snp_time_4 <- loci_snp_time_filter[300001:400000,]
#loci_snp_base_4 <- loci_snp_base_filter[300001:400000,]
#loci_snp_time_5 <- loci_snp_time_filter[400001:500000,]
#loci_snp_base_5 <- loci_snp_base_filter[400001:500000,]
#loci_snp_time_6 <- loci_snp_time_filter[500001:600000,]
#loci_snp_base_6 <- loci_snp_base_filter[500001:600000,]
#loci_snp_time_7 <- loci_snp_time_filter[600001:700000,]
#loci_snp_base_7 <- loci_snp_base_filter[600001:700000,]
#loci_snp_time_8 <- loci_snp_time_filter[700001:800000,]
#loci_snp_base_8 <- loci_snp_base_filter[700001:800000,]
#loci_snp_time_9 <- loci_snp_time_filter[800001:900000,]
#loci_snp_base_9 <- loci_snp_base_filter[800001:900000,]
#loci_snp_time_10 <- loci_snp_time_filter[900001:1000000,]
#loci_snp_base_10 <- loci_snp_base_filter[900001:1000000,]
#loci_snp_time_11 <- loci_snp_time_filter[1000001:1100000,]
#loci_snp_base_11 <- loci_snp_base_filter[1000001:1100000,]
#loci_snp_time_12 <- loci_snp_time_filter[1100001:1200000,]
#loci_snp_base_12 <- loci_snp_base_filter[1100001:1200000,]
#loci_snp_time_13 <- loci_snp_time_filter[1200001:1300000,]
#loci_snp_base_13 <- loci_snp_base_filter[1200001:1300000,]
#loci_snp_time_14 <- loci_snp_time_filter[1300001:1400000,]
#loci_snp_base_14 <- loci_snp_base_filter[1300001:1400000,]

#loci_snp_time_15 <- loci_snp_time_filter[1400001:1500000,]
#loci_snp_base_15 <- loci_snp_base_filter[1400001:1500000,]
#loci_snp_time_16 <- loci_snp_time_filter[1500001:1600000,]
#loci_snp_base_16 <- loci_snp_base_filter[1500001:1600000,]
#loci_snp_time_17 <- loci_snp_time_filter[1600001:1700000,]
#loci_snp_base_17 <- loci_snp_base_filter[1600001:1700000,]

#loci_snp_time_18 <- loci_snp_time_filter[1700001:dim(loci_snp_time_filter)[1],]
#loci_snp_base_18 <- loci_snp_base_filter[1700001:dim(loci_snp_time_filter)[1],]

########################################################################################################


## Make table with climate change associated SNPs for timeseries using baseline climate
#freq_MAT_1 <- FAT_p(loci_snp_base_1,loci_snp_time_1,climate,"MAT")
#write_csv(freq_MAT_1, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_1.csv")
#freq_MAT_2 <- FAT_p(loci_snp_base_2,loci_snp_time_2,climate,"MAT")
#write_csv(freq_MAT_2, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_2.csv")
#freq_MAT_3 <- FAT_p(loci_snp_base_3,loci_snp_time_3,climate,"MAT")
#write_csv(freq_MAT_3, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_3.csv")
#freq_MAT_4 <- FAT_p(loci_snp_base_4,loci_snp_time_4,climate,"MAT")
#write_csv(freq_MAT_4, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_4.csv")
#freq_MAT_5 <- FAT_p(loci_snp_base_5,loci_snp_time_5,climate,"MAT")
#write_csv(freq_MAT_5, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_5.csv")
#freq_MAT_6 <- FAT_p(loci_snp_base_6,loci_snp_time_6,climate,"MAT")
#write_csv(freq_MAT_6, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_6.csv")
#freq_MAT_7 <- FAT_p(loci_snp_base_7,loci_snp_time_7,climate,"MAT")
#write_csv(freq_MAT_7, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_7.csv")
#freq_MAT_8 <- FAT_p(loci_snp_base_8,loci_snp_time_8,climate,"MAT")
#write_csv(freq_MAT_8, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_8.csv")
#freq_MAT_9 <- FAT_p(loci_snp_base_9,loci_snp_time_9,climate,"MAT")
#write_csv(freq_MAT_9, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_9.csv")
#freq_MAT_10 <- FAT_p(loci_snp_base_10,loci_snp_time_10,climate,"MAT")
#write_csv(freq_MAT_10, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_10.csv")
#freq_MAT_11 <- FAT_p(loci_snp_base_11,loci_snp_time_11,climate,"MAT")
#write_csv(freq_MAT_11, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_11.csv")
#freq_MAT_12 <- FAT_p(loci_snp_base_12,loci_snp_time_12,climate,"MAT")
#write_csv(freq_MAT_12, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_12.csv")
#freq_MAT_13 <- FAT_p(loci_snp_base_13,loci_snp_time_13,climate,"MAT")
#write_csv(freq_MAT_13, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_13.csv")
#freq_MAT_14 <- FAT_p(loci_snp_base_14,loci_snp_time_14,climate,"MAT")
#write_csv(freq_MAT_14, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_14.csv")
#freq_MAT_15 <- FAT_p(loci_snp_base_15,loci_snp_time_15,climate,"MAT")
#write_csv(freq_MAT_15, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_15.csv")
#freq_MAT_16 <- FAT_p(loci_snp_base_16,loci_snp_time_16,climate,"MAT")
#write_csv(freq_MAT_16, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_16.csv")
#freq_MAT_17 <- FAT_p(loci_snp_base_17,loci_snp_time_17,climate,"MAT")
#write_csv(freq_MAT_17, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_17.csv")

#freq_MAT_18 <- FAT_p(loci_snp_base_18,loci_snp_time_18,climate,"MAT")
#write_csv(freq_MAT_18, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_18.csv")



#freq_MAP_1 <- FAT_n(loci_snp_base,loci_snp_time,climate,"MAP")
#freq_CMD_1 <- FAT_p(loci_snp_base,loci_snp_time,climate,"MAP")

