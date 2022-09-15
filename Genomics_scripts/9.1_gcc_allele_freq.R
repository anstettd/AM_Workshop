##################################################################################
## Generate SNP proportion matrix per region
## For peak window SNPS (WZA) with BF (BayPass) > 2
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

#Import peak bf>2 Timeseries data
snp1_time <- read_csv("Genomics_scripts/Data/peak_bf2_time_mat.csv")
snp2_time <- read_csv("Genomics_scripts/Data/peak_bf2_time_map.csv")
snp5_time <- read_csv("Genomics_scripts/Data/peak_bf2_time_cmd.csv")

###################################################################################
#Import snp Baseline Data BF>10, in order to get basetime frequencies
#basetime = Spacetime combinations (population at moment in time) shared in common between the baseline and timeseries data
env1_base <- read_csv("Genomics_scripts/Data/peak_bf2_base_mat.csv")
env2_base <- read_csv("Genomics_scripts/Data/peak_bf2_base_map.csv")
env5_base <- read_csv("Genomics_scripts/Data/peak_bf2_base_cmd.csv")

#Filter baseline to BF>20, remove 4 variables, 
#ensure baseline has same SNPs as timeseries and remove any not in the timeseries 
env1_base <- env1_base %>% filter (chr_snp %in% as.character(snp1_time$chr_snp))
env2_base <- env2_base %>% filter (chr_snp %in% as.character(snp2_time$chr_snp))
env5_base <- env5_base %>% filter (chr_snp %in% as.character(snp5_time$chr_snp))


#Import Baseline Climate 
climate <- read_csv("Donor_selection/Data/climate_pop.csv")
climate <- climate %>% select(Site_Name:MAP,CMD) #Select relevant climate variables

#Import pop names (site_year names)
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")

#Modify pop_order to includ regional and SNP information
pop_order_wide <- pop_order %>% separate(V1,c("Pop","Year"),sep="_") %>%
  mutate(Region = ifelse(Pop == 8 | Pop == 9 | Pop == 10 | Pop == 11, "North",
                         ifelse(Pop == 1 | Pop == 12 | Pop == 2, "South","Centre")))

pop_order_V <- data.frame()
counter <- 1
for(i in 1:dim(pop_order_wide)[1]){
pop_order_V[counter,1] <- pop_order_wide$Region[i]
pop_order_V[counter,2] <- pop_order_wide$Pop[i]
pop_order_V[counter,3] <- pop_order_wide$Year[i]
pop_order_V[counter,4] <- paste("V",counter, sep="")
pop_order_V[counter,5] <- "A"

pop_order_V[counter+1,1] <- pop_order_wide$Region[i]
pop_order_V[counter+1,2] <- pop_order_wide$Pop[i]
pop_order_V[counter+1,3] <- pop_order_wide$Year[i]
pop_order_V[counter+1,4] <- paste("V",counter+1, sep="")
pop_order_V[counter+1,5] <- "B"
counter <- counter+2
}  
colnames(pop_order_V) <- c("Region","Pop","Year","V_ID","SNP")

#write_csv(pop_order_V, "Genomics_scripts/Data/pop_order_V.csv")

#Make region_year
region_year <- c("North_2010","North_2011","North_2012","North_2013","North_2014","North_2015","North_2016",
                 "Centre_2010","Centre_2011","Centre_2012","Centre_2013","Centre_2014","Centre_2015","Centre_2016",
                 "South_2010","South_2011","South_2012","South_2013","South_2014","South_2015","South_2016")


###################################################################################
###################################################################################
##setup functions
#Sum site-year comibinations to get region data
make_region_snp <- function(tabulation,data_in,region,year,snp){
  region_filter <- tabulation %>% filter(Region==region,Year==year,SNP==snp)
  time_region <- data_in %>% select(chr_snp,all_of(region_filter$V_ID)) %>% 
    mutate(region_sum= rowSums(across(where(is.numeric)))) %>%
    select(region_sum)
  return(time_region)
}

#Generates the AB table per region/year
weave <- function(tabulation,data_in){
  Regions<-c("North", "Centre", "South")
  Year<-c(2010, 2011, 2012, 2013, 2014, 2015, 2016)
  SNP<-c("A","B")
  table_out <- data_in$chr_snp
  for (i in 1:3){
    for(j in 1:7){
      for(k in 1:2){
        temp1 <- make_region_snp(tabulation,data_in,Regions[i],Year[j],SNP[k])
        table_out <- cbind(table_out, temp1)
      }
    }
  }
for(i in 2:dim(table_out)[2]){
  colnames(table_out)[i]<-paste("V",i-1, sep="")
}
colnames(table_out)[1] <- "chr_snp"
return(table_out)  
}

###################################################################################
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
## Make region,time,snp input table
env1_region_time <- weave(pop_order_V,snp1_time)
env2_region_time <- weave(pop_order_V,snp2_time)
env5_region_time <- weave(pop_order_V,snp5_time)

## Make table with climate change associated SNPs for timeseries using baseline climate
freq_MAT_1 <- FAT_p(env1_base,env1_region_time,climate,"MAT")
freq_MAP_1 <- FAT_n(env2_base,env2_region_time,climate,"MAP")
freq_CMD_1 <- FAT_p(env5_base,env5_region_time,climate,"MAP")

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
colnames(freq_MAT)<- region_year #name each pop/time combination
colnames(freq_MAP)<- region_year #name each pop/time combination
colnames(freq_CMD)<- region_year #name each pop/time combination

#Transpose and split up region_year
freq_MAT_T <- as.data.frame(t(freq_MAT)) %>% rownames_to_column ("region_year") %>% separate(region_year, c("Site","Year"))
freq_MAP_T <- as.data.frame(t(freq_MAP)) %>% rownames_to_column ("region_year") %>% separate(region_year, c("Site","Year"))
freq_CMD_T <- as.data.frame(t(freq_CMD)) %>% rownames_to_column ("region_year") %>% separate(region_year, c("Site","Year"))

# Write out climate change correlated timeseries frequency table 
write_csv(freq_MAT_T, "Genomics_scripts/Data/freq_MAT_peakbf2_region.csv")
write_csv(freq_MAP_T, "Genomics_scripts/Data/freq_MAP_peakbf2_region.csv")
write_csv(freq_CMD_T, "Genomics_scripts/Data/freq_CMD_peakbf2_region.csv")




  

