##################################################################################
## Generate stratified distribution of random non-climate associated slopes
## Author Daniel Anstett
## 
## 
## Last Modified April 15, 2022
###################################################################################
#Import libraries
library(tidyverse)

###################################################################################
#Functions for generating stratigied random distribution
###################################################################################
#Large
#pop_snp = p1A
#freq_count = freq_count_MAT_IDX
large <- function(pop_snp,freq_count){
  rand_pop <- data.frame()
  bin_fraction<-1/length(freq_count)
  bin_step<-seq(0,1,bin_fraction)
  
  for (i in 1: length(freq_count)){
    if(i==length(freq_count)){
      p1_n<-pop_snp %>% filter(snpA>=bin_step[i] & snpA <= bin_step[i+1])
      list_env_p1_n<- sample.int(dim(p1_n)[1],freq_count[i],replace = FALSE)
      
    }else{
      p1_n<-pop_snp %>% filter(snpA>=bin_step[i] & snpA < bin_step[i+1])
      list_env_p1_n<- sample.int(dim(p1_n)[1],freq_count[i],replace = FALSE)
      
    }
    
    for (j in 1:length(list_env_p1_n)){
      rand_pop <- rbind(rand_pop ,p1_n %>% 
                          filter(chr_snp==as.character(p1_n$chr_snp[list_env_p1_n[j]])))
    }
    
  }
  
  return(rand_pop)
}

#Freq Count function
freq_bins <- function(basetime){
  freq_count_calc<-data.frame()
  bin_size<-4
  Regions<-c("North", "Centre", "South")
  for (i in 1:3){
    
    #   if(i==1 || i==2 || i==5 || i==7 || i==10 || i==12){
    #    }else if (i==3 || i==4 || i==6 || i==11){
    #      bin_size<-4
    #    }else if (i==8){
    #      bin_size<-2
    #    }else if(i==9){
    #      bin_size<-5
    #    }
    
    bin_fraction<-1/bin_size
    bin_step<-seq(0,1,bin_fraction)
    
    for(j in 1:bin_size){
      test_ENV <- basetime %>% filter(Site==Regions[i]) %>% select (-Site,-Year)
      test_ENV <- as.data.frame(test_ENV)
      test_ENV <- as.numeric(test_ENV[1,])
      if (j==bin_size){
        freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV <= bin_step[j+1],na.rm=T )
      }else{
        freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV < bin_step[j+1],na.rm=T )
      }  
    }
  }
  return(freq_count_calc)
}

#Generate frequency matrix for prop A 
prop_A <- function(snp_table,pop_ID) {
  snp_prop_A<- snp_table %>% select (chr_snp)
  counter=2
  pop_num=1
  for (i in seq(1,dim(snp_table)[2]-1,2)) {
    
    snpA<-paste("V",i, sep="") #sets up string for snpA for pop num
    snpB<-paste("V",i+1, sep="") #sets up string for snpA for pop num
    P<-paste("P", pop_num, sep="") #sets up paper_ID (population ID)
    tmp<-snp_table %>% select(snpA,snpB) #graps all snps per paper_ID
    
    colnames(tmp)<-c("A", "B") #renames column headers
    
    snp_prop_A[,counter]<-as.numeric(tmp$A/(tmp$A + tmp$B)) #calc proportion for A, assign to output
    colnames (snp_prop_A)[counter]<-P 
    
    counter<-counter+1
    pop_num<-pop_num+1
  }
  
  colnames(snp_prop_A)<- pop_ID[,1] #name each pop/time combination
  rownames(snp_prop_A)<- snp_prop_A$chr_snp
  snp1A_T <- as.data.frame(t(snp_prop_A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
  colnames(snp1A_T) <- snp1A_T[1,]
  snp1A_T <- snp1A_T[-1,]
  colnames(snp1A_T)[1]<- "Site"
  colnames(snp1A_T)[2]<- "Year"
  
  return(snp1A_T)
}

glm_slopes<-function(snp_popX){
  
  snp_popX_slope<-data.frame()
  counter<-1
  
  for (i in 3:dim(snp_popX)[2]){
    chr<-colnames(snp_popX)[i]
    popSNP <- snp_popX %>% select(Site,Year,all_of(chr))
    colnames(popSNP)[3]<-"snp_ID"
    
    if(all(is.na(popSNP$snp_ID))==FALSE && length(unique(as.numeric(popSNP$snp_ID)))>1 
       && sum(!is.na( as.numeric((popSNP$snp_ID))))>1){
      popSNP <- na.omit(popSNP)
      rSNP <- glm(as.numeric(snp_ID) ~ Year, family = binomial, data = popSNP)
      snp_popX_slope[counter,1]<-unique(popSNP$Site)
      snp_popX_slope[counter,2]<-chr
      snp_popX_slope[counter,3]<-rSNP$coefficients[2]
      counter<-counter+1
    } else if((all(is.na(popSNP$snp_ID))==FALSE && length(unique(as.numeric(popSNP$snp_ID)))==1)){
      snp_popX_slope[counter,1]<-unique(popSNP$Site)
      snp_popX_slope[counter,2]<-chr
      snp_popX_slope[counter,3]<-0
      counter<-counter+1
    }else {
      snp_popX_slope[counter,1]<-unique(popSNP$Site)
      snp_popX_slope[counter,2]<-chr
      snp_popX_slope[counter,3]<-NA
      counter<-counter+1
    }
  }
  colnames(snp_popX_slope)<-c("Site","snp_ID","Slope")
  return(snp_popX_slope) 
}


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

## Manupulate entire SNP datatset for timeseries

#Import timeseries counts for A and B
#snp1_time <- read_csv("Genomics_scripts/Data/snp1_filter.csv")
#snp2_time <- read_csv("Genomics_scripts/Data/snp2_filter.csv")
#snp5_time <- read_csv("Genomics_scripts/Data/snp5_filter.csv")

#Import BF>0 baseline SNPs
MAT_bf0 <- read_csv("Genomics_scripts/Data/env1_BF0.csv")
MAP_bf0 <- read_csv("Genomics_scripts/Data/env2_BF0.csv")
CMD_bf0 <- read_csv("Genomics_scripts/Data/env5_BF0.csv")

#Import full snp table for timeseries
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
colnames(loci) <- c("Chromosome","SNP")
loci_united <- loci %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp <-cbind(loci_united,snp) #add snp lables to rows

#Import pop_order table with regional and V information
pop_order_V <- read_csv("Genomics_scripts/Data/pop_order_V.csv")

#Add Chromsome ID to pop order
pop_order_2 <- data.frame()
pop_order_2 [1,1] <- "chr_shp"
pop_order_2 <- rbind(pop_order_2,pop_order)

#Filter full snp table to remove climate associated SNPs
snp_swiss <-loci_snp %>% filter (!chr_snp %in% as.character(MAT_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(MAP_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(CMD_bf0$chr_snp))

###################################################################################
#Import transformed timeseries frequencies selected from BF>20 SNPs (i.e. observed data)
freq_MAT <- read_csv("Genomics_scripts/Data/freq_MAT_peakbf2_region.csv")
freq_MAP <- read_csv("Genomics_scripts/Data/freq_MAP_peakbf2_region.csv")
freq_CMD <- read_csv("Genomics_scripts/Data/freq_CMD_peakbf2_region.csv")
###################################################################################

#Filter frequency table
freq_MAT_1011 <- freq_MAT %>% filter(Site == "North" & Year==2010 |
                                       Site == "Centre" & Year==2010 |
                                       Site == "South" & Year==2010)

freq_MAP_1011 <- freq_MAP %>% filter(Site == "North" & Year==2010 |
                                       Site == "Centre" & Year==2010 |
                                       Site == "South" & Year==2010)

freq_CMD_1011 <- freq_CMD %>% filter(Site == "North" & Year==2010 |
                                       Site == "Centre" & Year==2010 |
                                       Site == "South" & Year==2010)


#Determine distribution of frequncy count table
freq_MAT_1011_T <- freq_MAT_1011 %>% select(-Year,-Site)
freq_MAT_1011_T <- as.data.frame(t(freq_MAT_1011_T))

freq_MAP_1011_T <- freq_MAP_1011 %>% select(-Year,-Site)
freq_MAP_1011_T <- as.data.frame(t(freq_MAP_1011_T))

freq_CMD_1011_T <- freq_CMD_1011 %>% select(-Year,-Site)
freq_CMD_1011_T <- as.data.frame(t(freq_CMD_1011_T))

#Display histograms
hist(freq_MAT_1011_T$V1)
hist(freq_MAT_1011_T$V2)
hist(freq_MAT_1011_T$V3)

hist(freq_MAP_1011_T$V1)
hist(freq_MAP_1011_T$V2)
hist(freq_MAP_1011_T$V3)

hist(freq_CMD_1011_T$V1)
hist(freq_CMD_1011_T$V2)
hist(freq_CMD_1011_T$V3)


#Make frequency count table
freq_count_MAT <- freq_bins(freq_MAT_1011)
freq_count_MAP <- freq_bins(freq_MAP_1011)
freq_count_CMD <- freq_bins(freq_CMD_1011)

#Break up freq count into individual vectors
freq_count_MAT_1 <- as.vector(freq_count_MAT[,1])
freq_count_MAT_2 <- as.vector(freq_count_MAT[,2])
freq_count_MAT_3 <- as.vector(freq_count_MAT[,3])

freq_count_MAP_1 <- as.vector(freq_count_MAP[,1])
freq_count_MAP_2 <- as.vector(freq_count_MAP[,2])
freq_count_MAP_3 <- as.vector(freq_count_MAP[,3])

freq_count_CMD_1 <- as.vector(freq_count_CMD[,1])
freq_count_CMD_2 <- as.vector(freq_count_CMD[,2])
freq_count_CMD_3 <- as.vector(freq_count_CMD[,3])



#######################################################################################################
#Get all basetime frequencies per region

#Select 2011 values for pop 8 and 10
#For the reset select 2010 values

#Filter pop_order_V to call for 2010 regional basetime V column names
pop_order_North <- pop_order_V %>% filter(Year==2010 & Region=="North")
pop_order_Centre <- pop_order_V %>% filter(Year==2010 & Region=="Centre")
pop_order_South <- pop_order_V %>% filter(Year==2010 & Region=="South")


#Filter Baseline A and B numbers for 3 basetime regions
loci_base_North <- snp_swiss %>% select(chr_snp,all_of(pop_order_North$V_ID)) %>% 
  mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_Centre <- snp_swiss %>% select(chr_snp,all_of(pop_order_Centre$V_ID)) %>% 
  mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_South <- snp_swiss %>% select(chr_snp,all_of(pop_order_South$V_ID)) %>% 
  mutate(region_sum= rowSums(across(where(is.numeric))))

#Calculate frequency for SNP A for 12 basetime pops
#Gives all basetime frequencies per region
p1A <- loci_base_North %>% mutate(snpA=(V11+V113)/(region_sum)) %>% select(chr_snp,snpA)
p2A <- loci_base_Centre %>% mutate(snpA=(V51+V65+V75+V85+V95)/(region_sum)) %>% select(chr_snp,snpA)
p3A <- loci_base_South %>% mutate(snpA=(V21+V31+V41)/(region_sum)) %>% select(chr_snp,snpA)


#Get region order as carried out in weave function
region_order <- data.frame()
region_order[1,1] <- "chr_snp"
counter <- 2
Regions<-c("North", "Centre", "South")
Year<-c(2010, 2011, 2012, 2013, 2014, 2015, 2016)

for (i in 1:3){
  for(j in 1:7){
    temp1 <- paste(Regions[i],Year[j],sep="_")
    region_order[counter,1] <- temp1
    counter <- counter + 1
  }
}

#######################################################################################################
# Stratified random sampling SNPs for a given population

set.seed(1)
#mat
rand_mat_pop1 <- large(p1A,freq_count_MAT_1)
rand_mat_pop2 <- large(p2A,freq_count_MAT_2)
rand_mat_pop3 <- large(p3A,freq_count_MAT_3)

rand_map_pop1 <- large(p1A,freq_count_MAP_1)
rand_map_pop2 <- large(p2A,freq_count_MAP_2)
rand_map_pop3 <- large(p3A,freq_count_MAP_3)

rand_cmd_pop1 <- large(p1A,freq_count_CMD_1)
rand_cmd_pop2 <- large(p2A,freq_count_CMD_2)
rand_cmd_pop3 <- large(p3A,freq_count_CMD_3)

#Get Full timeseires for each randomly selected neutral location
rand_time_AB_mat_1 <-loci_snp %>% filter (chr_snp %in% as.character(rand_mat_pop1$chr_snp)) #mat
rand_time_AB_mat_2 <-loci_snp %>% filter (chr_snp %in% as.character(rand_mat_pop2$chr_snp)) #mat
rand_time_AB_mat_3 <-loci_snp %>% filter (chr_snp %in% as.character(rand_mat_pop3$chr_snp)) #mat

rand_time_AB_map_1 <-loci_snp %>% filter (chr_snp %in% as.character(rand_map_pop1$chr_snp)) #map
rand_time_AB_map_2 <-loci_snp %>% filter (chr_snp %in% as.character(rand_map_pop2$chr_snp)) #map
rand_time_AB_map_3 <-loci_snp %>% filter (chr_snp %in% as.character(rand_map_pop3$chr_snp)) #map

rand_time_AB_cmd_1 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop1$chr_snp)) #cmd
rand_time_AB_cmd_2 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop2$chr_snp)) #cmd
rand_time_AB_cmd_3 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop3$chr_snp)) #cmd

#Convert filtered timeresires data from population to region
weave_mat_1 <- weave(pop_order_V,rand_time_AB_mat_1)
weave_mat_2 <- weave(pop_order_V,rand_time_AB_mat_2)
weave_mat_3 <- weave(pop_order_V,rand_time_AB_mat_3)

weave_map_1 <- weave(pop_order_V,rand_time_AB_map_1)
weave_map_2 <- weave(pop_order_V,rand_time_AB_map_2)
weave_map_3 <- weave(pop_order_V,rand_time_AB_map_3)

weave_cmd_1 <- weave(pop_order_V,rand_time_AB_cmd_1)
weave_cmd_2 <- weave(pop_order_V,rand_time_AB_cmd_2)
weave_cmd_3 <- weave(pop_order_V,rand_time_AB_cmd_3)

# Calc SNP A and transpose
rand_mat_pop1 <- prop_A(weave_mat_1,region_order)
rand_mat_pop2 <- prop_A(weave_mat_2,region_order)
rand_mat_pop3 <- prop_A(weave_mat_3,region_order)

rand_map_pop1 <- prop_A(weave_map_1,region_order)
rand_map_pop2 <- prop_A(weave_map_2,region_order)
rand_map_pop3 <- prop_A(weave_map_3,region_order)

rand_cmd_pop1 <- prop_A(weave_cmd_1,region_order)
rand_cmd_pop2 <- prop_A(weave_cmd_2,region_order)
rand_cmd_pop3 <- prop_A(weave_cmd_3,region_order)

# Filter for site (its currently duplicated since timeseries filtering was on 12 pop datatset)
rand_mat_pop1 <- rand_mat_pop1 %>% filter(Site=="North")
rand_mat_pop2 <- rand_mat_pop2 %>% filter(Site=="Centre")
rand_mat_pop3 <- rand_mat_pop3 %>% filter(Site=="South")

rand_map_pop1 <- rand_map_pop1 %>% filter(Site=="North")
rand_map_pop2 <- rand_map_pop2 %>% filter(Site=="Centre")
rand_map_pop3 <- rand_map_pop3 %>% filter(Site=="South")

rand_cmd_pop1 <- rand_cmd_pop1 %>% filter(Site=="North")
rand_cmd_pop2 <- rand_cmd_pop2 %>% filter(Site=="Centre")
rand_cmd_pop3 <- rand_cmd_pop3 %>% filter(Site=="South")


## Gather and merge rand for each env
gather_mat_pop1 <- rand_mat_pop1 %>% gather(SNP_ID,SNP_Freq,3:dim(rand_mat_pop1)[2])
gather_mat_pop2 <- rand_mat_pop2 %>% gather(SNP_ID,SNP_Freq,3:dim(rand_mat_pop2)[2])
gather_mat_pop3 <- rand_mat_pop3 %>% gather(SNP_ID,SNP_Freq,3:dim(rand_mat_pop3)[2])

gather_map_pop1 <- rand_map_pop1 %>% gather(SNP_ID,SNP_Freq,3:dim(rand_map_pop1)[2])
gather_map_pop2 <- rand_map_pop2 %>% gather(SNP_ID,SNP_Freq,3:dim(rand_map_pop2)[2])
gather_map_pop3 <- rand_map_pop3 %>% gather(SNP_ID,SNP_Freq,3:dim(rand_map_pop3)[2])

gather_cmd_pop1 <- rand_cmd_pop1 %>% gather(SNP_ID,SNP_Freq,3:dim(rand_cmd_pop1)[2])
gather_cmd_pop2 <- rand_cmd_pop2 %>% gather(SNP_ID,SNP_Freq,3:dim(rand_cmd_pop2)[2])
gather_cmd_pop3 <- rand_cmd_pop3 %>% gather(SNP_ID,SNP_Freq,3:dim(rand_cmd_pop3)[2])

rand_gathered_MAT <- rbind(gather_mat_pop1,gather_mat_pop2,gather_mat_pop3)
rand_gathered_MAP <- rbind(gather_map_pop1,gather_map_pop2,gather_map_pop3)
rand_gathered_CMD <- rbind(gather_cmd_pop1,gather_cmd_pop2,gather_cmd_pop3)

#Export gathered rand for each env
write_csv(rand_gathered_MAT, "Genomics_scripts/Data/rand_gathered_MAT_peakbf2_region.csv")
write_csv(rand_gathered_MAP, "Genomics_scripts/Data/rand_gathered_MAP_peakbf2_region.csv")
write_csv(rand_gathered_CMD, "Genomics_scripts/Data/rand_gathered_CMD_peakbf2_region.csv")





