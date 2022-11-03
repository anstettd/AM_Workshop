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

##Freq Count functions

#MAT
freq_bins_mat <- function(basetime){
  freq_count_calc<-data.frame()
  for (i in 1:12){
    
    if(i==1 || i==7){
      bin_size<-10
    }else if (i==3 || i==4 || i==6 || i==11){
      bin_size<-4
    }else if (i==8){
      bin_size<-2
    }else if(i==2 || i==5 || i==9 || i==10 || i==12){
      bin_size<-5
    }
    
    bin_fraction<-1/bin_size
    bin_step<-seq(0,1,bin_fraction)
    
    for(j in 1:bin_size){
      test_ENV <- basetime %>% filter(Site==i) %>% select (-Site,-Year)
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

#MAP
freq_bins_map <- function(basetime){
  freq_count_calc<-data.frame()
  for (i in 1:12){
    
    if(i==1 || i==2 || i==5 || i==7 || i==10){
      bin_size<-10
    }else if (i==3 || i==4 || i==6 || i==11){
      bin_size<-4
    }else if (i==8){
      bin_size<-2
    }else if(i==9 || i==12){
      bin_size<-5
    }
    
    bin_fraction<-1/bin_size
    bin_step<-seq(0,1,bin_fraction)
    
    for(j in 1:bin_size){
      test_ENV <- basetime %>% filter(Site==i) %>% select (-Site,-Year)
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

#CMD
freq_bins_cmd <- function(basetime){
  freq_count_calc<-data.frame()
  for (i in 1:12){
    
    if(i==1 ||i==5 || i==7 || i==10){
      bin_size<-10
    }else if (i==3 || i==4 || i==6 || i==11){
      bin_size<-4
    }else if (i==8){
      bin_size<-2
    }else if( i==2 || i==9 || i==12){
      bin_size<-5
    }
    
    bin_fraction<-1/bin_size
    bin_step<-seq(0,1,bin_fraction)
    
    for(j in 1:bin_size){
      test_ENV <- basetime %>% filter(Site==i) %>% select (-Site,-Year)
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

#Add Chromsome ID to pop order
pop_order_2 <- data.frame()
pop_order_2 [1,1] <- "chr_shp"
pop_order_2 <- rbind(pop_order_2,pop_order)

#Import pop_order table with regional and V information
pop_order_V <- read_csv("Genomics_scripts/Data/pop_order_V.csv")

#Filter full snp table to remove climate associated SNPs
snp_swiss <-loci_snp %>% filter (!chr_snp %in% as.character(MAT_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(MAP_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(CMD_bf0$chr_snp))

###################################################################################
#Import transformed timeseries frequencies
freq_MAT <- read_csv("Genomics_scripts/Data/freq_MAT_peakbf5.csv")
freq_MAP <- read_csv("Genomics_scripts/Data/freq_MAP_peakbf5.csv")
freq_CMD <- read_csv("Genomics_scripts/Data/freq_CMD_peakbf5.csv")
###################################################################################

#Filter frequency table
freq_MAT_1011 <- freq_MAT %>% filter(Site == 1 & Year==2010 |
                                       Site == 2 & Year==2010 |
                                       Site == 3 & Year==2010 |
                                       Site == 4 & Year==2010 |
                                       Site == 5 & Year==2010 |
                                       Site == 6 & Year==2010 |
                                       Site == 7 & Year==2010 |
                                       Site == 8 & Year==2011 |
                                       Site == 9 & Year==2010 |
                                       Site == 10 & Year==2011 |
                                       Site == 11 & Year==2010 |
                                       Site == 12 & Year==2010)

freq_MAP_1011 <- freq_MAP %>% filter(Site == 1 & Year==2010 |
                                       Site == 2 & Year==2010 |
                                       Site == 3 & Year==2010 |
                                       Site == 4 & Year==2010 |
                                       Site == 5 & Year==2010 |
                                       Site == 6 & Year==2010 |
                                       Site == 7 & Year==2010 |
                                       Site == 8 & Year==2011 |
                                       Site == 9 & Year==2010 |
                                       Site == 10 & Year==2011 |
                                       Site == 11 & Year==2010 |
                                       Site == 12 & Year==2010)

freq_CMD_1011 <- freq_CMD %>% filter(Site == 1 & Year==2010 |
                                       Site == 2 & Year==2010 |
                                       Site == 3 & Year==2010 |
                                       Site == 4 & Year==2010 |
                                       Site == 5 & Year==2010 |
                                       Site == 6 & Year==2010 |
                                       Site == 7 & Year==2010 |
                                       Site == 8 & Year==2011 |
                                       Site == 9 & Year==2010 |
                                       Site == 10 & Year==2011 |
                                       Site == 11 & Year==2010 |
                                       Site == 12 & Year==2010)

#Make frequency count table
freq_count_MAT <- freq_bins_mat(freq_MAT_1011)
freq_count_MAP <- freq_bins_map(freq_MAP_1011)
freq_count_CMD <- freq_bins_cmd(freq_CMD_1011)


######################################################################################################
#Break up freq count into individual vectors
######################################################################################################
freq_count_MAT_1 <- as.vector(freq_count_MAT[,1])
freq_count_MAT_2 <- as.vector(freq_count_MAT[,2])
freq_count_MAT_3 <- as.vector(freq_count_MAT[,3])
freq_count_MAT_4 <- as.vector(freq_count_MAT[,4])
freq_count_MAT_5 <- as.vector(freq_count_MAT[,5])
freq_count_MAT_6 <- as.vector(freq_count_MAT[,6])
freq_count_MAT_7 <- as.vector(freq_count_MAT[,7])
freq_count_MAT_8 <- as.vector(freq_count_MAT[,8])
freq_count_MAT_9 <- as.vector(freq_count_MAT[,9])
freq_count_MAT_10 <- as.vector(freq_count_MAT[,10])
freq_count_MAT_11 <- as.vector(freq_count_MAT[,11])
freq_count_MAT_12 <- as.vector(freq_count_MAT[,12])


freq_count_MAT_1<-freq_count_MAT_1[!is.na(freq_count_MAT_1)]
freq_count_MAT_2<-freq_count_MAT_2[!is.na(freq_count_MAT_2)]
freq_count_MAT_3<-freq_count_MAT_3[!is.na(freq_count_MAT_3)]
freq_count_MAT_4<-freq_count_MAT_4[!is.na(freq_count_MAT_4)]
freq_count_MAT_5<-freq_count_MAT_5[!is.na(freq_count_MAT_5)]
freq_count_MAT_6<-freq_count_MAT_6[!is.na(freq_count_MAT_6)]
freq_count_MAT_7<-freq_count_MAT_7[!is.na(freq_count_MAT_7)]
freq_count_MAT_8<-freq_count_MAT_8[!is.na(freq_count_MAT_8)]
freq_count_MAT_9<-freq_count_MAT_9[!is.na(freq_count_MAT_9)]
freq_count_MAT_10<-freq_count_MAT_10[!is.na(freq_count_MAT_10)]
freq_count_MAT_11<-freq_count_MAT_11[!is.na(freq_count_MAT_11)]
freq_count_MAT_12<-freq_count_MAT_12[!is.na(freq_count_MAT_12)]

freq_count_MAP_1 <- as.vector(freq_count_MAP[,1])
freq_count_MAP_2 <- as.vector(freq_count_MAP[,2])
freq_count_MAP_3 <- as.vector(freq_count_MAP[,3])
freq_count_MAP_4 <- as.vector(freq_count_MAP[,4])
freq_count_MAP_5 <- as.vector(freq_count_MAP[,5])
freq_count_MAP_6 <- as.vector(freq_count_MAP[,6])
freq_count_MAP_7 <- as.vector(freq_count_MAP[,7])
freq_count_MAP_8 <- as.vector(freq_count_MAP[,8])
freq_count_MAP_9 <- as.vector(freq_count_MAP[,9])
freq_count_MAP_10 <- as.vector(freq_count_MAP[,10])
freq_count_MAP_11 <- as.vector(freq_count_MAP[,11])
freq_count_MAP_12 <- as.vector(freq_count_MAP[,12])

freq_count_MAP_1<-freq_count_MAP_1[!is.na(freq_count_MAP_1)]
freq_count_MAP_2<-freq_count_MAP_2[!is.na(freq_count_MAP_2)]
freq_count_MAP_3<-freq_count_MAP_3[!is.na(freq_count_MAP_3)]
freq_count_MAP_4<-freq_count_MAP_4[!is.na(freq_count_MAP_4)]
freq_count_MAP_5<-freq_count_MAP_5[!is.na(freq_count_MAP_5)]
freq_count_MAP_6<-freq_count_MAP_6[!is.na(freq_count_MAP_6)]
freq_count_MAP_7<-freq_count_MAP_7[!is.na(freq_count_MAP_7)]
freq_count_MAP_8<-freq_count_MAP_8[!is.na(freq_count_MAP_8)]
freq_count_MAP_9<-freq_count_MAP_9[!is.na(freq_count_MAP_9)]
freq_count_MAP_10<-freq_count_MAP_10[!is.na(freq_count_MAP_10)]
freq_count_MAP_11<-freq_count_MAP_11[!is.na(freq_count_MAP_11)]
freq_count_MAP_12<-freq_count_MAP_12[!is.na(freq_count_MAP_12)]


freq_count_CMD_1 <- as.vector(freq_count_CMD[,1])
freq_count_CMD_2 <- as.vector(freq_count_CMD[,2])
freq_count_CMD_3 <- as.vector(freq_count_CMD[,3])
freq_count_CMD_4 <- as.vector(freq_count_CMD[,4])
freq_count_CMD_5 <- as.vector(freq_count_CMD[,5])
freq_count_CMD_6 <- as.vector(freq_count_CMD[,6])
freq_count_CMD_7 <- as.vector(freq_count_CMD[,7])
freq_count_CMD_8 <- as.vector(freq_count_CMD[,8])
freq_count_CMD_9 <- as.vector(freq_count_CMD[,9])
freq_count_CMD_10 <- as.vector(freq_count_CMD[,10])
freq_count_CMD_11 <- as.vector(freq_count_CMD[,11])
freq_count_CMD_12 <- as.vector(freq_count_CMD[,12])

freq_count_CMD_1<-freq_count_CMD_1[!is.na(freq_count_CMD_1)]
freq_count_CMD_2<-freq_count_CMD_2[!is.na(freq_count_CMD_2)]
freq_count_CMD_3<-freq_count_CMD_3[!is.na(freq_count_CMD_3)]
freq_count_CMD_4<-freq_count_CMD_4[!is.na(freq_count_CMD_4)]
freq_count_CMD_5<-freq_count_CMD_5[!is.na(freq_count_CMD_5)]
freq_count_CMD_6<-freq_count_CMD_6[!is.na(freq_count_CMD_6)]
freq_count_CMD_7<-freq_count_CMD_7[!is.na(freq_count_CMD_7)]
freq_count_CMD_8<-freq_count_CMD_8[!is.na(freq_count_CMD_8)]
freq_count_CMD_9<-freq_count_CMD_9[!is.na(freq_count_CMD_9)]
freq_count_CMD_10<-freq_count_CMD_10[!is.na(freq_count_CMD_10)]
freq_count_CMD_11<-freq_count_CMD_11[!is.na(freq_count_CMD_11)]
freq_count_CMD_12<-freq_count_CMD_12[!is.na(freq_count_CMD_12)]


#######################################################################################################
#Make SNP table for neutral loci

#Filter pop_order_V to call for 2010 regional basetime V column names
pop_V1 <- pop_order_V %>% filter(Year==2010 & Pop==1)
pop_V2 <- pop_order_V %>% filter(Year==2010 & Pop==2)
pop_V3 <- pop_order_V %>% filter(Year==2010 & Pop==3)
pop_V4 <- pop_order_V %>% filter(Year==2010 & Pop==4)
pop_V5 <- pop_order_V %>% filter(Year==2010 & Pop==5)
pop_V6 <- pop_order_V %>% filter(Year==2010 & Pop==6)

pop_V7 <- pop_order_V %>% filter(Year==2010 & Pop==7)
pop_V8 <- pop_order_V %>% filter(Year==2011 & Pop==8)
pop_V9 <- pop_order_V %>% filter(Year==2010 & Pop==9)
pop_V10 <- pop_order_V %>% filter(Year==2011 & Pop==10)
pop_V11 <- pop_order_V %>% filter(Year==2010 & Pop==11)
pop_V12 <- pop_order_V %>% filter(Year==2010 & Pop==11)


#Filter Baseline A and B numbers for 3 basetime regions
loci_base_1 <- snp_swiss %>% select(chr_snp,all_of(pop_V1$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_2 <- snp_swiss %>% select(chr_snp,all_of(pop_V2$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_3 <- snp_swiss %>% select(chr_snp,all_of(pop_V3$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_4 <- snp_swiss %>% select(chr_snp,all_of(pop_V4$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_5 <- snp_swiss %>% select(chr_snp,all_of(pop_V5$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_6 <- snp_swiss %>% select(chr_snp,all_of(pop_V6$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))

loci_base_7 <- snp_swiss %>% select(chr_snp,all_of(pop_V7$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_8 <- snp_swiss %>% select(chr_snp,all_of(pop_V8$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_9 <- snp_swiss %>% select(chr_snp,all_of(pop_V9$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_10 <- snp_swiss %>% select(chr_snp,all_of(pop_V10$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_11 <- snp_swiss %>% select(chr_snp,all_of(pop_V11$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_12 <- snp_swiss %>% select(chr_snp,all_of(pop_V12$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))


#Calculate frequency for SNP A for 12 basetime pops
#Gives all basetime frequencies per region
p1A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(region_sum)) %>% select(chr_snp,snpA)
p2A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(region_sum)) %>% select(chr_snp,snpA)
p3A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(region_sum)) %>% select(chr_snp,snpA)
p4A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(region_sum)) %>% select(chr_snp,snpA)
p5A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(region_sum)) %>% select(chr_snp,snpA)
p6A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(region_sum)) %>% select(chr_snp,snpA)

p7A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(region_sum)) %>% select(chr_snp,snpA)
p8A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(region_sum)) %>% select(chr_snp,snpA)
p9A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(region_sum)) %>% select(chr_snp,snpA)
p10A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(region_sum)) %>% select(chr_snp,snpA)
p11A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(region_sum)) %>% select(chr_snp,snpA)
p12A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(region_sum)) %>% select(chr_snp,snpA)



#Filter full snp table to remove climate associated SNPs
rand_slope_mat_out<-data.frame()
rand_slope_map_out<-data.frame()
rand_slope_cmd_out<-data.frame()


#######################################################################################################
# Stratified random sampling SNPs for a given population
for (seed_num in 1:1000){
set.seed(seed_num)
#mat
rand_mat_pop1 <- large(p1A,freq_count_MAT_1)
rand_mat_pop2 <- large(p2A,freq_count_MAT_2)
rand_mat_pop3 <- large(p3A,freq_count_MAT_3)
rand_mat_pop4 <- large(p4A,freq_count_MAT_4)
rand_mat_pop5 <- large(p5A,freq_count_MAT_5)
rand_mat_pop6 <- large(p6A,freq_count_MAT_6)
rand_mat_pop7 <- large(p7A,freq_count_MAT_7)
rand_mat_pop8 <- large(p8A,freq_count_MAT_8)
rand_mat_pop9 <- large(p9A,freq_count_MAT_9)
rand_mat_pop10 <- large(p10A,freq_count_MAT_10)
rand_mat_pop11 <- large(p11A,freq_count_MAT_11)
rand_mat_pop12 <- large(p12A,freq_count_MAT_12)

rand_map_pop1 <- large(p1A,freq_count_MAP_1)
rand_map_pop2 <- large(p2A,freq_count_MAP_2)
rand_map_pop3 <- large(p3A,freq_count_MAP_3)
rand_map_pop4 <- large(p4A,freq_count_MAP_4)
rand_map_pop5 <- large(p5A,freq_count_MAP_5)
rand_map_pop6 <- large(p6A,freq_count_MAP_6)
rand_map_pop7 <- large(p7A,freq_count_MAP_7)
rand_map_pop8 <- large(p8A,freq_count_MAP_8)
rand_map_pop9 <- large(p9A,freq_count_MAP_9)
rand_map_pop10 <- large(p10A,freq_count_MAP_10)
rand_map_pop11 <- large(p11A,freq_count_MAP_11)
rand_map_pop12 <- large(p12A,freq_count_MAP_12)

rand_cmd_pop1 <- large(p1A,freq_count_CMD_1)
rand_cmd_pop2 <- large(p2A,freq_count_CMD_2)
rand_cmd_pop3 <- large(p3A,freq_count_CMD_3)
rand_cmd_pop4 <- large(p4A,freq_count_CMD_4)
rand_cmd_pop5 <- large(p5A,freq_count_CMD_5)
rand_cmd_pop6 <- large(p6A,freq_count_CMD_6)
rand_cmd_pop7 <- large(p7A,freq_count_CMD_7)
rand_cmd_pop8 <- large(p8A,freq_count_CMD_8)
rand_cmd_pop9 <- large(p9A,freq_count_CMD_9)
rand_cmd_pop10 <- large(p10A,freq_count_CMD_10)
rand_cmd_pop11 <- large(p11A,freq_count_CMD_11)
rand_cmd_pop12 <- large(p12A,freq_count_CMD_12)

#Get Full timeseires for each randomly selected neutral location
rand_time_AB_mat_1 <-loci_snp %>% filter (chr_snp %in% as.character(rand_mat_pop1$chr_snp)) #mat
rand_time_AB_mat_2 <-loci_snp %>% filter (chr_snp %in% as.character(rand_mat_pop2$chr_snp)) #mat
rand_time_AB_mat_3 <-loci_snp %>% filter (chr_snp %in% as.character(rand_mat_pop3$chr_snp)) #mat
rand_time_AB_mat_4 <-loci_snp %>% filter (chr_snp %in% as.character(rand_mat_pop4$chr_snp)) #mat
rand_time_AB_mat_5 <-loci_snp %>% filter (chr_snp %in% as.character(rand_mat_pop5$chr_snp)) #mat
rand_time_AB_mat_6 <-loci_snp %>% filter (chr_snp %in% as.character(rand_mat_pop6$chr_snp)) #mat
rand_time_AB_mat_7 <-loci_snp %>% filter (chr_snp %in% as.character(rand_mat_pop7$chr_snp)) #mat
rand_time_AB_mat_8 <-loci_snp %>% filter (chr_snp %in% as.character(rand_mat_pop8$chr_snp)) #mat
rand_time_AB_mat_9 <-loci_snp %>% filter (chr_snp %in% as.character(rand_mat_pop9$chr_snp)) #mat
rand_time_AB_mat_10 <-loci_snp %>% filter (chr_snp %in% as.character(rand_mat_pop10$chr_snp)) #mat
rand_time_AB_mat_11 <-loci_snp %>% filter (chr_snp %in% as.character(rand_mat_pop11$chr_snp)) #mat
rand_time_AB_mat_12 <-loci_snp %>% filter (chr_snp %in% as.character(rand_mat_pop12$chr_snp)) #mat

rand_time_AB_map_1 <-loci_snp %>% filter (chr_snp %in% as.character(rand_map_pop1$chr_snp)) #map
rand_time_AB_map_2 <-loci_snp %>% filter (chr_snp %in% as.character(rand_map_pop2$chr_snp)) #map
rand_time_AB_map_3 <-loci_snp %>% filter (chr_snp %in% as.character(rand_map_pop3$chr_snp)) #map
rand_time_AB_map_4 <-loci_snp %>% filter (chr_snp %in% as.character(rand_map_pop4$chr_snp)) #map
rand_time_AB_map_5 <-loci_snp %>% filter (chr_snp %in% as.character(rand_map_pop5$chr_snp)) #map
rand_time_AB_map_6 <-loci_snp %>% filter (chr_snp %in% as.character(rand_map_pop6$chr_snp)) #map
rand_time_AB_map_7 <-loci_snp %>% filter (chr_snp %in% as.character(rand_map_pop7$chr_snp)) #map
rand_time_AB_map_8 <-loci_snp %>% filter (chr_snp %in% as.character(rand_map_pop8$chr_snp)) #map
rand_time_AB_map_9 <-loci_snp %>% filter (chr_snp %in% as.character(rand_map_pop9$chr_snp)) #map
rand_time_AB_map_10 <-loci_snp %>% filter (chr_snp %in% as.character(rand_map_pop10$chr_snp)) #map
rand_time_AB_map_11 <-loci_snp %>% filter (chr_snp %in% as.character(rand_map_pop11$chr_snp)) #map
rand_time_AB_map_12 <-loci_snp %>% filter (chr_snp %in% as.character(rand_map_pop12$chr_snp)) #map

rand_time_AB_cmd_1 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop1$chr_snp)) #cmd
rand_time_AB_cmd_2 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop2$chr_snp)) #cmd
rand_time_AB_cmd_3 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop3$chr_snp)) #cmd
rand_time_AB_cmd_4 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop4$chr_snp)) #cmd
rand_time_AB_cmd_5 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop5$chr_snp)) #cmd
rand_time_AB_cmd_6 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop6$chr_snp)) #cmd
rand_time_AB_cmd_7 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop7$chr_snp)) #cmd
rand_time_AB_cmd_8 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop8$chr_snp)) #cmd
rand_time_AB_cmd_9 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop9$chr_snp)) #cmd
rand_time_AB_cmd_10 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop10$chr_snp)) #cmd
rand_time_AB_cmd_11 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop11$chr_snp)) #cmd
rand_time_AB_cmd_12 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop12$chr_snp)) #cmd

# Calc SNP A and transpose
rand_mat_pop1 <- prop_A(rand_time_AB_mat_1,pop_order_2)
rand_mat_pop2 <- prop_A(rand_time_AB_mat_2,pop_order_2)
rand_mat_pop3 <- prop_A(rand_time_AB_mat_3,pop_order_2)
rand_mat_pop4 <- prop_A(rand_time_AB_mat_4,pop_order_2)
rand_mat_pop5 <- prop_A(rand_time_AB_mat_5,pop_order_2)
rand_mat_pop6 <- prop_A(rand_time_AB_mat_6,pop_order_2)
rand_mat_pop7 <- prop_A(rand_time_AB_mat_7,pop_order_2)
rand_mat_pop8 <- prop_A(rand_time_AB_mat_8,pop_order_2)
rand_mat_pop9 <- prop_A(rand_time_AB_mat_9,pop_order_2)
rand_mat_pop10 <- prop_A(rand_time_AB_mat_10,pop_order_2)
rand_mat_pop11 <- prop_A(rand_time_AB_mat_11,pop_order_2)
rand_mat_pop12 <- prop_A(rand_time_AB_mat_12,pop_order_2)

rand_map_pop1 <- prop_A(rand_time_AB_map_1,pop_order_2)
rand_map_pop2 <- prop_A(rand_time_AB_map_2,pop_order_2)
rand_map_pop3 <- prop_A(rand_time_AB_map_3,pop_order_2)
rand_map_pop4 <- prop_A(rand_time_AB_map_4,pop_order_2)
rand_map_pop5 <- prop_A(rand_time_AB_map_5,pop_order_2)
rand_map_pop6 <- prop_A(rand_time_AB_map_6,pop_order_2)
rand_map_pop7 <- prop_A(rand_time_AB_map_7,pop_order_2)
rand_map_pop8 <- prop_A(rand_time_AB_map_8,pop_order_2)
rand_map_pop9 <- prop_A(rand_time_AB_map_9,pop_order_2)
rand_map_pop10 <- prop_A(rand_time_AB_map_10,pop_order_2)
rand_map_pop11 <- prop_A(rand_time_AB_map_11,pop_order_2)
rand_map_pop12 <- prop_A(rand_time_AB_map_12,pop_order_2)

rand_cmd_pop1 <- prop_A(rand_time_AB_cmd_1,pop_order_2)
rand_cmd_pop2 <- prop_A(rand_time_AB_cmd_2,pop_order_2)
rand_cmd_pop3 <- prop_A(rand_time_AB_cmd_3,pop_order_2)
rand_cmd_pop4 <- prop_A(rand_time_AB_cmd_4,pop_order_2)
rand_cmd_pop5 <- prop_A(rand_time_AB_cmd_5,pop_order_2)
rand_cmd_pop6 <- prop_A(rand_time_AB_cmd_6,pop_order_2)
rand_cmd_pop7 <- prop_A(rand_time_AB_cmd_7,pop_order_2)
rand_cmd_pop8 <- prop_A(rand_time_AB_cmd_8,pop_order_2)
rand_cmd_pop9 <- prop_A(rand_time_AB_cmd_9,pop_order_2)
rand_cmd_pop10 <- prop_A(rand_time_AB_cmd_10,pop_order_2)
rand_cmd_pop11 <- prop_A(rand_time_AB_cmd_11,pop_order_2)
rand_cmd_pop12 <- prop_A(rand_time_AB_cmd_12,pop_order_2)

# Filter for site (its currently duplicated since timeseries filtering was on 12 pop datatset)
rand_mat_pop1 <- rand_mat_pop1 %>% filter(Site==1)
rand_mat_pop2 <- rand_mat_pop2 %>% filter(Site==2)
rand_mat_pop3 <- rand_mat_pop3 %>% filter(Site==3)
rand_mat_pop4 <- rand_mat_pop4 %>% filter(Site==4)
rand_mat_pop5 <- rand_mat_pop5 %>% filter(Site==5)
rand_mat_pop6 <- rand_mat_pop6 %>% filter(Site==6)
rand_mat_pop7 <- rand_mat_pop7 %>% filter(Site==7)
rand_mat_pop8 <- rand_mat_pop8 %>% filter(Site==8)
rand_mat_pop9 <- rand_mat_pop9 %>% filter(Site==9)
rand_mat_pop10 <- rand_mat_pop10 %>% filter(Site==10)
rand_mat_pop11 <- rand_mat_pop11 %>% filter(Site==11)
rand_mat_pop12 <- rand_mat_pop12 %>% filter(Site==12)

rand_map_pop1 <- rand_map_pop1 %>% filter(Site==1)
rand_map_pop2 <- rand_map_pop2 %>% filter(Site==2)
rand_map_pop3 <- rand_map_pop3 %>% filter(Site==3)
rand_map_pop4 <- rand_map_pop4 %>% filter(Site==4)
rand_map_pop5 <- rand_map_pop5 %>% filter(Site==5)
rand_map_pop6 <- rand_map_pop6 %>% filter(Site==6)
rand_map_pop7 <- rand_map_pop7 %>% filter(Site==7)
rand_map_pop8 <- rand_map_pop8 %>% filter(Site==8)
rand_map_pop9 <- rand_map_pop9 %>% filter(Site==9)
rand_map_pop10 <- rand_map_pop10 %>% filter(Site==10)
rand_map_pop11 <- rand_map_pop11 %>% filter(Site==11)
rand_map_pop12 <- rand_map_pop12 %>% filter(Site==12)

rand_cmd_pop1 <- rand_cmd_pop1 %>% filter(Site==1)
rand_cmd_pop2 <- rand_cmd_pop2 %>% filter(Site==2)
rand_cmd_pop3 <- rand_cmd_pop3 %>% filter(Site==3)
rand_cmd_pop4 <- rand_cmd_pop4 %>% filter(Site==4)
rand_cmd_pop5 <- rand_cmd_pop5 %>% filter(Site==5)
rand_cmd_pop6 <- rand_cmd_pop6 %>% filter(Site==6)
rand_cmd_pop7 <- rand_cmd_pop7 %>% filter(Site==7)
rand_cmd_pop8 <- rand_cmd_pop8 %>% filter(Site==8)
rand_cmd_pop9 <- rand_cmd_pop9 %>% filter(Site==9)
rand_cmd_pop10 <- rand_cmd_pop10 %>% filter(Site==10)
rand_cmd_pop11 <- rand_cmd_pop11 %>% filter(Site==11)
rand_cmd_pop12 <- rand_cmd_pop12 %>% filter(Site==12)

#Get slopes
rand_slope_mat_pop1 <- glm_slopes(rand_mat_pop1)
rand_slope_mat_pop2 <- glm_slopes(rand_mat_pop2)
rand_slope_mat_pop3 <- glm_slopes(rand_mat_pop3)
rand_slope_mat_pop4 <- glm_slopes(rand_mat_pop4)
rand_slope_mat_pop5 <- glm_slopes(rand_mat_pop5)
rand_slope_mat_pop6 <- glm_slopes(rand_mat_pop6)
rand_slope_mat_pop7 <- glm_slopes(rand_mat_pop7)
rand_slope_mat_pop8 <- glm_slopes(rand_mat_pop8)
rand_slope_mat_pop9 <- glm_slopes(rand_mat_pop9)
rand_slope_mat_pop10 <- glm_slopes(rand_mat_pop10)
rand_slope_mat_pop11 <- glm_slopes(rand_mat_pop11)
rand_slope_mat_pop12 <- glm_slopes(rand_mat_pop12)

rand_slope_map_pop1 <- glm_slopes(rand_map_pop1)
rand_slope_map_pop2 <- glm_slopes(rand_map_pop2)
rand_slope_map_pop3 <- glm_slopes(rand_map_pop3)
rand_slope_map_pop4 <- glm_slopes(rand_map_pop4)
rand_slope_map_pop5 <- glm_slopes(rand_map_pop5)
rand_slope_map_pop6 <- glm_slopes(rand_map_pop6)
rand_slope_map_pop7 <- glm_slopes(rand_map_pop7)
rand_slope_map_pop8 <- glm_slopes(rand_map_pop8)
rand_slope_map_pop9 <- glm_slopes(rand_map_pop9)
rand_slope_map_pop10 <- glm_slopes(rand_map_pop10)
rand_slope_map_pop11 <- glm_slopes(rand_map_pop11)
rand_slope_map_pop12 <- glm_slopes(rand_map_pop12)

rand_slope_cmd_pop1 <- glm_slopes(rand_cmd_pop1)
rand_slope_cmd_pop2 <- glm_slopes(rand_cmd_pop2)
rand_slope_cmd_pop3 <- glm_slopes(rand_cmd_pop3)
rand_slope_cmd_pop4 <- glm_slopes(rand_cmd_pop4)
rand_slope_cmd_pop5 <- glm_slopes(rand_cmd_pop5)
rand_slope_cmd_pop6 <- glm_slopes(rand_cmd_pop6)
rand_slope_cmd_pop7 <- glm_slopes(rand_cmd_pop7)
rand_slope_cmd_pop8 <- glm_slopes(rand_cmd_pop8)
rand_slope_cmd_pop9 <- glm_slopes(rand_cmd_pop9)
rand_slope_cmd_pop10 <- glm_slopes(rand_cmd_pop10)
rand_slope_cmd_pop11 <- glm_slopes(rand_cmd_pop11)
rand_slope_cmd_pop12 <- glm_slopes(rand_cmd_pop12)


#Bind populations for each env
rand_slope_mat <- rbind(rand_slope_mat_pop1,
                      rand_slope_mat_pop2,
                      rand_slope_mat_pop3,
                      rand_slope_mat_pop4,
                      rand_slope_mat_pop5,
                      rand_slope_mat_pop6,
                      rand_slope_mat_pop7,
                      rand_slope_mat_pop8,
                      rand_slope_mat_pop9,
                      rand_slope_mat_pop10,
                      rand_slope_mat_pop11,
                      rand_slope_mat_pop12)

rand_slope_map <- rbind(rand_slope_map_pop1,
                      rand_slope_map_pop2,
                      rand_slope_map_pop3,
                      rand_slope_map_pop4,
                      rand_slope_map_pop5,
                      rand_slope_map_pop6,
                      rand_slope_map_pop7,
                      rand_slope_map_pop8,
                      rand_slope_map_pop9,
                      rand_slope_map_pop10,
                      rand_slope_map_pop11,
                      rand_slope_map_pop12)

rand_slope_cmd <- rbind(rand_slope_cmd_pop1,
                      rand_slope_cmd_pop2,
                      rand_slope_cmd_pop3,
                      rand_slope_cmd_pop4,
                      rand_slope_cmd_pop5,
                      rand_slope_cmd_pop6,
                      rand_slope_cmd_pop7,
                      rand_slope_cmd_pop8,
                      rand_slope_cmd_pop9,
                      rand_slope_cmd_pop10,
                      rand_slope_cmd_pop11,
                      rand_slope_cmd_pop12)
rand_slope_mat <- rand_slope_mat %>% mutate (Seed_ID = seed_num)
rand_slope_map <- rand_slope_map %>% mutate (Seed_ID = seed_num)
rand_slope_cmd <- rand_slope_cmd %>% mutate (Seed_ID = seed_num)
#Export each joint df

rand_slope_mat_out<-rbind(rand_slope_mat_out,rand_slope_mat)
rand_slope_map_out<-rbind(rand_slope_map_out,rand_slope_map)
rand_slope_cmd_out<-rbind(rand_slope_cmd_out,rand_slope_cmd)

print(seed_num)

}

#Save large files in folder outside of github
setwd("~/Dropbox/AM_Workshop/Large_files")

write_csv(rand_slope_mat_out, "rand_slope_mat_multi_peakbf5.csv")
write_csv(rand_slope_map_out, "rand_slope_map_multi_peakbf5.csv")
write_csv(rand_slope_cmd_out, "rand_slope_cmd_multi_peakbf5.csv")

setwd("~/Dropbox/AM_Workshop/AM_Workshop/")








rand_slope_cmd_out<-data.frame()

#CMD ONLY

#######################################################################################################
# Stratified random sampling SNPs for a given population
for (seed_num in 1:1000){
  set.seed(seed_num)
  
  rand_cmd_pop1 <- large(p1A,freq_count_CMD_1)
  rand_cmd_pop2 <- large(p2A,freq_count_CMD_2)
  rand_cmd_pop3 <- large(p3A,freq_count_CMD_3)
  rand_cmd_pop4 <- large(p4A,freq_count_CMD_4)
  rand_cmd_pop5 <- large(p5A,freq_count_CMD_5)
  rand_cmd_pop6 <- large(p6A,freq_count_CMD_6)
  rand_cmd_pop7 <- large(p7A,freq_count_CMD_7)
  rand_cmd_pop8 <- large(p8A,freq_count_CMD_8)
  rand_cmd_pop9 <- large(p9A,freq_count_CMD_9)
  rand_cmd_pop10 <- large(p10A,freq_count_CMD_10)
  rand_cmd_pop11 <- large(p11A,freq_count_CMD_11)
  rand_cmd_pop12 <- large(p12A,freq_count_CMD_12)
  
  #Get Full timeseires for each randomly selected neutral location
  rand_time_AB_cmd_1 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop1$chr_snp)) #cmd
  rand_time_AB_cmd_2 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop2$chr_snp)) #cmd
  rand_time_AB_cmd_3 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop3$chr_snp)) #cmd
  rand_time_AB_cmd_4 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop4$chr_snp)) #cmd
  rand_time_AB_cmd_5 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop5$chr_snp)) #cmd
  rand_time_AB_cmd_6 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop6$chr_snp)) #cmd
  rand_time_AB_cmd_7 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop7$chr_snp)) #cmd
  rand_time_AB_cmd_8 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop8$chr_snp)) #cmd
  rand_time_AB_cmd_9 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop9$chr_snp)) #cmd
  rand_time_AB_cmd_10 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop10$chr_snp)) #cmd
  rand_time_AB_cmd_11 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop11$chr_snp)) #cmd
  rand_time_AB_cmd_12 <-loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_pop12$chr_snp)) #cmd
  
  # Calc SNP A and transpose
  rand_cmd_pop1 <- prop_A(rand_time_AB_cmd_1,pop_order_2)
  rand_cmd_pop2 <- prop_A(rand_time_AB_cmd_2,pop_order_2)
  rand_cmd_pop3 <- prop_A(rand_time_AB_cmd_3,pop_order_2)
  rand_cmd_pop4 <- prop_A(rand_time_AB_cmd_4,pop_order_2)
  rand_cmd_pop5 <- prop_A(rand_time_AB_cmd_5,pop_order_2)
  rand_cmd_pop6 <- prop_A(rand_time_AB_cmd_6,pop_order_2)
  rand_cmd_pop7 <- prop_A(rand_time_AB_cmd_7,pop_order_2)
  rand_cmd_pop8 <- prop_A(rand_time_AB_cmd_8,pop_order_2)
  rand_cmd_pop9 <- prop_A(rand_time_AB_cmd_9,pop_order_2)
  rand_cmd_pop10 <- prop_A(rand_time_AB_cmd_10,pop_order_2)
  rand_cmd_pop11 <- prop_A(rand_time_AB_cmd_11,pop_order_2)
  rand_cmd_pop12 <- prop_A(rand_time_AB_cmd_12,pop_order_2)
  
  # Filter for site (its currently duplicated since timeseries filtering was on 12 pop datatset)
  rand_cmd_pop1 <- rand_cmd_pop1 %>% filter(Site==1)
  rand_cmd_pop2 <- rand_cmd_pop2 %>% filter(Site==2)
  rand_cmd_pop3 <- rand_cmd_pop3 %>% filter(Site==3)
  rand_cmd_pop4 <- rand_cmd_pop4 %>% filter(Site==4)
  rand_cmd_pop5 <- rand_cmd_pop5 %>% filter(Site==5)
  rand_cmd_pop6 <- rand_cmd_pop6 %>% filter(Site==6)
  rand_cmd_pop7 <- rand_cmd_pop7 %>% filter(Site==7)
  rand_cmd_pop8 <- rand_cmd_pop8 %>% filter(Site==8)
  rand_cmd_pop9 <- rand_cmd_pop9 %>% filter(Site==9)
  rand_cmd_pop10 <- rand_cmd_pop10 %>% filter(Site==10)
  rand_cmd_pop11 <- rand_cmd_pop11 %>% filter(Site==11)
  rand_cmd_pop12 <- rand_cmd_pop12 %>% filter(Site==12)
  
  #Get slopes
  rand_slope_cmd_pop1 <- glm_slopes(rand_cmd_pop1)
  rand_slope_cmd_pop2 <- glm_slopes(rand_cmd_pop2)
  rand_slope_cmd_pop3 <- glm_slopes(rand_cmd_pop3)
  rand_slope_cmd_pop4 <- glm_slopes(rand_cmd_pop4)
  rand_slope_cmd_pop5 <- glm_slopes(rand_cmd_pop5)
  rand_slope_cmd_pop6 <- glm_slopes(rand_cmd_pop6)
  rand_slope_cmd_pop7 <- glm_slopes(rand_cmd_pop7)
  rand_slope_cmd_pop8 <- glm_slopes(rand_cmd_pop8)
  rand_slope_cmd_pop9 <- glm_slopes(rand_cmd_pop9)
  rand_slope_cmd_pop10 <- glm_slopes(rand_cmd_pop10)
  rand_slope_cmd_pop11 <- glm_slopes(rand_cmd_pop11)
  rand_slope_cmd_pop12 <- glm_slopes(rand_cmd_pop12)
  
  
  #Bind populations for each env
  rand_slope_cmd <- rbind(rand_slope_cmd_pop1,
                          rand_slope_cmd_pop2,
                          rand_slope_cmd_pop3,
                          rand_slope_cmd_pop4,
                          rand_slope_cmd_pop5,
                          rand_slope_cmd_pop6,
                          rand_slope_cmd_pop7,
                          rand_slope_cmd_pop8,
                          rand_slope_cmd_pop9,
                          rand_slope_cmd_pop10,
                          rand_slope_cmd_pop11,
                          rand_slope_cmd_pop12)
  rand_slope_cmd <- rand_slope_cmd %>% mutate (Seed_ID = seed_num)
  #Export each joint df
  rand_slope_cmd_out<-rbind(rand_slope_cmd_out,rand_slope_cmd)
  
  print(seed_num)
  
}

