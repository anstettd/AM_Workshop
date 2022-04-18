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
      list_env_p1_n<- sample.int(dim(p1_n)[1],freq_count[i])
      
    }else{
      p1_n<-pop_snp %>% filter(snpA>=bin_step[i] & snpA < bin_step[i+1])
      list_env_p1_n<- sample.int(dim(p1_n)[1],freq_count[i])
      
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
  for (i in 1:12){
    
    if(i==1 || i==2 || i==5 || i==7 || i==10 || i==12){
      bin_size<-10
    }else if (i==3 || i==4 || i==6 || i==11){
      bin_size<-4
    }else if (i==8){
      bin_size<-2
    }else if(i==9){
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



###################################################################################

## Manupulate entire SNP datatset for timeseries

#Import timeseries counts for A and B
#snp1_time <- read_csv("Genomics_scripts/Data/snp1_filter.csv")
#snp2_time <- read_csv("Genomics_scripts/Data/snp2_filter.csv")
#snp5_time <- read_csv("Genomics_scripts/Data/snp5_filter.csv")

#Import BF>3 baseline SNPs
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

#Filter full snp table to remove climate associated SNPs
snp_swiss <-loci_snp %>% filter (!chr_snp %in% as.character(MAT_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(MAP_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(CMD_bf0$chr_snp))

###################################################################################
#Import transformed timeseries frequencies
freq_MAT <- read_csv("Genomics_scripts/Data/freq_MAT.csv")
freq_MAP <- read_csv("Genomics_scripts/Data/freq_MAP.csv")
freq_CMD <- read_csv("Genomics_scripts/Data/freq_CMD.csv")
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
freq_count_MAT <- freq_bins(freq_MAT_1011)
freq_count_MAP <- freq_bins(freq_MAP_1011)
freq_count_CMD <- freq_bins(freq_CMD_1011)


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

#Select 2011 values for pop 8 and 10
#For the reset select 2010 values

#Filter Baseline A and B numbers for 12 timebase pops
loci_base_p1 <- snp_swiss %>% select(chr_snp,V33,V34)
loci_base_p2 <- snp_swiss %>% select(chr_snp,V43,V44)
loci_base_p3 <- snp_swiss %>% select(chr_snp,V53,V54)
loci_base_p4 <- snp_swiss %>% select(chr_snp,V67,V68)
loci_base_p5 <- snp_swiss %>% select(chr_snp,V77,V78)
loci_base_p6 <- snp_swiss %>% select(chr_snp,V87,V88)

loci_base_p7 <- snp_swiss %>% select(chr_snp,V97,V98)
loci_base_p8 <- snp_swiss %>% select(chr_snp,V107,V108)
loci_base_p9 <- snp_swiss %>% select(chr_snp,V115,V116)
loci_base_p10 <- snp_swiss %>% select(chr_snp,V3,V4)
loci_base_p11 <- snp_swiss %>% select(chr_snp,V13,V14)
loci_base_p12 <- snp_swiss %>% select(chr_snp,V23,V24)

#Calculate frequency for SNP A for 12 timebase pops
p1A <- loci_base_p1 %>% mutate(snpA=V33/(V33+V34)) %>% select(-V33,-V34)
p2A <- loci_base_p2 %>% mutate(snpA=V43/(V43+V44)) %>% select(-V43,-V44)
p3A <- loci_base_p3 %>% mutate(snpA=V53/(V53+V54)) %>% select(-V53,-V54)
p4A <- loci_base_p4 %>% mutate(snpA=V67/(V67+V68)) %>% select(-V67,-V68)
p5A <- loci_base_p5 %>% mutate(snpA=V77/(V77+V78)) %>% select(-V77,-V78)
p6A <- loci_base_p6 %>% mutate(snpA=V87/(V87+V88)) %>% select(-V87,-V88)

p7A <- loci_base_p7 %>% mutate(snpA=V97/(V97+V98)) %>% select(-V97,-V98)
p8A <- loci_base_p8 %>% mutate(snpA=V107/(V107+V108)) %>% select(-V107,-V108)
p9A <- loci_base_p9 %>% mutate(snpA=V115/(V115+V116)) %>% select(-V115,-V116)
p10A <- loci_base_p10 %>% mutate(snpA=V3/(V3+V4)) %>% select(-V3,-V4)
p11A <- loci_base_p11 %>% mutate(snpA=V13/(V13+V14)) %>% select(-V13,-V14)
p12A <- loci_base_p12 %>% mutate(snpA=V23/(V23+V24)) %>% select(-V23,-V24)


#######################################################################################################
# Stratified random sampling SNPs for a given population

set.seed(1)
#MAT
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

rand_MAP_pop1 <- large(p1A,freq_count_MAP_1)
rand_MAP_pop2 <- large(p2A,freq_count_MAP_2)
rand_MAP_pop3 <- large(p3A,freq_count_MAP_3)
rand_MAP_pop4 <- large(p4A,freq_count_MAP_4)
rand_MAP_pop5 <- large(p5A,freq_count_MAP_5)
rand_MAP_pop6 <- large(p6A,freq_count_MAP_6)
rand_MAP_pop7 <- large(p7A,freq_count_MAP_7)
rand_MAP_pop8 <- large(p8A,freq_count_MAP_8)
rand_MAP_pop9 <- large(p9A,freq_count_MAP_9)
rand_MAP_pop10 <- large(p10A,freq_count_MAP_10)
rand_MAP_pop11 <- large(p11A,freq_count_MAP_11)
rand_MAP_pop12 <- large(p12A,freq_count_MAP_12)

rand_CMD_pop1 <- large(p1A,freq_count_CMD_1)
rand_CMD_pop2 <- large(p2A,freq_count_CMD_2)
rand_CMD_pop3 <- large(p3A,freq_count_CMD_3)
rand_CMD_pop4 <- large(p4A,freq_count_CMD_4)
rand_CMD_pop5 <- large(p5A,freq_count_CMD_5)
rand_CMD_pop6 <- large(p6A,freq_count_CMD_6)
rand_CMD_pop7 <- large(p7A,freq_count_CMD_7)
rand_CMD_pop8 <- large(p8A,freq_count_CMD_8)
rand_CMD_pop9 <- large(p9A,freq_count_CMD_9)
rand_CMD_pop10 <- large(p10A,freq_count_CMD_10)
rand_CMD_pop11 <- large(p11A,freq_count_CMD_11)
rand_CMD_pop12 <- large(p12A,freq_count_CMD_12)

#Add pop ID
rand_mat_pop1 <- rand_mat_pop1 %>% mutate(Site=1)
rand_mat_pop2 <- rand_mat_pop2 %>% mutate(Site=2)
rand_mat_pop3 <- rand_mat_pop3 %>% mutate(Site=3)
rand_mat_pop4 <- rand_mat_pop4 %>% mutate(Site=4)
rand_mat_pop5 <- rand_mat_pop5 %>% mutate(Site=5)
rand_mat_pop6 <- rand_mat_pop6 %>% mutate(Site=6)
rand_mat_pop7 <- rand_mat_pop7 %>% mutate(Site=7)
rand_mat_pop8 <- rand_mat_pop8 %>% mutate(Site=8)
rand_mat_pop9 <- rand_mat_pop9 %>% mutate(Site=9)
rand_mat_pop10 <- rand_mat_pop10 %>% mutate(Site=10)
rand_mat_pop11 <- rand_mat_pop11 %>% mutate(Site=11)
rand_mat_pop12 <- rand_mat_pop12 %>% mutate(Site=12)

rand_MAP_pop1 <- rand_MAP_pop1 %>% mutate(Site=1)
rand_MAP_pop2 <- rand_MAP_pop2 %>% mutate(Site=2)
rand_MAP_pop3 <- rand_MAP_pop3 %>% mutate(Site=3)
rand_MAP_pop4 <- rand_MAP_pop4 %>% mutate(Site=4)
rand_MAP_pop5 <- rand_MAP_pop5 %>% mutate(Site=5)
rand_MAP_pop6 <- rand_MAP_pop6 %>% mutate(Site=6)
rand_MAP_pop7 <- rand_MAP_pop7 %>% mutate(Site=7)
rand_MAP_pop8 <- rand_MAP_pop8 %>% mutate(Site=8)
rand_MAP_pop9 <- rand_MAP_pop9 %>% mutate(Site=9)
rand_MAP_pop10 <- rand_MAP_pop10 %>% mutate(Site=10)
rand_MAP_pop11 <- rand_MAP_pop11 %>% mutate(Site=11)
rand_MAP_pop12 <- rand_MAP_pop12 %>% mutate(Site=12)

rand_CMD_pop1 <- rand_CMD_pop1 %>% mutate(Site=1)
rand_CMD_pop2 <- rand_CMD_pop2 %>% mutate(Site=2)
rand_CMD_pop3 <- rand_CMD_pop3 %>% mutate(Site=3)
rand_CMD_pop4 <- rand_CMD_pop4 %>% mutate(Site=4)
rand_CMD_pop5 <- rand_CMD_pop5 %>% mutate(Site=5)
rand_CMD_pop6 <- rand_CMD_pop6 %>% mutate(Site=6)
rand_CMD_pop7 <- rand_CMD_pop7 %>% mutate(Site=7)
rand_CMD_pop8 <- rand_CMD_pop8 %>% mutate(Site=8)
rand_CMD_pop9 <- rand_CMD_pop9 %>% mutate(Site=9)
rand_CMD_pop10 <- rand_CMD_pop10 %>% mutate(Site=10)
rand_CMD_pop11 <- rand_CMD_pop11 %>% mutate(Site=11)
rand_CMD_pop12 <- rand_CMD_pop12 %>% mutate(Site=12)

#Bind populations for each env
rand_mat_key <- rbind(rand_mat_pop1,
                  rand_mat_pop2,
                  rand_mat_pop3,
                  rand_mat_pop4,
                  rand_mat_pop5,
                  rand_mat_pop6,
                  rand_mat_pop7,
                  rand_mat_pop8,
                  rand_mat_pop9,
                  rand_mat_pop10,
                  rand_mat_pop11,
                  rand_mat_pop12)

rand_map_key <- rbind(rand_MAP_pop1,
                  rand_MAP_pop2,
                  rand_MAP_pop3,
                  rand_MAP_pop4,
                  rand_MAP_pop5,
                  rand_MAP_pop6,
                  rand_MAP_pop7,
                  rand_MAP_pop8,
                  rand_MAP_pop9,
                  rand_MAP_pop10,
                  rand_MAP_pop11,
                  rand_MAP_pop12)

rand_cmd_key <- rbind(rand_CMD_pop1,
                  rand_CMD_pop2,
                  rand_CMD_pop3,
                  rand_CMD_pop4,
                  rand_CMD_pop5,
                  rand_CMD_pop6,
                  rand_CMD_pop7,
                  rand_CMD_pop8,
                  rand_CMD_pop9,
                  rand_CMD_pop10,
                  rand_CMD_pop11,
                  rand_CMD_pop12)

#Sort 
rand_mat <- loci_snp %>% filter (chr_snp %in% as.character(rand_mat_key$chr_snp))
rand_map <- loci_snp %>% filter (chr_snp %in% as.character(rand_map_key$chr_snp))
rand_cmd <- loci_snp %>% filter (chr_snp %in% as.character(rand_cmd_key$chr_snp))



#Export each joint df
write_csv(rand_mat, "Genomics_scripts/Data/rand_mat.csv")
write_csv(rand_map, "Genomics_scripts/Data/rand_map.csv")
write_csv(rand_cmd, "Genomics_scripts/Data/rand_cmd.csv")

