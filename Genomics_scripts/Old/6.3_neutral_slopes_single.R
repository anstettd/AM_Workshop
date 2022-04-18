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

#######################################################################################################
## Setup timeseries frequencies
#######################################################################################################

#Import transformed timeseries frequencies
freq_MAT <- read_csv("Genomics_scripts/Data/freq_MAT.csv")
freq_MAP <- read_csv("Genomics_scripts/Data/freq_MAP.csv")
freq_CMD <- read_csv("Genomics_scripts/Data/freq_CMD.csv")

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
freq_count_MAT <- data.frame()
freq_count_MAP <- data.frame()
freq_count_CMD <- data.frame()

for (i in 1:12){
  test_MAT <- freq_MAT_1011 %>% filter(Site==i) %>% select (-Site,-Year)
  test_MAT <- as.data.frame(test_MAT)
  test_MAT <- as.numeric(test_MAT[1,])
  freq_count_MAT[1,i] <- sum(test_MAT >= 0 & test_MAT <0.1, na.rm = T)
  freq_count_MAT[2,i] <- sum(test_MAT >= 0.1 & test_MAT <0.2, na.rm = T)
  freq_count_MAT[3,i] <- sum(test_MAT >= 0.2 & test_MAT <0.3, na.rm = T)
  freq_count_MAT[4,i] <- sum(test_MAT >= 0.3 & test_MAT <0.4, na.rm = T)
  freq_count_MAT[5,i] <- sum(test_MAT >= 0.4 & test_MAT <0.5, na.rm = T)
  
  freq_count_MAT[6,i] <- sum(test_MAT >= 0.5 & test_MAT <0.6, na.rm = T)
  freq_count_MAT[7,i] <- sum(test_MAT >= 0.6 & test_MAT <0.7, na.rm = T)
  freq_count_MAT[8,i] <- sum(test_MAT >= 0.7 & test_MAT <0.8, na.rm = T)
  freq_count_MAT[9,i] <- sum(test_MAT >= 0.8 & test_MAT <0.9, na.rm = T)
  freq_count_MAT[10,i] <- sum(test_MAT >= 0.9 & test_MAT <=1, na.rm = T)
}

for (i in 1:12){
  test_MAP <- freq_MAP_1011 %>% filter(Site==i) %>% select (-Site,-Year)
  test_MAP <- as.data.frame(test_MAP)
  test_MAP <- as.numeric(test_MAP[1,])
  freq_count_MAP[1,i] <- sum(test_MAP >= 0 & test_MAP <0.1, na.rm = T)
  freq_count_MAP[2,i] <- sum(test_MAP >= 0.1 & test_MAP <0.2, na.rm = T)
  freq_count_MAP[3,i] <- sum(test_MAP >= 0.2 & test_MAP <0.3, na.rm = T)
  freq_count_MAP[4,i] <- sum(test_MAP >= 0.3 & test_MAP <0.4, na.rm = T)
  freq_count_MAP[5,i] <- sum(test_MAP >= 0.4 & test_MAP <0.5, na.rm = T)
  
  freq_count_MAP[6,i] <- sum(test_MAP >= 0.5 & test_MAP <0.6, na.rm = T)
  freq_count_MAP[7,i] <- sum(test_MAP >= 0.6 & test_MAP <0.7, na.rm = T)
  freq_count_MAP[8,i] <- sum(test_MAP >= 0.7 & test_MAP <0.8, na.rm = T)
  freq_count_MAP[9,i] <- sum(test_MAP >= 0.8 & test_MAP <0.9, na.rm = T)
  freq_count_MAP[10,i] <- sum(test_MAP >= 0.9 & test_MAP <=1, na.rm = T)
}

for (i in 1:12){
  test_CMD <- freq_CMD_1011 %>% filter(Site==i) %>% select (-Site,-Year)
  test_CMD <- as.data.frame(test_CMD)
  test_CMD <- as.numeric(test_CMD[1,])
  freq_count_CMD[1,i] <- sum(test_CMD >= 0 & test_CMD <0.1, na.rm = T)
  freq_count_CMD[2,i] <- sum(test_CMD >= 0.1 & test_CMD <0.2, na.rm = T)
  freq_count_CMD[3,i] <- sum(test_CMD >= 0.2 & test_CMD <0.3, na.rm = T)
  freq_count_CMD[4,i] <- sum(test_CMD >= 0.3 & test_CMD <0.4, na.rm = T)
  freq_count_CMD[5,i] <- sum(test_CMD >= 0.4 & test_CMD <0.5, na.rm = T)
  
  freq_count_CMD[6,i] <- sum(test_CMD >= 0.5 & test_CMD <0.6, na.rm = T)
  freq_count_CMD[7,i] <- sum(test_CMD >= 0.6 & test_CMD <0.7, na.rm = T)
  freq_count_CMD[8,i] <- sum(test_CMD >= 0.7 & test_CMD <0.8, na.rm = T)
  freq_count_CMD[9,i] <- sum(test_CMD >= 0.8 & test_CMD <0.9, na.rm = T)
  freq_count_CMD[10,i] <- sum(test_CMD >= 0.9 & test_CMD <=1, na.rm = T)
}

colnames(freq_count_MAT) <- c("ID_1","ID_2","ID_3","ID_4","ID_5","ID_6",
                              "ID_7","ID_8","ID_9","ID_10","ID_11","ID_12")
colnames(freq_count_MAP) <- c("ID_1","ID_2","ID_3","ID_4","ID_5","ID_6",
                              "ID_7","ID_8","ID_9","ID_10","ID_11","ID_12")
colnames(freq_count_CMD) <- c("ID_1","ID_2","ID_3","ID_4","ID_5","ID_6",
                              "ID_7","ID_8","ID_9","ID_10","ID_11","ID_12")

#write_csv(freq_count_MAT, "Genomics_scripts/Data/freq_count_MAT.csv")
#write_csv(freq_count_MAP, "Genomics_scripts/Data/freq_count_MAP.csv")
#write_csv(freq_count_CMD, "Genomics_scripts/Data/freq_count_CMD.csv")

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
#Filter per_frequency

#Pop 1
p1_1 <- p1A %>% filter(snpA >= 0 & snpA <0.1)
p1_2 <- p1A %>% filter(snpA >= 0.1 & snpA <0.2)
p1_3 <- p1A %>% filter(snpA >= 0.2 & snpA <0.3)
p1_4 <- p1A %>% filter(snpA >= 0.3 & snpA <0.4)
p1_5 <- p1A %>% filter(snpA >= 0.4 & snpA <0.5)

p1_6 <- p1A %>% filter(snpA >= 0.5 & snpA <0.6)
p1_7 <- p1A %>% filter(snpA >= 0.6 & snpA <0.7)
p1_8 <- p1A %>% filter(snpA >= 0.7 & snpA <0.8)
p1_9 <- p1A %>% filter(snpA >= 0.8 & snpA <0.9)
p1_10 <- p1A %>% filter(snpA >= 0.9 & snpA <=1)

#######################################################################################################
##Sample each frequency per pop 

set.seed(1)

#Pop1
#Generate random number list
list_mat_p1_1 <- sample.int(dim(p1_1)[1],freq_count_MAT$ID_1[1])
list_mat_p1_2 <- sample.int(dim(p1_2)[1],freq_count_MAT$ID_1[2])
list_mat_p1_3 <- sample.int(dim(p1_3)[1],freq_count_MAT$ID_1[3])
list_mat_p1_4 <- sample.int(dim(p1_4)[1],freq_count_MAT$ID_1[4])
list_mat_p1_5 <- sample.int(dim(p1_5)[1],freq_count_MAT$ID_1[5])
list_mat_p1_6 <- sample.int(dim(p1_6)[1],freq_count_MAT$ID_1[6])
list_mat_p1_7 <- sample.int(dim(p1_7)[1],freq_count_MAT$ID_1[7])
list_mat_p1_8 <- sample.int(dim(p1_8)[1],freq_count_MAT$ID_1[8])
list_mat_p1_9 <- sample.int(dim(p1_9)[1],freq_count_MAT$ID_1[9])
list_mat_p1_10 <- sample.int(dim(p1_10)[1],freq_count_MAT$ID_1[10])

#Make data frame
rand_mat_p1_1 <- data.frame()


#Select randomized columns
for (i in 1:length(list_mat_p1_1)){
    rand_mat_p1_1<- rbind(rand_mat_p1_1 ,p1_1 %>% 
                            filter(chr_snp==as.character(p1_1$chr_snp[list_mat_p1_1[i]])))
}


for (i in 1:length(list_mat_p1_2)){
  rand_mat_p1_2[i,] <- p1_2 %>% filter(chr_snp==as.character(p1_2$chr_snp[list_mat_p1_2[i]] ))
}
for (i in 1:length(list_mat_p1_3)){
  rand_mat_p1_3[i,] <- p1_3 %>% filter(chr_snp==as.character(p1_3$chr_snp[list_mat_p1_3[i]] ))
}
for (i in 1:length(list_mat_p1_4)){
  rand_mat_p1_4[i,] <- p1_4 %>% filter(chr_snp==as.character(p1_4$chr_snp[list_mat_p1_4[i]] ))
}
for (i in 1:length(list_mat_p1_5)){
  rand_mat_p1_5[i,] <- p1_5 %>% filter(chr_snp==as.character(p1_5$chr_snp[list_mat_p1_5[i]] ))
}

for (i in 1:length(list_mat_p1_6)){
  rand_mat_p1_6[i,] <- p1_6 %>% filter(chr_snp==as.character(p1_6$chr_snp[list_mat_p1_6[i]] ))
}
for (i in 1:length(list_mat_p1_7)){
  rand_mat_p1_7[i,] <- p1_7 %>% filter(chr_snp==as.character(p1_7$chr_snp[list_mat_p1_7[i]] ))
}
for (i in 1:length(list_mat_p1_8)){
  rand_mat_p1_8[i,] <- p1_8 %>% filter(chr_snp==as.character(p1_8$chr_snp[list_mat_p1_8[i]] ))
}
for (i in 1:length(list_mat_p1_9)){
  rand_mat_p1_9[i,] <- p1_9 %>% filter(chr_snp==as.character(p1_9$chr_snp[list_mat_p1_9[i]] ))
}
for (i in 1:length(list_mat_p1_0)){
  rand_mat_p1_10[i,] <- p1_10 %>% filter(chr_snp==as.character(p1_10$chr_snp[list_mat_p1_10[i]] ))
}

#Bind lists of snps together

mat_rand_snp_pop1 <- rbind(as.data.frame(rand_mat_p1_1),as.data.frame(rand_mat_p1_2))













for(i in 1:10){
  nam <- paste("list_mat_p1", i, sep = "_")
  rrange <- paste("p1",i,sep="_")
  nam <- sample.int(dim(rrange)[i],freq_count_MAT$ID_1[i]) 
}










#loci_base_p<-head(loci_base_p1)
snpA_p1 <- data.frame()
for(i in 1:dim(loci_base_p1)[1]){
    tmp_total<-as.numeric(loci_base_p1[i,2]) + as.numeric(loci_base_p1[i,3])
    snpA_p1[i,1]<-as.numeric(loci_base_p1[i,2])/tmp_total
  }


#Caution significant computational time. Will use ~3 GB of ram
snpA_p1 <- data.frame()
counter<-1
for (i in seq (2,dim(loci_base_p1)[2]-1,2)){
  for(j in 1:dim(loci_base_p1)[1]){
    tmp_total<-as.numeric(loci_base_p1[j,i]) + as.numeric(loci_base_p1[j,i+1])
    snpA_p1[j,counter]<-as.numeric(loci_base_p1[j,i])/tmp_total
  }
  counter<-counter+1
}

colnames(loci_base) <- c("chr_snp","Pop1A","Pop1B",
                         "Pop2A","Pop2B",
                         "Pop3A","Pop3B",
                         "Pop4A","Pop4B",
                         "Pop5A","Pop5B",
                         "Pop6A","Pop6B",
                         "Pop7A","Pop7B",
                         "Pop8A","Pop8B",
                         "Pop9A","Pop9B",
                         "Pop10A","Pop10B",
                         "Pop11A","Pop11B",
                         "Pop12A","Pop12B")


colnames(snpA)<- pop_order[,1] #name each pop/time combination
rownames(snpA)<- snpA$chr_snp








