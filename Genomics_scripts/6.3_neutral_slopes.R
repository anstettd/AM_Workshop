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
snp1_time <- read_csv("Genomics_scripts/Data/snp1_filter.csv")
snp2_time <- read_csv("Genomics_scripts/Data/snp2_filter.csv")
snp5_time <- read_csv("Genomics_scripts/Data/snp5_filter.csv")

#Import full snp table for timeseries
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
colnames(loci) <- c("Chromosome","SNP")
loci_united <- loci %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp <-cbind(loci_united,snp) #add snp lables to rows

#Filter full snp table to remove climate associated SNPs
snp_swiss <-loci_snp %>% filter (!chr_snp %in% as.character(snp1_time$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(snp1_time$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(snp1_time$chr_snp))

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

write_csv(freq_count_MAT, "Genomics_scripts/Data/freq_count_MAT.csv")
write_csv(freq_count_MAP, "Genomics_scripts/Data/freq_count_MAP.csv")
write_csv(freq_count_CMD, "Genomics_scripts/Data/freq_count_CMD.csv")

#######################################################################################################
#Make SNP table for neutral loci
#Caution significant computational time. Will use ~3 GB of ram

snpA <- data.frame()
counter<-1
for (i in seq (2,dim(snp_swiss)[2]-1,2)){
  for(j in 1:dim(snp_swiss)[1]){
    tmp_total<-as.numeric(snp_swiss[j,i]) + as.numeric(snp_swiss[j,i+1])
    snpA[j,counter]<-as.numeric(snp_swiss[j,i])/tmp_total
  }
  counter<-counter+1
}
colnames(snpA)<- pop_order[,1] #name each pop/time combination
rownames(snpA)<- snpA$chr_snp

snpA_base <- snpA %>% select(chr_snp,P1_2010,P2_2010,P3_2010,P4_2010,P5_2010,P6_2010,
                             P7_2010,P8_2011,P9_2010,P10_2011,P11_2010,P12_2010)






