#############################################################################################################
## Get SNPs names for 1, 5, and 10th percentile of WZA windows
## Author Daniel Anstett
## 
## 
## 
## Last Modified Sept 2, 2022
#############################################################################################################
#Import libraries

library(tidyverse)
library(Kendall)
library(plotly)

#############################################################################################################
#Functions

#Tabulate SNPs in each chatergory

snp_count <- function(df1,df5,df10){
  
  #Filter by BF >0
  peak_bf0_1 <- df1 %>% filter(BF>=0)
  peak_bf0_5 <- df5 %>% filter(BF>=0)
  peak_bf0_10 <- df10 %>% filter(BF>=0)
  
  #Filter by BF >0.5
  peak_bf05_1 <- df1 %>% filter(BF>=0.5)
  peak_bf05_5 <- df5 %>% filter(BF>=0.5)
  peak_bf05_10 <- df10 %>% filter(BF>=0.5)
  
  #Filter by BF >1.3
  peak_bf13_1 <- df1 %>% filter(BF>=1.3)
  peak_bf13_5 <- df5 %>% filter(BF>=1.3)
  peak_bf13_10 <- df10 %>% filter(BF>=1.3)
  
  #Filter by BF >2
  peak_bf2_1 <- df1 %>% filter(BF>=2)
  peak_bf2_5 <- df5 %>% filter(BF>=2)
  peak_bf2_10 <- df10 %>% filter(BF>=2)
  
  #Filter by BF >5
  peak_bf5_1 <- df1 %>% filter(BF>=5)
  peak_bf5_5 <- df5 %>% filter(BF>=5)
  peak_bf5_10 <- df10 %>% filter(BF>=5)
  
  #Filter by BF >10
  peak_bf10_1 <- df1 %>% filter(BF>=10)
  peak_bf10_5 <- df5 %>% filter(BF>=10)
  peak_bf10_10 <- df10 %>% filter(BF>=10)
  
  #Filter by BF >20
  peak_bf20_1 <- df1 %>% filter(BF>=20)
  peak_bf20_5 <- df5 %>% filter(BF>=20)
  peak_bf20_10 <- df10 %>% filter(BF>=20)
  
  #Filter by BF >20
  peak_bf30_1 <- df1 %>% filter(BF>=30)
  peak_bf30_5 <- df5 %>% filter(BF>=30)
  peak_bf30_10 <- df10 %>% filter(BF>=30)
  
  
  snps_available <- data.frame()
  
  snps_available[1,1] <-  dim(df1)[1]
  snps_available[1,2] <-  dim(df5)[1]
  snps_available[1,3] <-  dim(df10)[1]
  
  snps_available[2,1] <-  dim(peak_bf0_1)[1]
  snps_available[2,2] <-  dim(peak_bf0_5)[1]
  snps_available[2,3] <-  dim(peak_bf0_10)[1]
  
  snps_available[3,1] <-  dim(peak_bf05_1)[1]
  snps_available[3,2] <-  dim(peak_bf05_5)[1]
  snps_available[3,3] <-  dim(peak_bf05_10)[1]
  
  snps_available[4,1] <-  dim(peak_bf13_1)[1]
  snps_available[4,2] <-  dim(peak_bf13_5)[1]
  snps_available[4,3] <-  dim(peak_bf13_10)[1]
  
  snps_available[5,1] <-  dim(peak_bf2_1)[1]
  snps_available[5,2] <-  dim(peak_bf2_5)[1]
  snps_available[5,3] <-  dim(peak_bf2_10)[1]
  
  snps_available[6,1] <-  dim(peak_bf5_1)[1]
  snps_available[6,2] <-  dim(peak_bf5_5)[1]
  snps_available[6,3] <-  dim(peak_bf5_10)[1]
  
  snps_available[7,1] <-  dim(peak_bf10_1)[1]
  snps_available[7,2] <-  dim(peak_bf10_5)[1]
  snps_available[7,3] <-  dim(peak_bf10_10)[1]
  
  snps_available[8,1] <-  dim(peak_bf20_1)[1]
  snps_available[8,2] <-  dim(peak_bf20_5)[1]
  snps_available[8,3] <-  dim(peak_bf20_10)[1]
  
  snps_available[9,1] <-  dim(peak_bf30_1)[1]
  snps_available[9,2] <-  dim(peak_bf30_5)[1]
  snps_available[9,3] <-  dim(peak_bf30_10)[1]
  
  
  return(snps_available)
}



#Calculate windows lost after filter by BF

win_lost <- function(df1,df5,df10){
  
  #Identify max BF for each window
  max_bf_mat <- df1 %>% group_by(V4) %>% summarise(max = max(BF, na.rm=TRUE))
  max_bf_map <- df5 %>% group_by(V4) %>% summarise(max = max(BF, na.rm=TRUE))
  max_bf_cmd <- df10 %>% group_by(V4) %>% summarise(max = max(BF, na.rm=TRUE))
  

  windows_lost <- data.frame()

  windows_lost[1,1] <- dim(max_bf_mat)[1] 
  windows_lost[1,2] <- dim(max_bf_map)[1] 
  windows_lost[1,3] <- dim(max_bf_cmd)[1] 
  
  windows_lost[2,1] <- dim(max_bf_mat)[1] - dim(max_bf_mat %>% filter(max>0))[1]
  windows_lost[2,2] <- dim(max_bf_map)[1] - dim(max_bf_map %>% filter(max>0))[1]
  windows_lost[2,3] <- dim(max_bf_cmd)[1] - dim(max_bf_cmd %>% filter(max>0))[1]
  
  windows_lost[3,1] <- dim(max_bf_mat)[1] - dim(max_bf_mat %>% filter(max>0.5))[1]
  windows_lost[3,2] <- dim(max_bf_map)[1] - dim(max_bf_map %>% filter(max>0.5))[1]
  windows_lost[3,3] <- dim(max_bf_cmd)[1] - dim(max_bf_cmd %>% filter(max>0.5))[1]
  
  windows_lost[4,1] <- dim(max_bf_mat)[1] - dim(max_bf_mat %>% filter(max>1.3))[1]
  windows_lost[4,2] <- dim(max_bf_map)[1] - dim(max_bf_map %>% filter(max>1.3))[1]
  windows_lost[4,3] <- dim(max_bf_cmd)[1] - dim(max_bf_cmd %>% filter(max>1.3))[1]
  
  windows_lost[5,1] <- dim(max_bf_mat)[1] - dim(max_bf_mat %>% filter(max>2))[1]
  windows_lost[5,2] <- dim(max_bf_map)[1] - dim(max_bf_map %>% filter(max>2))[1]
  windows_lost[5,3] <- dim(max_bf_cmd)[1] - dim(max_bf_cmd %>% filter(max>2))[1]
  
  windows_lost[6,1] <- dim(max_bf_mat)[1] - dim(max_bf_mat %>% filter(max>5))[1]
  windows_lost[6,2] <- dim(max_bf_map)[1] - dim(max_bf_map %>% filter(max>5))[1]
  windows_lost[6,3] <- dim(max_bf_cmd)[1] - dim(max_bf_cmd %>% filter(max>5))[1]
  
  windows_lost[7,1] <- dim(max_bf_mat)[1] - dim(max_bf_mat %>% filter(max>10))[1]
  windows_lost[7,2] <- dim(max_bf_map)[1] - dim(max_bf_map %>% filter(max>10))[1]
  windows_lost[7,3] <- dim(max_bf_cmd)[1] - dim(max_bf_cmd %>% filter(max>10))[1]
  
  windows_lost[8,1] <- dim(max_bf_mat)[1] - dim(max_bf_mat %>% filter(max>20))[1]
  windows_lost[8,2] <- dim(max_bf_map)[1] - dim(max_bf_map %>% filter(max>20))[1]
  windows_lost[8,3] <- dim(max_bf_cmd)[1] - dim(max_bf_cmd %>% filter(max>20))[1]
  
  windows_lost[9,1] <- dim(max_bf_mat)[1] - dim(max_bf_mat %>% filter(max>30))[1]
  windows_lost[9,2] <- dim(max_bf_map)[1] - dim(max_bf_map %>% filter(max>30))[1]
  windows_lost[9,3] <- dim(max_bf_cmd)[1] - dim(max_bf_cmd %>% filter(max>30))[1]
  
  colnames(windows_lost) <- c("MAT","MAP","CMD")
  rownames(windows_lost) <- c("all_BF","BF>0","BF>0.5","BF>1.3","BF>2","BF>5","BF>10","BF>20","BF>30")
  
  
  return(windows_lost)
}


#############################################################################################################


#Import files
snps_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_mat.csv")
snps_map <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_map.csv")
snps_cmd <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_cmd.csv")

wza_win_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_win_mat.csv")
wza_win_map <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_win_map.csv")
wza_win_cmd <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_win_cmd.csv")

#Import Large Loci Win
loci_win <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/loci_win.csv")                                                                                 
loci_win <- loci_win %>% unite(col="chr_snp", c("chr","snp"), sep="_")

#Get Quantile Lists
wza_win_mat_1 <- wza_win_mat %>% filter(WZA > quantile(wza_win_mat$WZA, probs = 0.99))
wza_win_mat_5 <- wza_win_mat %>% filter(WZA > quantile(wza_win_mat$WZA, probs = 0.95))
wza_win_mat_10 <- wza_win_mat %>% filter(WZA > quantile(wza_win_mat$WZA, probs = 0.90))

wza_win_map_1 <- wza_win_map %>% filter(WZA > quantile(wza_win_map$WZA, probs = 0.99))
wza_win_map_5 <- wza_win_map %>% filter(WZA > quantile(wza_win_map$WZA, probs = 0.95))
wza_win_map_10 <- wza_win_map %>% filter(WZA > quantile(wza_win_map$WZA, probs = 0.90))

wza_win_cmd_1 <- wza_win_cmd %>% filter(WZA > quantile(wza_win_cmd$WZA, probs = 0.99))
wza_win_cmd_5 <- wza_win_cmd %>% filter(WZA > quantile(wza_win_cmd$WZA, probs = 0.95))
wza_win_cmd_10 <- wza_win_cmd %>% filter(WZA > quantile(wza_win_cmd$WZA, probs = 0.90))


#Filter SNPs for precentile
snps_mat_peak_1 <- snps_mat %>% filter(win %in% wza_win_mat_1$win)
snps_mat_peak_5 <- snps_mat %>% filter(win %in% wza_win_mat_5$win)
snps_mat_peak_10 <- snps_mat %>% filter(win %in% wza_win_mat_10$win)

snps_map_peak_1 <- snps_map %>% filter(win %in% wza_win_map_1$win)
snps_map_peak_5 <- snps_map %>% filter(win %in% wza_win_map_5$win)
snps_map_peak_10 <- snps_map %>% filter(win %in% wza_win_map_10$win)

snps_cmd_peak_1 <- snps_cmd %>% filter(win %in% wza_win_cmd_1$win)
snps_cmd_peak_5 <- snps_cmd %>% filter(win %in% wza_win_cmd_5$win)
snps_cmd_peak_10 <- snps_cmd %>% filter(win %in% wza_win_cmd_10$win)


#Import snp env associations (Baseline)
env1 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_1_trim.tsv",header=F, sep=" ")
env2 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_2_trim.tsv",header=F, sep=" ")
env5 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_5_trim.tsv",header=F, sep=" ")

#Name Columns
colnames(env1) <- c("Chromosome","SNP","Env","BF")
colnames(env2) <- c("Chromosome","SNP","Env","BF")
colnames(env5) <- c("Chromosome","SNP","Env","BF")

env1_united <- env1 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env2_united <- env2 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env5_united <- env5 %>% unite(chr_snp,"Chromosome","SNP",sep="_")

#############################################################################################################
#Filter baseline by BF=>30

##Filter Bayes factor by peak windows
env1_united_bf30 <- env1_united %>% filter(BF>30)
env2_united_bf30 <- env2_united %>% filter(BF>30)
env5_united_bf30 <- env5_united %>% filter(BF>30)


#Filter Bayes factor by peak windows
snp_mat_1 <- env1_united %>% filter(BF<=30) %>% filter(chr_snp %in% as.character(snps_mat_peak_1$chr_snp))
snp_mat_5 <- env1_united %>% filter(BF<=30) %>% filter(chr_snp %in% as.character(snps_mat_peak_5$chr_snp))
snp_mat_10 <- env1_united%>% filter(BF<=30)  %>% filter(chr_snp %in% as.character(snps_mat_peak_10$chr_snp))

snp_map_1 <- env2_united %>% filter(BF<=30) %>% filter(chr_snp %in% as.character(snps_map_peak_1$chr_snp))
snp_map_5 <- env2_united %>% filter(BF<=30) %>% filter(chr_snp %in% as.character(snps_map_peak_5$chr_snp))
snp_map_10 <- env2_united %>% filter(BF<=30) %>% filter(chr_snp %in% as.character(snps_map_peak_10$chr_snp))

snp_cmd_1 <- env5_united %>% filter(BF<=30) %>% filter(chr_snp %in% as.character(snps_cmd_peak_1$chr_snp))
snp_cmd_5 <- env5_united %>% filter(BF<=30) %>% filter(chr_snp %in% as.character(snps_cmd_peak_5$chr_snp))
snp_cmd_10 <- env5_united %>% filter(BF<=30) %>% filter(chr_snp %in% as.character(snps_cmd_peak_10$chr_snp))

#Merge
snp_mat_1 <- rbind(env1_united_bf30,snp_mat_1)
snp_mat_5 <- rbind(env1_united_bf30,snp_mat_5)
snp_mat_10 <- rbind(env1_united_bf30,snp_mat_10)

snp_map_1 <- rbind(env1_united_bf30,snp_map_1)
snp_map_5 <- rbind(env2_united_bf30,snp_map_5)
snp_map_10 <- rbind(env2_united_bf30,snp_map_10)

snp_cmd_1 <- rbind(env5_united_bf30,snp_cmd_1)
snp_cmd_5 <- rbind(env5_united_bf30,snp_cmd_5)
snp_cmd_10 <- rbind(env5_united_bf30,snp_cmd_10)

#Window calculator
snp_mat_1_win <- left_join(snp_mat_1,loci_win, by="chr_snp")
snp_mat_5_win <- left_join(snp_mat_5,loci_win, by="chr_snp")
snp_mat_10_win <- left_join(snp_mat_10,loci_win, by="chr_snp")

snp_map_1_win <- left_join(snp_map_1,loci_win, by="chr_snp")
snp_map_5_win <- left_join(snp_map_5,loci_win, by="chr_snp")
snp_map_10_win <- left_join(snp_map_10,loci_win, by="chr_snp")

snp_cmd_1_win <- left_join(snp_cmd_1,loci_win, by="chr_snp")
snp_cmd_5_win <- left_join(snp_cmd_5,loci_win, by="chr_snp")
snp_cmd_10_win <- left_join(snp_cmd_10,loci_win, by="chr_snp")





#################################################################################################

snps_avail_mat <- snp_count(snp_mat_1,snp_mat_5,snp_mat_10)
snps_avail_map <- snp_count(snp_map_1,snp_map_5,snp_map_10)
snps_avail_cmd <- snp_count(snp_cmd_1,snp_cmd_5,snp_cmd_10)
snps_avail_precentiles <- cbind(snps_avail_mat,snps_avail_map,snps_avail_cmd)

write_csv(snps_avail_precentiles,"Genomics_scripts/Data/snps_avail_precentiles.csv")

win_lost_mat <- win_lost(snp_mat_1_win,snp_mat_5_win,snp_mat_10_win)
win_lost_map <- win_lost(snp_map_1_win,snp_map_5_win,snp_map_10_win)
win_lost_cmd <- win_lost(snp_cmd_1_win,snp_cmd_5_win,snp_cmd_10_win)
win_lost_percentiles <- cbind(win_lost_mat,win_lost_map,win_lost_cmd)

write_csv(win_lost_percentiles,"Genomics_scripts/Data/win_lost_percentile.csv")







