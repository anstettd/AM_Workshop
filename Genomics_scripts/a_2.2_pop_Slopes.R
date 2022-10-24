##################################################################################
## Get BF and Window info for all Slopes
## 
## Author Daniel Anstett
## 
## Last Modified October 12, 2022
###################################################################################
library(tidyverse)
library(fs)

#Impot Loci Windows
loci_win <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/loci_win.csv")
colnames(loci_win) <- c("Chromosome","SNP","snp_c","Window")
loci_win <- loci_win %>% unite(chr_snp,"Chromosome","SNP",sep="_")

###################################################################################
##MAT

#Import BF
env1 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_1_trim.tsv",header=F, sep=" ")
colnames(env1) <- c("Chromosome","SNP","Env","BF")
env1 <- env1 %>% unite(chr_snp,"Chromosome","SNP",sep="_")

#Import WZA
wza_win_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_win_mat.csv")
wza_win_mat <- wza_win_mat %>% select(win,WZA,approx_p)

#Batch Import Slopes
path_files <- "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT"
csv_files <- fs::dir_ls(path=path_files, pattern = "\\all_slopes_MAT*.csv$")
mat_all <- csv_files %>% map_dfr(read_csv)
colnames(mat_all)[2] <- "chr_snp"
mat_all <- na.omit(mat_all)

#Left Join BF, window ID and window p-value to slope df
mat_BF <- left_join(mat_all,env1,by="chr_snp")
mat_BF_win <- left_join(mat_BF,loci_win,by="chr_snp")
slope_mat <- left_join(mat_BF_win,wza_win_mat,by=c("Window"="win"))
write_csv(slope_mat, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/slope_mat.csv")


###################################################################################
##MAP

#Import BF
env2 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_1_trim.tsv",header=F, sep=" ")
colnames(env2) <- c("Chromosome","SNP","Env","BF")
env2 <- env2 %>% unite(chr_snp,"Chromosome","SNP",sep="_")

#Import WZA
wza_win_map <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_win_map.csv")
wza_win_map <- wza_win_map %>% select(win,WZA,approx_p)

#Batch Import Slopes
path_files <- "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP"
csv_files <- fs::dir_ls(path=path_files, pattern = "\\all_slopes_MAP*.csv$")
map_all <- csv_files %>% map_dfr(read_csv)
colnames(map_all)[2] <- "chr_snp"
map_all <- na.omit(map_all)

#Left Join BF, window ID and window p-value to slope df
map_BF <- left_join(map_all,env2,by="chr_snp")
map_BF_win <- left_join(map_BF,loci_win,by="chr_snp")
slope_map <- left_join(map_BF_win,wza_win_map,by=c("Window"="win"))
write_csv(slope_map, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/slope_map.csv")


###################################################################################
##CMD

#Import BF
env5 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_1_trim.tsv",header=F, sep=" ")
colnames(env5) <- c("Chromosome","SNP","Env","BF")
env5 <- env5 %>% unite(chr_snp,"Chromosome","SNP",sep="_")

#Import WZA
wza_win_cmd <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_win_cmd.csv")
wza_win_cmd <- wza_win_cmd %>% select(win,WZA,approx_p)

#Batch Import Slopes
path_files <- "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD"
csv_files <- fs::dir_ls(path=path_files, pattern = "\\all_slopes_CMD*.csv$")
cmd_all <- csv_files %>% map_dfr(read_csv)
colnames(cmd_all)[2] <- "chr_snp"
cmd_all <- na.omit(cmd_all)

#Left Join BF, window ID and window p-value to slope df
cmd_BF <- left_join(cmd_all,env5,by="chr_snp")
cmd_BF_win <- left_join(cmd_BF,loci_win,by="chr_snp")
slope_cmd <- left_join(cmd_BF_win,wza_win_cmd,by=c("Window"="win"))
write_csv(slope_cmd, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/slope_cmd.csv")




