##################################################################################
## Get positive slope info per population
## Author Daniel Anstett
## 
## 
## Last Modified September 19, 2022
###################################################################################

#Library install and import
library(tidyverse) # for data manipulation

#Bring in dataframe with offset information
offset_pop <- read_csv("Genomics_scripts/Data/offset_pop.csv")

#Bring in peak bf>2 SNPs per env variabile
snp1_time <- read_csv("Genomics_scripts/Data/peak_bf2_time_mat.csv")
snp2_time <- read_csv("Genomics_scripts/Data/peak_bf2_time_map.csv")
snp5_time <- read_csv("Genomics_scripts/Data/peak_bf2_time_cmd.csv")
snp_time <- rbind(snp1_time[,1],snp2_time[,1],snp5_time[,1])
snp_unique <- unique(snp_time)


#Import observed slopes per population
mat_slope <- read_csv("Genomics_scripts/Data/freq_MAT_slope_peakbf2.csv")
map_slope <- read_csv("Genomics_scripts/Data/freq_MAP_slope_peakbf2.csv")
cmd_slope <- read_csv("Genomics_scripts/Data/freq_cmd_slope_peakbf2.csv")
env_slope <- rbind(mat_slope,map_slope,cmd_slope)

#Filter env slope to only have unique values per population
env_slope_1_duplicate <- env_slope %>% filter(Site==1)
env_slope_1_unique <- unique(env_slope_1_duplicate$snp_ID)
env_slope_1<-env_slope_1_duplicate[match(env_slope_1_unique,env_slope_1_duplicate$snp_ID),]

env_slope_2_duplicate <- env_slope %>% filter(Site==2)
env_slope_2_unique <- unique(env_slope_2_duplicate$snp_ID)
env_slope_2<-env_slope_2_duplicate[match(env_slope_2_unique,env_slope_2_duplicate$snp_ID),]

env_slope_3_duplicate <- env_slope %>% filter(Site==3)
env_slope_3_unique <- unique(env_slope_3_duplicate$snp_ID)
env_slope_3<-env_slope_3_duplicate[match(env_slope_3_unique,env_slope_3_duplicate$snp_ID),]

env_slope_4_duplicate <- env_slope %>% filter(Site==4)
env_slope_4_unique <- unique(env_slope_4_duplicate$snp_ID)
env_slope_4<-env_slope_4_duplicate[match(env_slope_4_unique,env_slope_4_duplicate$snp_ID),]

env_slope_5_duplicate <- env_slope %>% filter(Site==5)
env_slope_5_unique <- unique(env_slope_5_duplicate$snp_ID)
env_slope_5<-env_slope_5_duplicate[match(env_slope_5_unique,env_slope_5_duplicate$snp_ID),]

env_slope_6_duplicate <- env_slope %>% filter(Site==6)
env_slope_6_unique <- unique(env_slope_6_duplicate$snp_ID)
env_slope_6<-env_slope_6_duplicate[match(env_slope_6_unique,env_slope_6_duplicate$snp_ID),]

env_slope_7_duplicate <- env_slope %>% filter(Site==7)
env_slope_7_unique <- unique(env_slope_7_duplicate$snp_ID)
env_slope_7<-env_slope_7_duplicate[match(env_slope_7_unique,env_slope_7_duplicate$snp_ID),]

env_slope_8_duplicate <- env_slope %>% filter(Site==8)
env_slope_8_unique <- unique(env_slope_8_duplicate$snp_ID)
env_slope_8<-env_slope_8_duplicate[match(env_slope_8_unique,env_slope_8_duplicate$snp_ID),]

env_slope_9_duplicate <- env_slope %>% filter(Site==9)
env_slope_9_unique <- unique(env_slope_9_duplicate$snp_ID)
env_slope_9<-env_slope_9_duplicate[match(env_slope_9_unique,env_slope_9_duplicate$snp_ID),]

env_slope_10_duplicate <- env_slope %>% filter(Site==10)
env_slope_10_unique <- unique(env_slope_10_duplicate$snp_ID)
env_slope_10<-env_slope_10_duplicate[match(env_slope_10_unique,env_slope_10_duplicate$snp_ID),]

env_slope_11_duplicate <- env_slope %>% filter(Site==11)
env_slope_11_unique <- unique(env_slope_11_duplicate$snp_ID)
env_slope_11<-env_slope_11_duplicate[match(env_slope_11_unique,env_slope_11_duplicate$snp_ID),]

env_slope_12_duplicate <- env_slope %>% filter(Site==12)
env_slope_12_unique <- unique(env_slope_12_duplicate$snp_ID)
env_slope_12<-env_slope_12_duplicate[match(env_slope_12_unique,env_slope_12_duplicate$snp_ID),]

#Filter for slopes > 0.2
env_slope_1 <- env_slope_1 %>% filter(Slope>=0.2)
env_slope_2 <- env_slope_2 %>% filter(Slope>=0.2)
env_slope_3 <- env_slope_3 %>% filter(Slope>=0.2)
env_slope_4 <- env_slope_4 %>% filter(Slope>=0.2)
env_slope_5 <- env_slope_5 %>% filter(Slope>=0.2)
env_slope_6 <- env_slope_6 %>% filter(Slope>=0.2)
env_slope_7 <- env_slope_7 %>% filter(Slope>=0.2)
env_slope_8 <- env_slope_8 %>% filter(Slope>=0.2)
env_slope_9 <- env_slope_9 %>% filter(Slope>=0.2)
env_slope_10 <- env_slope_10 %>% filter(Slope>=0.2)
env_slope_11 <- env_slope_11 %>% filter(Slope>=0.2)
env_slope_12 <- env_slope_12 %>% filter(Slope>=0.2)



#Input number of slopes >0.2 into offset data frame
offset_slope <- data.frame()
offset_slope[1,1] <- dim(env_slope_1)[1]
offset_slope[2,1] <- dim(env_slope_2)[1]
offset_slope[3,1] <- dim(env_slope_3)[1]
offset_slope[4,1] <- dim(env_slope_4)[1]
offset_slope[5,1] <- dim(env_slope_5)[1]
offset_slope[6,1] <- dim(env_slope_6)[1]
offset_slope[7,1] <- dim(env_slope_7)[1]
offset_slope[8,1] <- dim(env_slope_8)[1]
offset_slope[9,1] <- dim(env_slope_9)[1]
offset_slope[10,1] <- dim(env_slope_10)[1]
offset_slope[11,1] <- dim(env_slope_11)[1]
offset_slope[12,1] <- dim(env_slope_12)[1]
colnames(offset_slope) <- "Moderate_Selection"

offset_pop_s <- cbind(offset_pop,offset_slope)

write_csv(offset_pop_s,"Genomics_scripts/Data/offset_pop_s.csv")

