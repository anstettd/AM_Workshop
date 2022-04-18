##################################################################################
## Setup timeseries frequencies
## Author Daniel Anstett
## 
## 
## Last Modified April 15, 2022
###################################################################################
#Import libraries
library(tidyverse)
###################################################################################
#Make Function
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
      
      freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV>=bin_step[j+1],na.rm=T )
    }
  }
  return(freq_count_calc)
}


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