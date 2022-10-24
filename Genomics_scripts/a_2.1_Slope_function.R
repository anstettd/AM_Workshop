##################################################################################
## Get slopes (strength of selection) for all BF
## 
## Author Jose Lazaro
## 
## Last Modified October 6, 2022
###################################################################################

library(reshape2)
library(dplyr)
library(tidyverse)
###################################################################################
#Slope Function

slope_melt <- function(df) {
  freq_1.melted = reshape2::melt(df, id.vars = c("Site", "Year"))
  colnames(freq_1.melted)[3] <- "snp_ID"
  colnames(freq_1.melted)[4] <- "freq"
  
  freq_1.melted_mod <- na.omit(freq_1.melted)
  
  freq_1.slope <- group_by(freq_1.melted_mod, Site, snp_ID) %>%
    arrange(Site, Year) %>%
    summarize(Slope = glm(freq ~ Year, family = binomial)$coefficients[2])
  return(freq_1.slope)
}



###################################################################################
#MAT

freq_MAT_1 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_1.csv")
mat_slope_1 <- slope_melt(freq_MAT_1)
#write_csv(mat_slope_1, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_1.csv")
rm(freq_MAT_1)
rm(mat_slope_1)

freq_MAT_2 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_2.csv")
mat_slope_2 <- slope_melt(freq_MAT_2)
#write_csv(mat_slope_2, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_2.csv")
rm(freq_MAT_2)
rm(mat_slope_2)

freq_MAT_3 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_3.csv")
mat_slope_3 <- slope_melt(freq_MAT_3)
#write_csv(mat_slope_3, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_3.csv")
rm(freq_MAT_3)
rm(mat_slope_3)

freq_MAT_4 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_4.csv")
mat_slope_4 <- slope_melt(freq_MAT_4)
#write_csv(mat_slope_4, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_4.csv")
rm(freq_MAT_4)
rm(mat_slope_4)

freq_MAT_5 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_5.csv")
mat_slope_5 <- slope_melt(freq_MAT_5)
#write_csv(mat_slope_5, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_5.csv")
rm(freq_MAT_5)
rm(mat_slope_5)

freq_MAT_6 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_6.csv")
mat_slope_6 <- slope_melt(freq_MAT_6)
#write_csv(mat_slope_6, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_6.csv")
rm(freq_MAT_6)
rm(mat_slope_6)

freq_MAT_7 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_7.csv")
mat_slope_7 <- slope_melt(freq_MAT_7)
#write_csv(mat_slope_7, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_7.csv")
rm(freq_MAT_7)
rm(mat_slope_7)

freq_MAT_8 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_8.csv")
mat_slope_8 <- slope_melt(freq_MAT_8)
#write_csv(mat_slope_8, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_8.csv")
rm(freq_MAT_8)
rm(mat_slope_8)

freq_MAT_9 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_9.csv")
mat_slope_9 <- slope_melt(freq_MAT_9)
#write_csv(mat_slope_9, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_9.csv")
rm(freq_MAT_9)
rm(mat_slope_9)

freq_MAT_10 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_10.csv")
mat_slope_10 <- slope_melt(freq_MAT_10)
#write_csv(mat_slope_10, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_10.csv")
rm(freq_MAT_10)
rm(mat_slope_10)

freq_MAT_11 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_11.csv")
mat_slope_11 <- slope_melt(freq_MAT_11)
#write_csv(mat_slope_11, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_11.csv")
rm(freq_MAT_11)
rm(mat_slope_11)

freq_MAT_12 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_12.csv")
mat_slope_12 <- slope_melt(freq_MAT_12)
#write_csv(mat_slope_12, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_12.csv")
rm(freq_MAT_12)
rm(mat_slope_12)

freq_MAT_13 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_13.csv")
mat_slope_13 <- slope_melt(freq_MAT_13)
#write_csv(mat_slope_13, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_13.csv")
rm(freq_MAT_13)
rm(mat_slope_13)

freq_MAT_14 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_14.csv")
mat_slope_14 <- slope_melt(freq_MAT_14)
#write_csv(mat_slope_14, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_14.csv")
rm(freq_MAT_14)
rm(mat_slope_14)

freq_MAT_15 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_15.csv")
mat_slope_15 <- slope_melt(freq_MAT_15)
#write_csv(mat_slope_15, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_15.csv")
rm(freq_MAT_15)
rm(mat_slope_15)

freq_MAT_16 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_16.csv")
mat_slope_16 <- slope_melt(freq_MAT_16)
#write_csv(mat_slope_16, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_16.csv")
rm(freq_MAT_16)
rm(mat_slope_16)

freq_MAT_17 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_17.csv")
mat_slope_17 <- slope_melt(freq_MAT_17)
#write_csv(mat_slope_17, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_17.csv")
rm(freq_MAT_17)
rm(mat_slope_17)

freq_MAT_18 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAT_all_slopes_18.csv")
mat_slope_18 <- slope_melt(freq_MAT_18)
#write_csv(mat_slope_18, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAT_18.csv")
rm(freq_MAT_18)
rm(mat_slope_18)



###################################################################################
#MAP

freq_MAP_1 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_1.csv")
MAP_slope_1 <- slope_melt(freq_MAP_1)
write_csv(MAP_slope_1, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_01.csv")
rm(freq_MAP_1)
rm(MAP_slope_1)

freq_MAP_2 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_2.csv")
MAP_slope_2 <- slope_melt(freq_MAP_2)
write_csv(MAP_slope_2, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_02.csv")
rm(freq_MAP_2)
rm(MAP_slope_2)

freq_MAP_3 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_3.csv")
MAP_slope_3 <- slope_melt(freq_MAP_3)
write_csv(MAP_slope_3, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_03.csv")
rm(freq_MAP_3)
rm(MAP_slope_3)

freq_MAP_4 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_4.csv")
MAP_slope_4 <- slope_melt(freq_MAP_4)
write_csv(MAP_slope_4, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_04.csv")
rm(freq_MAP_4)
rm(MAP_slope_4)

freq_MAP_5 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_5.csv")
MAP_slope_5 <- slope_melt(freq_MAP_5)
write_csv(MAP_slope_5, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_05.csv")
rm(freq_MAP_5)
rm(MAP_slope_5)

freq_MAP_6 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_6.csv")
MAP_slope_6 <- slope_melt(freq_MAP_6)
write_csv(MAP_slope_6, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_06.csv")
rm(freq_MAP_6)
rm(MAP_slope_6)

freq_MAP_7 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_7.csv")
MAP_slope_7 <- slope_melt(freq_MAP_7)
write_csv(MAP_slope_7, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_07.csv")
rm(freq_MAP_7)
rm(MAP_slope_7)

freq_MAP_8 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_8.csv")
MAP_slope_8 <- slope_melt(freq_MAP_8)
write_csv(MAP_slope_8, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_08.csv")
rm(freq_MAP_8)
rm(MAP_slope_8)

freq_MAP_9 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_9.csv")
MAP_slope_9 <- slope_melt(freq_MAP_9)
write_csv(MAP_slope_9, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_09.csv")
rm(freq_MAP_9)
rm(MAP_slope_9)

freq_MAP_10 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_10.csv")
MAP_slope_10 <- slope_melt(freq_MAP_10)
write_csv(MAP_slope_10, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_10.csv")
rm(freq_MAP_10)
rm(MAP_slope_10)

freq_MAP_11 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_11.csv")
MAP_slope_11 <- slope_melt(freq_MAP_11)
write_csv(MAP_slope_11, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_11.csv")
rm(freq_MAP_11)
rm(MAP_slope_11)

freq_MAP_12 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_12.csv")
MAP_slope_12 <- slope_melt(freq_MAP_12)
write_csv(MAP_slope_12, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_12.csv")
rm(freq_MAP_12)
rm(MAP_slope_12)

freq_MAP_13 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_13.csv")
MAP_slope_13 <- slope_melt(freq_MAP_13)
write_csv(MAP_slope_13, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_13.csv")
rm(freq_MAP_13)
rm(MAP_slope_13)

freq_MAP_14 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_14.csv")
MAP_slope_14 <- slope_melt(freq_MAP_14)
write_csv(MAP_slope_14, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_14.csv")
rm(freq_MAP_14)
rm(MAP_slope_14)

freq_MAP_15 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_15.csv")
MAP_slope_15 <- slope_melt(freq_MAP_15)
write_csv(MAP_slope_15, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_15.csv")
rm(freq_MAP_15)
rm(MAP_slope_15)

freq_MAP_16 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_16.csv")
MAP_slope_16 <- slope_melt(freq_MAP_16)
write_csv(MAP_slope_16, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_16.csv")
rm(freq_MAP_16)
rm(MAP_slope_16)

freq_MAP_17 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_17.csv")
MAP_slope_17 <- slope_melt(freq_MAP_17)
write_csv(MAP_slope_17, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_17.csv")
rm(freq_MAP_17)
rm(MAP_slope_17)

freq_MAP_18 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_MAP_all_slopes_18.csv")
MAP_slope_18 <- slope_melt(freq_MAP_18)
write_csv(MAP_slope_18, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_MAP/all_slopes_MAP_18.csv")
rm(freq_MAP_18)
rm(MAP_slope_18)



###################################################################################

#CMD
freq_CMD_1 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_1.csv")
CMD_slope_1 <- slope_melt(freq_CMD_1)
write_csv(CMD_slope_1, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_01.csv")
rm(freq_CMD_1)
rm(CMD_slope_1)

freq_CMD_2 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_2.csv")
CMD_slope_2 <- slope_melt(freq_CMD_2)
write_csv(CMD_slope_2, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_02.csv")
rm(freq_CMD_2)
rm(CMD_slope_2)

freq_CMD_3 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_3.csv")
CMD_slope_3 <- slope_melt(freq_CMD_3)
write_csv(CMD_slope_3, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_03.csv")
rm(freq_CMD_3)
rm(CMD_slope_3)

freq_CMD_4 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_4.csv")
CMD_slope_4 <- slope_melt(freq_CMD_4)
write_csv(CMD_slope_4, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_04.csv")
rm(freq_CMD_4)
rm(CMD_slope_4)

freq_CMD_5 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_5.csv")
CMD_slope_5 <- slope_melt(freq_CMD_5)
write_csv(CMD_slope_5, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_05.csv")
rm(freq_CMD_5)
rm(CMD_slope_5)

freq_CMD_6 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_6.csv")
CMD_slope_6 <- slope_melt(freq_CMD_6)
write_csv(CMD_slope_6, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_06.csv")
rm(freq_CMD_6)
rm(CMD_slope_6)

freq_CMD_7 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_7.csv")
CMD_slope_7 <- slope_melt(freq_CMD_7)
write_csv(CMD_slope_7, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_07.csv")
rm(freq_CMD_7)
rm(CMD_slope_7)

freq_CMD_8 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_8.csv")
CMD_slope_8 <- slope_melt(freq_CMD_8)
write_csv(CMD_slope_8, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_08.csv")
rm(freq_CMD_8)
rm(CMD_slope_8)

freq_CMD_9 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_9.csv")
CMD_slope_9 <- slope_melt(freq_CMD_9)
write_csv(CMD_slope_9, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_09.csv")
rm(freq_CMD_9)
rm(CMD_slope_9)

freq_CMD_10 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_10.csv")
CMD_slope_10 <- slope_melt(freq_CMD_10)
write_csv(CMD_slope_10, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_10.csv")
rm(freq_CMD_10)
rm(CMD_slope_10)

freq_CMD_11 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_11.csv")
CMD_slope_11 <- slope_melt(freq_CMD_11)
write_csv(CMD_slope_11, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_11.csv")
rm(freq_CMD_11)
rm(CMD_slope_11)

freq_CMD_12 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_12.csv")
CMD_slope_12 <- slope_melt(freq_CMD_12)
write_csv(CMD_slope_12, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_12.csv")
rm(freq_CMD_12)
rm(CMD_slope_12)

freq_CMD_13 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_13.csv")
CMD_slope_13 <- slope_melt(freq_CMD_13)
write_csv(CMD_slope_13, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_13.csv")
rm(freq_CMD_13)
rm(CMD_slope_13)

freq_CMD_14 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_14.csv")
CMD_slope_14 <- slope_melt(freq_CMD_14)
write_csv(CMD_slope_14, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_14.csv")
rm(freq_CMD_14)
rm(CMD_slope_14)

freq_CMD_15 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_15.csv")
CMD_slope_15 <- slope_melt(freq_CMD_15)
write_csv(CMD_slope_15, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_15.csv")
rm(freq_CMD_15)
rm(CMD_slope_15)

freq_CMD_16 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_16.csv")
CMD_slope_16 <- slope_melt(freq_CMD_16)
write_csv(CMD_slope_16, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_16.csv")
rm(freq_CMD_16)
rm(CMD_slope_16)

freq_CMD_17 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_17.csv")
CMD_slope_17 <- slope_melt(freq_CMD_17)
write_csv(CMD_slope_17, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_17.csv")
rm(freq_CMD_17)
rm(CMD_slope_17)

freq_CMD_18 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/freq_CMD_all_slopes_18.csv")
CMD_slope_18 <- slope_melt(freq_CMD_18)
write_csv(CMD_slope_18, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/all_slopes_CMD/all_slopes_CMD_18.csv")
rm(freq_CMD_18)
rm(CMD_slope_18)




