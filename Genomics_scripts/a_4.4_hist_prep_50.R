##################################################################################
## Generate input for strength of selection graphs
## Processing for both observation and permuted CI
## Author Daniel Anstett
## 
## 
## Last Modified May 17, 2022
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

###################################################################################
#Function for generating stratigied random distribution
###################################################################################

#get_range
get_range <- function(df,site){
  # env Example
  env_slope <- df %>% filter(Site==site)
  
  #Make dataframe
  env_obs <- data.frame()

  #Get ranges  
  env_range1 <- env_slope %>% filter(Slope >= -55 & Slope< -45)
  env_range2 <- env_slope %>% filter(Slope >= -45 & Slope < -35)
  env_range3 <- env_slope %>% filter(Slope >= -35 & Slope < -25)
  env_range4 <- env_slope %>% filter(Slope >= -25 & Slope< -15)
  env_range5 <- env_slope %>% filter(Slope >= -15 & Slope < -5)
  env_range6 <- env_slope %>% filter(Slope >= -5 & Slope < 5)
  env_range7 <- env_slope %>% filter(Slope >= 5 & Slope < 15)
  env_range8 <- env_slope %>% filter(Slope >= 15 & Slope < 25)
  env_range9 <- env_slope %>% filter(Slope >= 25 & Slope < -35)
  env_range10 <- env_slope %>% filter(Slope >= 35 & Slope < 45)
  env_range11 <- env_slope %>% filter(Slope >= 45 & Slope < 55)

  
  #get dim and put into dataframe
  env_obs[1,1] <- dim(env_range1)[1]
  env_obs[2,1] <- dim(env_range2)[1]
  env_obs[3,1] <- dim(env_range3)[1]
  env_obs[4,1] <- dim(env_range4)[1]
  env_obs[5,1] <- dim(env_range5)[1]
  env_obs[6,1] <- dim(env_range6)[1]
  env_obs[7,1] <- dim(env_range7)[1]
  env_obs[8,1] <- dim(env_range8)[1]
  env_obs[9,1] <- dim(env_range9)[1]
  env_obs[10,1] <- dim(env_range10)[1]
  env_obs[11,1] <- dim(env_range11)[1]
  return(env_obs)
}


###################################################################################
#Import observed slopes
mat_slope2 <- read_csv("Genomics_scripts/Data/freq_MAT_slope_peakbf5.csv")
map_slope2 <- read_csv("Genomics_scripts/Data/freq_MAP_slope_peakbf5.csv")
cmd_slope2 <- read_csv("Genomics_scripts/Data/freq_cmd_slope_peakbf5.csv")

setwd("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files")

#Import large files with 1000 X random slopes
mat_bars2 <- read_csv("rand_slope_mat_multi_peakbf5.csv")
map_bars2 <- read_csv("rand_slope_map_multi_peakbf5.csv")
cmd_bars2 <- read_csv("rand_slope_cmd_multi_peakbf5.csv")

setwd("/Users/daniel_anstett/Dropbox/AM_Workshop/AM_Workshop")

#Setup data frame
env_obs <- as.data.frame(c(-50,-40,-30,-20,-10,0,10,20,30,40,50))

###################################################################################
#Run get range function
#mat
mat_obs1 <- get_range(mat_slope2,1)
mat_obs2 <- get_range(mat_slope2,2)
mat_obs3 <- get_range(mat_slope2,3)
mat_obs4 <- get_range(mat_slope2,4)
mat_obs5 <- get_range(mat_slope2,5)
mat_obs6 <- get_range(mat_slope2,6)
mat_obs7 <- get_range(mat_slope2,7)
mat_obs8 <- get_range(mat_slope2,8)
mat_obs9 <- get_range(mat_slope2,9)
mat_obs10 <- get_range(mat_slope2,10)
mat_obs11 <- get_range(mat_slope2,11)
mat_obs12 <- get_range(mat_slope2,12)
mat_obs <- cbind(env_obs,mat_obs1,mat_obs2,mat_obs3,mat_obs4,mat_obs5,mat_obs6,
                  mat_obs7,mat_obs8,mat_obs9,mat_obs10,mat_obs11,mat_obs12)
colnames(mat_obs) = c("S","p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12")

#map
map_obs1 <- get_range(map_slope2,1)
map_obs2 <- get_range(map_slope2,2)
map_obs3 <- get_range(map_slope2,3)
map_obs4 <- get_range(map_slope2,4)
map_obs5 <- get_range(map_slope2,5)
map_obs6 <- get_range(map_slope2,6)
map_obs7 <- get_range(map_slope2,7)
map_obs8 <- get_range(map_slope2,8)
map_obs9 <- get_range(map_slope2,9)
map_obs10 <- get_range(map_slope2,10)
map_obs11 <- get_range(map_slope2,11)
map_obs12 <- get_range(map_slope2,12)
map_obs <- cbind(env_obs,map_obs1,map_obs2,map_obs3,map_obs4,map_obs5,map_obs6,
                 map_obs7,map_obs8,map_obs9,map_obs10,map_obs11,map_obs12)
colnames(map_obs) = c("S","p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12")

#cmd
cmd_obs1 <- get_range(cmd_slope2,1)
cmd_obs2 <- get_range(cmd_slope2,2)
cmd_obs3 <- get_range(cmd_slope2,3)
cmd_obs4 <- get_range(cmd_slope2,4)
cmd_obs5 <- get_range(cmd_slope2,5)
cmd_obs6 <- get_range(cmd_slope2,6)
cmd_obs7 <- get_range(cmd_slope2,7)
cmd_obs8 <- get_range(cmd_slope2,8)
cmd_obs9 <- get_range(cmd_slope2,9)
cmd_obs10 <- get_range(cmd_slope2,10)
cmd_obs11 <- get_range(cmd_slope2,11)
cmd_obs12 <- get_range(cmd_slope2,12)
cmd_obs <- cbind(env_obs,cmd_obs1,cmd_obs2,cmd_obs3,cmd_obs4,cmd_obs5,cmd_obs6,
                 cmd_obs7,cmd_obs8,cmd_obs9,cmd_obs10,cmd_obs11,cmd_obs12)
colnames(cmd_obs) = c("S","p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12")

###############################################################################################
# Get CI
###############################################################################################
#Example of getting ci 

mat_ci_prep <- data.frame()
for(j in 1:12){
  for(i in 0:10){
    mat_bars <- mat_bars2 %>% filter(Site==j)
    mat_test <- mat_bars %>% filter(Slope >= (-55+10*i) & Slope< (-45+10*i))
    mat_count <- mat_test %>% count(Seed_ID)
    mat_ci_prep[1+i,j*2-1] <- as.numeric(quantile(mat_count$n, probs = c(0.025, 0.975))[1])
    mat_ci_prep[1+i,j*2] <- as.numeric(quantile(mat_count$n, probs = c(0.025, 0.975))[2])
  }
}

map_ci_prep <- data.frame()
for(j in 1:12){
  for(i in 0:10){
    map_bars <- map_bars2 %>% filter(Site==j)
    map_test <- map_bars %>% filter(Slope >= (-55+10*i) & Slope< (-45+10*i))
    map_count <- map_test %>% count(Seed_ID)
    map_ci_prep[1+i,j*2-1] <- as.numeric(quantile(map_count$n, probs = c(0.025, 0.975))[1])
    map_ci_prep[1+i,j*2] <- as.numeric(quantile(map_count$n, probs = c(0.025, 0.975))[2])
  }
}

cmd_ci_prep <- data.frame()
for(j in 1:12){
  for(i in 0:10){
    cmd_bars <- cmd_bars2 %>% filter(Site==j)
    cmd_test <- cmd_bars %>% filter(Slope >= (-55+10*i) & Slope< (-45+10*i))
    cmd_count <- cmd_test %>% count(Seed_ID)
    cmd_ci_prep[1+i,j*2-1] <- as.numeric(quantile(cmd_count$n, probs = c(0.025, 0.975))[1])
    cmd_ci_prep[1+i,j*2] <- as.numeric(quantile(cmd_count$n, probs = c(0.025, 0.975))[2])
  }
}

colnames(mat_ci_prep) <- c("p1_low","p1_up","p2_low","p2_up","p3_low","p3_up","p4_low","p4_up"
                           ,"p5_low","p5_up","p6_low","p6_up","p7_low","p7_up","p8_low","p8_up"
                           ,"p9_low","p9_up","p10_low","p10_up","p11_low","p11_up","p12_low","p12_up")
colnames(map_ci_prep) <- c("p1_low","p1_up","p2_low","p2_up","p3_low","p3_up","p4_low","p4_up"
                           ,"p5_low","p5_up","p6_low","p6_up","p7_low","p7_up","p8_low","p8_up"
                           ,"p9_low","p9_up","p10_low","p10_up","p11_low","p11_up","p12_low","p12_up")
colnames(cmd_ci_prep) <- c("p1_low","p1_up","p2_low","p2_up","p3_low","p3_up","p4_low","p4_up"
                           ,"p5_low","p5_up","p6_low","p6_up","p7_low","p7_up","p8_low","p8_up"
                           ,"p9_low","p9_up","p10_low","p10_up","p11_low","p11_up","p12_low","p12_up")
                           

###############################################################################################
# Merge data frames and export
mat_obs_ci <- cbind(mat_obs,mat_ci_prep)
map_obs_ci <- cbind(map_obs,map_ci_prep)
cmd_obs_ci <- cbind(cmd_obs,cmd_ci_prep)

#Name each dataframe as a variable
mat_obs_ci$env <- "MAT"
map_obs_ci$env <- "MAP"
cmd_obs_ci$env <- "CMD"

#Bind into single data frame
env_obs_ci <- rbind(mat_obs_ci,map_obs_ci,cmd_obs_ci)

#Make NAs zeros
env_obs_ci[is.na(env_obs_ci)] <- 0
mat_obs_ci[is.na(mat_obs_ci)] <- 0
map_obs_ci[is.na(map_obs_ci)] <- 0
cmd_obs_ci[is.na(cmd_obs_ci)] <- 0



#Export
write_csv(mat_obs_ci, "Genomics_scripts/Data/mat_obs_ci_peakbf5_slope50.csv")
write_csv(map_obs_ci, "Genomics_scripts/Data/map_obs_ci_peakbf5_slope50.csv")
write_csv(cmd_obs_ci, "Genomics_scripts/Data/cmd_obs_ci_peakbf5_slope50.csv")
write_csv(env_obs_ci, "Genomics_scripts/Data/env_obs_ci_peakbf5_slope50.csv")



