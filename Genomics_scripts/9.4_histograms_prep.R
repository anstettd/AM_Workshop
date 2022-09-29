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
  env_range1 <- env_slope %>% filter(Slope >= -2 & Slope< -1.8)
  env_range2 <- env_slope %>% filter(Slope >= -1.8 & Slope < -1.6)
  env_range3 <- env_slope %>% filter(Slope >= -1.6 & Slope < -1.4)
  env_range4 <- env_slope %>% filter(Slope >= -1.4 & Slope < -1.2)
  env_range5 <- env_slope %>% filter(Slope >= -1.2 & Slope < -1)
  env_range6 <- env_slope %>% filter(Slope >= -1 & Slope < -0.8)
  env_range7 <- env_slope %>% filter(Slope >= -0.8 & Slope < -0.6)
  env_range8 <- env_slope %>% filter(Slope >= -0.6 & Slope < -0.4)
  env_range9 <- env_slope %>% filter(Slope >= -0.4 & Slope < -0.2)
  env_range10 <- env_slope %>% filter(Slope >= -0.2 & Slope < 0)
  env_range11 <- env_slope %>% filter(Slope >= 0 & Slope < 0.2)
  env_range12 <- env_slope %>% filter(Slope >= 0.2 & Slope < 0.4)
  env_range13 <- env_slope %>% filter(Slope >= 0.3 & Slope < 0.6)
  env_range14 <- env_slope %>% filter(Slope >= 0.4 & Slope < 0.8)
  env_range15 <- env_slope %>% filter(Slope >= 0.8 & Slope < 1)
  env_range16 <- env_slope %>% filter(Slope >= 1 & Slope < 1.2)
  env_range17 <- env_slope %>% filter(Slope >= 1.2 & Slope < 1.4)
  env_range18 <- env_slope %>% filter(Slope >= 1.4 & Slope < 1.6)
  env_range19 <- env_slope %>% filter(Slope >= 1.6 & Slope < 1.8)
  env_range20 <- env_slope %>% filter(Slope >= 1.8 & Slope < 2)
  
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
  env_obs[12,1] <- dim(env_range12)[1]
  env_obs[13,1] <- dim(env_range13)[1]
  env_obs[14,1] <- dim(env_range14)[1]
  env_obs[15,1] <- dim(env_range15)[1]
  env_obs[16,1] <- dim(env_range16)[1]
  env_obs[17,1] <- dim(env_range17)[1]
  env_obs[18,1] <- dim(env_range18)[1]
  env_obs[19,1] <- dim(env_range19)[1]
  env_obs[20,1] <- dim(env_range20)[1]
  return(env_obs)
}





###################################################################################
#Import observed slopes
mat_slope <- read_csv("Genomics_scripts/Data/freq_MAT_slope_peakbf2_region.csv")
map_slope <- read_csv("Genomics_scripts/Data/freq_MAP_slope_peakbf2_region.csv")
cmd_slope <- read_csv("Genomics_scripts/Data/freq_cmd_slope_peakbf2_region.csv")

setwd("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files")

#Import large files with 1000 X random slopes
mat_bars <- read_csv("rand_slope_mat_multi_peakbf2_region.csv")
map_bars <- read_csv("rand_slope_map_multi_peakbf2_region.csv")
cmd_bars <- read_csv("rand_slope_cmd_multi_peakbf2_region.csv")

setwd("/Users/daniel_anstett/Dropbox/AM_Workshop/AM_Workshop")


###################################################################################
## Organize slope data
###################################################################################
#Filter data within s = -2 to 2
mat_slope2 <- mat_slope %>% filter(Slope<2 & Slope>-2) 
map_slope2 <- map_slope %>% filter(Slope<2 & Slope>-2) 
cmd_slope2 <- cmd_slope %>% filter(Slope<2 & Slope>-2) 

mat_bars2 <- mat_bars %>% filter(Slope<2 & Slope>-2) 
map_bars2 <- map_bars %>% filter(Slope<2 & Slope>-2) 
cmd_bars2 <- cmd_bars %>% filter(Slope<2 & Slope>-2) 

#Setup data frame
env_obs <- as.data.frame(c(-1.9,-1.7,-1.5,-1.3,-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,
                              0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9))

###################################################################################
#Run get range function
#MAT
mat_obs1 <- get_range(mat_slope2,"North")
mat_obs2 <- get_range(mat_slope2,"Centre")
mat_obs3 <- get_range(mat_slope2,"South")
mat_obs <- cbind(env_obs,mat_obs1,mat_obs2,mat_obs3)
colnames(mat_obs) = c("S","North","Centre","South")

#map
map_obs1 <- get_range(map_slope2,"North")
map_obs2 <- get_range(map_slope2,"Centre")
map_obs3 <- get_range(map_slope2,"South")
map_obs <- cbind(env_obs,map_obs1,map_obs2,map_obs3)
colnames(map_obs) = c("S","North","Centre","South")

#cmd
cmd_obs1 <- get_range(cmd_slope2,"North")
cmd_obs2 <- get_range(cmd_slope2,"Centre")
cmd_obs3 <- get_range(cmd_slope2,"South")
cmd_obs <- cbind(env_obs,cmd_obs1,cmd_obs2,cmd_obs3)
colnames(cmd_obs) = c("S","North","Centre","South")

###############################################################################################
# Get CI
###############################################################################################
#Example of getting ci 
regions <- c("North","Centre","South")
mat_ci_prep <- data.frame()
for(j in 1:3){
  for(i in 0:19){
    mat_bars <- mat_bars2 %>% filter(Site==regions[j])
    mat_test <- mat_bars %>% filter(Slope >= (-2+0.2*i) & Slope< (-1.8+0.2*i))
    mat_count <- mat_test %>% count(Seed_ID)
    mat_ci_prep[1+i,j*2-1] <- as.numeric(quantile(mat_count$n, probs = c(0.025, 0.975))[1])
    mat_ci_prep[1+i,j*2] <- as.numeric(quantile(mat_count$n, probs = c(0.025, 0.975))[2])
  }
}

map_ci_prep <- data.frame()
for(j in 1:3){
  for(i in 0:19){
    map_bars <- map_bars2 %>% filter(Site==regions[j])
    map_test <- map_bars %>% filter(Slope >= (-2+0.2*i) & Slope< (-1.8+0.2*i))
    map_count <- map_test %>% count(Seed_ID)
    map_ci_prep[1+i,j*2-1] <- as.numeric(quantile(map_count$n, probs = c(0.025, 0.975))[1])
    map_ci_prep[1+i,j*2] <- as.numeric(quantile(map_count$n, probs = c(0.025, 0.975))[2])
  }
}

cmd_ci_prep <- data.frame()
for(j in 1:3){
  for(i in 0:19){
    cmd_bars <- cmd_bars2 %>% filter(Site==regions[j])
    cmd_test <- cmd_bars %>% filter(Slope >= (-2+0.2*i) & Slope< (-1.8+0.2*i))
    cmd_count <- cmd_test %>% count(Seed_ID)
    cmd_ci_prep[1+i,j*2-1] <- as.numeric(quantile(cmd_count$n, probs = c(0.025, 0.975))[1])
    cmd_ci_prep[1+i,j*2] <- as.numeric(quantile(cmd_count$n, probs = c(0.025, 0.975))[2])
  }
}

colnames(mat_ci_prep) <- c("p1_low","p1_up","p2_low","p2_up","p3_low","p3_up")
colnames(map_ci_prep) <- c("p1_low","p1_up","p2_low","p2_up","p3_low","p3_up")
colnames(cmd_ci_prep) <- c("p1_low","p1_up","p2_low","p2_up","p3_low","p3_up")
                           

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


#Export
write_csv(mat_obs_ci, "mat_obs_ci_peakbf2_region.csv")
write_csv(map_obs_ci, "map_obs_ci_peakbf2_region.csv")
write_csv(cmd_obs_ci, "cmd_obs_ci_peakbf2_region.csv")
write_csv(env_obs_ci, "env_obs_ci_peakbf2_region.csv")



