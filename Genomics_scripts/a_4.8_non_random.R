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
#Import observed slopes
mat_slope <- read_csv("Genomics_scripts/Data/freq_MAT_slope_peakbf5.csv") %>% mutate(env="MAT")
map_slope <- read_csv("Genomics_scripts/Data/freq_MAP_slope_peakbf5.csv") %>% mutate(env="MAP")
cmd_slope <- read_csv("Genomics_scripts/Data/freq_cmd_slope_peakbf5.csv") %>% mutate(env="CMD")


#Import files that give slopes
mat_obs_ci <- read_csv("Genomics_scripts/Data/mat_obs_ci_peakbf5.csv")
map_obs_ci <- read_csv("Genomics_scripts/Data/map_obs_ci_peakbf5.csv")
cmd_obs_ci <- read_csv("Genomics_scripts/Data/cmd_obs_ci_peakbf5.csv")



#Filter for positive slope bin MAT
mat_pos <- mat_obs_ci %>% filter(S>0)
mat_obs <- mat_pos %>% select(p1:p12) 
mat_rand <-mat_pos %>% select(p1_up,p2_up,p3_up,p4_up,p5_up,p6_up,p7_up,p8_up,p9_up,p10_up,p11_up,p12_up) 
instruction_mat <- mat_obs - mat_rand 
instruction_mat <- cbind (instruction_mat,mat_pos$S)
colnames(instruction_mat)[13] <- "S"

#Filter for positive slope bin MAP
map_pos <- map_obs_ci %>% filter(S>0)
map_obs <- map_pos %>% select(p1:p12) 
map_rand <-map_pos %>% select(p1_up,p2_up,p3_up,p4_up,p5_up,p6_up,p7_up,p8_up,p9_up,p10_up,p11_up,p12_up) 
instruction_map <- map_obs - map_rand 
instruction_map <- cbind (instruction_map,map_pos$S)
colnames(instruction_map)[13] <- "S"

#Filter for positive slope bin CMD
cmd_pos <- cmd_obs_ci %>% filter(S>0)
cmd_obs <- cmd_pos %>% select(p1:p12) 
cmd_rand <-cmd_pos %>% select(p1_up,p2_up,p3_up,p4_up,p5_up,p6_up,p7_up,p8_up,p9_up,p10_up,p11_up,p12_up) 
instruction_cmd <- cmd_obs - cmd_rand 
instruction_cmd <- cbind (instruction_cmd,cmd_pos$S)
colnames(instruction_cmd)[13] <- "S"



#####################################################
#Get non-random slopes

#mat
slope.out.mat <- data.frame()
for (j in 1:12){
  instruction_env <- instruction_mat %>% select(eval(paste("p",j, sep="")),S)
  instruction_env <- instruction_env %>% filter(instruction_env[,1]>0)
  env_slope_pop <- mat_slope %>% filter(Site==j)
  
  if(dim(instruction_env)[1]!=0){
  instruction_env <- instruction_env %>%  mutate(low_S=S-0.1,high_S=S+0.1)
  
  for (i in 1:2){
    data.temp <- instruction_env[i,]
    slope.temp <- env_slope_pop %>% filter(Slope >= data.temp$low_S & Slope < data.temp$high_S) %>% 
      mutate(Site=j)
    slope.out.mat <- rbind(slope.out.mat,slope.temp)
  }
}
}

#map
slope.out.map <- data.frame()
for (j in 1:12){
  instruction_env <- instruction_map %>% select(eval(paste("p",j, sep="")),S)
  instruction_env <- instruction_env %>% filter(instruction_env[,1]>0)
  env_slope_pop <- map_slope %>% filter(Site==j)

  if(dim(instruction_env)[1]!=0){
    instruction_env <- instruction_env %>%  mutate(low_S=S-0.1,high_S=S+0.1)
    
    for (i in 1:dim(instruction_env)[1]){
      data.temp <- instruction_env[i,]
      slope.temp <- env_slope_pop %>% filter(Slope >= data.temp$low_S & Slope < data.temp$high_S) %>% 
        mutate(Site=j)
      slope.out.map <- rbind(slope.out.map,slope.temp)
    }
  }  
}

#cmd
slope.out.cmd <- data.frame()
for (j in 1:12){
  instruction_env <- instruction_cmd %>% select(eval(paste("p",j, sep="")),S)
  instruction_env <- instruction_env %>% filter(instruction_env[,1]>0)
  env_slope_pop <- cmd_slope %>% filter(Site==j)
  
  if(dim(instruction_env)[1]!=0){
    instruction_env <- instruction_env %>%  mutate(low_S=S-0.1,high_S=S+0.1)
    
    for (i in 1:dim(instruction_env)[1]){
      data.temp <- instruction_env[i,]
      slope.temp <- env_slope_pop %>% filter(Slope >= data.temp$low_S & Slope < data.temp$high_S) %>% 
        mutate(Site=j)
      slope.out.cmd <- rbind(slope.out.cmd,slope.temp)
    }
  }  
}

snp_list <- rbind(slope.out.mat,slope.out.map,slope.out.cmd)



############################################################################################################
#Make unique snp_list across env

#Diagnostics
snp_list_p1 <- snp_list %>% filter(Site==1)
#Not unique SNP ID
unique(snp_list_p1$snp_ID[duplicated(snp_list_p1$snp_ID)]) 

#Remove not unique SNPs
snp_list_filter <- data.frame()
for(i in 1:12){
  snp_list_p <- snp_list %>% filter(Site==i)
  snp_list_p_filtered <- snp_list_p %>% filter(duplicated(snp_ID) == FALSE)
  snp_list_filter <- rbind(snp_list_filter,snp_list_p_filtered)
}




write_csv(snp_list,"Genomics_scripts/Data/snp_list.csv")
write_csv(snp_list_filter,"Genomics_scripts/Data/snp_list_unique.csv")




############################################################################################################


