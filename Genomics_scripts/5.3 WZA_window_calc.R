#############################################################################################################
## Calc WZA for windows
## Author Daniel Anstett
## 
## 
## Modified from Tom Booker WZA Vignette
## Last Modified August 3, 2022
#############################################################################################################
#Import libraries

library(tidyverse)
library(Kendall)

#Import files

snps_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_mat.csv")
snps_map <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_map.csv")
snps_cmd <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_cmd.csv")


###########################################################################################################





###########################################################################################################
#Calc WZA


# Now convert the empirical p-values to z scores
snps_mat$z_score <- qnorm(snps_mat$empirical_p, lower.tail = F)
snps_map$z_score <- qnorm(snps_map$empirical_p, lower.tail = F)
snps_cmd$z_score <- qnorm(snps_cmd$empirical_p, lower.tail = F)


# Calculate p_bar*q_bar
snps_mat$pbar_qbar <- snps_mat$p_bar*snps_mat$q_bar
snps_map$pbar_qbar <- snps_map$p_bar*snps_mat$q_bar
snps_cmd$pbar_qbar <- snps_cmd$p_bar*snps_mat$q_bar


# Calculate the WZA for each window
#MAT
WZA_numerator_mat <- tapply(snps_mat$z_score*snps_mat$pbar_qbar, snps_mat$win, sum)
WZA_denominator_mat <- sqrt( tapply((snps_mat$pbar_qbar*snps_mat$pbar_qbar), snps_mat$win, sum) )
WZA_mat <- WZA_numerator_mat/WZA_denominator_mat
pos_mat <- tapply(snps_mat$snp_c, snps_mat$win, mean)
win_mat <- snps_mat %>% distinct(chr, win) #Get window and chromosome info into WZA data frame
WZA_df_mat_df <- data.frame(WZA = WZA_mat, pos = pos_mat)
WZA_df_mat <- cbind(win_mat,WZA_df_mat_df)

#MAP
WZA_numerator_map <- tapply(snps_map$z_score*snps_map$pbar_qbar, snps_map$win, sum)
WZA_denominator_map <- sqrt( tapply((snps_map$pbar_qbar*snps_map$pbar_qbar), snps_map$win, sum) )
WZA_map <- WZA_numerator_map/WZA_denominator_map
pos_map <- tapply(snps_map$snp_c, snps_map$win, mean)
win_map <- snps_map %>% distinct(chr, win) #Get window and chromosome info into WZA data frame
WZA_df_map_df <- data.frame(WZA = WZA_map, pos = pos_map)
WZA_df_map <- cbind(win_map,WZA_df_map_df)

#CMD
WZA_numerator_cmd <- tapply(snps_cmd$z_score*snps_cmd$pbar_qbar, snps_cmd$win, sum)
WZA_denominator_cmd <- sqrt( tapply((snps_cmd$pbar_qbar*snps_cmd$pbar_qbar), snps_cmd$win, sum) )
WZA_cmd <- WZA_numerator_cmd/WZA_denominator_cmd
pos_cmd <- tapply(snps_cmd$snp_c, snps_cmd$win, mean)
win_cmd <- snps_cmd %>% distinct(chr, win) #Get window and chromosome info into WZA data frame
WZA_df_cmd_df <- data.frame(WZA = WZA_cmd, pos = pos_cmd)
WZA_df_cmd <- cbind(win_cmd,WZA_df_cmd_df)

###########################################################################################################
#Calc empirical p-value based on WZA scores

#MAT
WZA_mean_mat <- mean(WZA_df_mat$WZA) ## Calculate the mean and sd of the WZA distribution
WZA_sd_mat <- sd(WZA_df_mat$WZA)
WZA_df_mat$approx_p <- 2*pnorm(-abs(WZA_df_mat$WZA), mean = WZA_mean_mat, sd= WZA_sd_mat) ## Calculate an approximate p-value based on the assumption of normality

#MAP
WZA_mean_map <- mean(WZA_df_map$WZA) ## Calculate the mean and sd of the WZA distribution
WZA_sd_map <- sd(WZA_df_map$WZA)
WZA_df_map$approx_p <- 2*pnorm(-abs(WZA_df_map$WZA), mean = WZA_mean_map, sd= WZA_sd_map) ## Calculate an approximate p-value based on the assumption of normality

#CMD
WZA_mean_cmd <- mean(WZA_df_cmd$WZA) ## Calculate the mean and sd of the WZA distribution
WZA_sd_cmd <- sd(WZA_df_cmd$WZA)
WZA_df_cmd$approx_p <- 2*pnorm(-abs(WZA_df_cmd$WZA), mean = WZA_mean_cmd, sd= WZA_sd_cmd) ## Calculate an approximate p-value based on the assumption of normality


write_csv(WZA_df_mat, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_win_mat.csv")
write_csv(WZA_df_map, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_win_map.csv")
write_csv(WZA_df_cmd, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_win_cmd.csv")


