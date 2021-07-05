##################################################################################
## Merging of data from sequencing centre (NanuQ) with locality data
## Author Daniel Anstett
## 
##
## Last Modified July 4, 2021
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)
library(ape)

#Import data
omega_1 <- read_csv("Genomics_scripts/Data/core_model_rand_1_mat_omega.csv",col_names = FALSE)
omega_2 <- read_csv("Genomics_scripts/Data/core_model_rand_2_mat_omega.csv",col_names = FALSE)
omega_3 <- read_csv("Genomics_scripts/Data/core_model_rand_3_mat_omega.csv",col_names = FALSE)

#Change dataframe into matrix
omega_1 <- as.matrix(omega_1)
omega_2 <- as.matrix(omega_2)
omega_3 <- as.matrix(omega_3)

#Matrix correlation
cor.test(omega_1,omega_2)
cor.test(omega_1,omega_3)
cor.test(omega_2,omega_3)


#Mantel Tests
mantel.test(omega_1,omega_2,nperm=999)
mantel.test(omega_2,omega_3,nperm=999)
mantel.test(omega_3,omega_1,nperm=999)




