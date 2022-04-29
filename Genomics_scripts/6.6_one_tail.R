##################################################################################
## Make SNP slopes histograms
## Author Daniel Anstett
## 
## Currently done just for ENV 1,2 and 5 (MAT, MAP and CMD)
## Last Modified Apr 29, 2022
###################################################################################
#Import libraries
library(tidyverse)

###################################################################################

snp1A_slope <- read_csv("Genomics_scripts/Data/freq_MAT_slope.csv") %>% 
  mutate(Clim="MAT",
         Type="Observed") %>% 
  filter(Slope<2) %>% 
  filter(Slope>-2)
snp2A_slope <- read_csv("Genomics_scripts/Data/freq_MAP_slope.csv") %>% 
  mutate(Clim="MAP",
         Type="Observed") %>% 
  filter(Slope<2) %>% 
  filter(Slope>-2)
snp5A_slope <- read_csv("Genomics_scripts/Data/freq_CMD_slope.csv") %>% 
  mutate(Clim="CMD",
         Type="Observed") %>% 
  filter(Slope<2) %>% 
  filter(Slope>-2)

snp1A_rslope <- read_csv("Genomics_scripts/Data/rand_slope_mat.csv") %>% 
  mutate(Clim="MAT",
         Type="Random") %>% 
  filter(Slope<2) %>% 
  filter(Slope>-2)
snp2A_rslope <- read_csv("Genomics_scripts/Data/rand_slope_map.csv") %>% 
  mutate(Clim="MAP",
         Type="Random") %>% 
  filter(Slope<2) %>% 
  filter(Slope>-2)
snp5A_rslope <- read_csv("Genomics_scripts/Data/rand_slope_cmd.csv") %>% 
  mutate(Clim="CMD",
         Type="Random") %>% 
  filter(Slope<2) %>% 
  filter(Slope>-2)

dat <- bind_rows(snp1A_slope, snp2A_slope, snp5A_slope, snp1A_rslope, snp2A_rslope, snp5A_rslope)

clim_vec <- c("MAT", "MAP", "CMD")

# Singleton example
#snp1A_slope_1 <- snp1A_slope %>% filter(Site==1)
#wilcox.test(snp1A_slope_1$Slope, mu = 0, alternative = "greater")

# Loop to test whether median of climate-associated loci is >0 (one-tailed test; null is that median is <= 0)
wil_p <- data.frame()
for (i in 1:12){
  dat_i <- dat %>% filter(Site==i) %>% filter(Type=="Observed") 
  for (j in 1:3){
    dat_ij <- dat_i %>% filter(Clim==clim_vec[j])
    wil_test <- wilcox.test(dat_ij$Slope, mu = 0, alternative = "greater")
    wil_p[i,j] <-wil_test$p.value
  }
}

colnames(wil_p) <- clim_vec

#write_csv(wil_p, "Genomics_scripts/Data/wilcox_slope_p-values.csv")
write_csv(wil_p, "Genomics_scripts/Data/wilcox_slope_p-values_trim2.csv")

# Loop to test whether mean of random loci is different from 0 (two-tailed test; null is that mean is = 0)
wil_p_rand <- data.frame()

for (i in 1:12){
  dat_i <- dat %>% filter(Site==i) %>% filter(Type=="Random") 
  for (j in 1:3){
    dat_ij <- dat_i %>% filter(Clim==clim_vec[j])
    wil_rtest <- wilcox.test(dat_ij$Slope, mu = 0)
    wil_p_rand[i,j] <-wil_rtest$p.value
 }
}
colnames(wil_p_rand) <- clim_vec

#write_csv(wil_p_rand, "Genomics_scripts/Data/wilcox_slope_random_p-values.csv")
write_csv(wil_p_rand, "Genomics_scripts/Data/wilcox_slope_random_p-values_trim2.csv")

# Loops to test for random loci, when mean =/= 0, which way is it going?
wil_p_rand <- data.frame()

for (i in 1:12){
  dat_i <- dat %>% filter(Site==i) %>% filter(Type=="Random") 
  for (j in 1:3){
    dat_ij <- dat_i %>% filter(Clim==clim_vec[j])
    wil_rtest <- wilcox.test(dat_ij$Slope, mu = 0, alternative="greater")
    wil_p_rand[i,j] <-wil_rtest$p.value
  }
}

for (i in 1:12){
  dat_i <- dat %>% filter(Site==i) %>% filter(Type=="Random") 
  for (j in 1:3){
    dat_ij <- dat_i %>% filter(Clim==clim_vec[j])
    wil_rtest <- wilcox.test(dat_ij$Slope, mu = 0, alternative="less")
    wil_p_rand[i,j] <-wil_rtest$p.value
  }
}

# Two-sample tests for whether observed slopes are greater than random
wil_p_twosamp <- data.frame()

for (i in 1:12){
  dat_i <- dat %>% filter(Site==i) 
  for (j in 1:3){
    dat_ij <- dat_i %>% filter(Clim==clim_vec[j])
    wil_2test <- wilcox.test(dat_ij$Slope[dat_ij$Type=="Observed"], 
                             dat_ij$Slope[dat_ij$Type=="Random"], 
                             mu = 0, alternative="greater")
    wil_p_twosamp[i,j] <-wil_2test$p.value
  }
}

colnames(wil_p_twosamp) <- clim_vec

#write_csv(wil_p_twosamp, "Genomics_scripts/Data/wilcox_slope_twosamp_p-values.csv")
write_csv(wil_p_twosamp, "Genomics_scripts/Data/wilcox_slope_twosamp_p-values_trim2.csv")
