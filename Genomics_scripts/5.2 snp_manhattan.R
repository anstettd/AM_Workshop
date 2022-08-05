#############################################################################################################
## Make manhattan plots for WZA snps
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
#Plots on SNP data
# Make histogram
ggplot(data = snps_mat[snps_mat$snp_c >800000,],
       aes(x = p_val))+
  geom_histogram( bins = 50)+
  scale_y_continuous("Count")+
  scale_x_continuous(expression(italic(p)*"-value"))+
  theme_bw()

ggplot(data = snps_mat[snps_map$snp_c >800000,],
       aes(x = p_val))+
  geom_histogram( bins = 50)+
  scale_y_continuous("Count")+
  scale_x_continuous(expression(italic(p)*"-value"))+
  theme_bw()

ggplot(data = snps_mat[snps_cmd$snp_c >800000,],
       aes(x = p_val))+
  geom_histogram( bins = 50)+
  scale_y_continuous("Count")+
  scale_x_continuous(expression(italic(p)*"-value"))+
  theme_bw()


#Manhattan
ggplot(data = snps_mat, 
       aes(x = snp_c/1e6, 
           y = -log10(p_val),
           fill = p_bar*q_bar))+
  geom_point(shape = 21)+
  geom_hline(aes(yintercept = quantile(-log10(p_val), 0.99)), col = "red", lty = 2, lwd = 1)+
  scale_fill_gradient(expression(italic(bar(p)*bar(q))), low = "white", high = "black")+
  scale_y_continuous(expression(-log[10]*"(p-value)"), limits = c(0,10))+
  scale_x_continuous("Position (Mbp)")+
  theme_classic()
ggsave("myplot.png")

ggplot(data = snps_map, 
       aes(x = snp_c/1e6, 
           y = -log10(p_val),
           fill = p_bar*q_bar))+
  geom_point(shape = 21)+
  geom_hline(aes(yintercept = quantile(-log10(p_val), 0.99)), col = "red", lty = 2, lwd = 1)+
  scale_fill_gradient(expression(italic(bar(p)*bar(q))), low = "white", high = "black")+
  scale_y_continuous(expression(-log[10]*"(p-value)"), limits = c(0,10))+
  scale_x_continuous("Position (Mbp)")+
  theme_classic()

ggplot(data = snps_cmd, 
       aes(x = snp_c/1e6, 
           y = -log10(p_val),
           fill = p_bar*q_bar))+
  geom_point(shape = 21)+
  geom_hline(aes(yintercept = quantile(-log10(p_val), 0.99)), col = "red", lty = 2, lwd = 1)+
  scale_fill_gradient(expression(italic(bar(p)*bar(q))), low = "white", high = "black")+
  scale_y_continuous(expression(-log[10]*"(p-value)"), limits = c(0,10))+
  scale_x_continuous("Position (Mbp)")+
  theme_classic()

