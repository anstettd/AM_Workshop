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

geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=0.5) +

#Manhattan
man_snp_mat <- ggplot(data = snps_mat, 
       aes(x = snp_c/1e6, 
           y = -log10(p_val),
           fill = p_bar*q_bar))+
  geom_point(aes(color=as.factor(chr), alpha=0.8, size=0.5))+
  scale_color_manual(values = rep(c("black", "darkgoldenrod"), 22 )) +
  geom_hline(aes(yintercept = quantile(-log10(p_val), 0.99)), col = "red", lty = 2, lwd = 1)+
#  scale_fill_gradient(expression(italic(bar(p)*bar(q))), low = "white", high = "black")+
  scale_y_continuous(expression(-log[10]*"(p-value)"), limits = c(0,10))+
  scale_x_continuous("Position (Mbp)")+
  theme_classic()+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=15, face="bold"),
    axis.text.y = element_text(size=15,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 0.8, face="bold",hjust=0.5)
  )
ggsave("/Users/daniel_anstett/Dropbox/a_Papers/Genomics_paper/Graphs/WZA_Graphs/man_snp_mat.png",
       man_snp_mat, width=10, height = 5, units = "in")

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

