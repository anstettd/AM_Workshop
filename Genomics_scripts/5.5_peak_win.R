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
library(plotly)

#Import files
snps_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_mat.csv")
snps_map <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_map.csv")
snps_cmd <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_cmd.csv")

wza_win_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_win_mat.csv")
wza_win_map <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_win_map.csv")
wza_win_cmd <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_win_cmd.csv")

#Filter by Bonferonnii correction alpha critical = 1.29423e-06, aka 5.887988 sigma
mat_bon <- wza_win_mat %>% filter(approx_p<1.29423e-06)
map_bon <- wza_win_map %>% filter(approx_p<1.29423e-06)
cmd_bon <- wza_win_cmd %>% filter(approx_p<1.29423e-06)

#Get peak windows only
mat_bon_peak <- mat_bon %>% filter(!win %in% c(9933:9936,9939,9940))
map_bon_peak <- map_bon %>% filter(!win %in% c(20883:20885,20887))
cmd_bon_peak <- cmd_bon %>% filter(!win %in% c(9933,9935,9937,9939,9940,20883,20884,20886,20888))

#Filter SNPs for peak
snps_mat_peak <- snps_mat %>% filter(win %in% mat_bon_peak$win)
snps_map_peak <- snps_mat %>% filter(win %in% map_bon_peak$win)
snps_cmd_peak <- snps_mat %>% filter(win %in% cmd_bon_peak$win)





#View empirical p-value plots using plotly

#CMD
wza_empri_cmd <- ggplot(data = wza_win_cmd, aes( x = pos/1e6, y = -log10(approx_p)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(397,398))+
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_cmd)[1])), col = "red", lty = 2, lwd = 1)+
  scale_color_manual(values = rep(c("black", "magenta3"), 22 )) +
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
    axis.title.y = element_text(color="black", size=18,vjust = 1.4, face="bold",hjust=0.5)
  ) 
#wza_empri_cmd
ggplotly(wza_empri_cmd)

#MAP
wza_empri_map <- ggplot(data = wza_win_map, aes( x = pos/1e6, y = -log10(approx_p)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(397,399))+
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_map)[1])), col = "red", lty = 2, lwd = 1)+
  scale_color_manual(values = rep(c("black", "deepskyblue"), 22 )) +
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
    axis.title.y = element_text(color="black", size=18,vjust = 1.4, face="bold",hjust=0.5)
  ) 
#wza_empri_map
ggplotly(wza_empri_map)

#MAT
wza_empri_mat <- ggplot(data = wza_win_mat, aes( x = pos/1e6, y = -log10(approx_p)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(98,100))+
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_mat)[1])), col = "red", lty = 2, lwd = 1)+
  scale_color_manual(values = rep(c("black", "darkgoldenrod"), 22 )) +
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
    axis.title.y = element_text(color="black", size=18,vjust = 1.4, face="bold",hjust=0.5)
  ) 
#wza_empri_mat
ggplotly(wza_empri_mat)


#Add chr-snp column
snps_mat_peak <- snps_mat_peak %>% unite(chr_snp,"chr","snp",sep="_")
snps_map_peak <- snps_map_peak %>% unite(chr_snp,"chr","snp",sep="_")
snps_cmd_peak <- snps_cmd_peak %>% unite(chr_snp,"chr","snp",sep="_")


#Export peak snps
write_csv(snps_mat_peak,"Genomics_scripts/Data/snps_peak_mat.csv")
write_csv(snps_map_peak,"Genomics_scripts/Data/snps_peak_map.csv")
write_csv(snps_cmd_peak,"Genomics_scripts/Data/snps_peak_cmd.csv")

#Export peak windows
write_csv(mat_bon_peak,"Genomics_scripts/Data/peak_window_mat.csv")
write_csv(map_bon_peak,"Genomics_scripts/Data/peak_window_map.csv")
write_csv(cmd_bon_peak,"Genomics_scripts/Data/peak_window_cmd.csv")




