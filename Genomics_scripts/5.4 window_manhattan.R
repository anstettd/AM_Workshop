#############################################################################################################
## Plot Manhattan for WZA windows
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

wza_win_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_win_mat.csv")
wza_win_map <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_win_map.csv")
wza_win_cmd <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_win_cmd.csv")


###########################################################################################################


#Plot Distribution
ggplot(data = wza_win_mat, aes( x = WZA))+
  geom_histogram( bins = 50)+
  scale_y_continuous("Count")+
  theme_classic() #Approx normal

ggplot(data = wza_win_map, aes( x = WZA))+
  geom_histogram( bins = 50)+
  scale_y_continuous("Count")+
  theme_classic() #Approx normal

ggplot(data = wza_win_cmd, aes( x = WZA))+
  geom_histogram( bins = 50)+
  scale_y_continuous("Count")+
  theme_classic() #Approx normal


# Plot WZA score per window

#MAT
wza_mat_man <- ggplot(data = wza_win_mat, aes( x = pos/1e6, y = WZA))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("WZA Score", limits=c(-16,20))+
  scale_x_continuous("Position (Mbp)")+
  geom_hline(aes(yintercept = quantile(WZA, 0.99)), col = "red", lty = 2, lwd = 1)+
  geom_hline(aes(yintercept = quantile(WZA, 0.95)), col = "orange", lty = 2, lwd = 1)+
  geom_hline(aes(yintercept = quantile(WZA, 0.90)), col = "skyblue", lty = 2, lwd = 1)+
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
    axis.title.y = element_text(color="black", size=20,vjust = 0.8, face="bold",hjust=0.5)
  )
#wza_mat_man 
ggsave("/Users/daniel_anstett/Dropbox/a_Papers/Genomics_paper/Graphs/WZA_Graphs/wza_mat_man.png",
       wza_mat_man, width=10, height = 5, units = "in")

#MAP
wza_map_man <- ggplot(data = wza_win_map, aes( x = pos/1e6, y = WZA))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("WZA Score", limits=c(-16,20))+
  scale_x_continuous("Position (Mbp)")+
  geom_hline(aes(yintercept = quantile(WZA, 0.99)), col = "red", lty = 2, lwd = 1)+
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
    axis.title.y = element_text(color="black", size=20,vjust = 0.8, face="bold",hjust=0.5)
  )
#wza_map_man 
ggsave("/Users/daniel_anstett/Dropbox/a_Papers/Genomics_paper/Graphs/WZA_Graphs/wza_map_man.png",
       wza_map_man, width=10, height = 5, units = "in")


#CMD
wza_cmd_man <- ggplot(data = wza_win_cmd, aes( x = pos/1e6, y = WZA))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("WZA Score", limits=c(-16,20))+
  scale_x_continuous("Position (Mbp)")+
  geom_hline(aes(yintercept = quantile(WZA, 0.99)), col = "red", lty = 2, lwd = 1)+
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
    axis.title.y = element_text(color="black", size=20,vjust = 0.8, face="bold",hjust=0.5)
  )
#wza_cmd_man 
ggsave("/Users/daniel_anstett/Dropbox/a_Papers/Genomics_paper/Graphs/WZA_Graphs/wza_cmd_man.png",
       wza_cmd_man, width=10, height = 5, units = "in")



#Get 1% of windows as a list
#hit_set <- row.names(wza_win_mat[wza_win_mat$WZA > quantile(wza_win_mat$WZA, 0.99),])
#hit_set

#######################################################################################################


# WZA Empirical p-value

#MAT
wza_empri_mat <- ggplot(data = wza_win_mat, aes( x = pos/1e6, y = -log10(approx_p)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)")+
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
wza_empri_mat
ggsave("/Users/daniel_anstett/Dropbox/a_Papers/Genomics_paper/Graphs/WZA_Graphs/wza_empri_mat.png",
       wza_empri_mat, width=10, height = 5, units = "in")


#MAP
wza_empri_map <- ggplot(data = wza_win_map, aes( x = pos/1e6, y = -log10(approx_p)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)")+
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
wza_empri_map
ggsave("/Users/daniel_anstett/Dropbox/a_Papers/Genomics_paper/Graphs/WZA_Graphs/wza_empri_map.png",
       wza_empri_map, width=10, height = 5, units = "in")


#CMD
wza_empri_cmd <- ggplot(data = wza_win_cmd, aes( x = pos/1e6, y = -log10(approx_p)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)")+
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
wza_empri_cmd
ggsave("/Users/daniel_anstett/Dropbox/a_Papers/Genomics_paper/Graphs/WZA_Graphs/wza_empri_cmd.png",
       wza_empri_cmd, width=10, height = 5, units = "in")





# FDR Correction

#MAT
wza_empri_fdr_mat <- ggplot(data = wza_win_mat, aes( x = pos/1e6, y = -log10(p.adjust(approx_p,method = "fdr"))))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(WZA q-value)", limits = c(0,15))+
  scale_x_continuous("Position (Mbp)")+
  geom_hline(aes(yintercept = -log10(0.05)), col = "red", lty = 2, lwd = 1)+
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
    axis.title.y = element_text(color="black", size=20,vjust = 0.8, face="bold",hjust=0.5)
  )
wza_empri_fdr_mat
ggsave("/Users/daniel_anstett/Dropbox/a_Papers/Genomics_paper/Graphs/WZA_Graphs/wza_empri_fdr_mat.png",
       wza_empri_fdr_mat, width=10, height = 5, units = "in")


#MAP
wza_empri_fdr_map <- ggplot(data = wza_win_map, aes( x = pos/1e6, y = -log10(p.adjust(approx_p,method = "fdr"))))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(WZA q-value)", limits = c(0,15))+
  scale_x_continuous("Position (Mbp)")+
  geom_hline(aes(yintercept = -log10(0.05)), col = "red", lty = 2, lwd = 1)+
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
    axis.title.y = element_text(color="black", size=20,vjust = 0.8, face="bold",hjust=0.5)
  )
wza_empri_fdr_map
ggsave("/Users/daniel_anstett/Dropbox/a_Papers/Genomics_paper/Graphs/WZA_Graphs/wza_empri_fdr_map.png",
       wza_empri_fdr_map, width=10, height = 5, units = "in")


#CMD
wza_empri_fdr_cmd <- ggplot(data = wza_win_cmd, aes( x = pos/1e6, y = -log10(p.adjust(approx_p,method = "fdr"))))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(WZA q-value)", limits = c(0,15))+
  scale_x_continuous("Position (Mbp)")+
  geom_hline(aes(yintercept = -log10(0.05)), col = "red", lty = 2, lwd = 1)+
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
    axis.title.y = element_text(color="black", size=20,vjust = 0.8, face="bold",hjust=0.5)
  )
wza_empri_fdr_cmd
ggsave("/Users/daniel_anstett/Dropbox/a_Papers/Genomics_paper/Graphs/WZA_Graphs/wza_empri_fdr_cmd.png",
       wza_empri_fdr_cmd, width=10, height = 5, units = "in")

