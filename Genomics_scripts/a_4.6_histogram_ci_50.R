##################################################################################
## Generate input for strength of selection graphs
## Processing for both observation and permuted CI
## Author Daniel Anstett
## 
## Modified from https://jkzorz.github.io/2020/05/17/Error-bars.html
## Last Modified May 17, 2022
###################################################################################
#Function
#hist_ci <- function(df,site,low_CI,high_CI){
#  emv_hist <- ggplot(df,aes(x=S,y=site,ymin=low_CI,ymax=p2_high_CI))+
#    geom_bar(colour = "black", stat = "identity", width = 0.2, fill = "lightblue1")+
#    geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.12) +
#    labs(x = "Strength of Selection", y = "Count") +
#    theme_classic()
#  emv_hist <- emv_hist +
#    theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
#          axis.title = element_text(size = 20, face = "bold"), 
#          axis.text.y = element_text(size = 16, face = "bold"))
#  #  scale_y_continuous(expand = c(0,0))
#  return(emv_hist)
#}

###################################################################################
#Import libraries
library(tidyverse)

#Import files
mat_obs_ci <- read_csv("Genomics_scripts/Data/mat_obs_ci_peakbf5_slope50.csv")
map_obs_ci <- read_csv("Genomics_scripts/Data/map_obs_ci_peakbf5_slope50.csv")
cmd_obs_ci <- read_csv("Genomics_scripts/Data/cmd_obs_ci_peakbf5_slope50.csv")
env_obs_ci <- read_csv("Genomics_scripts/Data/env_obs_ci_peakbf5_slope50.csv")


###################################################################################
## Histogram graphs with CI
###################################################################################

# Site 1
mat_p1_hist <- ggplot(env_obs_ci,aes(x=S,y=p1,ymin=p1_low,ymax=p1_up))+
  geom_bar(colour = "black", stat = "identity", width = 10, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 6) +
  labs(x = "Strength of Selection (P1)", y = "# of SNPs") +
  scale_x_continuous(breaks=seq(-50,50,25))+
  scale_y_continuous(limits=c(0,130),breaks=seq(0,125,25))+
  theme_classic()
mat_p1_hist <- mat_p1_hist +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
mat_p1_hist <- mat_p1_hist + facet_wrap(.~env)
mat_p1_hist 
ggsave("Graphs_Oct_22/3_01env_hist_ci_p1.pdf",width=10, height = 6, units = "in")

# Site 12
mat_p12_hist <- ggplot(env_obs_ci,aes(x=S,y=p12,ymin=p12_low,ymax=p12_up))+
  geom_bar(colour = "black", stat = "identity", width = 10, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 6) +
  labs(x = "Strength of Selection (P12)", y = "# of SNPs") +
 scale_x_continuous(breaks=seq(-50,50,25))+
  scale_y_continuous(limits=c(0,130),breaks=seq(0,125,25))+
  theme_classic()
mat_p12_hist <- mat_p12_hist +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
mat_p12_hist <- mat_p12_hist + facet_wrap(.~env)
mat_p12_hist 
ggsave("Graphs_Oct_22/3_02env_hist_ci_p12.pdf",width=10, height = 6, units = "in")

# Site 2
mat_p2_hist <- ggplot(env_obs_ci,aes(x=S,y=p2,ymin=p2_low,ymax=p2_up))+
  geom_bar(colour = "black", stat = "identity", width = 10, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 6) +
  labs(x = "Strength of Selection (P2)", y = "# of SNPs") +
 scale_x_continuous(breaks=seq(-50,50,25))+
  scale_y_continuous(limits=c(0,130),breaks=seq(0,125,25))+
  theme_classic()
mat_p2_hist <- mat_p2_hist +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
mat_p2_hist <- mat_p2_hist + facet_wrap(.~env)
mat_p2_hist 
ggsave("Graphs_Oct_22/3_03env_hist_ci_p2.pdf",width=10, height = 6, units = "in")

# Site 3
mat_p3_hist <- ggplot(env_obs_ci,aes(x=S,y=p3,ymin=p3_low,ymax=p3_up))+
  geom_bar(colour = "black", stat = "identity", width = 10, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 6) +
  labs(x = "Strength of Selection (P3)", y = "# of SNPs") +
 scale_x_continuous(breaks=seq(-50,50,25))+
  scale_y_continuous(limits=c(0,130),breaks=seq(0,125,25))+
  theme_classic()
mat_p3_hist <- mat_p3_hist +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
mat_p3_hist <- mat_p3_hist + facet_wrap(.~env)
mat_p3_hist 
ggsave("Graphs_Oct_22/3_04env_hist_ci_p3.pdf",width=10, height = 6, units = "in")

# Site 4
mat_p4_hist <- ggplot(env_obs_ci,aes(x=S,y=p4,ymin=p4_low,ymax=p4_up))+
  geom_bar(colour = "black", stat = "identity", width = 10, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 6) +
  labs(x = "Strength of Selection (P4)", y = "# of SNPs") +
 scale_x_continuous(breaks=seq(-50,50,25))+
  scale_y_continuous(limits=c(0,130),breaks=seq(0,125,25))+
  theme_classic()
mat_p4_hist <- mat_p4_hist +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
mat_p4_hist <- mat_p4_hist + facet_wrap(.~env)
mat_p4_hist 
ggsave("Graphs_Oct_22/3_05env_hist_ci_p4.pdf",width=10, height = 6, units = "in")

# Site 5
mat_p5_hist <- ggplot(env_obs_ci,aes(x=S,y=p5,ymin=p5_low,ymax=p5_up))+
  geom_bar(colour = "black", stat = "identity", width = 10, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 6) +
  labs(x = "Strength of Selection (P5)", y = "# of SNPs") +
 scale_x_continuous(breaks=seq(-50,50,25))+
  scale_y_continuous(limits=c(0,130))+
  theme_classic()
mat_p5_hist <- mat_p5_hist +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
mat_p5_hist <- mat_p5_hist + facet_wrap(.~env)
mat_p5_hist 
ggsave("Graphs_Oct_22/3_06env_hist_ci_p5.pdf",width=10, height = 6, units = "in")

# Site 6
mat_p6_hist <- ggplot(env_obs_ci,aes(x=S,y=p6,ymin=p6_low,ymax=p6_up))+
  geom_bar(colour = "black", stat = "identity", width = 10, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 6) +
  labs(x = "Strength of Selection (P6)", y = "# of SNPs") +
  scale_x_continuous(breaks=seq(-50,50,25))+
  scale_y_continuous(limits=c(0,130),breaks=seq(0,125,25))+
  theme_classic()
mat_p6_hist <- mat_p6_hist +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
mat_p6_hist <- mat_p6_hist + facet_wrap(.~env)
mat_p6_hist 
ggsave("Graphs_Oct_22/3_07env_hist_ci_p6.pdf",width=10, height = 6, units = "in")

# Site 7
mat_p7_hist <- ggplot(env_obs_ci,aes(x=S,y=p7,ymin=p7_low,ymax=p7_up))+
  geom_bar(colour = "black", stat = "identity", width = 10, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 6) +
  labs(x = "Strength of Selection (P7)", y = "# of SNPs") +
 scale_x_continuous(breaks=seq(-50,50,25))+
  scale_y_continuous(limits=c(0,130),breaks=seq(0,125,25))+
  theme_classic()
mat_p7_hist <- mat_p7_hist +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
mat_p7_hist <- mat_p7_hist + facet_wrap(.~env)
mat_p7_hist 
ggsave("Graphs_Oct_22/3_08env_hist_ci_p7.pdf",width=10, height = 6, units = "in")

# Site 8
mat_p8_hist <- ggplot(env_obs_ci,aes(x=S,y=p8,ymin=p8_low,ymax=p8_up))+
  geom_bar(colour = "black", stat = "identity", width = 10, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 6) +
  labs(x = "Strength of Selection (P8)", y = "# of SNPs") +
 scale_x_continuous(breaks=seq(-50,50,25))+
  scale_y_continuous(limits=c(0,130),breaks=seq(0,125,25))+
  theme_classic()
mat_p8_hist <- mat_p8_hist +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
mat_p8_hist <- mat_p8_hist + facet_wrap(.~env)
mat_p8_hist 
ggsave("Graphs_Oct_22/3_09env_hist_ci_p8.pdf",width=10, height = 6, units = "in")

# Site 9
mat_p9_hist <- ggplot(env_obs_ci,aes(x=S,y=p9,ymin=p9_low,ymax=p9_up))+
  geom_bar(colour = "black", stat = "identity", width = 10, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 6) +
  labs(x = "Strength of Selection (P9)", y = "# of SNPs") +
 scale_x_continuous(breaks=seq(-50,50,25))+
  scale_y_continuous(limits=c(0,130),breaks=seq(0,125,25))+
  theme_classic()
mat_p9_hist <- mat_p9_hist +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
mat_p9_hist <- mat_p9_hist + facet_wrap(.~env)
mat_p9_hist 
ggsave("Graphs_Oct_22/3_10env_hist_ci_p9.pdf",width=10, height = 6, units = "in")

# Site 10
mat_p10_hist <- ggplot(env_obs_ci,aes(x=S,y=p10,ymin=p10_low,ymax=p10_up))+
  geom_bar(colour = "black", stat = "identity", width = 10, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 6) +
  labs(x = "Strength of Selection (P10)", y = "# of SNPs") +
 scale_x_continuous(breaks=seq(-50,50,25))+
  scale_y_continuous(limits=c(0,130),breaks=seq(0,125,25))+
  theme_classic()
mat_p10_hist <- mat_p10_hist +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
mat_p10_hist <- mat_p10_hist + facet_wrap(.~env)
mat_p10_hist 
ggsave("Graphs_Oct_22/3_11env_hist_ci_p10.pdf",width=10, height = 6, units = "in")

# Site 11
mat_p11_hist <- ggplot(env_obs_ci,aes(x=S,y=p11,ymin=p11_low,ymax=p11_up))+
  geom_bar(colour = "black", stat = "identity", width = 10, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 6) +
  labs(x = "Strength of Selection (P11)", y = "# of SNPs") +
 scale_x_continuous(breaks=seq(-50,50,25))+
  scale_y_continuous(limits=c(0,130),breaks=seq(0,125,25))+
  theme_classic()
mat_p11_hist <- mat_p11_hist +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
mat_p11_hist <- mat_p11_hist + facet_wrap(.~env)
mat_p11_hist 
ggsave("Graphs_Oct_22/3_12env_hist_ci_p11.pdf",width=10, height = 6, units = "in")




