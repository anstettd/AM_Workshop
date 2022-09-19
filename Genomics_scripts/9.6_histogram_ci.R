##################################################################################
## Generate input for strength of selection graphs
## Processing for both observation and permuted CI
## Author Daniel Anstett
## 
## Modified from https://jkzorz.github.io/2020/05/17/Error-bars.html
## Last Modified May 17, 2022
###################################################################################
###################################################################################
#Import libraries
library(tidyverse)

#Import files
mat_obs_ci <- read_csv("mat_obs_ci_peakbf2_region.csv")
map_obs_ci <- read_csv("map_obs_ci_peakbf2_region.csv")
cmd_obs_ci <- read_csv("cmd_obs_ci_peakbf2_region.csv")
env_obs_ci <- read_csv("env_obs_ci_peakbf2_region.csv")


###################################################################################
## Histogram graphs with CI
###################################################################################

# North
mat_p1_hist <- ggplot(env_obs_ci,aes(x=S,y=North,ymin=p1_low,ymax=p1_up))+
  geom_bar(colour = "black", stat = "identity", width = 0.2, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.12) +
  labs(x = "Strength of Selection North", y = "# of SNPs") +
#  scale_y_continuous(limits=c(0,300))+
  theme_classic()
mat_p1_hist <- mat_p1_hist +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
mat_p1_hist <- mat_p1_hist + facet_wrap(.~env)
mat_p1_hist 
ggsave("Graphs_CI_peak/1_1env_hist_ci_peakbf2_North.pdf",width=10, height = 6, units = "in")

# Centre
mat_p2_hist <- ggplot(env_obs_ci,aes(x=S,y=Centre,ymin=p2_low,ymax=p2_up))+
  geom_bar(colour = "black", stat = "identity", width = 0.2, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.12) +
  labs(x = "Strength of Selection Center", y = "# of SNPs") +
  #  scale_y_continuous(limits=c(0,300))+
  theme_classic()
mat_p2_hist <- mat_p2_hist +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
mat_p2_hist <- mat_p2_hist + facet_wrap(.~env)
mat_p2_hist 
ggsave("Graphs_CI_peak/1_2env_hist_ci_peakbf2_Center.pdf",width=10, height = 6, units = "in")

# South
mat_p3_hist <- ggplot(env_obs_ci,aes(x=S,y=South,ymin=p3_low,ymax=p3_up))+
  geom_bar(colour = "black", stat = "identity", width = 0.2, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.12) +
  labs(x = "Strength of Selection South", y = "# of SNPs") +
  #  scale_y_continuous(limits=c(0,300))+
  theme_classic()
mat_p3_hist <- mat_p3_hist +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
mat_p3_hist <- mat_p3_hist + facet_wrap(.~env)
mat_p3_hist 
ggsave("Graphs_CI_peak/1_3env_hist_ci_peakbf2_South.pdf",width=10, height = 6, units = "in")




