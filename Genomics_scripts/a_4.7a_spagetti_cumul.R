################################################################################################################
## Regression plots of BF>20 SNPs for random and observed
## Author Daniel Anstett
## 
## 
## Last Modified May 19, 2021
################################################################################################################
#Import libraries
library(tidyverse)

################################################################################################################
#Functions

get_cumul <- function(freq_df,cumul_df,pop_num){
  freq_env_pop <- freq_df %>% filter(Site==pop_num)
  snp_list_pop <- cumul_df %>% filter(Site==pop_num)
  freq_mat_pop <- freq_env_pop %>% filter(env=="MAT")
  freq_map_pop <- freq_env_pop %>% filter(env=="MAP")
  freq_cmd_pop <- freq_env_pop %>% filter(env=="CMD")
  snp_list_pop_mat <- snp_list_pop %>% filter(env=="MAT")
  snp_list_pop_map <- snp_list_pop %>% filter(env=="MAP")
  snp_list_pop_cmd <- snp_list_pop %>% filter(env=="CMD")
  freq_mat_cumul <- freq_mat_pop %>% filter(SNP_ID %in% snp_list_pop_mat$snp_ID)
  freq_map_cumul <- freq_map_pop %>% filter(SNP_ID %in% snp_list_pop_map$snp_ID)
  freq_cmd_cumul <- freq_cmd_pop %>% filter(SNP_ID %in% snp_list_pop_cmd$snp_ID)
  freq_env_cumul <- rbind(freq_mat_cumul,freq_map_cumul,freq_cmd_cumul)
  return(freq_env_cumul)
}


#get_cumul <- function(freq_df,cumul_df,pop_num){
#  freq_env_pop <- freq_df %>% filter(Site==pop_num)
#  snp_list_pop <- cumul_df %>% filter(Site==pop_num)
#  freq_mat_cumul <- freq_env_pop %>% filter(SNP_ID %in% snp_list_pop$snp_ID)
#  return(freq_mat_cumul)
#}

#Define Theme

theme_spaghetti <- function(){ 
  
  theme_classic() %+replace%    #replace elements we want to change
    theme(legend.position = "none",
          axis.text.x = element_text(size = 14, face = "bold", angle = 45,hjust = 1, vjust = 1), 
          axis.title = element_text(size =16, face = "bold"), 
          axis.text.y = element_text(size = 14, face = "bold"),legend.title = element_blank(),
          legend.text = element_text(size=12,face="bold"),
          strip.background = element_blank(), 
          strip.text.x=element_text(size=14,face="bold",hjust=0,vjust=-1.2))
}

theme( 
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0, face="bold"),
  axis.title.y = element_text(color="black", size=20, vjust = 1.6,face="bold",hjust=0.5 ,angle=90)
)
################################################################################################################

#Import cumul SNPs
snp_list <- read.csv("Genomics_scripts/Data/snp_list.csv")

#Import timeseries frequencies
freq_mat <- read_csv("Genomics_scripts/Data/freq_MAT_peakbf5.csv")
freq_map <- read_csv("Genomics_scripts/Data/freq_MAP_peakbf5.csv")
freq_cmd <- read_csv("Genomics_scripts/Data/freq_CMD_peakbf5.csv")

#Gather data frames
freq_mat <- freq_mat %>% gather(SNP_ID,SNP_Freq,3:dim(freq_mat)[2]) %>% mutate(env="MAT")
freq_map <- freq_map %>% gather(SNP_ID,SNP_Freq,3:dim(freq_map)[2]) %>% mutate(env="MAP")
freq_cmd <- freq_cmd %>% gather(SNP_ID,SNP_Freq,3:dim(freq_cmd)[2]) %>% mutate(env="CMD")
freq_env <- rbind(freq_mat,freq_map,freq_cmd)

################################################################################################################
#Single Case
freq_env_pop <- freq_env %>% filter(Site==1)
snp_list_pop <- snp_list %>% filter(Site==1)
freq_mat_cumul <- freq_env_pop %>% filter(SNP_ID %in% snp_list_pop$snp_ID)
################################################################################################################

#Generate timeseires frequency data frames filter for cumul SNP list 
#Done for all 12 populations
for(i in 1:12){
  assign(paste("freq_cumul_pop", i, sep="_"), get_cumul(freq_env,snp_list,i))
}

################################################################################################################
#Make graphs
#p1
spag_cumul_1 <- ggplot(data=freq_cumul_pop_1,aes(Year,SNP_Freq,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.25,cex=0.6,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~env) + theme_spaghetti()
ggsave("Graphs_Oct_22/2_01env_both_freqchange_paper_p1.pdf",spag_cumul_1,width=12, height = 6, units = "in")

#p12
spag_cumul_12 <- ggplot(data=freq_cumul_pop_12,aes(Year,SNP_Freq,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.25,cex=0.6,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~env) + theme_spaghetti()
ggsave("Graphs_Oct_22/2_02env_both_freqchange_paper_p12.pdf",spag_cumul_12,width=9, height = 6, units = "in")

#p2
spag_cumul_2 <- ggplot(data=freq_cumul_pop_2,aes(Year,SNP_Freq,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.25,cex=0.6,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~env) + theme_spaghetti()
ggsave("Graphs_Oct_22/2_03env_both_freqchange_paper_p2.pdf",width=12, height = 6, units = "in")

#p3
spag_cumul_3 <- ggplot(data=freq_cumul_pop_3,aes(Year,SNP_Freq,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.25,cex=0.6,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~env) + theme_spaghetti()
ggsave("Graphs_Oct_22/2_04env_both_freqchange_paper_p3.pdf",width=12, height = 6, units = "in")

#p4
spag_cumul_4 <- ggplot(data=freq_cumul_pop_4,aes(Year,SNP_Freq,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.25,cex=0.6,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~env) + theme_spaghetti()
ggsave("Graphs_Oct_22/2_05env_both_freqchange_paper_p4.pdf",width=12, height = 6, units = "in")

#p5
spag_cumul_5 <- ggplot(data=freq_cumul_pop_5,aes(Year,SNP_Freq,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.25,cex=0.6,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~env) + theme_spaghetti()
ggsave("Graphs_Oct_22/2_06env_both_freqchange_paper_p5.pdf",width=12, height = 6, units = "in")

#p6
spag_cumul_6 <- ggplot(data=freq_cumul_pop_6,aes(Year,SNP_Freq,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.25,cex=0.6,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~env) + theme_spaghetti()
ggsave("Graphs_Oct_22/2_07env_both_freqchange_paper_p6.pdf",width=12, height = 6, units = "in")

#p7
spag_cumul_7 <- ggplot(data=freq_cumul_pop_7,aes(Year,SNP_Freq,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.25,cex=0.6,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~env) + theme_spaghetti()
ggsave("Graphs_Oct_22/2_08env_both_freqchange_paper_p7.pdf",width=5, height = 6, units = "in")

#p8
spag_cumul_8 <- ggplot(data=freq_cumul_pop_8,aes(Year,SNP_Freq,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.25,cex=0.6,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~env) + theme_spaghetti()
ggsave("Graphs_Oct_22/2_09env_both_freqchange_paper_p8.pdf",width=9, height = 6, units = "in")

#p9
spag_cumul_9 <- ggplot(data=freq_cumul_pop_9,aes(Year,SNP_Freq,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.25,cex=0.6,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~env) + theme_spaghetti()
ggsave("Graphs_Oct_22/2_10env_both_freqchange_paper_9.pdf",width=9, height = 6, units = "in")

#p10
spag_cumul_10 <- ggplot(data=freq_cumul_pop_10,aes(Year,SNP_Freq,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.25,cex=0.6,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~env) + theme_spaghetti()
ggsave("Graphs_Oct_22/2_11env_both_freqchange_paper_p10.pdf",width=9, height = 6, units = "in")

#p11
spag_cumul_11 <- ggplot(data=freq_cumul_pop_11,aes(Year,SNP_Freq,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.25,cex=0.6,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~env) + theme_spaghetti()
ggsave("Graphs_Oct_22/2_12env_both_freqchange_paper_p11.pdf",width=12, height = 6, units = "in")










