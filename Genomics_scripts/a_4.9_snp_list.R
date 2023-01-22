##################################################################################
## Arrange SNP list per population
## 
## Author Daniel Anstett
## 
## 
## Last Modified May 17, 2022
##################################################################################################################
#Import libraries
library(tidyverse)

##################################################################################################################
#Import files
snp_list_unique <- read.csv("Genomics_scripts/Data/snp_list_unique.csv")
snp_list <- read.csv("Genomics_scripts/Data/snp_list.csv")
snp_list_mat <- snp_list  %>% filter(env=="MAT")
snp_list_map <- snp_list  %>% filter(env=="MAP")
snp_list_cmd <- snp_list  %>% filter(env=="CMD")

snp_list_mat_pop <- snp_list_mat %>% group_by(snp_ID) %>% summarise(SNP_count = n()) %>% mutate(env="MAT")
snp_list_map_pop <- snp_list_map %>% group_by(snp_ID) %>% summarise(SNP_count = n()) %>% mutate(env="MAP")
snp_list_cmd_pop <- snp_list_cmd %>% group_by(snp_ID) %>% summarise(SNP_count = n()) %>% mutate(env="CMD")
snp_list_pop <- snp_list_unique %>% group_by(snp_ID) %>% summarise(SNP_count = n())

snp_list_env_pop <- rbind(snp_list_mat_pop,snp_list_map_pop,snp_list_cmd_pop)


## Filter SNPS represented at least 4 times
snp_list_pop_4 <- snp_list_pop %>% filter(SNP_count>3)
snp_list_unique_4 <- snp_list_unique %>% filter(snp_ID %in% as.character(snp_list_pop_4$snp_ID))

#unique_4_spread <- snp_list_unique_4 %>% select(Site,snp_ID) %>% spread(snp_ID,Site)

#snp_pop_overlap <- data.frame()
#snp_rec <- unique(unique_4_spread$snp_ID)
#for(i in 1:length(snp_rec)) {
#  snp_input <- unique_4_spread %>% filter(snp_ID==as.character(snp_rec[i]))
#  snp_df <- data.frame(t(snp_input$Site))
#  snp_pop_overlap[i,1] <-  as.character(snp_rec[i])
#  snp_pop_overlap[i,2:dim(snp_df)[2]+1] <- snp_df
#}


  
  


##################################################################################################################
# Plot histograms of # of times SNP is in each population

mat_hist <- ggplot(snp_list_pop,aes(x=SNP_count))+
  geom_histogram(binwidth = 0.5)+
  labs(x = "Number of Populations", y = "Number of SNPs") +
  scale_x_continuous(n.breaks=6)+
  theme_classic() +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
#,
#strip.background = element_blank(),
#strip.text.x = element_blank())
#mat_p1_hist <- mat_p1_hist + facet_wrap(.~env)
mat_hist 
ggsave("Graphs_Oct_22/parallel_snp_unique.pdf",width=10, height = 6, units = "in")



mat_p1_hist <- ggplot(snp_list_env_pop,aes(x=SNP_count))+
  geom_histogram(binwidth = 0.5)+
  labs(x = "Number of Populations", y = "Number of SNPs") +
  scale_x_continuous(n.breaks=6)+
  theme_classic()
mat_p1_hist <- mat_p1_hist +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
        #,
        #strip.background = element_blank(),
        #strip.text.x = element_blank())
mat_p1_hist <- mat_p1_hist + facet_wrap(.~env)
mat_p1_hist 
ggsave("Graphs_Oct_22/parallel_snp.pdf",width=10, height = 6, units = "in")




##################################################################################################################
# Plot # of times population is repeated 4 or more times

pop_rep4 <- ggplot(snp_list_unique_4,aes(x=Site))+
  geom_histogram(binwidth = 0.5)+
  labs(x = "Populations", y = "SNPs Repeated 4 or More Times") +
  scale_x_continuous(n.breaks=12)+
  theme_classic() +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
#,
#strip.background = element_blank(),
#strip.text.x = element_blank())
#mat_p1_hist <- mat_p1_hist + facet_wrap(.~env)
pop_rep4
ggsave("Graphs_Oct_22/pop_repeated.pdf",width=10, height = 6, units = "in")

