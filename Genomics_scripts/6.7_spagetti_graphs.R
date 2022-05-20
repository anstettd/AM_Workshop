##################################################################################
## Regression plots of BF>20 SNPs for random and observed
## Author Daniel Anstett
## 
## 
## Last Modified May 19, 2021
###################################################################################
#Import libraries
library(tidyverse)

#Import timeseries frequencies
freq_MAT <- read_csv("Genomics_scripts/Data/freq_MAT.csv")
freq_MAP <- read_csv("Genomics_scripts/Data/freq_MAP.csv")
freq_CMD <- read_csv("Genomics_scripts/Data/freq_CMD.csv")

#Gather data frames
freq_MAT <- freq_MAT %>% gather(SNP_ID,SNP_Freq,3:275)
freq_MAP <- freq_MAP %>% gather(SNP_ID,SNP_Freq,3:571)
freq_CMD <- freq_CMD %>% gather(SNP_ID,SNP_Freq,3:304)

#Import gathered rand (seed=1) data frame from 6.3
rand_MAT <- read_csv("Genomics_scripts/Data/rand_gathered_MAT.csv")
rand_MAP <- read_csv("Genomics_scripts/Data/rand_gathered_MAP.csv")
rand_CMD <- read_csv("Genomics_scripts/Data/rand_gathered_CMD.csv")


#Add Type
freq_MAT$Type <- "Climate Associated"
freq_MAP$Type <- "Climate Associated"
freq_CMD$Type <- "Climate Associated"

rand_MAT$Type <- "Random"
rand_MAP$Type <- "Random"
rand_CMD$Type <- "Random"

#Merge Climate Associated and Random
joint_MAT <- rbind(freq_MAT,rand_MAT)
joint_MAP <- rbind(freq_MAP,rand_MAP)
joint_CMD <- rbind(freq_CMD,rand_CMD)


###################################################################################
###Test plot
##MAT

#Site 1
freq_MAT_p1 <- freq_MAT %>% filter(Site==4)
rand_MAT_p1 <- rand_MAT %>% filter(Site==4)
joint_MAT_p1 <- joint_MAT %>% filter(Site==4)



ggplot(data=rand_MAT_p1,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F, alpha=.4,cex=0.6) + theme_classic() +
  labs(y="SNP Frequency",x="Year") + scale_color_manual(values=c("blue")) + 
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))

ggplot(data=freq_MAT_p1,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F, alpha=.4,cex=0.6) + theme_classic() +
  labs(y="SNP Frequency",x="Year") + scale_color_manual(values=c("red")) + 
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))

ggplot(data=joint_MAT_p1,aes(Year,SNP_Freq,group=SNP_ID,color=Type)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F, alpha=.2,cex=0.6) + theme_classic() +
  labs(y="SNP Frequency",x="Year") + scale_color_manual(values=c("red", "grey")) + theme(
    axis.text.x = element_text(size = 16, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
    axis.title = element_text(size = 20, face = "bold"), 
    axis.text.y = element_text(size = 16, face = "bold"))












ggplot(data=long_012_merged,aes(year,geno,group=environment,color=quantile)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F, alpha=.7,cex=1.5) +
  theme_bw() +
  labs(y="Frequency",x="Year")