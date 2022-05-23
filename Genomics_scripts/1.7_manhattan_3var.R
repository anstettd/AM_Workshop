##################################################################################
## Filter out strong evidence for snps associated with climate
## Author Daniel Anstett & Julia Anstett
## 
##
## Last Modified July 6, 2021
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)
library(qqman)

#Import chromosome size
chr_size <- read_csv("Genomics_scripts/Data/chr_size.csv")
chr_size[,3] <- cumsum(chr_size$size) #get cumulative chromosome position
colnames(chr_size)[3] <- "poz"


#Import snp env associations
env1 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_1_trim.tsv",header=F, sep=" ")
env2 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_2_trim.tsv",header=F, sep=" ")
env5 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_5_trim.tsv",header=F, sep=" ")


#Name Columns
colnames(env1) <- c("Chromosome","SNP","Env","BF")
colnames(env2) <- c("Chromosome","SNP","Env","BF")
colnames(env3) <- c("Chromosome","SNP","Env","BF")
colnames(env4) <- c("Chromosome","SNP","Env","BF")
colnames(env5) <- c("Chromosome","SNP","Env","BF")
colnames(env6) <- c("Chromosome","SNP","Env","BF")
colnames(env7) <- c("Chromosome","SNP","Env","BF")
colnames(env8) <- c("Chromosome","SNP","Env","BF")
colnames(env9) <- c("Chromosome","SNP","Env","BF")

#Make Chromosome SNP variable for easier left_joining later
#Make a cumulative basepair value
env1_united <- env1 %>% unite(chr_snp,Chromosome,SNP)
env1_united <- env1_united %>% select(chr_snp)
env1_united <- cbind(env1_united,env1)
env1_united <- env1_united%>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env2_united <- env2 %>% unite(chr_snp,Chromosome,SNP)
env2_united <- env2_united %>% select(chr_snp)
env2_united <- cbind(env2_united,env2)
env2_united <- env2_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env5_united <- env5 %>% unite(chr_snp,Chromosome,SNP)
env5_united <- env5_united %>% select(chr_snp)
env5_united <- cbind(env5_united,env5)
env5_united <- env5_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))


###
#Make Manhattan Plots

#ENV 1 - MAT 
axisdf_1 <- env1_united %>% group_by(CHR) %>% summarize(center=( max(BP) + min(BP) ) / 2 )

ggplot(env1_united, aes(x=BP, y=BF)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=0.5) +
  scale_color_manual(values = rep(c("black", "darkgoldenrod"), 22 )) +
  geom_hline(yintercept=20, linetype="dashed",color = "red", size=2) +
  # custom X axis:
  scale_x_continuous(expand = c(0, 0), label = axisdf_1$CHR, breaks= axisdf_1$center) +
  scale_y_continuous(breaks = c(-20,-10,0,10,20,30,40)) +
  theme_classic() +
  labs(
    y = "Log Bayes Factor",
    x = "Position")+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=25, face="bold"),
    axis.text.y = element_text(size=25,face="bold"),
    axis.title.x = element_text(color="black", size=0, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=0,vjust = 2, face="bold",hjust=0.5))

#Export. File is very large and somewhat hard to hande/open. 
ggsave("/Users/daniel_anstett/Dropbox/AM_Workshop/large_graphs/graph_manhattan_1.pdf",width=9, height=6)



#ENV 2 - MAP 
axisdf_2 <- env2_united %>% group_by(CHR) %>% summarize(center=( max(BP) + min(BP) ) / 2 )

ggplot(env2_united, aes(x=BP, y=BF)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=0.5) +
  scale_color_manual(values = rep(c("black", "deepskyblue"), 22 )) +
  geom_hline(yintercept=20, linetype="dashed",color = "red", size=2) +
  # custom X axis:
  scale_x_continuous(expand = c(0, 0), label = axisdf_2$CHR, breaks= axisdf_2$center) +
  scale_y_continuous(breaks = c(-20,-10,0,10,20,30,40)) +
  theme_classic() +
  labs(
    y = "Bayes Factor",
    x = "Position")+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=25, face="bold"),
    axis.text.y = element_text(size=25,face="bold"),
    axis.title.x = element_text(color="black", size=0, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=0,vjust = 2, face="bold",hjust=0.5)
  )

#Export. File is very large and somewhat hard to hande/open. 
ggsave("/Users/daniel_anstett/Dropbox/AM_Workshop/large_graphs/graph_manhattan_2.pdf",width=9, height=6)




#ENV 5 - CMD 
axisdf_5 <- env5_united %>% group_by(CHR) %>% summarize(center=( max(BP) + min(BP) ) / 2 )

ggplot(env5_united, aes(x=BP, y=BF)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=0.5) +
  scale_color_manual(values = rep(c("black", "magenta3"), 22 )) +
  geom_hline(yintercept=20, linetype="dashed",color = "blue", size=2) +
  # custom X axis:
  scale_x_continuous(expand = c(0, 0), label = axisdf_5$CHR, breaks= axisdf_5$center) +
  scale_y_continuous(breaks = c(-20,-10,0,10,20,30,40)) +
  theme_classic() +
  labs(
    y = "Bayes Factor",
    x = "Position")+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=25, face="bold"),
    axis.text.y = element_text(size=25,face="bold"),
    axis.title.x = element_text(color="black", size=0, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=0,vjust = 2, face="bold",hjust=0.5)
  )

#Export. File is very large and somewhat hard to hande/open. 
ggsave("/Users/daniel_anstett/Dropbox/AM_Workshop/large_graphs/graph_manhattan_5.pdf",width=9, height=6)




################
#Test plot with only two chromosomes
env1_united_chr1_2 <- env1_united %>% filter(CHR<3)

axisdf_1 <- env1_united_chr1_2 %>% group_by(CHR) %>% summarize(center=( max(BP) + min(BP) ) / 2 )

ggplot(env1_united_chr1_2, aes(x=BP, y=BF)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=0.5) +
  scale_color_manual(values = rep(c("black", "magenta2"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "deepskyblue1", size=0.9) +
  # custom X axis:
  scale_x_continuous(expand = c(0, 0), label = axisdf_1$CHR, breaks= axisdf_1$center) +
  scale_y_continuous(breaks = c(-20,-10,0,10,20,30,40)) +
  theme_classic() +
    labs(
      y = "Bayes Factor",
      x = "Position")+
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
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5)
  )

  












