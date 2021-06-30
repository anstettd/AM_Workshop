##########################################################################################################
## Plotting of hard filtering paramters for goldset
## Author Daniel Anstett
## 
##
## Last Modified May 17, 2021
##########################################################################################################


##########################################################################################################
#Import libraries
library(tidyverse)
library(cowplot)

##########################################################################################################
#Import Data
indel <- read.table("Genomics_scripts/Data/indel_gold_var.tsv", sep = "\t", fill = TRUE , header = TRUE)
snp <- read.table("Genomics_scripts/Data/snp_gold_var.tsv", sep = "\t", fill = TRUE , header = TRUE)
dp_values_indel <- read.csv("Genomics_scripts/Data/indel_cov_summary.csv", header=T)
dp_values_indel <- dp_values_indel %>% mutate(lower=d1-d2,upper=d1+d2)
dp_values_snp <- read.csv("Genomics_scripts/Data/snp_cov_summary.csv", header=T)
dp_values_snp <- dp_values_snp %>% mutate(lower=d1-d2,upper=d1+d2)




##########################################################################################################
#Rearrange data

#AF Rearrangement
indel_af <- indel %>% select(CHROM,POS,TYPE,AF)
indel_af <- indel_af %>% filter(!grepl(",",AF))

snp_af <- snp %>% select(CHROM,POS,TYPE,AF)
snp_af <- snp_af %>% filter(!grepl(",",AF))

#DP per chromosome
#indel
indel_dp <- indel %>% select(CHROM,DP)
indel_dp_chr1 <-indel_dp %>% filter(CHROM == "CE10_chr1")
indel_dp_chr2 <-indel_dp %>% filter(CHROM == "CE10_chr2")
indel_dp_chr3 <-indel_dp %>% filter(CHROM == "CE10_chr3")
indel_dp_chr4 <-indel_dp %>% filter(CHROM == "CE10_chr4")
indel_dp_chr5 <-indel_dp %>% filter(CHROM == "CE10_chr5")
indel_dp_chr6 <-indel_dp %>% filter(CHROM == "CE10_chr6")
indel_dp_chr7 <-indel_dp %>% filter(CHROM == "CE10_chr7")
indel_dp_chr8 <-indel_dp %>% filter(CHROM == "CE10_chr8")

#snp
snp_dp <- snp %>% select(CHROM,DP)
snp_dp_chr1 <-snp_dp %>% filter(CHROM == "CE10_chr1")
snp_dp_chr2 <-snp_dp %>% filter(CHROM == "CE10_chr2")
snp_dp_chr3 <-snp_dp %>% filter(CHROM == "CE10_chr3")
snp_dp_chr4 <-snp_dp %>% filter(CHROM == "CE10_chr4")
snp_dp_chr5 <-snp_dp %>% filter(CHROM == "CE10_chr5")
snp_dp_chr6 <-snp_dp %>% filter(CHROM == "CE10_chr6")
snp_dp_chr7 <-snp_dp %>% filter(CHROM == "CE10_chr7")
snp_dp_chr8 <-snp_dp %>% filter(CHROM == "CE10_chr8")

##########################################################################################################
#Make Functions

## Establish function to make histogram plot with one vertical line
hist_line <- function(dat,variable,hard,title,name) {
  hist.hard <- ggplot(dat, aes(X=variable))+
    geom_histogram(aes(variable))+
    geom_vline(xintercept=c(hard),color="#0000CC")+
    scale_y_continuous(name="Count")+
    scale_x_continuous(name=name)+
    theme_classic()
  hist.hard <- hist.hard  + theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))+
    ggtitle(title)
  return(hist.hard)
}

## Establish function to make histogram plot with two vertical lines
hist_line2 <- function(dat,variable,hard1,hard2,title,name) {
  hist.hard <- ggplot(dat, aes(X=variable))+
    geom_histogram(aes(variable))+
    geom_vline(xintercept=c(hard1),color="#0000CC")+
    geom_vline(xintercept=c(hard2),color="#0000CC")+
    scale_y_continuous(name="Count")+
    scale_x_continuous(name=name)+
    theme_classic()
  hist.hard <- hist.hard  + theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))+
    ggtitle(title)
  return(hist.hard)
}


## Establish function to make histogram plot with two vertical lines, range
hist_line2_range <- function(dat,variable,hard1,hard2,title,name) {
  hist.hard <- ggplot(dat, aes(X=variable))+
    geom_histogram(aes(variable))+
    geom_vline(xintercept=c(hard1),color="#0000CC")+
    geom_vline(xintercept=c(hard2),color="#0000CC")+
    scale_y_continuous(name="Count")+
    scale_x_continuous(name=name,limits=c(-50,500))+
    theme_classic()
  hist.hard <- hist.hard  + theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))+
    ggtitle(title)
  return(hist.hard)
}


##########################################################################################################
##########################################################################################################
#Based on MQ > 50.0 && AN > 24*0.9 && SOR > -1.0 && SOR < 1.0 && AF > 0.25 && ExcessHet < 10 
# && ExcessHet>-4.5 && BaseQRankSum > -1.0 && BaseQRankSum < 1.0 && DP

#Make histogram plots of each value for
#MQ
MQ_1 <- hist_line(indel,indel$MQ,50,"(indel) keep MQ > 50.0","MQ")
MQ_1
MQ_2 <- hist_line(snp,snp$MQ,50,"(SNP) keep MQ > 50.0","MQ")
MQ_2

#AN
AN_1 <- hist_line(indel,indel$AN,16*0.9,"(indel) keep AN > 14.4","AN") #AN
AN_1
AN_2 <- hist_line(snp,snp$AN,16*0.9,"(SNP) keep AN > 14.4","AN") #AN
AN_2

#AF
AF_1 <- hist_line(indel_af,as.numeric(as.character(indel_af$AF)),0.25,"(indel) keep AF>0.25", "AF")
AF_1
AF_2 <- hist_line(snp_af,as.numeric(as.character(snp_af$AF)),0.25, "(snp) keep AF>0.25", "AF")
AF_2

#SOR
SOR_1 <- hist_line2(indel,indel$SOR,-1,1,"(indel) keep SOR > -1.0, SOR <1.0","SOR") #SOR
SOR_1
SOR_2 <- hist_line2(snp,snp$SOR,-1,1,"(SNP) keep SOR > -1.0, SOR <1.0","SOR") #SOR
SOR_2

#ExcessHet
EH_1 <- hist_line2(indel,indel$ExcessHet,-4.5,10,"(indel) keep >-4.5, <10","ExcessHet") #ExcessHet
EH_1
EH_2 <- hist_line2(snp,snp$ExcessHet,-4.5,10,"(SNP) keep >-4.5, <10","ExcessHet") #ExcessHet
EH_2




#DP
#indel
dp_indel_chr1 <- hist_line2_range(indel_dp_chr1,indel_dp_chr1$DP,dp_values_indel$lower[1], 
                  dp_values_indel$upper[1], "(chr1 indel) keep DP between lines", "DP")
dp_indel_chr2 <- hist_line2_range(indel_dp_chr2,indel_dp_chr2$DP,dp_values_indel$lower[2], 
                dp_values_indel$upper[2], "(chr2 indel) keep DP between lines", "DP")
dp_indel_chr3 <- hist_line2_range(indel_dp_chr3,indel_dp_chr3$DP,dp_values_indel$lower[3], 
                 dp_values_indel$upper[3], "(chr3 indel) keep DP between lines", "DP")
dp_indel_chr4 <- hist_line2_range(indel_dp_chr4,indel_dp_chr4$DP,dp_values_indel$lower[4], 
                 dp_values_indel$upper[4], "(chr4 indel) keep DP between lines", "DP")
dp_indel_chr5 <- hist_line2_range(indel_dp_chr5,indel_dp_chr5$DP,dp_values_indel$lower[5], 
                 dp_values_indel$upper[5], "(chr5 indel) keep DP between lines", "DP")
dp_indel_chr6 <- hist_line2_range(indel_dp_chr6,indel_dp_chr6$DP,dp_values_indel$lower[6], 
                 dp_values_indel$upper[6], "(chr6 indel) keep DP between lines", "DP")
dp_indel_chr7 <- hist_line2_range(indel_dp_chr7,indel_dp_chr7$DP,dp_values_indel$lower[7], 
                 dp_values_indel$upper[7], "(chr7 indel) keep DP between lines", "DP")
dp_indel_chr8 <- hist_line2_range(indel_dp_chr8,indel_dp_chr8$DP,dp_values_indel$lower[8], 
                 dp_values_indel$upper[8], "(chr8 indel) keep DP between lines", "DP")

#DP
#snp
dp_snp_chr1 <- hist_line2_range(snp_dp_chr1,snp_dp_chr1$DP,dp_values_snp$lower[1], 
                 dp_values_snp$upper[1], "(chr1 snp) keep DP between lines", "DP")
dp_snp_chr2 <- hist_line2_range(snp_dp_chr2,snp_dp_chr2$DP,dp_values_snp$lower[2], 
                 dp_values_snp$upper[2], "(chr2 snp) keep DP between lines", "DP")
dp_snp_chr3 <- hist_line2_range(snp_dp_chr3,snp_dp_chr3$DP,dp_values_snp$lower[3], 
                 dp_values_snp$upper[3], "(chr3 snp) keep DP between lines", "DP")
dp_snp_chr4 <- hist_line2_range(snp_dp_chr4,snp_dp_chr4$DP,dp_values_snp$lower[4], 
                 dp_values_snp$upper[4], "(chr4 snp) keep DP between lines", "DP")
dp_snp_chr5 <- hist_line2_range(snp_dp_chr5,snp_dp_chr5$DP,dp_values_snp$lower[5], 
                 dp_values_snp$upper[5], "(chr5 snp) keep DP between lines", "DP")
dp_snp_chr6 <- hist_line2_range(snp_dp_chr6,snp_dp_chr6$DP,dp_values_snp$lower[6], 
                 dp_values_snp$upper[6], "(chr6 snp) keep DP between lines", "DP")
dp_snp_chr7 <- hist_line2_range(snp_dp_chr7,snp_dp_chr7$DP,dp_values_snp$lower[7], 
                 dp_values_snp$upper[7], "(chr7 snp) keep DP between lines", "DP")
dp_snp_chr8 <- hist_line2_range(snp_dp_chr8,snp_dp_chr8$DP,dp_values_snp$lower[8], 
                 dp_values_snp$upper[8], "(chr8 snp) keep DP between lines", "DP")

#Export Cowplots

## Cowplot most export at 11 X 6 inches
plot_grid(MQ_1, MQ_2,AN_1,AN_2,AF_1,AF_2,SOR_1,SOR_2,EH_1,EH_2, ncol = 2)
## Cowplot dp indel export at 8 X 6 inches
plot_grid(dp_indel_chr1,dp_indel_chr2,dp_indel_chr3,dp_indel_chr4,
          dp_indel_chr5,dp_indel_chr6,dp_indel_chr7,dp_indel_chr8, ncol = 2)
## Cowplot dp snp export at 8 X 6 inches
plot_grid(dp_snp_chr1,dp_snp_chr2,dp_snp_chr3,dp_snp_chr4,
          dp_snp_chr5,dp_snp_chr6,dp_snp_chr7,dp_snp_chr8, ncol = 2)







#Single example
# plot 
hist.DF.SLA <- ggplot(indel, aes(X=MQ))+
  geom_histogram(aes(MQ))+
  geom_vline(xintercept=c(50),color="#0000CC")+
#  scale_y_continuous(name="Count")+
#  scale_x_continuous(name="Angle")+
  theme_classic()
hist.DF.SLA <- hist.DF.SLA  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
hist1.N <- hist.DF.SLA  
hist1.N



