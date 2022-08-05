##################################################################################
## WZA Make window
## Author Daniel Anstett
## 
## 
## Last Modified August 2, 2022
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)
library(Kendall)

#Import full snp table for baseline
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
colnames(loci) <- c("chr","snp")
loci_united <- loci %>% unite(chr_snp,"chr","snp",sep="_")

#Make cumulative SNP ID
length_c <- read_csv("Genomics_scripts/Data/chr_length.csv")
loci_win<-loci %>% mutate(snp_c = ifelse(chr == "CE10_chr1", snp,
                                       ifelse(chr == "CE10_chr2", snp+length_c$length_cumu[1],
                                              ifelse(chr == "CE10_chr3", snp+length_c$length_cumu[2],
                                                     ifelse(chr == "CE10_chr4", snp+length_c$length_cumu[3],
                                                            ifelse(chr == "CE10_chr5", snp+length_c$length_cumu[4],
                                                                   ifelse(chr == "CE10_chr6", snp+length_c$length_cumu[5],
                                                                          ifelse(chr == "CE10_chr7", snp+length_c$length_cumu[6], snp+length_c$length_cumu[7]))))))))

#Calculate 10k window
chr_names <- unique(loci_win$chr)
window_size <- 10000
counter <- 0
loci_win_out<-c()
for (j in 1:8){
  loci_chr<-loci_win %>% filter(chr==as.character(chr_names[j]))
  for (i in 1:dim(loci_chr)[1]){
  temp <- loci_chr$snp_c[i]
  if(temp/window_size <=1){
    loci_chr[i,4] <- counter
  } else {
    window_size <- window_size + 10000
    counter <- counter + 1
    loci_chr[i,4] <- counter  
  }
  print(i)
  }
  loci_win_out<-append(loci_win_out, loci_chr[,4])
  counter<-counter+1
  print(chr_names[j])
}

loci_win$V4 <-loci_win_out

#Too large to store on github. Store locally
write_csv(loci_win, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/loci_win.csv")                                                                                 

                                                                                 
###################################################################################








