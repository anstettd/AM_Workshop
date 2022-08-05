#############################################################################################################
## Calculate windows for WZA
## Author Daniel Anstett
## 
## Modified from Tom Booker WZA Vignette
## Last Modified August 3, 2022
#############################################################################################################
#Functions
#Calculate frequency
prop_A <- function(snp_table) {
  snp_prop_A<- snp_table %>% select (chr_snp)
  counter=2
  pop_num=1
  for (i in seq(1,dim(snp_table)[2]-1,2)) {
    
    snpA<-paste("V",i, sep="") #sets up string for snpA for pop num
    snpB<-paste("V",i+1, sep="") #sets up string for snpA for pop num
    P<-paste("P", pop_num, sep="") #sets up paper_ID (population ID)
    tmp<-snp_table %>% select(snpA,snpB) #graps all snps per paper_ID
    
    colnames(tmp)<-c("A", "B") #renames column headers
    
    snp_prop_A[,counter]<-tmp$A/(tmp$A + tmp$B) #calc proportion for A, assign to output
    colnames (snp_prop_A)[counter]<-P 
    
    counter<-counter+1
    pop_num<-pop_num+1
  }
  #snp_prop_A[is.na(snp_prop_A)] <- 0
  return(snp_prop_A)
}

#############################################################################################################
#Import libraries
library(tidyverse)
library(Kendall)

#Import Climate Data
climate <- read_csv("Genomics_scripts/Data/env_baseline.csv",col_names=FALSE)
climate2 <- filter(climate[1:2,])
climate5 <- filter(climate[5,])
climate_wza <- rbind(climate2,climate5)
mat_pop <- as.numeric(climate_wza[1,])
map_pop <- as.numeric(climate_wza[2,])
cmd_pop <- as.numeric(climate_wza[3,])


optima <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/WZA_vignette/environments.1_0.5_192.alleleFreqs.csv", header = F)

#Import full snp table for baseline
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
loci_win<-read_csv("Genomics_scripts/Data/loci_win.csv")

#Calculate frequency using prop_A
colnames(loci) <- c("chr","snp")
loci_united <- loci %>% unite(chr_snp,"chr","snp",sep="_")
snp_chr_snp <- cbind(loci_united,snp)
freq <- prop_A(snp_chr_snp)
all_data <-cbind(loci_win,freq) #add snp lables to rows
colnames(all_data)[4] <- "win"


#############################################################################################################
#############################################################################################################
#Tutorial

#Filter for minor allele frequency
# Let's split the data so that we have an allele freq. dataframe and a dataframe with SNP info
snps_all <- all_data[,1:5]

# let's take a look at the data, there should be SNPs in gene0 on chromosome 1 
str(snps_all)

# extract frequency data for each SNP
freqs_all <- all_data[,5:dim(all_data)[2]] # There are 40 demes in the dataset, so columns 4 to 43 represent all the demes

# take a look at the freuqency for the first SNP
as.numeric(freqs_all[1,])

# As you can see, most SNPs are at a low minor allele frequency - we remove the low frequency ones before applying the WZA.

# The first thing we'll do is calculate the average allele frequency across populations for each SNP and add it to the SNP info dataframe
#Calc row means
all_data<-all_data %>% mutate(p_bar = rowMeans(select(all_data, P1:P55), na.rm = TRUE))
all_data$q_bar <- 1 - all_data$p_bar

# take the minimum to get the MAF
all_data$MAF <- pmin(all_data$p_bar, all_data$q_bar)

# now we'll apply the MAF filter to the data
#separate
snps <- snps_all[all_data$MAF >= 0.05, ]
freqs <- freqs_all[all_data$MAF >= 0.05, ]

#together
all_data_filter <- all_data %>% filter(MAF>=0.05)

#############################################################################################################
#Calc WZA

#Make a little wrapper function for cor.test to return the summary stat and the p-value as a vector
cor_test_wrapper <- function(p_vec, env_vector){
  correlation_result <- cor.test(p_vec, env_vector, method = "kendall", exact = F)
  return(c(correlation_result$estimate,
           correlation_result$p.val))
}

## Use the apply function to use the wrapper function on each line of the DF - this step takes a while...
cor_results_1 <- apply(as.matrix(freqs[,-1]),
                      1,
                      function(x) cor_test_wrapper(x, mat_pop))
cor_results_2 <- apply(as.matrix(freqs[,-1]),
                       1,
                       function(x) cor_test_wrapper(x, map_pop))
cor_results_3 <- apply(as.matrix(freqs[,-1]),
                       1,
                       function(x) cor_test_wrapper(x, cmd_pop))

# transpose the result and store as a dataframe
cor_results_1_t <- as.data.frame(t(cor_results_1))
cor_results_2_t <- as.data.frame(t(cor_results_2))
cor_results_3_t <- as.data.frame(t(cor_results_3))

# give the column informative names
names(cor_results_1_t) <- c("Kendall", "p_val")
names(cor_results_2_t) <- c("Kendall", "p_val")
names(cor_results_3_t) <- c("Kendall", "p_val")

# combine the result with the SNP table
snps_mat <- cbind(snps , all_data_filter[61:63] , cor_results_1_t)
snps_map <- cbind(snps , all_data_filter[61:63] , cor_results_2_t)
snps_cmd <- cbind(snps , all_data_filter[61:63] , cor_results_3_t)

snps_mat <- snps_mat %>% na.omit(snps_mat)
snps_map <- snps_map %>% na.omit(snps_map)
snps_cmd <- snps_cmd %>% na.omit(snps_cmd)

min(snps_mat$p_val)
min(snps_map$p_val)
min(snps_cmd$p_val)

# Make empirical p-value
snps_mat$empirical_p <- rank(snps_mat$p_val)/length(snps_mat$p_val)
snps_map$empirical_p <- rank(snps_mat$p_val)/length(snps_map$p_val)
snps_cmd$empirical_p <- rank(snps_mat$p_val)/length(snps_cmd$p_val)

#Too large to store on github. Store locally
write_csv(snps_mat, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_mat.csv")
write_csv(snps_map, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_map.csv")     
write_csv(snps_cmd, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_cmd.csv")     


