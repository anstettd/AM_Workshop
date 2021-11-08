library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly = TRUE)



pop_sample <- args[1] #a file with two columns,first column: sample name, second column: popname. No header
VCF <- args[2] #index file in the same directory
SAVE_DIR <- args[3] #directory to save the analysis. ends with /
PERL_SCRIPT <-  args[4] #directory to vcf2baypass.pl

# pop_sample <- "/home/mjahani/scratch/daniel/baseline_pop_id.txt" #a file with two columns,first column: sample name, second column: popname
# VCF <- "/home/mjahani/scratch/daniel/baseline_filtered_variants.vcf.gz"
# SAVE_DIR <- "/home/mjahani/scratch/daniel/" #ends with /
# PERL_SCRIPT <- "/home/mjahani/scratch/daniel/vcf2baypass.pl"


fread(pop_sample,
      header = F) %>%
  distinct(V1) %>%
  nrow() -> N_SAMPLES


system( #hard filter VCF for quality flags
  paste0("gatk VariantsToTable -V ",
         VCF,
         " -F DP -O ",
         paste0(SAVE_DIR,gsub(".vcf.gz","",gsub(".*/","",VCF)),".DP.table")
  )
)

#extract DP values
fread(paste0(SAVE_DIR,gsub(".vcf.gz","",gsub(".*/","",VCF)),".DP.table"), 
      header = F,
      skip = 1) %>%
  mutate(V1=as.numeric(V1)) %>% 
  pull(V1) -> DPs

mean(DPs) -> mean
sd(DPs) -> SD

mean+SD -> HIGH
mean-SD -> LOW


system( #hard filter VCF for quality flags
  paste0("gatk SelectVariants -V ",
         VCF,
         " -select 'QUAL >= 20.0 && MQ >= 40.0 && AF > 0.03 && AF < 0.97 && AN >= ",
         floor((as.numeric(N_SAMPLES)*2)*0.8), #cutoff for allele Number %80
         " && DP>= ",
         LOW,
         " && DP<= ",
         HIGH,
         "' -O ",
         paste0(SAVE_DIR,gsub(".vcf.gz","",gsub(".*/","",VCF)),".QUAL20_MQ40_AN80_MAF0.03_DP1SD.vcf")
  )
)


system( #using perl script to convert VCF to Baypass format+snp Ids and population order
  paste("cat",
         paste0(SAVE_DIR,gsub(".vcf.gz","",gsub(".*/","",VCF)),".QUAL20_MQ40_AN80_MAF0.03_DP1SD.vcf"),
         "| perl",
         PERL_SCRIPT,
         pop_sample,
         paste0(SAVE_DIR,gsub(".vcf.gz","",gsub(".*/","",VCF)),".QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table"),
         sep = " "
  )
)

system(
  paste0("rm ",
         paste0(SAVE_DIR,gsub(".vcf.gz","",gsub(".*/","",VCF)),".DP.table " ),
         paste0(SAVE_DIR,gsub(".vcf.gz","",gsub(".*/","",VCF)),".QUAL20_MQ40_AN80_MAF0.03_DP1SD.vcf "),
         paste0(SAVE_DIR,gsub(".vcf.gz","",gsub(".*/","",VCF)),".QUAL20_MQ40_AN80_MAF0.03_DP1SD.vcf.tbi")
  )
)

fread(pop_sample,
      header = F) %>%
  group_by(V2) %>%
  tally() %>% 
  mutate(cut_off=as.numeric(n)*2*0.5) %>% #cutoff for each pop allele numer %50
  left_join(
    fread(paste0(SAVE_DIR,gsub(".vcf.gz","",gsub(".*/","",VCF)),".QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order"),
          sep = "\t",
          header=F),
    ., by = c("V1" = "V2")) -> pop_cut_off # calculating a cut of for pops with less than %50 call rate 
  

fread(paste0(SAVE_DIR,gsub(".vcf.gz","",gsub(".*/","",VCF)),".QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table"),
      sep = " ",
      header=F) %>% 
  rowwise() %>%
  mutate(sum_all = sum(across(everything()))) %>% 
  cbind(
    unite(
      fread(paste0(SAVE_DIR,gsub(".vcf.gz","",gsub(".*/","",VCF)),".QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci"),
            sep = "\t",
            header=F),
      "ID",
      V1:V2,
      sep = ":",
      remove = T),
        .)  -> SNP_TABLE 

seq(1,
    by = 2,
    len = nrow(pop_cut_off)) -> POPS

table <- paste0(SAVE_DIR,gsub(".vcf.gz","",gsub(".*/","",VCF)),".QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table")

for (i in 1:length(POPS)) {
  SNP_TABLE %>%
    mutate(!!sym(paste0("pop",i)) := !!sym(paste0("V",POPS[i])) + !!sym(paste0("V",POPS[i]+1))) %>% #population column summary of both alleles
    mutate(!!sym(paste0("pop",i,"_1")) := ifelse(!!sym(paste0("pop",i)) >= as.numeric(pop_cut_off[i,3]),!!sym(paste0("V",POPS[i])),0)) %>% #if the first allele is count smaller than the the threshold replace it with 0
    mutate(!!sym(paste0("pop",i,"_2")) := ifelse(!!sym(paste0("pop",i)) >= as.numeric(pop_cut_off[i,3]),!!sym(paste0("V",POPS[i]+1)),0)) %>% #if the second allele is count smaller than the the threshold replace it with 0
    select(- !!sym(paste0("pop",i)),- !!sym(paste0("V",POPS[i])),- !!sym(paste0("V",POPS[i]+1))) -> SNP_TABLE
}

SNP_TABLE %>%
  rowwise() %>%
  mutate(sum_all_miss = sum(across(starts_with("pop")))) %>%
  mutate(contain_missing_pop=ifelse(sum_all > sum_all_miss,"YES","NO")) %>%
  mutate(missing_percentage=1-(sum_all_miss/(sum(pop_cut_off$n)*2))) -> SNP_TABLE

SNP_TABLE %>%
  ungroup() %>%
  separate(ID, into = c("CHR","POS"), sep = ":", remove = F) -> SNP_TABLE

system(
  paste0("mkdir -p ",
         SAVE_DIR,
         "10k_random_genotype/"
  )
)


SNP_TABLE %>%
  slice_min(as.numeric(missing_percentage), n = 30000) %>%
  sample_n(30000, replace = F) %>%
  distinct(ID) -> BEST_30


BEST_30 %>%
  sample_n(10000, replace = F) %>%
  distinct(ID) %>%
  mutate(random=1) -> random_IDS1

BEST_30 %>%
  filter(!ID %in% pull(random_IDS1,ID)) %>%
  sample_n(10000, replace = F) %>%
  distinct(ID) %>%
  mutate(random=2) -> random_IDS2

BEST_30 %>%
  filter(!ID %in% pull(random_IDS1,ID)) %>%
  filter(!ID %in% pull(random_IDS2,ID)) %>%
  sample_n(10000, replace = F) %>%
  distinct(ID) %>%
  mutate(random=3) -> random_IDS3

rbind(rbind(random_IDS1,random_IDS2),random_IDS3) -> random_IDS


rm(random_IDS1,random_IDS2,random_IDS3,BEST_30)


for (i in 1:3) {
  
  SNP_TABLE %>%
    filter(ID %in% pull(filter(random_IDS,random==i),ID)) %>%
    select(-ID,
           -CHR,
           -POS,
           -sum_all,
           -sum_all_miss,
           -contain_missing_pop,
           -missing_percentage) %>%
    fwrite(paste0(
      paste0(SAVE_DIR,
            "10k_random_genotype/"),
    gsub(".*/","",table),"10K_random",i),
           col.names = F,
           sep = " ",
           quote = F)
  gsub(".vcf.gz","",gsub(".*/","",VCF)) 

  SNP_TABLE %>%
    filter(ID %in% pull(filter(random_IDS,random==i),ID)) %>%
    select(CHR,
           POS) %>%
    fwrite(paste0(paste0(SAVE_DIR,
                         "10k_random_genotype/"
    ),
    gsub(".*/","",table),"10K_random_",i,"_loci"),
           col.names = F,
           sep = " ",
           quote = F)

}

rm(random_IDS)

system(
  paste0("mkdir -p ",
         SAVE_DIR,
         "5k_genotype/"
  )
)


seq(0,nrow(SNP_TABLE),5000)-> SLC
for (i in 1:length(SLC)) {

SNP_TABLE %>% 
  select(-ID,
         -CHR,
         -POS,
         -sum_all,
         -sum_all_miss,
         -contain_missing_pop,
         -missing_percentage) %>%
  slice((as.numeric(SLC[i])+1):(as.numeric(SLC[i])+5000)) %>%
    fwrite(paste0(paste0(SAVE_DIR,
                         "5k_genotype/"
    ),gsub(".*/","",table),"5K_genotype_",i),
           col.names = F,
           sep = " ",
           quote = F)
  
  SNP_TABLE %>% 
    select(CHR,
           POS) %>%
    slice((as.numeric(SLC[i])+1):(as.numeric(SLC[i])+5000)) %>%
    fwrite(paste0(paste0(SAVE_DIR,
                         "5k_genotype/"
    ),gsub(".*/","",table),"5K_genotype_",i,"_loci"),
           col.names = F,
           sep = " ",
           quote = F)
    
}

print(done)

