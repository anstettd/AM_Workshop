#Large
#pop_snp = p1A
#freq_count = freq_count_MAT_IDX
large <- function(pop_snp,freq_count){
  p1_1 <- pop_snp %>% filter(snpA >= 0 & snpA <0.1)
  p1_2 <- pop_snp %>% filter(snpA >= 0.1 & snpA <0.2)
  p1_3 <- pop_snp %>% filter(snpA >= 0.2 & snpA <0.3)
  p1_4 <- pop_snp %>% filter(snpA >= 0.3 & snpA <0.4)
  p1_5 <- pop_snp %>% filter(snpA >= 0.4 & snpA <0.5)
  
  p1_6 <- pop_snp %>% filter(snpA >= 0.5 & snpA <0.6)
  p1_7 <- pop_snp %>% filter(snpA >= 0.6 & snpA <0.7)
  p1_8 <- pop_snp %>% filter(snpA >= 0.7 & snpA <0.8)
  p1_9 <- pop_snp %>% filter(snpA >= 0.8 & snpA <0.9)
  p1_10 <- pop_snp %>% filter(snpA >= 0.9 & snpA <=1)
  
  list_mat_p1_1 <- sample.int(dim(p1_1)[1],freq_count[1])
  list_mat_p1_2 <- sample.int(dim(p1_2)[1],freq_count[2])
  list_mat_p1_3 <- sample.int(dim(p1_3)[1],freq_count[3])
  list_mat_p1_4 <- sample.int(dim(p1_4)[1],freq_count[4])
  list_mat_p1_5 <- sample.int(dim(p1_5)[1],freq_count[5])
  list_mat_p1_6 <- sample.int(dim(p1_6)[1],freq_count[6])
  list_mat_p1_7 <- sample.int(dim(p1_7)[1],freq_count[7])
  list_mat_p1_8 <- sample.int(dim(p1_8)[1],freq_count[8])
  list_mat_p1_9 <- sample.int(dim(p1_9)[1],freq_count[9])
  list_mat_p1_10 <- sample.int(dim(p1_10)[1],freq_count[10])
  
  rand_pop <- data.frame()
  
  for (i in 1:length(list_mat_p1_1)){
    rand_pop <- rbind(rand_pop ,p1_1 %>% 
                        filter(chr_snp==as.character(p1_1$chr_snp[list_mat_p1_1[i]])))
  }
  for (i in 1:length(list_mat_p1_2)){
    rand_pop <- rbind(rand_pop ,p1_2 %>% 
                        filter(chr_snp==as.character(p1_2$chr_snp[list_mat_p1_2[i]])))
  }
  for (i in 1:length(list_mat_p1_3)){
    rand_pop <- rbind(rand_pop ,p1_3 %>% 
                        filter(chr_snp==as.character(p1_3$chr_snp[list_mat_p1_3[i]])))
  }
  for (i in 1:length(list_mat_p1_4)){
    rand_pop <- rbind(rand_pop ,p1_4 %>% 
                        filter(chr_snp==as.character(p1_4$chr_snp[list_mat_p1_4[i]])))
  }
  for (i in 1:length(list_mat_p1_5)){
    rand_pop <- rbind(rand_pop ,p1_5 %>% 
                        filter(chr_snp==as.character(p1_5$chr_snp[list_mat_p1_5[i]])))
  }
  
  for (i in 1:length(list_mat_p1_6)){
    rand_pop <- rbind(rand_pop ,p1_6 %>% 
                        filter(chr_snp==as.character(p1_6$chr_snp[list_mat_p1_6[i]])))
  }
  for (i in 1:length(list_mat_p1_7)){
    rand_pop <- rbind(rand_pop ,p1_7 %>% 
                        filter(chr_snp==as.character(p1_7$chr_snp[list_mat_p1_7[i]])))
  }
  for (i in 1:length(list_mat_p1_8)){
    rand_pop <- rbind(rand_pop ,p1_8 %>% 
                        filter(chr_snp==as.character(p1_8$chr_snp[list_mat_p1_8[i]])))
  }
  for (i in 1:length(list_mat_p1_9)){
    rand_pop <- rbind(rand_pop ,p1_9 %>% 
                        filter(chr_snp==as.character(p1_9$chr_snp[list_mat_p1_9[i]])))
  }
  for (i in 1:length(list_mat_p1_10)){
    rand_pop <- rbind(rand_pop ,p1_10 %>% 
                        filter(chr_snp==as.character(p1_10$chr_snp[list_mat_p1_10[i]])))
  }
  
  return(rand_pop)
}