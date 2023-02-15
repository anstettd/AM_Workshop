##################################################################################
## Make venn diagram for MAT, MAP and CMD SNPs at BF>20
## Author Daniel Anstett 
## 
##
## Last Modified July 8, 2021
###################################################################################
library(tidyverse)
library(VennDiagram)

###################################################################################

#Import double filtered Baseline SNP data 
env1 <- read_csv("Genomics_scripts/Data/win_bf_mat30_5.csv")
env2 <- read_csv("Genomics_scripts/Data/win_bf_map30_5.csv")
env5 <- read_csv("Genomics_scripts/Data/win_bf_cmd30_5.csv")

env_all <- rbind(env1,env2,env5)

#write out unique SNP set
snp_set <- env_all %>% distinct(chr_snp)
#write_csv(snp_set, "Genomics_scripts/Data/snp_set.csv")



#Venn Diagram of SNP overlap
MAT_set <- pull(env1[,1])
MAP_set <- pull(env2[,1])
CMD_set <- pull(env5[,1])

ven_3 <- list(MAT_set,MAP_set,CMD_set)

VD_1<-venn.diagram(x=ven_3,
                   category.names = c("MAT" , "MAP" , "CMD"),
                   #fill = c("#F8766D", "#00BFC4","#C77CFF"),
                   fill = c("yellow", "blue","pink"),
                   #alpha = c(0.5, 0.5, 0.5),
                   main.cex = 5,
                   cat.cex = 0, cex=1.8,
                   fontface = "bold",
                   filename = NULL,
)

pdf("Graphs_overlap/venn_double_filter.pdf")
grid.newpage()
grid.draw(VD_1)
dev.off()
