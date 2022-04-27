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

#Import BF>10 Baseline SNP data 
env1 <- read_csv("Genomics_scripts/Data/env1_adapt.csv")
env2 <- read_csv("Genomics_scripts/Data/env2_adapt.csv")
env5 <- read_csv("Genomics_scripts/Data/env5_adapt.csv")

#Filter to BF>20 data
env120 <- env1 %>% filter(BF>20)
env220 <- env2 %>% filter(BF>20)
env520 <- env5 %>% filter(BF>20)

#Venn Diagram of SNP overlap
MAT_set <- pull(env120[,1])
MAP_set <- pull(env220[,1])
CMD_set <- pull(env520[,1])

ven_3 <- list(MAT_set,MAP_set,CMD_set)

VD_1<-venn.diagram(x=ven_3,
                   category.names = c("MAT" , "MAP" , "CMD"),
                   fill = c("#F8766D", "#00BFC4","#C77CFF"),
                   #alpha = c(0.5, 0.5, 0.5),
                   main.cex = 5,
                   cat.cex = 0, cex=1.8,
                   fontface = "bold",
                   filename = NULL,
)

pdf("Graphs_overlap/venn_20.pdf")
grid.newpage()
grid.draw(VD_1)
dev.off()
