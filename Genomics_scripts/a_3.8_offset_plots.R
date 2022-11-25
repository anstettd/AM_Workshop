##################################################################################
## Gradient forest plots
## Author Daniel Anstett
## 
## 
## Last Modified September 19, 2022
###################################################################################

#Library install and import
library(tidyverse) 
library(car)
library(ggrepel)

#Import data
offset_pop_nocumul <- read_csv("Genomics_scripts/Data/offset_pop_lambda.csv")
slope.cumul <- read_csv("Genomics_scripts/Data/slope.cumul.csv")
offset_pop <- cbind(offset_pop_nocumul,slope.cumul)
colnames(offset_pop)[18] <- "cumul_MAT" ; colnames(offset_pop)[19] <- "cumul_MAP"

offset_pop_mod <- offset_pop %>% filter(Paper_ID!=10)

#stats
lm.1215 <- lm(lambda_slope~offset_1215,data=offset_pop)
lm.1215_mod <- lm(lambda_slope~offset_1215,data=offset_pop_mod)

lm.dist <- lm(lambda_slope~clim_distance,data=offset_pop)
lm.dist_mod <- lm(lambda_slope~clim_distance,data=offset_pop_mod)


#Without pop 10
summary(lm.1215_mod)
Anova(lm.1215_mod,type="III")

summary(lm.dist_mod)
Anova(lm.dist_mod,type="III")


# With pop 10
summary(lm.1215)
Anova(lm.1215,type="III")

summary(lm.dist)
Anova(lm.dist,type="III")




#######################################################################################################3####

#2012-2015 offset plotted against lambda
ggplot(offset_pop, aes(x=offset_1215, y=lambda_slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
#  geom_label_repel(aes(label = Paper_ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2012-2015 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
#ggsave("Graphs_Oct_22/offset_lambda_slope_1215.pdf",width=8, height = 6, units = "in")



#2012-2015 offset plotted against lambda, rm pop 10
ggplot(offset_pop_mod, aes(x=offset_1215, y=lambda_slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = Paper_ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2012-2015 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
#ggsave("Graphs_Oct_22/offset_lambda_slope_1215_mod.pdf",width=7, height = 5, units = "in")

###########################

#climate distance plotted against lambda
ggplot(offset_pop, aes(x=clim_distance, y=lambda_slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = Paper_ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2012-2015 Climate Distance")+
  scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs_Oct_22/5_distance_lambda.pdf",width=8, height = 6, units = "in")

#climate distance plotted against lambda, rm pop 10
ggplot(offset_pop_mod, aes(x=clim_distance, y=lambda_slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = Paper_ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="2012-2015 Climate Distance")+
  scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs_Oct_22/5_distance_lambda_mod.pdf",width=7, height = 5, units = "in")




#######################################################################################################3####

#Cumulative slope plotted against lambda. No Site 10
#MAT
ggplot(offset_pop_mod, aes(x=cumul_MAT, y=lambda_slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  #  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = Paper_ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="Cumul. Slope MAT")+
  scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(1.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs_Oct_22/lamda_cumul_slope_mat.pdf",width=7, height = 5, units = "in")

#MAP
ggplot(offset_pop_mod, aes(x=cumul_MAP, y=lambda_slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  #  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = Paper_ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="Cumul. Slope MAP")+
  scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(1.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs_Oct_22/lamda_cumul_slope_map.pdf",width=7, height = 5, units = "in")


#CMD
ggplot(offset_pop_mod, aes(x=cumul_CMD, y=lambda_slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  #  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = Paper_ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="Cumul. Slope CMD")+
  scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(1.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("Graphs_Oct_22/lamda_cumul_slope_cmd.pdf",width=7, height = 5, units = "in")

