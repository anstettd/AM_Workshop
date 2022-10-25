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
offset_pop <- read_csv("Genomics_scripts/Data/offset_pop_lambda.csv")
offset_pop_mod <- offset_pop %>% filter(Paper_ID!=10)

#stats
lm.1215 <- lm(lambda_slope~offset_1215,data=offset_pop)
lm.1215_mod <- lm(lambda_slope~offset_1215,data=offset_pop_mod)
#lm.2014 <- lm(lambda_slope~offset_2014,data=offset_pop)
#lm.2015 <- lm(lambda_slope~offset_2015,data=offset_pop)
#lm.2016 <- lm(lambda_slope~offset_2016,data=offset_pop)

summary(lm.1215)
summary(lm.1215_mod)

Anova(lm.1215,type="III")
Anova(lm.1215_mod,type="III")
#Anova(lm.2014,type="III")
#Anova(lm.2015,type="III")
#Anova(lm.2016,type="III")




#plots
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



#plots
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








