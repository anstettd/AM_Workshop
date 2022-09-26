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
offset_pop_s <- read_csv("Genomics_scripts/Data/offset_pop_s.csv")
offset_pop_s <- offset_pop_s %>% mutate(Region = ifelse(Lat >= 40, "North", 
                                            ifelse((Lat >35) & (Lat <40), "Centre","South")))
#stats
lm.2016 <- lm(Moderate_Selection~offset_2016,data=offset_pop_s)
lm.2015 <- lm(Moderate_Selection~offset_2015,data=offset_pop_s)
lm.2014 <- lm(Moderate_Selection~offset_2014,data=offset_pop_s)

Anova(lm.2016,type="III")
Anova(lm.2015,type="III")
Anova(lm.2014,type="III")



#plots
ggplot(offset_pop_s, aes(x=offset_2016, y=Moderate_Selection, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
#  geom_label_repel(aes(label = Paper_ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Moderate Selection (S > 0.2)")+
  scale_x_continuous(name="Genetic Offset 2016")+
  scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
ggsave("Offset_slope/1_offset_slope_2016.pdf",width=8, height = 6, units = "in")



ggplot(offset_pop_s, aes(x=offset_2015, y=Moderate_Selection, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = Paper_ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Moderate Selection (S > 0.2)")+
  scale_x_continuous(name="Genetic Offset 2015")+
  scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
ggsave("Offset_slope/1_offset_slope_2015.pdf",width=8, height = 6, units = "in")

ggplot(offset_pop_s, aes(x=offset_2014, y=Moderate_Selection)) + 
  geom_point()+
  geom_smooth(method=lm)+
  scale_y_continuous(name="Moderate Selection")+
  scale_x_continuous(name="Genetic Offset 2014")+
  theme_classic() + theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
ggsave("Offset_slope/1_offset_slope_2014.pdf",width=10, height = 6, units = "in")

ggplot(offset_pop_s, aes(x=offset_2013, y=Moderate_Selection)) + 
  geom_point()+
  geom_smooth(method=lm)+
  scale_y_continuous(name="Moderate Selection")+
  scale_x_continuous(name="Genetic Offset 2013")+
  theme_classic() + theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
#ggsave("Offset_slope/1_offset_slope_2013.pdf",width=10, height = 6, units = "in")

ggplot(offset_pop_s, aes(x=offset_2012, y=Moderate_Selection)) + 
  geom_point()+
  geom_smooth(method=lm)+
  scale_y_continuous(name="Moderate Selection")+
  scale_x_continuous(name="Genetic Offset 2012")+
  theme_classic() + theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
#ggsave("Offset_slope/1_offset_slope_2012.pdf",width=10, height = 6, units = "in")

ggplot(offset_pop_s, aes(x=offset_2011, y=Moderate_Selection)) + 
  geom_point()+
  geom_smooth(method=lm)+
  scale_y_continuous(name="Moderate Selection")+
  scale_x_continuous(name="Genetic Offset 2011")+
  theme_classic() + theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
#ggsave("Offset_slope/1_offset_slope_2011.pdf",width=10, height = 6, units = "in")








#Plot offset against latitude
#4.5
ggplot(offset_pop_s, aes(x=Elevation, y=offset_4.5_peakbf2)) + 
  geom_point()+
  geom_smooth(method=lm)+
  scale_y_continuous(name="Genetic Offset RCP4.5")+
  scale_x_continuous(name="Elevation (m)")+
  theme_classic() + theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))

ggplot(offset_pop_s, aes(x=Lat, y=offset_4.5_peakbf2)) + 
  geom_point()+
  geom_smooth(method=lm)+
  scale_y_continuous(name="Genetic Offset RCP4.5")+
  scale_x_continuous(name="Latitude")+
  theme_classic() + theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))

#8.5
ggplot(offset_pop_s, aes(x=Elevation, y=offset_8.5_peakbf2)) + 
  geom_point()+
  geom_smooth(method=lm)+
  scale_y_continuous(name="Genetic Offset RCP8.5")+
  scale_x_continuous(name="Elevation (m)")+
  theme_classic() + theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))

ggplot(offset_pop_s, aes(x=Lat, y=offset_8.5_peakbf2)) + 
  geom_point()+
  geom_smooth(method=lm)+
  scale_y_continuous(name="Genetic Offset RCP8.5")+
  scale_x_continuous(name="Latitude")+
  theme_classic() + theme(
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=12,face="bold"),
    axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))







