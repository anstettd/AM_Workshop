##################################################################################
## Plot Slopes Against BF and Window
## 
## Author Daniel Anstett
## 
## Last Modified October 12, 2022
###################################################################################
library(tidyverse)

#Import Slopes
slope_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/slope_mat.csv")
#slope_map <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/slope_map.csv")
#slope_cmd <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/slope_cmd.csv")
slope_mat <- slope_mat %>% mutate(log_p=-log(approx_p,base=10))
slope_mat25 <- slope_mat %>% filter(Slope>-25) %>% filter(Slope<25)
slope_mat15 <- slope_mat %>% filter(Slope>-15) %>% filter(Slope<15)
slope_mat5 <- slope_mat %>% filter(Slope>-5) %>% filter(Slope<5)
slope_mat2 <- slope_mat %>% filter(Slope>-2) %>% filter(Slope<2)


#Make Sites
slope_mat_p1 <- slope_mat %>% filter(Site==1)

###################################################################################

#BF All range

slope_bf <- ggplot(slope_mat, aes(x=Slope, y=BF)) + 
  geom_point(aes(), size = 0.2)+
  #geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = Paper_ID))+
  scale_y_continuous(name="Log Bayes Factor",breaks=seq(-40,40,20))+
  scale_x_continuous(name="Slope")+
  #scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))+
  facet_wrap(~Site)
  

ggsave("/Users/daniel_anstett/Dropbox/AM_Workshop/large_graphs/slope_bf_mat_facet50.png",
       plot=slope_bf, width=8, height = 6, units = "in")

#BF -25 to 25
slope_bf25 <- ggplot(slope_mat25, aes(x=Slope, y=BF)) + 
  geom_point(aes(), size = 0.2)+
  #geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = Paper_ID))+
  scale_y_continuous(name="Log Bayes Factor",breaks=seq(-40,40,20))+
  scale_x_continuous(name="Slope")+
  #,limits=c(-25,25))+
  #scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))+
  facet_wrap(~Site)

ggsave("/Users/daniel_anstett/Dropbox/AM_Workshop/large_graphs/slope_bf_mat_facet25.png",
       plot=slope_bf25, width=8, height = 6, units = "in")

#BF -5 to 5
slope_bf5 <- ggplot(slope_mat5, aes(x=Slope, y=BF)) + 
  geom_point(aes(), size = 0.2)+
  #geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = Paper_ID))+
  scale_y_continuous(name="Log Bayes Factor",breaks=seq(-40,40,20))+
  scale_x_continuous(name="Slope")+
  #,limits=c(-25,25))+
  #scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))+
  facet_wrap(~Site)

ggsave("/Users/daniel_anstett/Dropbox/AM_Workshop/large_graphs/slope_bf_mat_facet5.png",
       plot=slope_bf5, width=8, height = 6, units = "in")

#BF -2 to 2
slope_bf2 <- ggplot(slope_mat2, aes(x=Slope, y=BF)) + 
  geom_point(aes(), size = 0.2)+
  #geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = Paper_ID))+
  scale_y_continuous(name="Log Bayes Factor",breaks=seq(-40,40,20))+
  scale_x_continuous(name="Slope")+
  #,limits=c(-25,25))+
  #scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))+
  facet_wrap(~Site)

ggsave("/Users/daniel_anstett/Dropbox/AM_Workshop/large_graphs/slope_bf_mat_facet2.png",
       plot=slope_bf2, width=8, height = 6, units = "in")



###################################################################################

#Window

#log_p All range

slope_log_p <- ggplot(slope_mat, aes(x=Slope, y=log_p)) + 
  geom_point(aes(), size = 0.2)+
  #geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = Paper_ID))+
#  scale_y_continuous(name="Log Bayes Factor",breaks=seq(-40,40,20))+
  scale_x_continuous(name="Slope")+
  #scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))+
  facet_wrap(~Site)
ggsave("/Users/daniel_anstett/Dropbox/AM_Workshop/large_graphs/slope_log_p_mat_facet50.png",
       plot=slope_log_p, width=8, height = 6, units = "in")

#log_p -25 to 25
slope_log_p25 <- ggplot(slope_mat25, aes(x=Slope, y=log_p)) + 
  geom_point(aes(), size = 0.2)+
  #geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = Paper_ID))+
#  scale_y_continuous(name="Log Bayes Factor",breaks=seq(-40,40,20))+
  scale_x_continuous(name="Slope")+
  #,limits=c(-25,25))+
  #scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))+
  facet_wrap(~Site)

ggsave("/Users/daniel_anstett/Dropbox/AM_Workshop/large_graphs/slope_log_p_mat_facet25.png",
       plot=slope_log_p25, width=8, height = 6, units = "in")

#log_p -5 to 5
slope_log_p5 <- ggplot(slope_mat5, aes(x=Slope, y=log_p)) + 
  geom_point(aes(), size = 0.2)+
  #geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = Paper_ID))+
  #scale_y_continuous(name="Log Bayes Factor",breaks=seq(-40,40,20))+
  scale_x_continuous(name="Slope")+
  #,limits=c(-25,25))+
  #scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))+
  facet_wrap(~Site)

ggsave("/Users/daniel_anstett/Dropbox/AM_Workshop/large_graphs/slope_log_p_mat_facet5.png",
       plot=slope_log_p5, width=8, height = 6, units = "in")

#log_p -2 to 2
slope_log_p2 <- ggplot(slope_mat2, aes(x=Slope, y=log_p)) + 
  geom_point(aes(), size = 0.2)+
  #geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = Paper_ID))+
  #scale_y_continuous(name="Log Bayes Factor",breaks=seq(-40,40,20))+
  scale_x_continuous(name="Slope")+
  #,limits=c(-25,25))+
  #scale_color_manual(values= c("North"="#3399FF", "Centre"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))+
  facet_wrap(~Site)

ggsave("/Users/daniel_anstett/Dropbox/AM_Workshop/large_graphs/slope_log_p_mat_facet2.png",
       plot=slope_log_p2, width=8, height = 6, units = "in")

