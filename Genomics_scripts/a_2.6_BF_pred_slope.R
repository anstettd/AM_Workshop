##################################################################################
## Plot Slopes Against BF and Window
## 
## Author Daniel Anstett
## 
## Last Modified October 12, 2022
###################################################################################
library(tidyverse)
library(car)

#Import Slopes
slope_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/slope_mat.csv")
#slope_map <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/slope_map.csv")
#slope_cmd <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/slope_cmd.csv")
slope_mat <- slope_mat %>% mutate(log_p=-log(approx_p,base=10))
slope_mat2 <- slope_mat %>% filter(Slope>-2.5) %>% filter(Slope<2.5)

#Test Case
slope_mat2_p1 <- slope_mat2 %>% filter(Site==1)
lm.mat_p1 <- lm(Slope~BF,data=slope_mat2_p1)
summary(lm.mat_p1)$adj.r.squared[1]

#Loop
clim_slope<-seq(1,12,1)
clim_slope <-as.data.frame(clim_slope)

for (i in 1:12){
  slope_pop <- slope_mat2 %>% filter(Site==i)
  lm.env <- lm(Slope~BF,data=slope_pop)
  clim_slope[i,2] <- summary(lm.env)$coefficients[2,4]
  clim_slope[i,3] <- summary(lm.env)$adj.r.squared[1]
}

###########################################################################################
#Graphs


#Population 1



#Pop 1
slope_mat2_p1 <- slope_mat2 %>% filter(Site==1)
slope_p1 <- ggplot(slope_mat2_p1, aes(x=BF, y=Slope)) + 
  geom_point(size=0.01)+
  geom_smooth(method=lm)+
  scale_y_continuous(name="Slope")+
  scale_x_continuous(name="Log Bayes Factor",breaks = seq(-25, 50, by = 25))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
#slope_p1
ggsave("/Users/daniel_anstett/Dropbox/AM_Workshop/large_graphs/bf_pred_slope_range2_p1.png",
       plot=slope_p1, width=8, height = 6, units = "in")


#full range
slope_p1 <- ggplot(slope_mat2_p1, aes(x=BF, y=Slope)) + 
  geom_point(size=0.01)+
  geom_smooth(method=lm)+
  scale_y_continuous(name="Slope")+
  scale_x_continuous(name="Log Bayes Factor")+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
#slope_p1
ggsave("/Users/daniel_anstett/Dropbox/AM_Workshop/large_graphs/bf_pred_slope_p1.png",
       plot=slope_p1, width=8, height = 6, units = "in")





#All pops
# slope -2 to 2
#slope_mat2_p1 <- slope_mat2 %>% filter(Site==3)
slope_p1 <- ggplot(slope_mat2, aes(x=BF, y=Slope)) + 
  geom_point(size=0.1)+
  geom_smooth(method=lm)+
  scale_y_continuous(name="Slope")+
  scale_x_continuous(name="Log Bayes Factor")+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))+
  facet_wrap(~Site)
#slope_p1
ggsave("/Users/daniel_anstett/Dropbox/AM_Workshop/large_graphs/bf_pred_slope_range2.png",
       plot=slope_p1, width=8, height = 6, units = "in")


#full range
#slope_mat2_p1 <- slope_mat %>% filter(Site==2)
slope_p1 <- ggplot(slope_mat, aes(x=BF, y=Slope)) + 
  geom_point(size=0.1)+
  geom_smooth(method=lm)+
  scale_y_continuous(name="Slope")+
  scale_x_continuous(name="Log Bayes Factor")+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))+
  facet_wrap(~Site)
#slope_p1
ggsave("/Users/daniel_anstett/Dropbox/AM_Workshop/large_graphs/bf_pred_slope_p2.png",
       plot=slope_p1, width=8, height = 6, units = "in")


