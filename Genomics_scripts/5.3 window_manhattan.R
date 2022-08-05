#############################################################################################################
## Calc WZA and make window manhattan plots
## Author Daniel Anstett
## 
## 
## Modified from Tom Booker WZA Vignette
## Last Modified August 3, 2022
#############################################################################################################
#Import libraries

library(tidyverse)
library(Kendall)

#Import files

snps_mat <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_mat.csv")
snps_map <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_map.csv")
snps_cmd <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_cmd.csv")


###########################################################################################################
# Make empirical p-value
snps_mat$empirical_p <- rank(snps_mat$p_val)/length(snps_mat$p_val)

## Plot those out p-values - it should look like a squashed version of the above Manhattan plot...

ggplot(data = snps_mat, 
       aes(x = snp_c/1e6, 
           y = -log10(empirical_p),
           fill = p_bar*q_bar))+
  geom_point(shape = 21)+
  geom_hline(aes(yintercept = quantile(-log10(empirical_p), 0.99)), col = "red", lty = 2, lwd = 1)+
  scale_fill_gradient(expression(italic(bar(p)*bar(q))), low = "white", high = "black")+
  scale_y_continuous(expression(-log[10]*"(p-value)"), limits = c(0,10))+
  scale_x_continuous("Position (Mbp)")+
  theme_bw()



###########################################################################################################
#Calc WZA


# Now convert the empirical p-values to z scores
snps_mat$z_score <- qnorm(snps_mat$empirical_p, lower.tail = F)
head(snps_mat)

# For convenience, let's calculate p_bar*q_bar and store it in the DF 
snps_mat$pbar_qbar <- snps_mat$p_bar*snps_mat$q_bar


# Now let's calculate the WZA for each gene
WZA_numerator_mat <- tapply(snps_mat$z_score*snps_mat$pbar_qbar,
                            snps_mat$win,
                            sum)
WZA_denominator_mat <- sqrt( tapply((snps_mat$pbar_qbar*snps_mat$pbar_qbar), 
                                    snps_mat$win, 
                                    sum) )
WZA_mat <- WZA_numerator_mat/WZA_denominator_mat


pos_mat <- tapply(snps_mat$snp_c,
                  snps_mat$win,
                  mean)
WZA_df_mat <- data.frame(WZA_mat = WZA_mat, 
                         pos_mat = pos_mat)

#Get window and chromosome info into WZA data frame
win_mat <- snps_mat %>% distinct(chr, win)




#Plot Distribution
ggplot(data = WZA_df_mat, aes( x = WZA_mat))+
  geom_histogram( bins = 50)+
  scale_y_continuous("Count")+
  theme_bw()

# Now let's plot out the result:
ggplot(data = WZA_df_mat, aes( x = pos/1e6, y = WZA_mat))+
  geom_point(size=0.5)+
  scale_y_continuous("WZA Score")+
  scale_x_continuous("Position (Mbp)")+
  theme_classic()+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=15, face="bold"),
    axis.text.y = element_text(size=15,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 0.8, face="bold",hjust=0.5)
  )



# 1% of genome
ggplot(data = WZA_df_mat, aes( x = pos/1e6, y = WZA_mat))+
  geom_point(size=0.5)+
  scale_y_continuous("WZA Score")+
  scale_x_continuous("Position (Mbp)")+
  geom_hline(aes(yintercept = quantile(WZA_mat, 0.99)), 
             col = "red",
             lty = 2,
             lwd = 1)+
  theme_classic()+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=15, face="bold"),
    axis.text.y = element_text(size=15,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 0.8, face="bold",hjust=0.5)
  )

#High windows
hit_set <- row.names(WZA_df[WZA_df$WZA > quantile(WZA_df$WZA, 0.99),])
hit_set

#######################################################################################################
#Calc empirical p-value based on WZA scores

## Calculate the mean and sd of the WZA distribution...
WZA_mean_mat <- mean(WZA_df_mat$WZA_mat)
WZA_sd_mat <- sd(WZA_df_mat$WZA_mat)

## Calculate an approximate p-value based on the assumption of normality...
WZA_df_mat$approx_p <- 2*pnorm(-abs(WZA_df_mat$WZA), mean = WZA_mean_mat, sd= WZA_sd_mat)

# WZA Empirical p-value
ggplot(data = WZA_df_mat, aes( x = pos/1e6, y = -log10(approx_p)))+
  geom_point(size=0.5)+
  scale_y_continuous("-log10(WZA Empirical p-value)")+
  geom_hline(aes(yintercept = -log10(0.05)), 
             col = "red",
             lty = 2,
             lwd = 1)+
  scale_x_continuous("Position (Mbp)")+
  theme_classic()+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=15, face="bold"),
    axis.text.y = element_text(size=15,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 0.8, face="bold",hjust=0.5)
  )

# FDR Correction
ggplot(data = WZA_df_mat, aes( x = pos/1e6, y = -log10(p.adjust(approx_p,method = "fdr"))))+
  geom_point( shape = 21, fill = "grey",size=0.5 )+
  scale_y_continuous("-log10(q-value)", limits = c(0,5))+
  geom_hline(aes(yintercept = -log10(0.05)), 
             col = "red",
             lty = 2,
             lwd = 1)+
  scale_x_continuous("Position (Mbp)")+
  theme_classic()+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=15, face="bold"),
    axis.text.y = element_text(size=15,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 0.8, face="bold",hjust=0.5)
  )


