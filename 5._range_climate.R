#Find range of 1981 to 2010 climate varaibles

library(tidyverse) # for data manipulation

normal.clim <- read_csv("Donor_selection/data/paper_ID_site_select_Normal_1981_2010Y.csv")

range(normal.clim$MAT)
range(normal.clim$MAP)