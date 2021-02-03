library(sf)
library(tmap)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)



calo <- states %>%
  filter(name_en=="Oregon" |
           name_en=="California")


#plot selected states, by wes bbox
tm_shape(calo) +
  tm_borders()
