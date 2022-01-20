## ----------------------------------------------------------------------- ##
                      # Lyon Tutorial - Random Forest
## ----------------------------------------------------------------------- ##
# Code written by: Nick J Lyon

# Clear environment
rm(list = ls())

# Get an object of the directory to this project file
myWD <- getwd()

# Get this package retrieving function
  ## This function will automatically load packages that you already have
  ## and will install packages you don't yet have then load them
ipak <- function(pkg){
  # Function written by Dr. Evan Fricke
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = T)
  sapply(pkg, require, character.only = T)
}

# Define the packages that the script needs
myPackages <- c("tidyverse", "randomForest", "permimp", "vegan")

# Load the packages
ipak(myPackages)

## --------------------------------------------------- ##
               # Part 1: Prepare Data ####
## --------------------------------------------------- ##
# We'll use data from the "vegan" package

# Load some lichen species data
data(varespec)

# And some chemical predictor variables to go along with those species
data(varechem)

# Get just one of the lichen species to use as a response variable
varespec.v2 <- dplyr::select(varespec, Callvulg)

# Bind the response and predictors together
data.obj <- cbind(varechem, varespec.v2)

# Check it out
head(data.obj)

# Want to learn more about the dataset?
?varespec

## --------------------------------------------------- ##
     # Part 2: Random Forest with randomForest ####
## --------------------------------------------------- ##
# Run the random forest
rf1 <- randomForest(Callvulg ~ .,
                    # The 'Y ~ .' format uses all other columns as predictors
                    # Makes formatting your data **crucial**
                    data = data.obj,
                    ntree = 1000,
                    # How many trees should be in the forest
                    mtry = 2,
                    # mtry is # variables / node in tree
                    na.action = na.omit,
                    keep.forest = T,
                    keep.inbag = T)

# Create a variable importance plot
randomForest::varImpPlot(x = rf1,
                         sort = T,
                         n.var = (ncol(data.obj) - 1),
                         main = "Variable Importance")
## Great!

## --------------------------------------------------- ##
 # Part 3: Conditional Permutation Importance (CPI) ####
## --------------------------------------------------- ##
# See vignette for details:
## cran.r-project.org/web/packages/permimp/vignettes/permimp-package.html

# Let's implement conditional permutation
  ## Should see a progress bar in the Console after you run the above
rf1.hiCond <- permimp::permimp(object = rf1,
                               conditional = T,
                               threshold = 0.95,
                               do_check = F)

# Plot distribution of importance
  ## Note this step may take a minute or two
plot(rf1.hiCond, type = "box", horizontal = T)

# Re-fit with a lower threshold
rf1.loCond <- permimp::permimp(object = rf1,
                               conditional = T,
                               threshold = 0.55,
                               do_check = F)

# Plot distribution of this importance
plot(rf1.loCond, type = "box", horizontal = T)

# Make objects of the importance values
rf1.loCond.vals <- data.frame(rf1.loCond$values)
rf1.hiCond.vals <- data.frame(rf1.hiCond$values)

# Side by side comparison of ranking
  ## Ignore warning message
rownames(rf1.loCond.vals)[order(-rf1.loCond.vals)]
rownames(rf1.hiCond.vals)[order(-rf1.hiCond.vals)]

## --------------------------------------------------- ##
          # Part 4: Exploratory Plotting ####
## --------------------------------------------------- ##
# This part is just for fun

# Set custom aesthetic
pref_theme <- theme_classic() + theme(axis.text = element_text(size = 13),
                                      axis.title = element_text(size = 15),
                                      legend.position = "none")

# Graph the top few metrics
ggplot(data.obj, aes(y = Callvulg, x = pH)) +
  geom_point(color = '#a6cee3') +
  geom_smooth(method = 'lm', color = 'black') +
  labs(x = "pH", y = "Callvulg Cover (%)") +
  pref_theme


ggplot(data.obj, aes(y = Callvulg, x = Al)) +
  geom_point(color = '#b2df8a') +
  geom_smooth(method = 'lm', color = 'black') +
  labs(x = "Aluminum", y = "Callvulg Cover (%)") +
  pref_theme


ggplot(data.obj, aes(y = Callvulg, x = P)) +
  geom_point(color = '#b2df8a') +
  geom_smooth(method = 'lm', color = 'black') +
  labs(x = "Phosphorous", y = "Callvulg Cover (%)") +
  pref_theme


ggplot(data.obj, aes(y = Callvulg, x = Baresoil)) +
  geom_point(color = '#33a02c') +
  geom_smooth(method = 'lm', color = 'black') +
  labs(x = "Baresoil (%)", y = "Callvulg Cover (%)") +
  pref_theme

# END ####

