## Initial data exploration - trifolium
# -- DATA IMPORT (NAT) ----
#1. set data pathway 
datpath <- "~/Desktop/USDA-compost/Data/Competition/"  ##NAT TO EDIT 

# list files in entered data folder
datfiles <- dats <- list.files(paste0(datpath, "Competition_CleanedData"), full.names = T)

# read in raw trifolium data
dat<- read.csv(paste0(datpath, "Competition_CleanedData/competition_trifolium_seeds_2021.csv"), strip.white = T)

## load packages
library(tidyverse)

## divide total seeds by total stems to standardize seed output
dat <- dat %>%
  mutate(seeds_per_stem = tot_seeds/tot_stems)


## VISUALLY EXPLORE DATA! 
## want to visualize seed output by precipitation condition
## how does this differ by competition background?

## seeds per stem vs. precipitation treatment; faceted by competition background
ggplot(dat, aes(x=ppt_trt, y=seeds_per_stem)) +
  geom_boxplot() +
  facet_wrap(~background)

## seeds per stem vs. competition background; faceted by precipitation treatment
ggplot(dat, aes(x=background, y=seeds_per_stem)) +
  geom_boxplot() +
  facet_wrap(~ppt_trt)


## ANOVAs 
  ## name the model something unique
  ## use function aov(response variable ~ predictor, data = your data frame)

ppt <- aov(seeds_per_stem~ppt_trt, data = dat)
tests <- summary(ppt)
str(tests)
tests # Pr(>F) = p-value 

## Tukey's Honest significant difference test; I like to use this when I have a lot of categories
## you can see the differences between each pair of categories
TukeyHSD(ppt) # p adj is p-values 


## A few potential next steps (and some pseudo-code just for an example): 

  ## we still have a wide range of nutrient treatments in this data set
  ## this could explain some of the variation we are seeing in seed output
  ## we might want to filter the data to only include control nutrient treatments (rather than N fertilization or compost)
  ## tidyverse is really good at filtering

no_nut <- dat %>%
  filter(nut_trt == "C") ## I moved the script to github adn can no longer see the data, 
  ## but go ahead and fill in the name of the nutrient column here and the specific value for the control plots 
  ## R will then filter out only control treatments

  ## we could also compare specific species -> also using the filtering function
  
  ## if you want to compare mean values between groups, I like to use tidyverse's group_by() and then summarize() functions 
mean <- dat %>%
  group_by(background, ppt_trt, any_other_variables_you_want_to_find_the_average_over) %>% 
  summarize(your_new_column_name = mean(seeds_per_stem))
## this block of code would take the mean seeds per stem of trifolium in each specific background-ppt treatment combo

mean <- dat %>%
  group_by(background, ppt_trt, nut_trt) %>% 
  summarize(b_ppt_nut = mean(seeds_per_stem))

## Nat Exploration A
## want to visualize seed output by nutrient treatment
## how does this differ by competition background?

## seeds per stem vs. nutrient treatment; faceted by competition background
ggplot(dat, aes(x=nut_trt, y=seeds_per_stem)) +
  geom_boxplot() +
  facet_wrap(~background)

## seeds per stem vs. competition background; faceted by nutrient treatment
ggplot(dat, aes(x=background, y=seeds_per_stem)) +
  geom_boxplot() +
  facet_wrap(~nut_trt)

nut_trt <- aov(seeds_per_stem~nut_trt, data = dat)
summary(nut_trt)
summary(nut_trt)[[1]][["Pr(>F)"]][1] 
tests <- summary(nut_trt) 
tests #easiest way to see p-value


## Nat Exploration B
## want to visualize total mass by precipitation treatment
## how does this differ by competition background?

## total mass vs. precipitation treatment; faceted by competition background
ggplot(dat, aes(x=ppt_trt, y=tot_mass)) +
  geom_boxplot() +
  facet_wrap(~background)

#total mass (and variance) seemed to be higher with irrigation (when comparing precip.treatments within one background)

## total mass vs. competition background; faceted by precipitation treatment
ggplot(dat, aes(x=background, y=tot_mass)) +
  geom_boxplot() +
  facet_wrap(~ppt_trt)

#total mass averages overall seemed to be higher in control precip.treatment 

## Nat Exploration C
## want to visualize total mass by nutrient treatment
## how does this differ by competition background?

## total mass vs. nutrient treatment; faceted by competition background
ggplot(dat, aes(x=nut_trt, y=tot_mass)) +
  geom_boxplot() +
  facet_wrap(~background)

## total mass vs. competition background; faceted by nutrient treatment
ggplot(dat, aes(x=background, y=tot_mass)) +
  geom_boxplot() +
  facet_wrap(~nut_trt)

## Nat Exploration D
## how do I compare seeds_per_stem to total mass? comparison of control vs. intraspecific (Trif)
## per precip.treatment
## per nutrient.treatment

### x-axis (1) Trif  vs. Control, (2) treatment (precip or nutrient), y-axis 0-2, legend for seed_per_stem (blue) and total mass (red) 
(stacked bar chart)

#example code: https://stackoverflow.com/questions/53454787/comparing-multiple-categorical-variables-in-r

# find if there is relationship between the variables (seed_per_stem vs. total mass); scatterplot
ggplot(dat, aes(x=tot_mass, y=seeds_per_stem, color = ppt_trt))+
  geom_point() +
  facet_wrap(~background)

#comparison of background 
ggplot(dat, aes(x=tot_mass, y=seeds_per_stem, color = background))+
  geom_point()+
  geom_smooth(method="lm")
#this is testing whether total mass is affecting seed per stem...background nor treatment is affecting relationships

ggplot(dat, aes(x=tot_mass, y=seeds_per_stem, color = ppt_trt))+
  geom_point()+ #when total mass increases, seeds per stem increase = positive correlation 
  geom_smooth(method="lm")
#similar relationship for all precip treatments 

ggplot(dat, aes(x=tot_mass, y=seeds_per_stem))+
  geom_point()+ #when total mass increases, seeds per stem increase = positive correlation 
  geom_smooth(method="lm") #shadowed area is estimate of error

# linear model for significance? 
seeds_totmass <- lm(seeds_per_stem~tot_mass, data =dat)
summary(seeds_totmass) 

# more helpful to used stacked bar chart to compare seed mass and stem mass as a percent of total mass


## Nat Exploration E
#use data simulation to compare the data I collected to simulated data 
#instead, compare to other NutNet sites?

##Step 1: simulate seeds and stems of TRIF based on treatment (only Control, no background)
##Step 2: simulate seeds and stems of TRIF based on treatment (in addition to background) 


tot_seed_mean <- mean(dat$tot_seeds,na.rm=TRUE)
tot_seed_sd <- sd(dat$tot_seeds,na.rm=TRUE)

tot_stem_mean <- mean(dat$tot_stems,na.rm=TRUE) #issues with -Inf in raw data?
tot_stem_sd <- sd(dat$tot_stems,na.rm=TRUE) #could probably set to 1 because we had min of 1 max of 3

Seeds <- rnorm(252, tot_seed_mean, tot_seed_sd)
Stems <- rnorm(252, tot_stem_mean, tot_stem_mean)

#I lost the reference source for this part...
simulated_dat <- data.frame(
  sub_condition = rep( c("Seeds", "Stems"), c(A_sub_n, B_sub_n) ),
  score = c(Seeds, Stems)
)

simulated_dat <- simulated_dat %>%
  mutate(Seeds_Per_Stem = Seeds/Stems)

#then compare seeds per stem of dat and simulated_dat using scatterplot? 

