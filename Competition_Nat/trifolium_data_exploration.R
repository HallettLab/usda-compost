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

## Tukey's Honest significant difference test; I like to use this when I have a lot of categories
## you can see the differences between each pair of categories
TukeyHSD(ppt)


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

nut_trt <- aov(seeds_per_stem~ppt_trt, data = dat)
summary(nut_trt)
summary(nut_trt)[[1]][["Pr(>F)"]][1] #It's the same p-value as ppt? 
tests <- summary(nut_trt) #Alternative way to see p-value...it's still the same  
str(tests)

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
## how do I compare seeds_per_stem to total mass? 

