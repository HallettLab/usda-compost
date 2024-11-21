# use hmsc to include nested error in hierarchical bayes model

# script purpose:
# native recruitment species comp is super patchy across blocks and plots. too patchy to model response in glmm, 
# > also too few obs per trt combo to model nested error.. bayes is an option to overcome these limitations of glmm approach
# for ESA 2023, ctw used gjam (clark et al. 2017) .. however gjam does not allow multiple random error terms
# compost has nested exp design, meaning precip trt subplots within soil plots are not independent.. and these are both nested in a block.

# possible problems:
# hmsc takes longer to run than gjam.. we shall see if it's an improvement or not.
# not sure if hmsc can accomodate data types as gjam can (e.g., fractional cover)

# script steps:
# read in prepped seedtrt dataset from 2021 and trt key
# prep (aggregate data collected over two months, possibly transform cover to count)
# determine best distribution to use
# possibly run simple glmer with nested structure to build expectations (also have gjam results for this)
# find a good model on just the species cover data
# incorporate traits: does it improve the model? (does bringing in trait info explain shed light on what patchy species id data can't?)
# > if traits do help, is this corroborated with another approach like ordination?


# -- SETUP -----
library(tidyverse)
library(Hmsc)
library(coda)
library(fitdistrplus)

# modify default settings
options(stringsAsFactors = F)
theme_set(theme_test())
na_vals <- c("" , " ","NA", NA)

# specify dropbox pathway
datpath <- "~/Dropbox/USDA-compost/Data/"
# specify path for writing out models and figs
#outpath <- "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/chapter/"

# read in cleaned native data
natlong <- read.csv(paste0(datpath, "Native/Native_CleanedData/Compost_Native_LongClean.csv"), strip.white = T, na.strings = na_vals)
natwide <- read.csv(paste0(datpath, "Native/Native_CleanedData/Compost_Native_WideClean.csv"), strip.white = T, na.strings = na_vals)
# read in treatment key
trtkey <- read.csv(paste0(datpath, "Compost_TreatmentKey.csv"),na.strings = na_vals, strip.white = T)
# read in master spp list (lookup table)
spplist <- read.csv(paste0(datpath, "Compost_SppList.csv"), na.strings = na_vals, strip.white = T)
# read in cleaned cover data for background comparison
coverlong <- read.csv(paste0(datpath, "Cover/Cover_CleanedData/Compost_Cover_LongClean.csv"), strip.white = T, na.strings = na_vals)

# check that all read in as expected
str(natlong)
str(natwide)
str(coverlong)
#fix date
natlong$date <- as.Date(natlong$date, format = "%Y-%m-%d")
natwide$date <- as.Date(natwide$date, format = "%Y-%m-%d")
coverlong$date <- as.Date(coverlong$date, format = "%Y-%m-%d")


# aggregate data over sampling season
# choose max pct_cov per species per plot per herbicide and seeding trt, unless grass (use mean), then sum fxnl groups
natcoarse_all <- natlong %>%
  dplyr::select(plot:mon, code4:pct_cover) %>%
  spread(code4, pct_cover, fill = 0) %>%
  # gather any spp codes [codenames won't have lowcase character]
  gather(code4, pct_cover, names(.)[!grepl("[a-z]{2,}", names(.))]) %>%
  left_join(distinct(natlong, code4, fxnl_grp, nativity, native_seeded)) %>%
  mutate(nativity = ifelse(nativity == "Unknown", "Exotic", nativity),
         native_seeded = ifelse(native_seeded == "No", "Background", "Seeded"),
         coarse_fxnl = paste(native_seeded, nativity, fxnl_grp,sep = " ")) %>%
  # by not including month as grouping factor, take max cover val of species over season
  group_by(plot, seedtrt, herbicide, fulltrt, block, ppt_trt, nut_trt, fxnl_grp, nativity, native_seeded, coarse_fxnl, code4) %>%
  reframe(nobs = length(mon),
          cover = ifelse(nobs == 2 & coarse_fxnl == "Background Exotic Grass", 
                         mean(pct_cover, na.rm = T), 
                         #max(pct_cover, na.rm = T), 
                         max(pct_cover, na.rm = T))) %>%
  distinct() %>%
  # sum by fxnl grp
  group_by(plot, seedtrt, herbicide, fulltrt, block, ppt_trt, nut_trt, fxnl_grp, nativity, native_seeded, coarse_fxnl) %>%
  reframe(totcov = sum(cover),
          spp = str_flatten(unique(code4[cover > 0]), collapse = ", "),
          S = length(unique(code4[cover > 0]))) %>%
  # change nut_trt control to XC
  mutate(nut_trt = dplyr::recode(nut_trt, N = "XC"),
         fulltrt = gsub("N", "XC", fulltrt))

# repeat for background unherbicided (i.e. main spp comp)
# > just for 2021
covcoarse <- subset(coverlong, yr == 2021) %>%
  dplyr::select(plot:sample_event, code4:pct_cover) %>%
  spread(code4, pct_cover, fill = 0) %>%
  # gather any spp codes [codenames won't have lowcase character]
  gather(code4, pct_cover, names(.)[!grepl("[a-z]{2,}", names(.))]) %>%
  left_join(distinct(coverlong, code4, fxnl_grp, nativity)) %>%
  group_by(plot, plotid, fulltrt, block, nut_trt, ppt_trt, yr, fxnl_grp, nativity, code4) %>%
  # take max cover for species between april and may
  reframe(cover = ifelse(grepl("Gra", fxnl_grp) & grepl("Exo|Unk", nativity), 
                         mean(pct_cover, na.rm = T), 
                         #max(pct_cover, na.rm = T), 
                         max(pct_cover, na.rm = T))) %>%
  distinct() %>%
  mutate(nativity = ifelse(nativity == "Unknown", "Exotic", nativity),
         native_seeded = "Background",
         seedtrt = "Unseeded",
         herbicide = "Non-herbicided",
         coarse_fxnl = paste(native_seeded, nativity, fxnl_grp,sep = " ")) %>%
  # by not including month as grouping factor, summing apr and may vals
  group_by(plot, yr, seedtrt, herbicide, fulltrt, block, ppt_trt, nut_trt, fxnl_grp, nativity, native_seeded, coarse_fxnl) %>%
  reframe(totcov = sum(cover),
          spp = str_flatten(unique(code4[cover > 0]), collapse = ", "),
          S = length(unique(code4[cover > 0]))) %>%
  ungroup() %>%
  # change nut_trt control to XC
  mutate(nut_trt = dplyr::recode(nut_trt, N = "XC"),
         fulltrt = gsub("N", "XC", fulltrt))

# create dataset with everything
allcov <- subset(covcoarse, yr == 2021, select = -yr) %>%
  rbind(natcoarse_all[names(.)])

# note: plot 33 not seeded because not herbicided (block 4, NW plot)


# -- DETERMINE BEST DISTRIBUTION -----
summary(natcoarse_all$totcov[grepl("Nat", natcoarse_all$seedtrt)]) # between 0 and 100
plotdist(natcoarse_all$totcov[grepl("Nat", natcoarse_all$seedtrt)], demp = T)
descdist(natcoarse_all$totcov[grepl("Nat", natcoarse_all$seedtrt)], boot = 10000)
descdist(natcoarse_ydata$IntroForb, boot = 10000)
descdist(natcoarse_ydata$IntroGrass, boot = 10000)
descdist(natcoarse_ydata$IntroNfixer, boot = 10000)
descdist(natcoarse_ydata$NatForb, boot = 10000)
descdist(natcoarse_ydata$SeededNatForb, boot = 10000)
descdist(natcoarse_ydata$SeededNatGrass, boot = 10000)
descdist(round(seedtrt_ydata$AVBA), boot = 10000)
descdist(round(seedtrt_ydata$ERBO), boot = 10000)
descdist(seedtrt_ydata$FEMI, boot = 10000)
fw <- fitdist(natcoarse_ydata$IntroForb[grepl("Nat", natcoarse_all$seedtrt)], method = "mle", "gamma")
summary(fw)
denscomp(list(fw), legendtext = "gamma")

# prep dataset with all neighbor abundance, and with neighbor forbs and grams split out
# can start with simplecov background exotics, sum across and retain functional groups
neighborhood <- subset(simplecov, grepl("Nat", seedtrt) & grepl("Back", fullgrp), select = -fullgrp) %>%
  mutate(fxnlgrp = casefold(fxnlgrp)) %>%
  spread(fxnlgrp, totcov, fill = 0) %>%
  mutate(neighbors = forb + grass)
nats <- c("TRCI", "BRCA", "FEMI", "NEMA", "ESCA")
seededspp <- dplyr::select(natwide, plot:date, all_of(nats)) %>%
  subset(!grepl("Unse", seedtrt)) %>%
  dplyr::select(plot, herbicide, mon, date, TRCI:ESCA) %>%
  gather(spp, cov, TRCI:ESCA) %>%
  group_by(plot, herbicide, spp) %>%
  reframe(cov = max(cov)) %>%
  spread(spp, cov)

# join to neighborhood
neighborhood <- merge(neighborhood, seededspp)
# neighhors have continuous cover, native seeded rounded for count data
neighborhood_counts <- mutate_at(neighborhood, .vars = nats, function(x) ifelse(x %in% c(1:9), x+1, ifelse(x > 0 & x<1, 1, x)))
# presence absence dataset
neighborhood_PA <- mutate_at(neighborhood, .vars = nats, function(x) ifelse(x>0, 1, 0))


# make datasets with unseeded plots as well
neighborhood_allplots <- subset(simplecov, grepl("Back", fullgrp), select = -fullgrp) %>%
  mutate(fxnlgrp = casefold(fxnlgrp)) %>%
  spread(fxnlgrp, totcov, fill = 0) %>%
  mutate(neighbors = forb + grass)
seededspp_allplots <- dplyr::select(natwide, plot:date, all_of(nats)) %>%
  dplyr::select(plot, seedtrt, herbicide, mon, date, TRCI:ESCA) %>%
  gather(spp, cov, TRCI:ESCA) %>%
  group_by(plot, seedtrt, herbicide, spp) %>%
  reframe(cov = max(cov)) %>%
  spread(spp, cov)
neighborhood_allplots <- left_join(neighborhood_allplots, seededspp_allplots) %>%
  # infill any NAs for nat cols with 0
  mutate_at(all_of(nats), function(x) ifelse(is.na(x), 0, x))


# -- DETERMINE DISTRIBUTIONS -----



