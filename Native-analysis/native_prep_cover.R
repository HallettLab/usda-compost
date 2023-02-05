# aggregate native datasets

# purpose:
# source script to run for cover by species and coarse functional groups
# preps species list for pairing (reduce to only species noted in native experiment)
# compiles trait data available for as many sub-experiment species as possible
# > trait data from: Brad (GH), ongoing Hallett-Larios, Julie (CU trait screening for ClimVar + Chapman), Lina + Ashley

# notes:
# plot 33 not seeded because not herbicided (block 4, NW plot)


# -- SETUP ------
# load needed libraries
library(tidyverse)

# modify default settings
options(stringsAsFactors = F)
theme_set(theme_test())
na_vals <- c("" , " ","NA", NA)


# path to compost dat specifically
datpath <- "~/Dropbox/USDA-compost/Data/"

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
#str(natlong); str(natwide); str(coverlong)

#fix date
natlong$date <- as.Date(natlong$date, format = "%Y-%m-%d")
natwide$date <- as.Date(natwide$date, format = "%Y-%m-%d")
coverlong$date <- as.Date(coverlong$date, format = "%Y-%m-%d")


  
# -- PREP COMPOST TRT KEY FOR NESTED ANALYSIS -----
# add whole plot ID to trtkey (block-nut_trt ID for nesting)
trtkey <- trtkey %>%
  # correct nut_trt control to XC
  mutate_at(.vars = c("nut_trt", "fulltrt", "plotid"), function(x) gsub("N", "XC", x)) %>%
  mutate(wholeplotID = paste0(block, nut_trt)) %>% 
  rename(subplotID = plotid) %>% # plotid = subplotID
  # rearrange cols
  subset(select = c(block, nut_trt, ppt_trt, fulltrt, wholeplotID, subplotID, plot)) %>%
  # add x,y order within blocks just in case
  mutate(ypos = rep(1:3, 12), # uphill to downhill
         # east to west
         xpos = rep(rep(1:3, each = 3), 4)) %>%
  # set factor levels
  mutate(nut_trt = factor(nut_trt, levels = c("XC", "F", "C")),
         ppt_trt = factor(ppt_trt, levels = c("XC", "D", "W")))
  


# -- AGGREGATE COVER -----
# 1. COARSE FXNL, ALL SUB-SUB-PLOTS -----
# ashley only updated cover vals for re-surveys for spp that had changed
# choose max pct_cov per species per plot per herbicide and seeding trt, then sum
coarsecov_nat <- natlong %>%
  mutate(nativity = ifelse(nativity == "Unknown", "Exotic", nativity),
         native_seeded = ifelse(native_seeded == "No", "Background", "Seeded"),
         coarse_fxnl = paste(native_seeded, nativity, fxnl_grp,sep = " ")) %>%
  # by not including month as grouping factor, summing apr and may vals
  group_by(plot, seedtrt, herbicide,fxnl_grp, nativity, native_seeded, code4) %>%
  mutate(maxcov = max(pct_cover)) %>%
  ungroup() %>%
  group_by(plot, seedtrt, herbicide, fxnl_grp, nativity, native_seeded, coarse_fxnl) %>%
  summarise(totcov = sum(pct_cover),
            spp = str_flatten(unique(code4), collapse = ", "),
            S = length(unique(code4))) %>%
  ungroup() %>%
  data.frame()

# repeat for background unherbicided (i.e. main spp comp)
# > for quick prelim, take average between apr and may each year
coarsecov_main_allyrs <- coverlong %>%
  dplyr::select(plot:sample_event, code4:pct_cover) %>%
  spread(code4, pct_cover, fill = 0) %>%
  gather(code4, pct_cover, AGSP:ncol(.)) %>%
  group_by(plot, yr, code4) %>%
  summarise(max_cover = max(pct_cover)) %>%
  ungroup() %>%
  left_join(distinct(spplist[c("code4", "fxnl_grp", "nativity")])) %>%
  mutate(nativity = ifelse(nativity == "Unknown", "Exotic", nativity),
         native_seeded = "Background",
         seedtrt = "Unseeded",
         herbicide = "Non-herbicided",
         coarse_fxnl = paste(native_seeded, nativity, fxnl_grp,sep = " ")) %>%
  # by not including month as grouping factor, summing apr and may vals
  group_by(plot, yr, seedtrt, herbicide, fxnl_grp, nativity, native_seeded, coarse_fxnl) %>%
  summarise(totcov = sum(max_cover),
            spp = str_flatten(unique(code4), collapse = ", "),
            S = length(unique(code4))) %>%
  ungroup() %>%
  data.frame()

# create dataset with everything
coarsecov_all_2021 <- subset(coarsecov_main_allyrs, yr == 2021, select = -yr) %>%
  rbind(coarsecov_nat[names(.)]) %>%
  # join additional trtkey nestedness info
  left_join(trtkey) %>%
  # set factor levels for seedtrt and herbicide
  mutate(seedtrt = factor(seedtrt, levels = c("Unseeded", "Native seeded")),
         herbicide = factor(herbicide, levels = c("Non-herbicided", "Herbicided"))
         ) 
                        

# 2. SIMPLIFY COARSE FXNL COV -----
# drop any background native plant that recruited.. just want to know if herbicide advantaged seeded natives against ambient non-natives
coarsecov_all_2021_simple <- dplyr::select(coarsecov_all_2021, -c(spp, S)) %>%
  #combine N-fixer and forb
  mutate(fxnl_grp2 = ifelse(grepl("Forb|N-fix", fxnl_grp), "Forb", fxnl_grp)) %>%
  group_by(plot, seedtrt, herbicide, fxnl_grp2, nativity, native_seeded) %>%
  summarise(totcov = sum(totcov), .groups = "drop_last") %>%
  # drop any background native
  subset(!(native_seeded == "Background" & nativity == "Native")) %>%
  unite(fullgrp, native_seeded, nativity, fxnl_grp2, sep = "_") %>%
  # make sure all possible groups repped
  spread(fullgrp, totcov, fill = 0) %>%
  # also drop plot 33 since no experiment there
  subset(plot != 33) %>%
  gather(fullgrp, totcov, Background_Exotic_Forb:ncol(.)) %>%
  # make simple fxnl grp
  mutate(fxnlgrp =  ifelse(grepl("Forb", fullgrp), "Forb", "Grass"),
         fxnlgrp = factor(fxnlgrp, levels = c("Forb", "Grass"))) %>%
  # rejoin trtment info
  left_join(trtkey) %>%
  data.frame()


# 3. PREP NEIGHBORHOOD DATA -----
# individual seeded spp in their own columns with sums of neighbors (grams, forbs, and alltog)
# prep dataset with all neighbor abundance, and with neighbor forbs and grams split out
# can start with simplecov background exotics, sum across and retain functional groups
neighbors <- subset(coarsecov_nat, grepl("Nat", seedtrt) & grepl("Back", native_seeded), select = c(plot:herbicide, coarse_fxnl, totcov)) %>%
  # simplify coarse_fxnl names
  mutate(coarse_fxnl = gsub ("Background| |-", "", coarse_fxnl)) %>%
  spread(coarse_fxnl, totcov, fill = 0) %>%
  mutate(ExoForbNfix = ExoticForb + ExoticNfixer,
         NatForbNfix = NativeForb + NativeNfixer,
         ExoticNeighbors = ExoticForb + ExoticGrass + ExoticNfixer,
         NativeNeighbors = NativeForb + NativeGrass + NativeNfixer,
         AllGrass = ExoticGrass + NativeGrass,
         AllForbNfix = ExoForbNfix + NatForbNfix,
         AllNeighbors = ExoticNeighbors + NativeNeighbors)

# list natives seeded
nats <- c("TRCI", "BRCA", "FEMI", "NEMA", "ESCA")
seededspp <- dplyr::select(natwide, plot:date, all_of(nats)) %>%
  subset(!grepl("Unse", seedtrt)) %>%
  dplyr::select(plot, herbicide, mon, date, TRCI:ESCA) %>%
  gather(spp, cov, TRCI:ESCA) %>%
  group_by(plot, herbicide, spp) %>%
  summarise(cov = max(cov), .groups = "drop_last") %>%
  spread(spp, cov)

# join to neighborhood
neighborhood <- merge(neighbors, seededspp) %>%
  #join trtkey
  left_join(trtkey) %>%
  #rearrange
  subset(select = c(plot, block:ypos, herbicide:TRCI)) %>%
  arrange(plot, herbicide)
# transform continuous abundance to count data for native only
##  names(neighbors)[!grepl("plot|seed|herb", names(neighbors))] <- if want to transform neighbors later
neighborhood_counts <- mutate_at(neighborhood, .vars = nats, function(x) ifelse(x %in% c(1:9), x+1, ifelse(x > 0 & x<1, 1, x)))
neighborhood_PA <- mutate_at(neighborhood, .vars = nats, function(x) ifelse(x>0, 1, 0))


