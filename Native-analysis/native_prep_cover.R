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
# str(natlong); str(natwide); str(coverlong)

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
         ppt_trt = factor(ppt_trt, levels = c("XC", "D", "W")),
         fulltrt = factor(fulltrt, levels = c("XCXC", "XCD", "XCW",
                                              "FXC", "FD", "FW",
                                              "CXC", "CD", "CW")),
         wholeplotID = factor(wholeplotID, 
                              levels = paste0(rep.int(c(1:4),3),rep(levels(nut_trt), each = 4))
                              ),
         subplotID = factor(subplotID,
                            levels = paste0(rep.int(c(1:4),4),rep(levels(fulltrt), each = 4))
                            ),
         block = factor(block))
# review
summary(trtkey)
levels(trtkey$wholeplotID)
levels(trtkey$subplotID) # good


# -- QA: SCREEN GRASS COVER CHANGES -----
# general note: sometimes May revisits reveal April grasses were mis-ID'd (e.g., if they didn't have infloresence)
# native seeded surveys: ashley only updated cover vals for re-surveys for spp that had changed
sort(unique(natlong$notes))
# focus comment search based on skim
with(natlong, sort(unique(notes[grepl("correct|update ", notes)])))
# screen for seeded spp downhill of plot (present, but seed migrated downslope)
# > plots 11 herb, 12 non herb, 2 [update AICA]
with(natlong, sort(unique(notes[grepl("down|outside", notes)])))
# look at entered cover for seeded spp outside plot
# > in manual review, anytime FEMI is noted outside plot it is also recorded with cover inside plot
# > BRCA noted outside plot is not present in one plot (assign trace cover) [plot 34, block4 CD, non-herbicided]

with(subset(coverlong, yr == 2021), sort(unique(notes)))
# ^ no notes in main cover dataset about mis-ID'd grasses

# rules for aggregating data across apr and may surveys: 
# 1. take max value of a species from both surveys UNLESS says to update (applies to 3 subplots)
# 2. if seeded spp found outside plot and not recorded as present in plot, assign trace value: it established, seed just migrated
# > if a plant was present in april and died by may, still count as it being able to grow: we are interested in whether treatments inhibit ability to establish,



# - 1. SPECIES LEVEL AGGREGATE COVER -----
# compile season cover per species dataset for all native seeded, and unseeded plots in herb/non-herb 2021
natsp <- subset(natlong, select = c(plot:mon, notes, code4, pct_cover, fxnl_grp:native_seeded)) %>%
  group_by(plot, seedtrt, herbicide, code4, fxnl_grp) %>%
  mutate(usemay = any(grepl("correct/update|update AICA",notes) & fxnl_grp == "Grass")) %>%
  ungroup() %>%
  reframe(pct_cover = ifelse(usemay, pct_cover[mon == 5], max(pct_cover)),.by = names(.)[!grepl("cover|usemay|notes|mon", names(.))]) %>%
  distinct()

# add row for BRCA in plot where outside frame and not recorded in
brca_row <- subset(natlong, grepl("BRCA outside", notes) & fxnl_grp == "Grass", select = names(natsp))[1,] %>%
  mutate(code4 = "BRCA", pct_cover = 0.01, nativity = "Native", native_seeded = "Yes")

# bind brca row to nat sp
natsp <- rbind(natsp, brca_row) %>%
  arrange(plot, seedtrt, herbicide, code4) %>%
  # set seedtrt and herbicide as factors
  mutate(seedtrt = factor(seedtrt, levels = c("Unseeded", "Native seeded")),
         herbicide = factor(herbicide, levels = c("Non-herbicided", "Herbicided")))

# be sure there aren't any duplicate entries
summary(duplicated(grouped_df(natsp, names(natsp)[!names(natsp) == "pct_cover"]))) # no duplicates
# review cover values
with(natsp, sort(unique(pct_cover))) # 1.01 looks like typo; some values above 1 have a 0.5 addition
natsp$pct_cover[natsp$pct_cover == 1.01] <- 1
# check all
summary(natsp) # ok

# prep ambient cover -- just select 2021 for comparison
nats <- unique(natsp$code4[natsp$native_seeded == "Yes"])
ambsp <- mutate(coverlong, native_seeded = ifelse(code4 %in% nats, "Yes", "No"),
                herbicide = "Non-herbicided", seedtrt = "Unseeded",
                # seed factor levels for row-binding
                seedtrt = factor(seedtrt, levels = c("Unseeded", "Native seeded")),
                herbicide = factor(herbicide, levels = c("Non-herbicided", "Herbicided"))) %>%
  subset(yr == 2021, select = c(plot:date, notes, code4, pct_cover, fxnl_grp, nativity, native_seeded, herbicide, seedtrt)) %>%
  reframe(pct_cover = max(pct_cover),.by = names(natsp)[!grepl("cover|usemay|notes|mon", names(natsp))]) %>%
  distinct()
# review cover vals
with(ambsp, sort(unique(pct_cover))) # looks okay
summary(ambsp)

# put both together, add column for count cover data
natambsp_cov <- rbind(natsp,ambsp) %>%
  arrange(plot, seedtrt, herbicide, code4) %>%
  mutate(nativity = ifelse(nativity == "Unknown", "Exotic", nativity),
         native_seeded = ifelse(native_seeded == "No", "Background", "Seeded"),
         coarse_fxnl = paste(native_seeded, nativity, fxnl_grp,sep = " "),
         # trace gets 1, all else gets rounded; add 0.1 sound 0.5 rounds (R doesn't round it always)
         count_cover = ifelse(pct_cover == 0.01, 1, round(pct_cover + 1.1))) %>%
  # drop trt cols to join trt key
  subset(select = -c(plotid:ppt_trt)) %>%
  left_join(trtkey) %>%
  # rearrance cols
  subset(select = c(plot, wholeplotID, subplotID, block:fulltrt, ypos, xpos, seedtrt, herbicide, yr, code4, pct_cover, count_cover, fxnl_grp:native_seeded, coarse_fxnl)) %>%
  data.frame()
  
# make sure no duplicates
summary(duplicated(grouped_df(natambsp_cov, names(natambsp_cov)[!grepl("cover", names(natambsp_cov))]))) # no duplicates
summary(natambsp_cov)
with(distinct(natambsp_cov, plot, fulltrt, seedtrt, herbicide), table(fulltrt[seedtrt == "Native seeded"])) # plot 33 does not have herbicided subplots



# - 2. FXNL GRP AGGREGATE COVER -----
# sum on coarse fxnl grp, make pct cov and 'count' cover
# need to calculate count differently otherwise some groups get inflated (e.g., 14.03 -> 22)
natambfxnl_cov <- grouped_df(natambsp_cov, names(natambsp_cov)[!grepl("cover|code4", names(natambsp_cov))]) %>%
  reframe(totcov_pct = sum(pct_cover),
          # insert placeholder col for count cov
          totcov_count = sum(pct_cover),
          spp = str_flatten(unique(code4), collapse = ", "),
          S = length(unique(code4))) %>%
  distinct() %>%
  #make cover count data (trace amounts --> 1)
  mutate(integer = floor(totcov_pct),
    decimals = str_extract(totcov_pct, "(?<=[.])[0-9]+"),
    # replace na in decimals with 0, and change .5 to .50
    decimals = ifelse(is.na(decimals), "00", 
                      ifelse(decimals == "5", "50", decimals)),
         tenths = as.numeric(substr(decimals, 1, 1)),
         hundredths = as.numeric(substr(decimals, 2, 2)),
         addct = ifelse(tenths == 5 & !is.na(hundredths), 1+ hundredths, 
                        ifelse(tenths == 5 & is.na(hundredths), 1, hundredths)),
    totcov_count = ifelse(is.na(addct), integer, integer + addct)) %>%
  # drop cols used to make count
  subset(select = -c(integer:addct)) %>%
  data.frame()



# - 3. WIDE NEIGHBORHOOD DATA -----
# individual seeded spp in their own columns with sums of neighbors (grams, forbs, and alltog)
# prep dataset with all neighbor abundance, and with neighbor forbs and grams split out
# can start with simplecov background exotics, sum across and retain functional groups
widefxnl_pct <- subset(natambfxnl_cov, select = -c(fxnl_grp:native_seeded, totcov_count:S)) %>%
  # simplify coarse_fxnl names
  mutate(coarse_fxnl = gsub ("Background| |-", "", coarse_fxnl)) %>%
  spread(coarse_fxnl, totcov_pct, fill = 0) %>%
  mutate(ExoForbNfix = ExoticForb + ExoticNfixer,
         NatForbNfix = NativeForb + NativeNfixer,
         ExoticNeighbors = ExoticForb + ExoticGrass + ExoticNfixer,
         NativeNeighbors = NativeForb + NativeGrass + NativeNfixer,
         AllGrass = ExoticGrass + NativeGrass,
         AllForbNfix = ExoForbNfix + NatForbNfix,
         AllNeighbors = ExoticNeighbors + NativeNeighbors)

# repeat for fxnl group count
widefxnl_count <- subset(natambfxnl_cov, select = -c(fxnl_grp:native_seeded, totcov_pct, spp:S)) %>%
  # simplify coarse_fxnl names
  mutate(coarse_fxnl = gsub ("Background| |-", "", coarse_fxnl)) %>%
  spread(coarse_fxnl, totcov_count, fill = 0) %>%
  mutate(ExoForbNfix = ExoticForb + ExoticNfixer,
         NatForbNfix = NativeForb + NativeNfixer,
         ExoticNeighbors = ExoticForb + ExoticGrass + ExoticNfixer,
         NativeNeighbors = NativeForb + NativeGrass + NativeNfixer,
         AllGrass = ExoticGrass + NativeGrass,
         AllForbNfix = ExoForbNfix + NatForbNfix,
         AllNeighbors = ExoticNeighbors + NativeNeighbors)

# repeat for species level
widesp_pct <- subset(natambsp_cov, select = -c(fxnl_grp:native_seeded, count_cover,fxnl_grp:coarse_fxnl)) %>%
  spread(code4, pct_cover, fill = 0)
widesp_count <- subset(natambsp_cov, select = -c(fxnl_grp:native_seeded, pct_cover,fxnl_grp:coarse_fxnl)) %>%
  spread(code4, count_cover, fill = 0)
# review
summary(widesp_pct) # ok


# create presence-absence forms
widefxnl_pa <- widefxnl_count %>%
  mutate_at(names(.)[grepl("Exo|Nat|All", names(.))], function(x) ifelse(x>0, 1, x))
# repeat for species level 
widesp_pa <- widesp_count %>%
  mutate_at(names(.)[grepl("^[A-Z]+", names(.))], function(x) ifelse(x>0, 1, x))
# review
summary(widesp_pa) # ok


