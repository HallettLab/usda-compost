# native recruitment GJAM analysis
# authors(s): ctw (caitlin.t.white@colorado.edu)

# made in prep for ESA abstract, to submit Feb 2023

# motivation:
# previously have tried mixed models, however not enough data points to model all species while accounting for treatments and controlling for quirks of each block
# Fall 2022 use gjamTime for KNS on a time-series monitoring dataset. That dataset was too messy for GJAM, *however* learned GJAM might be a good candidate for native recruitment dataset

# script purpose:
# read in all cleaned and prepped dat for native analyses
# specify gjam model
# run gjam model with species x treatment info only
# perhaps try bringing in GH trait data (it's either that or compare with climvar natdiv recruitment experiment)

# basics to address:
# were ppt and nut treatments effective? (can point to AS ms if out)
# was herbicide effective in reducing background flora?
# what was climate of the experiment year?
# did seeded species come up? 

# hypotheses to test:
# compost + watering may advantage exotics (all/any) over native flora (all/any)
# native forbs may persist with native or non-native grams, non-native grams should compete directly with native grams (non-native gram abundance will have strongest neg effect on native gram abundance)
# forbs generally, but native forbs in particular, may do relatively better (to control or wet plots) in drought conditions
# > if go trait route, conservative traits will enable persistence in unfertilized drought conditions.. exotic forbs may do relatively best in amended drought plots, exotic grams (or spp with aquisitive traits) in amended irrigated plots

# notes: 
# plot 33 not seeded because not herbicided (block 4, NW plot)



# -- SETUP -----
# load needed libraries
library(tidyverse)
library(cowplot)
# gjam packages and functions
library(gjam)
# for gjamTime
library(devtools)
# source supplemental functions for gjamTime
d <- "https://github.com/jimclarkatduke/gjam/blob/master/gjamTimeFunctions.R?raw=True"
source_url(d)

# modify default settings
options(stringsAsFactors = F)
theme_set(theme_test())
na_vals <- c("" , " ","NA", NA)

# specify dropbox pathway
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
str(natlong)
str(natwide)
str(coverlong)
#fix date
natlong$date <- as.Date(natlong$date, format = "%Y-%m-%d")
natwide$date <- as.Date(natwide$date, format = "%Y-%m-%d")
coverlong$date <- as.Date(coverlong$date, format = "%Y-%m-%d")




# -- DATA PREP ----
# want to know..
# 1) how did seeded natives do?
# 2) is there any difference between native cover in herbicide vs. not 
# 3) Did seeding natives matter at all for changing exotic cover (e.g. suppressing non-native cover?)

# ashley only updated cover vals for re-surveys for spp that had changed
# choose max pct_cov per species per plot per herbicide and seeding trt, then sum
natcoarse_all <- natlong %>%
  mutate(nativity = ifelse(nativity == "Unknown", "Exotic", nativity),
         native_seeded = ifelse(native_seeded == "No", "Background", "Seeded"),
         coarse_fxnl = paste(native_seeded, nativity, fxnl_grp,sep = " ")) %>%
  # by not including month as grouping factor, summing apr and may vals
  group_by(plot, seedtrt, herbicide, fulltrt, block, ppt_trt, nut_trt, fxnl_grp, nativity, native_seeded, code4) %>%
  mutate(maxcov = max(pct_cover)) %>%
  ungroup() %>%
  group_by(plot, seedtrt, herbicide, fulltrt, block, ppt_trt, nut_trt, fxnl_grp, nativity, native_seeded, coarse_fxnl) %>%
  summarise(totcov = sum(pct_cover),
            spp = str_flatten(unique(code4), collapse = ", "),
            S = length(unique(code4))) %>%
  ungroup() %>%
  # change nut_trt control to XC
  mutate(nut_trt = dplyr::recode(nut_trt, N = "XC"),
         fulltrt = gsub("N", "XC", fulltrt))

# repeat for background unherbicided (i.e. main spp comp)
# > for quick prelim, take average between apr and may each year
covcoarse <- coverlong %>%
  dplyr::select(plot:sample_event, code4:pct_cover) %>%
  spread(code4, pct_cover, fill = 0) %>%
  gather(code4, pct_cover, AGSP:ncol(.)) %>%
  group_by(plot, plotid, fulltrt, block, nut_trt, ppt_trt, yr, code4) %>%
  summarise(mean_cover = max(pct_cover)) %>%
  ungroup() %>%
  left_join(distinct(spplist[c("code4", "fxnl_grp", "nativity")])) %>%
  mutate(nativity = ifelse(nativity == "Unknown", "Exotic", nativity),
         native_seeded = "Background",
         seedtrt = "Unseeded",
         herbicide = "Non-herbicided",
         coarse_fxnl = paste(native_seeded, nativity, fxnl_grp,sep = " ")) %>%
  # by not including month as grouping factor, summing apr and may vals
  group_by(plot, yr, seedtrt, herbicide, fulltrt, block, ppt_trt, nut_trt, fxnl_grp, nativity, native_seeded, coarse_fxnl) %>%
  summarise(totcov = sum(mean_cover),
            spp = str_flatten(unique(code4), collapse = ", "),
            S = length(unique(code4))) %>%
  ungroup() %>%
  # change nut_trt control to XC
  mutate(nut_trt = dplyr::recode(nut_trt, N = "XC"),
         fulltrt = gsub("N", "XC", fulltrt))

# create dataset with everything
allcov <- subset(covcoarse, yr == 2021, select = -yr) %>%
  rbind(natcoarse_all[names(.)])

# note: plot 33 not seeded because not herbicided (block 4, NW plot)


# drop any background native plant that recruited.. just want to know if herbicide advantaged seeded natives against ambient non-natives
simplecov <- dplyr::select(allcov, -c(spp, S)) %>%
  #combine N-fixer and forb
  mutate(fxnl_grp2 = ifelse(grepl("Forb|N-fix", fxnl_grp), "Forb", fxnl_grp)) %>%
  group_by(plot, seedtrt, herbicide, fulltrt, block, ppt_trt, nut_trt, fxnl_grp2, nativity, native_seeded) %>%
  summarise(totcov = sum(totcov), .groups = "drop_last") %>%
  # drop any background native
  subset(!(native_seeded == "Background" & nativity == "Native")) %>%
  unite(fullgrp, native_seeded, nativity, fxnl_grp2, sep = "_") %>%
  # make sure all possible groups repped
  spread(fullgrp, totcov, fill = 0) %>%
  # also drop plot 33 since no experiment there
  subset(plot != 33) %>%
  gather(fullgrp, totcov, Background_Exotic_Forb:ncol(.))
# assign levels
simplecov$ppt_trt <- factor(simplecov$ppt_trt, levels = c("XC", "D", "W"))
simplecov$nut_trt <- factor(simplecov$nut_trt, levels = c("XC", "F", "C"))
simplecov$herbicide <- factor(simplecov$herbicide, levels = c("Non-herbicided", "Herbicided"))
simplecov$seedtrt <- factor(simplecov$seedtrt, levels = c("Unseeded", "Native seeded"))
simplecov$hillpos <- ifelse(simplecov$block < 3, "low", "high")
simplecov$hillpos <- factor(simplecov$hillpos, levels = c("low", "high"))
simplecov$fxnlgrp <- with(simplecov, ifelse(grepl("Forb", fullgrp), "Forb", "Grass"))
simplecov$fxnlgrp <- factor(simplecov$fxnlgrp, levels = c("Forb", "Grass"))
# compare native seeded in herb vs. non-herb to see if any differences
simplecov$nativity <- with(simplecov, ifelse(grepl("Nati", fullgrp), "Native", "Exotic"))
simplecov$nativity <- factor(simplecov$nativity)

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
  summarise(cov = max(cov), .groups = "drop_last") %>%
  spread(spp, cov)

# join to neighborhood
neighborhood <- merge(neighborhood, seededspp)
neighborhood_counts <- mutate_at(neighborhood, .vars = nats, function(x) ifelse(x %in% c(1:9), x+1, ifelse(x > 0 & x<1, 1, x)))
neighborhood_PA <- mutate_at(neighborhood, .vars = nats, function(x) ifelse(x>0, 1, 0))



# -- VIZ NATIVE RECRUITS ----
# check on abundance of each nat per treatment
left_join(seededspp, trtkey) %>%
  gather(code4, cover, BRCA:TRCI) %>%
  group_by(herbicide, block > 2) %>%
  ggplot(aes(code4, cover, group = paste(herbicide, block >2), col = block > 2)) +
  geom_jitter(aes(shape = herbicide), width = 0.2,height = 0, alpha = 0.5, size = 2) +
  stat_summary(aes(shape = herbicide), alpha = 0.7) +
  facet_grid(ppt_trt~nut_trt) +
  coord_flip()

# plot by species to see cover change across treatments relative for each spp
left_join(seededspp, trtkey) %>%
  gather(code4, cover, BRCA:TRCI) %>%
  subset(code4 != "TRCI") %>%
  mutate(hillpos = ifelse(block > 2, "Upslope", "Downslope")) %>%
  #group_by(herbicide, block > 2) %>%
  ggplot(aes(fulltrt, cover, group = paste(fulltrt, herbicide), col = herbicide)) +
  #geom_jitter(aes(shape = herbicide), width = 0.2,height = 0, alpha = 0.5, size = 2) +
  stat_summary(aes(shape = nut_trt), position = position_dodge(width = 0.4), alpha = 0.9) +
  scale_shape_manual(name = "Amendment", values = c(C = 16, F = 15, N = 4), labels = c("Compost", "Fert", "None")) +
  #scale_shape_manual(name = "Subplot", values = c(21,22)) +
  #scale_color_brewer(name = "Amendment", palette = "Set2", labels = c("Compost", "Fert", "None")) +
  scale_color_brewer(name = "Subplot", palette = "Set2") +
  #scale_fill_discrete(name = "Precip trt") +
  facet_grid(code4~hillpos, scales = "free")


# don't split by hillslope
left_join(seededspp, trtkey) %>%
  gather(code4, cover, BRCA:TRCI) %>%
  subset(code4 != "TRCI") %>%
  mutate(hillpos = ifelse(block > 2, "Upslope", "Downslope")) %>%
  #group_by(herbicide, block > 2) %>%
  ggplot(aes(fulltrt, cover, group = paste(fulltrt, herbicide), col = herbicide)) +
  #geom_jitter(aes(shape = herbicide), width = 0.2,height = 0, alpha = 0.5, size = 2) +
  stat_summary(aes(shape = nut_trt), position = position_dodge(width = 0.4), alpha = 0.9) +
  scale_shape_manual(name = "Amendment", values = c(C = 16, F = 15, N = 4), labels = c("Compost", "Fert", "None")) +
  #scale_shape_manual(name = "Subplot", values = c(21,22)) +
  #scale_color_brewer(name = "Amendment", palette = "Set2", labels = c("Compost", "Fert", "None")) +
  scale_color_brewer(name = "Subplot", palette = "Set2") +
  #scale_fill_discrete(name = "Precip trt") +
  facet_grid(code4~., scales = "free")


# overall experiment rank abundance (to assess dist of species)
sppcover <- group_by(natlong, yr, block, plot, fulltrt, nut_trt, ppt_trt, herbicide, seedtrt, code4) %>%
  summarise(max_cov = max(pct_cover)) %>%
  ungroup() %>%
  mutate(expgroup = paste(plot, fulltrt, herbicide, seedtrt, sep = "_")) %>%
  subset(seedtrt != "Unseeded")

# plot overall rank abundance
global_ranks <- group_by(sppcover, code4) %>%
  summarise(totcov = sum(max_cov),
            nplots = length(unique(expgroup)),
            prop_plots = nplots/length(unique(sppcover$expgroup))) %>%
  ungroup() %>%
  arrange(desc(totcov), desc(nplots)) %>%
  mutate(spprank = 1:nrow(.),
         code4_fac = factor(code4, levels = code4[order(spprank)])) %>%
  # join spp info
  left_join(unique(natlong[c("code4", "fxnl_grp", "nativity", "native_seeded")]))

gather(global_ranks, met, val, totcov:prop_plots) %>%
  subset(met != "nplots") %>%
  ggplot(aes(code4_fac, val, fill = fxnl_grp, col = native_seeded)) +
  geom_col(lwd = 2, alpha = 0.8) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_viridis_d() +
  scale_color_manual(values = c("transparent", "black")) +
  facet_grid(met~., scales = "free")

# check by treatment
trt_ranks <- group_by(sppcover, fulltrt) %>%
  mutate(allplots = length(unique(expgroup))) %>%
  group_by(fulltrt, code4) %>%
  summarise(totcov = sum(max_cov),
            nplots = length(unique(expgroup)),
            prop_plots = nplots/allplots) %>% # 105 = total # unique experiment subplots
  ungroup() %>%
  distinct() %>%
  arrange(fulltrt, desc(totcov), desc(nplots)) %>%
  group_by(fulltrt) %>%
  mutate(spprank = 1:length(code4)) %>%
         #code4_fac = factor(code4, levels = code4[order(spprank)])) %>%
  # join spp info
  left_join(unique(natlong[c("code4", "fxnl_grp", "nativity", "native_seeded")]))

ggplot(trt_ranks, aes(spprank, totcov, fill = fxnl_grp)) +
  geom_col(aes(col = native_seeded)) +
  geom_text(aes(spprank, totcov+1, label =code4),angle = 45, size = 2.5, hjust = 0, vjust = 0) + 
  scale_fill_viridis_d() +
  scale_color_manual(values = c("transparent", "black")) +
  facet_wrap(~fulltrt)

ggplot(trt_ranks, aes(spprank, prop_plots, fill = fxnl_grp)) +
  geom_col(aes(col = native_seeded)) +
  geom_text(aes(spprank, prop_plots, label =code4),angle = 90, size = 2.5, hjust = -0.05, vjust = 0.5 ) + 
  scale_y_continuous(limits = c(0,1.1)) +
  scale_fill_viridis_d() +
  scale_color_manual(values = c("transparent", "black")) +
  facet_wrap(~fulltrt)

# check by blocks
block_ranks <- group_by(sppcover, block) %>%
  mutate(allplots = length(unique(expgroup))) %>%
  group_by(block, code4) %>%
  summarise(totcov = sum(max_cov),
            nplots = length(unique(expgroup)),
            prop_plots = nplots/allplots) %>% # 105 = total # unique experiment subplots
  ungroup() %>%
  distinct() %>%
  arrange(block, desc(totcov), desc(nplots)) %>%
  group_by(block) %>%
  mutate(spprank = 1:length(code4)) %>%
  ungroup() %>%
  #code4_fac = factor(code4, levels = code4[order(spprank)])) %>%
  # join spp info
  left_join(unique(natlong[c("code4", "fxnl_grp", "nativity", "native_seeded")])) %>%
  data.frame()

ggplot(block_ranks, aes(spprank, totcov, fill = fxnl_grp)) +
  geom_col(aes(col = native_seeded)) +
  geom_text(aes(spprank, totcov+1, label =code4),angle = 45, size = 2.5, hjust = 0, vjust = 0) + 
  scale_fill_viridis_d() +
  scale_color_manual(values = c("transparent", "black")) +
  facet_wrap(~block)

ggplot(block_ranks, aes(spprank, prop_plots, fill = fxnl_grp)) +
  geom_col(aes(col = native_seeded), lwd = 2) +
  geom_text(aes(spprank, prop_plots, label =code4),angle = 90, size = 2.5, hjust = -0.1, vjust = 0.5 ) + 
  scale_y_continuous(limits = c(0,1.1)) +
  scale_fill_viridis_d() +
  scale_color_manual(values = c("transparent", "black")) +
  facet_wrap(~block)

# visualize interactions within blocks
ggplot(subset(simplecov, seedtrt == "Native seeded" & grepl("Exotic", fullgrp)), aes(nut_trt, totcov, col = ppt_trt, group = paste(ppt_trt, herbicide, fullgrp))) +
  geom_line(aes(lty = fullgrp)) +
  facet_wrap(herbicide~block, nrow = 2)
ggplot(subset(simplecov, seedtrt == "Native seeded" & grepl("Native", fullgrp)), aes(ppt_trt, totcov, col = nut_trt, group = paste(nut_trt, herbicide, fullgrp))) +
  geom_line(aes(lty = fullgrp)) +
  facet_wrap(herbicide~block, nrow = 2)

ggplot(subset(simplecov, seedtrt == "Native seeded" & ppt_trt == "XC"), aes(nut_trt, totcov, col = fullgrp, group = paste(ppt_trt, herbicide, fullgrp))) +
  geom_line(aes()) +
  facet_wrap(herbicide~block, nrow = 2)
ggplot(subset(simplecov, seedtrt == "Native seeded" & ppt_trt == "D"), aes(nut_trt, totcov, col = fullgrp, group = paste(ppt_trt, herbicide, fullgrp))) +
  geom_line(aes()) +
  facet_wrap(herbicide~block, nrow = 2)
ggplot(subset(simplecov, seedtrt == "Native seeded" & ppt_trt == "W"), aes(nut_trt, totcov, col = fullgrp, group = paste(ppt_trt, herbicide, fullgrp))) +
  geom_line(aes()) +
  facet_wrap(herbicide~block, nrow = 2)
subset(simplecov, seedtrt == "Native seeded") %>%
  mutate(ppt_trt = relevel(ppt_trt, "D")) %>%
ggplot(aes(block, totcov, col = fullgrp, lty = herbicide, group = paste(ppt_trt, fullgrp, herbicide))) +
  geom_point(aes()) +
  geom_smooth(se = F) +
  facet_wrap(nut_trt~ppt_trt, nrow = 3)
subset(neighborhood, select = c(plot:hillpos, BRCA:NEMA)) %>%
  gather(spp, totcov, BRCA:NEMA) %>%
  ggplot(aes(block, totcov, col = spp, lty = herbicide, group = paste(ppt_trt, spp, herbicide))) +
  geom_point(aes()) +
  geom_smooth(se = F) +
  facet_wrap(nut_trt~ppt_trt, nrow = 3)
subset(neighborhood, select = c(plot:hillpos, BRCA:NEMA)) %>%
  gather(spp, totcov, BRCA:NEMA) %>%
  mutate(ppt_trt = relevel(ppt_trt, "D")) %>%
  ggplot(aes(ppt_trt, totcov, col = spp, lty = herbicide, group = paste(spp, herbicide))) +
  geom_line(aes()) +
  scale_x_discrete(expand = c(0,0)) +
  facet_wrap(block~nut_trt, nrow = 4)




# -- GJAM ----
# -- Try coarse functional group to assess overall response to mgmt ----
# > note seeded natives are split from background natives here

# 4 types of cover groups per plot
natcoarse_ydata <- subset(natcoarse_all, select = c(plot, seedtrt, herbicide, coarse_fxnl, totcov)) %>%
  # clean up fxnl group names for spreading as colnames
  mutate(coarse_fxnl = gsub(" |-", "", coarse_fxnl)) %>%
  # make cover count data (trace amounts --> 1)
  # decimals = str_extract(totcov, "(?<=[.])[0-9]+"),
  # tenths = as.numeric(substr(decimals, 1, 1)),
  # hundredths = as.numeric(substr(decimals, 2, 2)),
  # addct = ifelse(tenths == 5 & !is.na(tenths), 1+ hundredths, hundredths))
  spread(coarse_fxnl, totcov, fill = 0) %>%
  mutate(abbrv_herb = ifelse(grepl("^Herbi", herbicide), "Herbi", "NoHerbi"),
         abbrv_seed2 = ifelse(grepl("Native", seedtrt), "NatSeed", "Unseed")) %>%
  unite(rowid, abbrv_herb, abbrv_seed2, sep = "", remove = T) %>%
  mutate(rowid = paste(rowid, plot, sep = "_")) %>%
  data.frame()

natcoarse_xdata <- subset(natcoarse_ydata, select = c(rowid, plot, seedtrt, herbicide)) %>%
  left_join(trtkey) %>%
  # convert env vars to factors
  mutate(seedtrt = factor(seedtrt, levels = c("Unseeded", "Native seeded")),
         herbicide = gsub("Non-herb", "NonHerb", herbicide),
         herbicide = factor(herbicide,ordered = T,  levels = c("NonHerbicided", "Herbicided")),
         nut_trt = factor(nut_trt, ordered = T, levels = c("N", "F", "C")),
         ppt_trt = factor(ppt_trt, ordered = T, levels = c("XC", "D", "W")),
         block = factor(block))
rownames(natcoarse_xdata) <- natcoarse_xdata$rowid
# clean up colnames
colnames(natcoarse_xdata) <- gsub("_trt", "Trt", colnames(natcoarse_xdata))

# assign rownames and check colsums on fxnl grps
rownames(natcoarse_ydata) <- natcoarse_ydata$rowid 
natcoarse_ydata <- subset(natcoarse_ydata, select = -c(plot, seedtrt, herbicide, rowid))  
sumcheck <- sapply(natcoarse_ydata, function(x) sum(x > 0))
sort(sumcheck)
natcoarse_ydata <- natcoarse_ydata[,sumcheck > 10]

# make effort
# effort is the sum of all % cover
natcoarse_effort_sum <- subset(natcoarse_all, select = c(plot, seedtrt, herbicide, coarse_fxnl, totcov)) %>%
  mutate(herbicide = gsub("Non-herb", "NonHerb", herbicide)) %>%
  group_by(plot, seedtrt, herbicide) %>%
  summarise(effort = sum(totcov)) %>%
  ungroup() %>%
  left_join(natcoarse_xdata[c("rowid", "plot", "seedtrt", "herbicide")]) %>%
  data.frame()
rownames(natcoarse_effort_sum) <- natcoarse_effort_sum$rowid
natcoarse_effort_sum <- natcoarse_effort_sum["effort"]
# assign total cover to all fxnl groups
natcoarse_effort <- natcoarse_ydata
natcoarse_effort[names(natcoarse_effort)] <- natcoarse_effort_sum$effort

# make effort and ydata matrices
natcoarse_ydata <- as.matrix(natcoarse_ydata)
natcoarse_effort <- as.matrix(natcoarse_effort)
# still not accepting controls, rename to alphabetize order
#natcoarse_xdata <- natcoarse()
# build model

# remake factors because still not working
natcoarse_xdata <- natcoarse_xdata %>%
  mutate(seedtrt = recode_factor(seedtrt, Unseeded = "0unseed", `Native seeded` = "1seed"),
         herbicide = recode_factor(herbicide, NonHerbicided = "0noHerb", Herbicided = "1herb"),
         pptTrt = recode_factor(pptTrt, XC = "0xc", D = "1d", W = "2w"),
         nutTrt = recode_factor(nutTrt, N = "0noNut", F = "1fert", C = "2comp")
  )
summary(natcoarse_xdata)
glimpse(natcoarse_xdata)
natcoarse_model <- as.formula(~ seedtrt + herbicide + pptTrt * nutTrt)
# build dimension reduction list (using same as in shown in tutorial)
#rlist <- list(r = 8, N = 10)

#make sure levels are set correctly (model output is alphabetizing ppt and nut trts)
#natcoarse_xdata$pptTrt <- relevel(natcoarse_xdata$pptTrt, ref = "XC")
#natcoarse_xdata$nutTrt <- relevel(natcoarse_xdata$nutTrt, ref = "N")
# build priors
natcoarse_prior <- gjamPriorTemplate(formula = natcoarse_model, 
                                     xdata =  natcoarse_xdata,
                                     ydata = natcoarse_ydata,
                                     lo = -Inf, 
                                     hi = Inf
) # not setting hi or lo limit, going with default setting ([-Inf, Inf]) for all spp

# build model list (same as in tutorial)
natcoarse_mlist <- list(ng=15000, burnin=5000, typeNames = 'CA', random = "block")
#betaPrior = prior, effort = elist <- these were used in tutorial, but looking at j clark's vignette doesn't look like I have to use them

# run model
natcoarse_gjam <- gjam(formula = natcoarse_model, 
                       xdata = natcoarse_xdata,
                       ydata = natcoarse_ydata,
                       modelList = natcoarse_mlist)
summary(natcoarse_gjam)
gjamPlot(natcoarse_gjam, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, SAVEPLOTS = T, 
                                         outFolder = "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/gjam_analyses/fxnl_group_pptxnut_blockrando/"))

str(natcoarse_xdata)

# try to plot interaction
xWC <- xWF <- xDC <- xDF <- natcoarse_gjam$inputs$xUnstand[1,]
xWC[c("pptTrt2w", "nutTrt1fert", "nutTrt2comp", "pptTrt2w:nutTrt2comp")]<- c(1,0,1,0)
xWC
xDC[c("pptTrt1d", "nutTrt1fert", "nutTrt2comp", "pptTrt1d:nutTrt2comp")] <- c(1,0,1,0)
testfit <- gjamIIE(natcoarse_gjam, xDC)
gjamIIEplot(testfit,
            #response = "BackgroundExoticForb", 
            #response = "SeededNativeForb",
            response = "SeededNativeGrass",
            effectMu = c("main", "direct", "ind"), effectSd = c("main", "direct", "ind"))
natcoarse_xdata


# -- Try with individual native-seeded cover -----
# 4 types of cover groups per plot
# drop unseeded trt
natcoarse_ydata <- subset(natcoarse_all, !grepl("Seeded", coarse_fxnl) & seedtrt != "Unseeded", select = c(plot, seedtrt, herbicide, coarse_fxnl, totcov)) %>%
  # clean up fxnl group names for spreading as colnames
  mutate(coarse_fxnl = gsub(" |-", "", coarse_fxnl)) %>%
  # add individual seeded spp
  rbind(gather(subset(neighborhood, select = c(plot, seedtrt, herbicide, BRCA:NEMA)), coarse_fxnl, totcov, BRCA:NEMA)) %>%
  # make cover count data (trace amounts --> 1)
  # decimals = str_extract(totcov, "(?<=[.])[0-9]+"),
  # tenths = as.numeric(substr(decimals, 1, 1)),
  # hundredths = as.numeric(substr(decimals, 2, 2)),
  # addct = ifelse(tenths == 5 & !is.na(tenths), 1+ hundredths, hundredths))
  spread(coarse_fxnl, totcov, fill = 0) %>%
  mutate(abbrv_herb = ifelse(grepl("^Herbi", herbicide), "Herbi", "NoHerbi"),
         abbrv_seed2 = ifelse(grepl("Native", seedtrt), "NatSeed", "Unseed")) %>%
  unite(rowid, abbrv_herb, abbrv_seed2, sep = "", remove = T) %>%
  mutate(rowid = paste(rowid, plot, sep = "_")) %>%
  data.frame()

natcoarse_xdata <- subset(natcoarse_ydata, select = c(rowid, plot, seedtrt, herbicide)) %>%
  left_join(trtkey) %>%
  # convert env vars to factors
  mutate(hillpos = ifelse(block <3, "Downslope", "Upslope"),
         hillpos = factor(hillpos),
         seedtrt = factor(seedtrt, levels = c("Unseeded", "Native seeded")),
         herbicide = gsub("Non-herb", "NonHerb", herbicide),
         herbicide = factor(herbicide,ordered = T,  levels = c("NonHerbicided", "Herbicided")),
         nut_trt = factor(nut_trt, ordered = T, levels = c("N", "F", "C")),
         ppt_trt = factor(ppt_trt, ordered = T, levels = c("XC", "D", "W")),
         block = factor(block))
rownames(natcoarse_xdata) <- natcoarse_xdata$rowid
# clean up colnames
colnames(natcoarse_xdata) <- gsub("_trt", "Trt", colnames(natcoarse_xdata))

# assign rownames and check colsums on fxnl grps
rownames(natcoarse_ydata) <- natcoarse_ydata$rowid 
natcoarse_ydata <- subset(natcoarse_ydata, select = -c(plot, seedtrt, herbicide, rowid))  
sumcheck <- sapply(natcoarse_ydata, function(x) sum(x > 0))
sort(sumcheck)
natcoarse_ydata <- natcoarse_ydata[,sumcheck > 10]

# make effort
# effort is the sum of all % cover
natcoarse_effort_sum <- subset(natcoarse_all, seedtrt != "Unseeded", select = c(plot, seedtrt, herbicide, coarse_fxnl, totcov)) %>%
  mutate(herbicide = gsub("Non-herb", "NonHerb", herbicide)) %>%
  group_by(plot, seedtrt, herbicide) %>%
  summarise(effort = sum(totcov)) %>%
  ungroup() %>%
  left_join(natcoarse_xdata[c("rowid", "plot", "seedtrt", "herbicide")]) %>%
  data.frame()
rownames(natcoarse_effort_sum) <- natcoarse_effort_sum$rowid
natcoarse_effort_sum <- natcoarse_effort_sum["effort"]
# assign total cover to all fxnl groups
natcoarse_effort <- natcoarse_ydata
natcoarse_effort[names(natcoarse_effort)] <- natcoarse_effort_sum$effort

# make effort and ydata matrices
natcoarse_ydata <- as.matrix(natcoarse_ydata)
natcoarse_effort <- as.matrix(natcoarse_effort)
# still not accepting controls, rename to alphabetize order
#natcoarse_xdata <- natcoarse()
# build model

# remake factors because still not working
natcoarse_xdata <- natcoarse_xdata %>%
  mutate(seedtrt = recode_factor(seedtrt, Unseeded = "0unseed", `Native seeded` = "1seed"),
         herbicide = recode_factor(herbicide, NonHerbicided = "0noHerb", Herbicided = "1herb"),
         pptTrt = recode_factor(pptTrt, XC = "0xc", D = "1d", W = "2w"),
         nutTrt = recode_factor(nutTrt, N = "0noNut", F = "1fert", C = "2comp")
  )
summary(natcoarse_xdata)
glimpse(natcoarse_xdata)
natcoarse_model <- as.formula(~ herbicide + pptTrt * nutTrt)
# build dimension reduction list (using same as in shown in tutorial)
#rlist <- list(r = 8, N = 10)

#make sure levels are set correctly (model output is alphabetizing ppt and nut trts)
#natcoarse_xdata$pptTrt <- relevel(natcoarse_xdata$pptTrt, ref = "XC")
#natcoarse_xdata$nutTrt <- relevel(natcoarse_xdata$nutTrt, ref = "N")
# build priors
natcoarse_prior <- gjamPriorTemplate(formula = natcoarse_model, 
                                     xdata =  natcoarse_xdata,
                                     ydata = natcoarse_ydata,
                                     lo = -Inf, 
                                     hi = Inf
) # not setting hi or lo limit, going with default setting ([-Inf, Inf]) for all spp

# build model list (same as in tutorial)
natcoarse_mlist <- list(ng=30000, burnin=5000, typeNames = 'CA', random = "block")
#betaPrior = prior, effort = elist <- these were used in tutorial, but looking at j clark's vignette doesn't look like I have to use them

# run model
natcoarse_gjam <- gjam(formula = natcoarse_model, 
                       xdata = natcoarse_xdata,
                       ydata = natcoarse_ydata,
                       modelList = natcoarse_mlist)
summary(natcoarse_gjam)
gjamPlot(natcoarse_gjam, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, SAVEPLOTS = T, 
                                         outFolder = "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/gjam_analyses/indiv_vfxnl_pptxnut_30kchains/"))

str(natcoarse_xdata)

# try to plot interaction
xWC <- xWF <- xDC <- xDF <- natcoarse_gjam$inputs$xUnstand[1,]
xWC[c("pptTrt2w", "nutTrt1fert", "nutTrt2comp", "pptTrt2w:nutTrt2comp")]<- c(1,0,1,0)
xWC
xDC[c("pptTrt1d", "nutTrt1fert", "nutTrt2comp", "pptTrt1d:nutTrt2comp")] <- c(1,0,1,0)
testfit <- gjamIIE(natcoarse_gjam, xDC)
gjamIIEplot(testfit,
            #response = "BackgroundExoticForb", 
            #response = "SeededNativeForb",
            response = "BackgroundExoticForb",
            effectMu = c("direct", "ind"), effectSd = c("direct", "ind"))


# try GJAM predict to compare model fit vs. observed
gjamPredict(natcoarse_gjam)



# -- Don't split hillslope, block as random effect -----
# 4 types of cover groups per plot
# drop unseeded trt
natcoarse_ydata <- subset(natcoarse_all, !grepl("Seeded", coarse_fxnl) & seedtrt != "Unseeded", select = c(plot, seedtrt, herbicide, coarse_fxnl, totcov)) %>%
  # clean up fxnl group names for spreading as colnames
  mutate(coarse_fxnl = gsub(" |-", "", coarse_fxnl)) %>%
  # add individual seeded spp
  rbind(gather(subset(neighborhood, select = c(plot, seedtrt, herbicide, BRCA:NEMA)), coarse_fxnl, totcov, BRCA:NEMA)) %>%
  mutate(herbicide = gsub("-h", "H", herbicide)) %>%
  spread(herbicide, totcov, fill = 0) %>%
  mutate(diffcov = NonHerbicided - Herbicided) %>%
  subset(select = c(plot, coarse_fxnl, diffcov)) %>%
  # make cover count data (trace amounts --> 1)
  # decimals = str_extract(totcov, "(?<=[.])[0-9]+"),
  # tenths = as.numeric(substr(decimals, 1, 1)),
  # hundredths = as.numeric(substr(decimals, 2, 2)),
  # addct = ifelse(tenths == 5 & !is.na(tenths), 1+ hundredths, hundredths))
  spread(coarse_fxnl, diffcov, fill = 0) %>%
  data.frame()

natcoarse_xdata <- subset(natcoarse_ydata, select = c(plot)) %>%
  left_join(trtkey) %>%
  # convert env vars to factors
  mutate(nut_trt = factor(nut_trt, ordered = T, levels = c("N", "F", "C")),
         ppt_trt = factor(ppt_trt, ordered = T, levels = c("XC", "D", "W")),
         block = factor(block))
# clean up colnames
colnames(natcoarse_xdata) <- gsub("_trt", "Trt", colnames(natcoarse_xdata))

# check colsums on fxnl grps
natcoarse_ydata <- subset(natcoarse_ydata, select = -c(plot))  
sumcheck <- sapply(natcoarse_ydata, function(x) sum(x != 0))
sort(sumcheck)
natcoarse_ydata <- natcoarse_ydata[,sumcheck > 10]

# make effort
# effort is the sum of all % cover
natcoarse_effort_sum <- subset(natcoarse_all, seedtrt != "Unseeded", select = c(plot, herbicide, coarse_fxnl, totcov)) %>%
  # sum cover across herbicide trtment subplots
  group_by(plot) %>%
  summarise(effort = sum(totcov)) %>%
  ungroup() %>%
  data.frame()
natcoarse_effort_sum <- natcoarse_effort_sum["effort"]
# assign total cover to all fxnl groups
natcoarse_effort <- natcoarse_ydata
natcoarse_effort[names(natcoarse_effort)] <- natcoarse_effort_sum$effort

# make effort and ydata matrices
natcoarse_ydata <- as.matrix(natcoarse_ydata)
natcoarse_effort <- as.matrix(natcoarse_effort)
# still not accepting controls, rename to alphabetize order
#natcoarse_xdata <- natcoarse()
# build model

# remake factors because still not working
natcoarse_xdata <- natcoarse_xdata %>%
  mutate(pptTrt = recode_factor(pptTrt, XC = "0xc", D = "1d", W = "2w"),
         nutTrt = recode_factor(nutTrt, N = "0noNut", F = "1fert", C = "2comp")
  )
summary(natcoarse_xdata)
glimpse(natcoarse_xdata)
natcoarse_model <- as.formula(~ pptTrt * nutTrt)
# build dimension reduction list (using same as in shown in tutorial)
#rlist <- list(r = 8, N = 10)

#make sure levels are set correctly (model output is alphabetizing ppt and nut trts)
#natcoarse_xdata$pptTrt <- relevel(natcoarse_xdata$pptTrt, ref = "XC")
#natcoarse_xdata$nutTrt <- relevel(natcoarse_xdata$nutTrt, ref = "N")
# build priors
natcoarse_prior <- gjamPriorTemplate(formula = natcoarse_model, 
                                     xdata =  natcoarse_xdata,
                                     ydata = natcoarse_ydata,
                                     lo = -Inf, 
                                     hi = Inf
) # not setting hi or lo limit, going with default setting ([-Inf, Inf]) for all spp

# build model list (same as in tutorial)
natcoarse_mlist <- list(ng=15000, burnin=6000, typeNames = 'CON', random = "block")
#betaPrior = prior, effort = elist <- these were used in tutorial, but looking at j clark's vignette doesn't look like I have to use them

# run model
natcoarse_gjam <- gjam(formula = natcoarse_model, 
                       xdata = natcoarse_xdata,
                       ydata = natcoarse_ydata,
                       modelList = natcoarse_mlist)
summary(natcoarse_gjam)
gjamPlot(natcoarse_gjam, plotPars = list(GRIDPLOTS = T))
gjamPlot(natcoarse_gjam, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, SAVEPLOTS = T, 
                                         outFolder = "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/gjam_analyses/indiv_vfxnl_pptxnut_blockrandom_5000x15000/"))

str(natcoarse_xdata)

# try to plot interaction
xWC <- xWF <- xDC <- xDF <- natcoarse_gjam$inputs$xUnstand[1,]
xWC[c("pptTrt2w", "nutTrt1fert", "nutTrt2comp", "pptTrt2w:nutTrt2comp")]<- c(1,0,1,0)
xWC
xDC[c("pptTrt1d", "nutTrt1fert", "nutTrt2comp", "pptTrt1d:nutTrt2comp")] <- c(1,0,1,0)
testfit <- gjamIIE(natcoarse_gjam, xDC)
gjamIIEplot(testfit,
            #response = "BackgroundExoticForb", 
            response = "BRCA",
            #response = "BackgroundExoticForb",
            #effectMu = c("int"), effectSd = c("int"))
            effectMu = c("direct", "ind"), effectSd = c("direct", "ind"))

# check model
natcoarse_gjam$fit$DIC
natcoarse_gjam$inputs$designTable
natcoarse_gjam$inputs$factorBeta$missFacSpec
natcoarse_gjam$parameters$sensTable
natcoarse_gjam$chains$bgibbs

# check sensitivity
cols  <- c( '#1f78b4', '#33a02c' )
natcoarse_sense <- gjamSensitivity(natcoarse_gjam)
ynames <- colnames(natcoarse_gjam$inputs$y)
natcoarse_sense_ca <- gjamSensitivity(natcoarse_gjam, ynames)
ylim = range(c(natcoarse_sense, natcoarse_sense_ca))
nt = ncol(natcoarse_sense)
boxplot(natcoarse_sense, boxwex = 0.2,  at = 1:nt - .15, col = cols[1], log='y',
         ylim = ylim, xaxt = 'n', xlab = 'Predictors', ylab='Sensitivity')
boxplot(natcoarse_sense_ca, boxwex = 0.2, at = 1:nt + .15, col = cols[2], add=T,
         xaxt = 'n')
axis( 1, at=1:nt,labels=colnames(natcoarse_sense))

data.frame(natcoarse_sense) %>%
  gather(predictor, val, 1:ncol(.)) %>%
  ggplot(aes(predictor, log(val))) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))

beta_df <- data.frame(natcoarse_gjam$parameters$betaMu) %>%
  mutate(predictor = rownames(.)) %>%
  gather(response, beta, 1:(ncol(.)-1))
betaSE_df <- data.frame(natcoarse_gjam$parameters$betaSe) %>%
  mutate(predictor = rownames(.)) %>%
  gather(response, betaSE, 1:(ncol(.)-1))
beta_df <- left_join(beta_df, betaSE_df)

ggplot(subset(beta_df, predictor != "intercept"), aes(predictor, beta)) +
  geom_hline(aes(yintercept = 0), lty = 2, col = "red") +
  geom_errorbar(aes(ymax = betaSE + beta, ymin = beta-betaSE)) +
  geom_point() +
  #theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(response~ .) +
  coord_flip()

natcoarse_gjam$parameters$betaStandXWmu
gjamOrdination(natcoarse_gjam, method = "NMDS")


# -- Try diffing cover in herb v nonherb ----
natcoarse_ydata <- subset(natcoarse_all, !grepl("Seeded", coarse_fxnl) & seedtrt != "Unseeded", select = c(plot, seedtrt, herbicide, coarse_fxnl, totcov)) %>%
  # clean up fxnl group names for spreading as colnames
  mutate(coarse_fxnl = gsub(" |-", "", coarse_fxnl)) %>%
  # add individual seeded spp
  rbind(gather(subset(neighborhood, select = c(plot, seedtrt, herbicide, BRCA:NEMA)), coarse_fxnl, totcov, BRCA:NEMA)) %>%
  # make cover count data (trace amounts --> 1)
  # decimals = str_extract(totcov, "(?<=[.])[0-9]+"),
  # tenths = as.numeric(substr(decimals, 1, 1)),
  # hundredths = as.numeric(substr(decimals, 2, 2)),
  # addct = ifelse(tenths == 5 & !is.na(tenths), 1+ hundredths, hundredths))
  spread(coarse_fxnl, totcov, fill = 0) %>%
  mutate(abbrv_herb = ifelse(grepl("^Herbi", herbicide), "Herbi", "NoHerbi"),
         abbrv_seed2 = ifelse(grepl("Native", seedtrt), "NatSeed", "Unseed")) %>%
  unite(rowid, abbrv_herb, abbrv_seed2, sep = "", remove = T) %>%
  mutate(rowid = paste(rowid, plot, sep = "_")) %>%
  data.frame()

natcoarse_xdata <- subset(natcoarse_ydata, select = c(rowid, plot, seedtrt, herbicide)) %>%
  left_join(trtkey) %>%
  # convert env vars to factors
  mutate(hillpos = ifelse(block <3, "Downslope", "Upslope"),
         hillpos = factor(hillpos),
         seedtrt = factor(seedtrt, levels = c("Unseeded", "Native seeded")),
         herbicide = gsub("Non-herb", "NonHerb", herbicide),
         herbicide = factor(herbicide,ordered = T,  levels = c("NonHerbicided", "Herbicided")),
         nut_trt = factor(nut_trt, ordered = T, levels = c("N", "F", "C")),
         ppt_trt = factor(ppt_trt, ordered = T, levels = c("XC", "D", "W")),
         block = factor(block))
rownames(natcoarse_xdata) <- natcoarse_xdata$rowid
# clean up colnames
colnames(natcoarse_xdata) <- gsub("_trt", "Trt", colnames(natcoarse_xdata))

# assign rownames and check colsums on fxnl grps
rownames(natcoarse_ydata) <- natcoarse_ydata$rowid 
natcoarse_ydata <- subset(natcoarse_ydata, select = -c(plot, seedtrt, herbicide, rowid))  
sumcheck <- sapply(natcoarse_ydata, function(x) sum(x > 0))
sort(sumcheck)
natcoarse_ydata <- natcoarse_ydata[,sumcheck > 10]

# make effort
# effort is the sum of all % cover
natcoarse_effort_sum <- subset(natcoarse_all, seedtrt != "Unseeded", select = c(plot, seedtrt, herbicide, coarse_fxnl, totcov)) %>%
  mutate(herbicide = gsub("Non-herb", "NonHerb", herbicide)) %>%
  group_by(plot, seedtrt, herbicide) %>%
  summarise(effort = sum(totcov)) %>%
  ungroup() %>%
  left_join(natcoarse_xdata[c("rowid", "plot", "seedtrt", "herbicide")]) %>%
  data.frame()
rownames(natcoarse_effort_sum) <- natcoarse_effort_sum$rowid
natcoarse_effort_sum <- natcoarse_effort_sum["effort"]
# assign total cover to all fxnl groups
natcoarse_effort <- natcoarse_ydata
natcoarse_effort[names(natcoarse_effort)] <- natcoarse_effort_sum$effort

# make effort and ydata matrices
natcoarse_ydata <- as.matrix(natcoarse_ydata)
natcoarse_effort <- as.matrix(natcoarse_effort)
# still not accepting controls, rename to alphabetize order
#natcoarse_xdata <- natcoarse()
# build model

# remake factors because still not working
natcoarse_xdata <- natcoarse_xdata %>%
  mutate(seedtrt = recode_factor(seedtrt, Unseeded = "0unseed", `Native seeded` = "1seed"),
         herbicide = recode_factor(herbicide, NonHerbicided = "0noHerb", Herbicided = "1herb"),
         pptTrt = recode_factor(pptTrt, XC = "0xc", D = "1d", W = "2w"),
         nutTrt = recode_factor(nutTrt, N = "0noNut", F = "1fert", C = "2comp")
  )
summary(natcoarse_xdata)
glimpse(natcoarse_xdata)
natcoarse_model <- as.formula(~ herbicide + pptTrt * nutTrt)
# build dimension reduction list (using same as in shown in tutorial)
#rlist <- list(r = 8, N = 10)

#make sure levels are set correctly (model output is alphabetizing ppt and nut trts)
#natcoarse_xdata$pptTrt <- relevel(natcoarse_xdata$pptTrt, ref = "XC")
#natcoarse_xdata$nutTrt <- relevel(natcoarse_xdata$nutTrt, ref = "N")
# build priors
natcoarse_prior <- gjamPriorTemplate(formula = natcoarse_model, 
                                     xdata =  natcoarse_xdata,
                                     ydata = natcoarse_ydata,
                                     lo = -Inf, 
                                     hi = Inf
) # not setting hi or lo limit, going with default setting ([-Inf, Inf]) for all spp

# build model list (same as in tutorial)
natcoarse_mlist <- list(ng=30000, burnin=5000, typeNames = 'CA', random = "block")
#betaPrior = prior, effort = elist <- these were used in tutorial, but looking at j clark's vignette doesn't look like I have to use them

# run model
natcoarse_gjam <- gjam(formula = natcoarse_model, 
                       xdata = natcoarse_xdata,
                       ydata = natcoarse_ydata,
                       modelList = natcoarse_mlist)
summary(natcoarse_gjam)
gjamPlot(natcoarse_gjam, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, SAVEPLOTS = T, 
                                         outFolder = "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/gjam_analyses/indiv_vfxnl_pptxnut_30kchains/"))





