# starter script to clean and format competition recruitment data from jan 2020 surveys

# notes: 
# 1/5/20: quick code so LMH and AS can preview data during CTW jan survey trip (1/2 - 1/6/2020), decide if need to do more in field while here
# CTW will clean up and expand code later
# 1/21/20: none of the code after density lookup table will work (updated data entry, need to now update code.. will get to it later this week)

# Important note about seed weights and seeding weights from repo wiki ("2019 10_07_seeding expt setup" page):
# "All seeds were pre-weighed to 8g/m2, which is 2g per half meter squared of raw seed weight (not including husk or awns)"
# > Question (1/22/20).. but CTW saw awns in ground in field during Jan 2020 trip.. when AS returns, confirm raw seed only weighed out and seeded, adjust max density possible in code accordingly
# > For now coding following what's written on GitHub wiki



# -- SETUP ----
rm(list = ls()) # clean enviro
# libraries needed
library(tidyverse)
# change default settings
na_vals <- c("", " ", NA, "NA")
options(stringsAsFactors = F)
theme_set(theme_bw())

# specify dropbox pathway (varies by user -- ctw can tweak this later)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/Competition/")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/Competition/"
}else{
  ## probs LMH and AS? (fix if not -- rproj is two levels down in github repo)
  datpath <- "~/Dropbox/USDA-compost/Data/Competition/"
}

# read in raw data
dats <- list.files(paste0(datpath, "Competition_EnteredData"), full.names = T)
phytos <- read.csv(dats[grep("phyto", dats, ignore.case = T)], na.strings = na_vals)
comp <- read.csv(dats[grep("competitor", dats)], na.strings = na_vals)
photokey <- read.csv(dats[grep("photo", dats)], na.strings = na_vals)

# load seed mass data to crunch qty seeded
# > reading in others dats for now (do we have our own measurements?)
# read in J. Larson dry seed mass to screen for overcounts in density (more likely spp present in background seed bank so density enhanced)
seed_mass <- read.csv("Data-cleaning/Larson_CA_dryseedmass.csv")
# read brad traits for HOMU seed mass
brad_traits <- read.csv("~/Dropbox/Cali\ Trait\ Data/Brad_Trait\ screening/Greenhouse_TraitScreening_10_SpeciesData.csv", skip = 11)


#quick check
glimpse(phytos)
glimpse(comp)
glimpse(seed_mass); sort(unique(seed_mass$species)) # has AVFA instead of AVBA but close enough for quick check
glimpse(brad_traits)
# note > CTW not sure if JL's dat include seed attachments or is for pure seed only.. repo wiki says AS weighed out seeds w/out attachments


# specify plotting cols
ppt_cols <- c(D = "brown", W = "dodgerblue", XC = "seagreen")
plant_cols <- c(AVBA = "darkgreen", HOMU = "lightgreen", TACA = "limegreen", LOMU = "blue", 
                ERBO  = "red", TRHI = "orchid", Control = "grey40")




# -- PREP LOOKUP TABLES ----
# make lookup table for subsample frame to scale to meter-square
scale_lt <- data.frame(subsample_cm = c("5x5", "10x10", "25x25", "50x50"),
                       area_cm2 = c(5*5, 10*10, 25*25, 50*50)) %>%
  # make half plot scale for reality check with density (i.e. does stems projected for 50x50cm exceed amount seeded?)
  mutate(scale_half = (50*50)/area_cm2,
         # full meter scale factor
         scale_m2 = (100*100)/area_cm2)

# join seed mass data to LUT to project max density possible per species
# > note: here is where to adjust based on whether seeds weighed out with or without awns/husks/attachments (looking at you ERBO)
seed_lt <- subset(seed_mass, grepl("avef|erob|lolm|taec|trih", species)) %>% #start with Julie's seed weights
  group_by(species) %>%
  summarise(Seed = mean(perseedwt)) %>%
  ungroup() %>%
  rename("ID" = "species") %>%
  # add jl for source
  mutate(source = "JL") %>%
  # append brad's seed mass data (BB has HOMU, JL does not. BB overlaps with JL on other spp we used)
  rbind(cbind(brad_traits[grepl("avef|erob|lolm|taec|trih|hormu", brad_traits$ID, ignore.case = T), c("ID", "Seed")], source = "BB")) %>%
  mutate(ID = casefold(ID, upper = T),
         ID = paste0(substr(ID, 1,2), substr(ID, 4,5)),
         # scale to half meter density -- seeded at 8g per m2 (2g per half m2)
         max_density_halfm2 = 2/Seed)  
# review
seed_lt # HOMU wgt is definitely without the seed attachment

# how do the two sources of info compare in seed weights?
ggplot(seed_lt, aes(ID, Seed, col = source)) +
  geom_point()  # all close but AVFA and ERBO..
# what does that mean for density?
ggplot(seed_lt, aes(ID, max_density_halfm2, col = source)) +
  geom_point() #hm..
# > Brad's allows more individuals per subplot (JL and BB agree for TRHI tho)
# press on with data prep, look at how affects recruitemnt rates further down in code..



# -- PREP COMPETITOR DENSITY ---
# tidy and standardize to density per m2
## want to gather all subsample rep counts (up to 4 per subplot) and remove NAs
## also have cols for QA resamples (either re-sampled subplot to confirm high densities, or subsampled control background)
## also have photo info to join..
## CTW also kept track of time start/stop per plot (to know, but also to help estimate how long this might take in spring)
##.. ultimately what's needed for analysis is tidy (long) format data table.. maybe add col for whether resample/QC sample then gather and project densities
## also maybe build other species table separately, then can add back in if needed

# start with add col to indicate qc samples
comp <- comp %>%
  #sort by date -- earlier first
  arrange(survey_date) %>%
  grouped_df(names(.)[(grep("block", names(.))):((grep("backgr", names(.))))]) %>%
  # number survey event
  mutate(survey_event = seq(1, length(survey_date),1)) %>%
  ungroup()

# compile experiment as designed (no extra spp, no background spp in control plot), only first samples (no resamples for QC)
# i.e. survey_event  == 1
comp_tidy <- filter(comp, survey_event == 1) %>%
  select(block:plot4) %>%
  # recode None in nut_trt to XC (to not confuse C with control)
  mutate(nut_trt = recode(nut_trt, "N" = "XC")) %>%
  gather(met, val, stems1:plot4) %>%
  mutate(rep = parse_number(met),
         met = gsub("[1-4]", "", met)) %>%
  spread(met, val) %>%
  filter(!is.na(stems)) %>%
  rename(subsample_cm = plot) %>%
  mutate(stems = as.numeric(stems)) %>%
  left_join(scale_lt) %>%
  mutate(density_half = stems*scale_half,
         density_m2 = stems * scale_m2) %>% # keep tidy df separate to look at density patchiness within plots 
  # check for any subsamples of different size within same subplot
  group_by(block, nut_trt, ppt_trt, subplot, background) %>%
  mutate(frame_check = length(unique(subsample_cm))>1) %>%
  ungroup()

# compile data for re-sampled plots and control background
comp_resample <- filter(comp, survey_event == 2) %>%
  select(block:plot4) %>%
  # recode None in nut_trt to XC (to not confuse C with control)
  mutate(nut_trt = recode(nut_trt, "N" = "XC")) %>%
  gather(met, val, stems1:plot4) %>%
  mutate(rep = parse_number(met),
         met = gsub("[1-4]", "", met)) %>%
  spread(met, val) %>%
  filter(!is.na(stems)) %>%
  rename(subsample_cm = plot) %>%
  mutate(stems = as.numeric(stems)) %>%
  left_join(scale_lt) %>%
  mutate(density_half = stems*scale_half,
         density_m2 = stems * scale_m2) %>% # keep tidy df separate to look at density patchiness within plots 
  # check for any subsamples of different size within same subplot
  group_by(block, nut_trt, ppt_trt, subplot, background) %>%
  mutate(frame_check = length(unique(subsample_cm))>1) %>%
  ungroup()

# compile data on extra spp sampled to qualify seeded recruitment (e.g. background ERBO density in control plot)
comp_extra <- filter(comp, !is.na(nonseed_spp1)) %>%
  # drop background competitor cols -- only interested in the extras
  select(-(c("time_start", "time_stop", names(.)[grep("^stems1", names(.))]: names(.)[grep("pict", names(.))]))) %>%
  gather(met,val, nonseed_spp1:nonseed_plot5) %>%
  filter(!is.na(val)) %>%
  #parse number
  mutate(rep = str_extract(met, "[:digit:]"),
         met = gsub("[0-9]", "", met)) %>%
  spread(met, val) %>%
  # rearrange cols a little
  dplyr::select(1:survey_event, nonseed_spp, rep:ncol(.)) %>%
  # reconvert stem count to numberic
  mutate(nonseed_stems = as.numeric(nonseed_stems)) %>%
  # join scale LUT to project densities
  left_join(scale_lt, by = c("nonseed_plot" = "subsample_cm")) %>%
  mutate(density_half = nonseed_stems*scale_half,
         density_m2 = nonseed_stems * scale_m2) %>% # keep tidy df separate to look at density patchiness within plots 
  # check for any subsamples of different size within same subplot
  group_by(block, nut_trt, ppt_trt, subplot, background, survey_event, nonseed_spp) %>%
  mutate(frame_check = length(unique(nonseed_plot))>1) %>%
  ungroup()
  
# calculate average density
# first remove any multi-frame size subsample (frame_check == T), defer to larger frame sample (there's only one case and it's 50x50 vs 25x25cm)
mean_density <- comp_tidy %>%
  filter(!(frame_check & subsample_cm == "25x25")) %>%
  group_by(block, nut_trt, ppt_trt, subplot, background) %>%
  summarise(mean_density_halfm2 = mean(density_half), # if took 4 reps at 25x25cm will average to correct 50x50cm tally 
            mean_density_1m2 = mean(density_m2),
            nobs = length(rep)) %>% # check
  ungroup()

# write out mean_density for ctw ebio grant proposal (temporary file until we troubleshoot seed masses/flagging)
write_csv(x = mean_density, path = paste0(datpath, "Competition_CleanedData/competitor_meandensity_jan2020.csv"))


# -- PREP PHYTOS ----
# recode N in nut_trt to XC
phytos <- mutate(phytos, nut_trt = recode(nut_trt, "N" = "XC"))


# -- VISUALIZE DATA ----
# competitor density so far..
# need to write this way for now because two spp only have 1 obs and stat_summ will NA those in mean and SE calc
group_by(comp_tidy, nut_trt, ppt_trt, background) %>%
  summarise(se_density = sd(density_m2)/sqrt(length(density_m2)),
            density_m2 = mean(density_m2)) %>%
  ggplot(aes(nut_trt, density_m2, col = background)) +
  geom_errorbar(aes(ymax = density_m2 + se_density, ymin = density_m2 - se_density), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  scale_color_manual(values = plant_cols) +
  ggtitle("Mean competitor density +- 1SE") +
  facet_grid(.~ppt_trt, scales = "free_y") # can also facet by background, trying alltog with colors to compare more directly

# plot data points
ggplot(mean_density, aes(nut_trt, mean_density_1m2)) +
  #geom_violin() +
  geom_jitter(alpha = 0.6, width = 0.1) +
  ggtitle("Block-level mean competitor density (B4, B2)") +
  facet_grid(background~ppt_trt, scales = "free_y")

# plot patchiness as 50x50cm density
ggplot(comp_tidy, aes(nut_trt,density_half, col = ppt_trt)) +
  #geom_violin() +
  geom_jitter(alpha = 0.6, width = 0.1) +
  scale_color_manual(values = ppt_cols) +
  labs(y = "Competitor density (0.5 m^2)", x = "Nutrient treatment", title = "Variability in 50x50cm competitor density by subsample size",
       subtitle = "QA?: Density greater the smaller the sampling frame used, but use small frame when plot looks hi dens") +
  facet_grid(background~subsample_cm, scales = "free_y")
# > thought/Q: I wonder how much greater density in 10x10 is results of using smaller frame, applying that to whole 50x50cm area vs. density numbers really are that much greater..


# -- VIEW PHYTOMETERS ----
# plot freq of 1, 2 or 3 stems present .. but at 3 different levels (ppt trt, nut add, competitor.. ctw tired, add better plots as needed to judge)
subset(phytos, !is.na(stems)) %>%
  group_by(nut_trt, ppt_trt, background, phyto, stems) %>%
  summarise(nobs = length(block)) %>%
  ungroup() %>%
  ggplot(aes(stems, nobs)) +
  geom_col() +
  labs(y = "Count", x = "# phyto recruits (3 max)",
  title = "Tallies of phyto recruits by phyto (x panel) and competitor (y panel)",
  subtitle = "Sums across ppt trts, nutrient trts, and blocks") +
  facet_grid(background~phyto)

subset(phytos, !is.na(stems)) %>%
  group_by(nut_trt, ppt_trt, background, phyto, stems) %>%
  summarise(nobs = length(block)) %>%
  ungroup() %>%
  ggplot(aes(stems, fill = ppt_trt)) +
  geom_bar(col = "grey10") +
  labs(x = "# phyto recruits (3 max)", title = "Tallies of phyto recruits by nutrient and ppt trts") +
  scale_fill_manual(values = ppt_cols) +
  facet_grid(nut_trt~phyto)
 
# recruitment trend 
ggplot(phytos, aes(phyto, stems)) +
  stat_summary() +
  labs(y = "Mean phyto recruits +- 1SE",
       title = "Phyto recruitment trend by nutrient (x panel) and ppt (y panel) trts",
       subtitle = "Average of all blocks and subplots pooled") +
  coord_flip() +
  facet_grid(ppt_trt~nut_trt) + 
  geom_hline(yintercept = 1, color = "grey")


# does phyto recruitment correlate with density (i.e. if AVBA phytos didn't recruit in XC drought, is AVBA comp density low there too?)
filter(phytos, !is.na(stems)) %>%
  group_by(block, nut_trt, ppt_trt, phyto) %>%
  summarise(se_stems = sd(stems)/sqrt(length(stems)),
            mean_stems = mean(stems)) %>%
  ungroup() %>%
  left_join(dplyr::select(mean_density, block:ppt_trt, background, mean_density_halfm2, mean_density_1m2), by = c("block", "nut_trt", "ppt_trt", "phyto" = "background")) %>%
  ggplot(aes(mean_density_halfm2, mean_stems, col = ppt_trt)) +
  geom_errorbar(aes(ymax = mean_stems + se_stems, ymin = mean_stems - se_stems), width = 0.2) +
  ## > not enough comp surveys yet for std errors
  #geom_errorbarh(aes(xmax = mean_density_halfm2 + se_density, xmin = mean_density_halfm2 + se_density), width = 0.2) +
  geom_point() +
  #geom_smooth(method = "lm", se =F, col = "black") +
  scale_color_manual(values = ppt_cols) +
  labs(y = "Species mean phyto recruitment", x = "Species mean 50x50cm density", 
       title = "Species mean comp density vs. mean recruitment by trts",
       subtitle = "i.e. if species came in well as background, did it recruit well in plot also?") +
  facet_grid(nut_trt~phyto, scales = "free_x")



# -- PREVIEW COMPETITION ---- 
# LOMU and TACA seemed to compete (impression while surveying) .. plot phyto recruitment vs competitor density for those two
## note: half of comp not done yet
subset(phytos, grepl("TACA|LOMU", background) & grepl("TACA|LOMU", phyto)) %>%
  # filter out same spp
  filter(!is.na(stems)) %>%
  #join background density
  left_join(mean_density) %>%
  ggplot(aes(mean_density_1m2, stems)) +
  geom_point(aes(col = phyto)) +
  labs(y = "# Recruits (3 max)", x = "Mean competitor density 1m2)",
       title = "TACA vs. LOMU: # recruits by competitor density") +
  scale_color_manual(values = plant_cols) +
  facet_wrap(~background, scales = "free")


# try all
# filter out same spp
filter(phytos, !is.na(stems) & background != "Control") %>%
  #join background density
  left_join(mean_density) %>%
  ggplot(aes(mean_density_1m2, stems, col = background)) +
  geom_jitter(height = 0.1) +
  labs(y = "# phyto recruits (3 max)", x = "Mean competitor density 1m2",
       title = "Phyto recruits by competitor density (panel = phyto)") +
  #scale_color_viridis_d() +
  scale_color_manual(values = plant_cols, name = "Competitor") +
  facet_wrap(~phyto, scales = "free")


# -- DATA/EXP QC -----
# 1) screen 0s for re-seeding
# how many 0s in comp density and in phytos per species?
filter(phytos, stems == 0) %>%
  group_by(phyto) %>%
  summarise(n_zeros = length(stems),
            # crunch effective non-recruitment rate (how many 0s at any plot [ignore 3 TPS -- just is there anything there])
            effect_nonrec_rate = (n_zeros/(35*6))*100) %>% # 35 comp plots * 6 subplots where seeded per comp plot
  # sort
  arrange(desc(n_zeros)) 
# > AVBA has most -- 31% of subplots have no AVBA phytos
# > 20% subplots have no ERBO *at phytos* but ERBO is often present in background
# others are generally there in background as well

# where didn't AVBA come in?
subset(phytos, stems == 0 & phyto == "AVBA") %>%
  # group_by(nut_trt, ppt_trt, background, phyto, stems) %>%
  # summarise(nobs = length(block)) %>%
  # ungroup() %>%
  ggplot(aes(block, fill = ppt_trt)) +
  geom_bar(col = "grey10") +
  labs(x = "Block", title = "# competition subplots where no AVBA phytos recruited",
       subtitle = "54 subplots per block for B1-3, B4 has 48 subplots") +
  scale_fill_manual(values = ppt_cols) #+
  #facet_grid(nut_trt~phyto)
# > definitely a recruitment-fail effect by hill position.. recruits better as move uphill-ish (b2 level with b1)


# how many 0s in comp dens?
nrow(subset(mean_density, mean_density_halfm2 == 0)) # no zeros in comp density

# 2) look for high counts in density
density_check <- dplyr::select(mean_density, block:mean_density_halfm2) %>%
  # create alt ID to join both BB and JL's projections for AVFA to AVBA
  mutate(ID = ifelse(background == "AVBA", "AVFA", background)) %>%
  left_join(dplyr::select(seed_lt, -Seed)) %>%
  filter(mean_density_halfm2 > max_density_halfm2)

# review
View(arrange(density_check, background, block, nut_trt))
# > AVBA weighs less than AVFA so 70s and 80s numbers may be fine.. but also B3 and B4 have a lot of background AVBA. e.g. 132 count is an AVBA-invaded plot
# > ERBO also abundant in background for blocks 1-3
# > TRHI high in background in lower blocks.. and is high in fertilized plot in B4, but not by much..


# 3) (after Jan 2020 trip) compare re-samples and QC samples to high counts identified above
# i.e. use comp_extra and comp_resample
rbind(cbind(comp_resample, survey_event = 2), cbind(comp_tidy, survey_event = 1)) %>%
  group_by(block, nut_trt, ppt_trt, subplot, background) %>%
  filter(length(unique(survey_event))>1) %>%
  #ungroup() %>%
  group_by(survey_event, subsample_cm, add = T) %>%
  summarise(mean_density_halfm2 = mean(density_half),
            se_halfm2 = sd(density_half)/sqrt(length(density_half))) %>%
  ungroup() %>%
  mutate(plotid = paste0(block, nut_trt, ppt_trt, subplot, background)) %>%
  left_join(seed_lt[seed_lt$source == "JL", c("ID", "max_density_halfm2")], by = c("background" = "ID")) %>%
  rename(JL_max = max_density_halfm2) %>%
  left_join(seed_lt[seed_lt$source == "BB", c("ID", "max_density_halfm2")], by = c("background" = "ID")) %>%
  ggplot(aes(as.factor(survey_event), mean_density_halfm2, group= plotid)) +
  geom_line(col = "grey50", position = position_dodge(width = 0.2), lty = 2) +
  geom_errorbar(aes(ymax = mean_density_halfm2 + se_halfm2, ymin = mean_density_halfm2 - se_halfm2), position = position_dodge(width = 0.2),  width = 0.1) +
  geom_point(aes(col = subsample_cm), size = 2, position = position_dodge(width = 0.2)) +
  #scale_color_manual(values = plant_cols)
  facet_wrap(~background, scales = "free_y") +
  labs(title = "Data QA: Competitor densities resurveyed for high values",
       subtitle = "Lesson: Don't use 10x10cm frames? Too noisy?",
       x = "Survey event",
       # I tried to plot JL and BB max densities within the panels but code didn't want to behave, so adding it as long string in caption
       caption = paste0("Max densities possible: ERBO, ", round(seed_lt$max_density_halfm2[seed_lt$source == "JL" & seed_lt$ID == "ERBO"], 0), "(JL), ", round(seed_lt$max_density_halfm2[seed_lt$source == "BB" & seed_lt$ID == "ERBO"],0), "(BB); TRHI: ", round(seed_lt$max_density_halfm2[seed_lt$source == "JL" & seed_lt$ID == "TRHI"], 0), "(JL), ", round(seed_lt$max_density_halfm2[seed_lt$source == "BB"& seed_lt$ID == "TRHI"],0), "(BB)")) +
  theme(plot.caption = element_text(hjust = 0))


# plot competitor densities vs. background species densities (comp_extra)
# a) visualize comp_extra as is
ggplot(comp_extra, aes(background, density_half, col = ppt_trt, shape = as.factor(block))) +
  geom_point(position = position_dodge(width = 0.3)) +
  scale_color_manual(values = ppt_cols) +
  facet_wrap(~nonseed_spp, scales = "free_x")

# b) Look at background mean densities vs. seeded competitor mean densities
group_by(comp_extra, block, nut_trt, ppt_trt, background, nonseed_spp) %>%
  summarise(mean_density_halfm2_nonseed = mean(density_half)) %>%
  ungroup() %>%
  mutate(nut_trt = recode(nut_trt, N= "XC")) %>%
  rename(subplot_background = background,
         background = nonseed_spp) %>%
  left_join(mean_density) %>%
  unite(grp, block, nut_trt, ppt_trt, sep = "", remove = T) %>%
  ggplot(aes(grp, mean_density_halfm2_nonseed)) +
  geom_point(aes(col = subplot_background), pch = 8) +
  geom_point(aes(grp, mean_density_halfm2)) +
  scale_color_manual(name = "Subplot\nsampled", values = plant_cols) +
  labs(y = "Mean density, half m^2", x = "Plot (block-nut-ppt)",
       title = "Seeded competitor density vs. background (unseeded) density",
       subtitle = "Black = seeded species (ours), colored = unseeded species (true background)") +
  facet_grid(.~background, scales = "free_x", space = "free_x")

