# starter script to clean and format competition recruitment data from jan 2020 surveys

# notes: 
# 1/5/20: quick code so LMH and AS can preview data during CTW jan survey trip (1/2 - 1/6/2020), decide if need to do more in field while here
# CTW will clean up and expand code later


# -- SETUP ----
library(tidyverse)
options(stringsAsFactors = F)
theme_set(theme_bw())

# specify dropbox pathway (varies by user -- ctw can tweak this later)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/Competition/Competition_EnteredData")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/Competition/Competition_EnteredData"
}else{
  ## probs LMH and AS? (fix if not -- rproj is two levels down in github repo)
  datpath <- "~/Dropbox/USDA-compost/Data/Competition/Competition_EnteredData"
}

dats <- list.files(datpath, full.names = T)
phytos <- read.csv(dats[grep("phyto", dats, ignore.case = T)])
comp <- read.csv(dats[grep("competitor", dats)])

# read in J. Larson dry seed mass to screen for overcounts in density (more likely spp present in background seed bank so density enhanced)
seed_mass <- read.csv("Data-cleaning/Larson_CA_dryseedmass.csv")
# read brad traits for HOMU seed mass
brad_traits <- read.csv("~/Dropbox/Cali\ Trait\ Data/Brad_Trait\ screening/Greenhouse_TraitScreening_10_SpeciesData.csv", skip = 11)


#quick check
glimpse(phytos)
glimpse(comp)
glimpse(seed_mass); sort(unique(seed_mass$species)) # has AVFA instead of AVBA but close enough for quick check
glimpse(brad_traits)
# note > CTW not sure if JL's dat include seed attachments or is for pure seed only.. guessing AS weighed w/ seed attachments on


# specify plotting cols
ppt_cols <- c(D = "brown", W = "dodgerblue", XC = "seagreen")
plant_cols <- c(AVBA = "darkgreen", HOMU = "lightgreen", TACA = "limegreen", LOMU = "blue", 
                ERBO  = "red", TRHI = "orchid")



# -- PREP LOOKUP TABLES ----
# make lookup table for subsample frame to scale to meter-square
scale_lt <- data.frame(subsample_cm = c("5x5", "10x10", "25x25", "50x50"),
                       area_cm2 = c(5*5, 10*10, 25*25, 50*50)) %>%
  # make half plot scale for reality check with density (i.e. does stems projected for 50x50cm exceed amount seeded?)
  mutate(scale_half = (50*50)/area_cm2,
         # full meter scale factor
         scale_m2 = (100*100)/area_cm2)

seed_lt <- subset(seed_mass, grepl("avef|erob|lolm|taec|trih", species)) %>%
  group_by(species) %>%
  summarise(Seed = mean(perseedwt)) %>%
  ungroup() %>%
  rename("ID" = "species") %>%
  rbind(brad_traits[brad_traits$ID == "HORMUR", c("ID", "Seed")]) %>%
  mutate(ID = casefold(ID, upper = T),
         ID = paste0(substr(ID, 1,2), substr(ID, 4,5)),
         # scale to half meter density -- seeded at 8g per m2 (2g per half m2)
         max_density_halfm2 = 2/Seed)  
# review
seed_lt # HOMU wgt is definitely without the seed attachment
  

# -- PREP COMPETITOR DENSITY ---
# tidy and standardize to density per m2
## want to gather all subsample rep counts (up to 4 per subplot) and remove NAs
comp_tidy <- select(comp, block:plot4) %>%
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
         density_m2 = stems * scale_m2) # keep tidy df separate to look at density patchiness within plots 

# calculate average density
mean_density <- group_by(comp_tidy, block, nut_trt, ppt_trt, subplot, background) %>%
  summarise(mean_density_halfm2 = mean(density_half), # if took 4 reps at 25x25cm will average to correct 50x50cm tally 
            mean_density_1m2 = mean(density_m2),
            nobs = length(rep)) %>% # check
  ungroup()


# -- PREP PHYTOS ----
# recode N in nut_trt to XC
phytos <- mutate(phytos, nut_trt = recode(nut_trt, "N" = "XC"))


# -- VISUALIZE DATA ----
# competitor density so far..
# need to write this way for now because two spp only have 1 obs and stat_summ will NA those in mean and SE calc
group_by(comp_tidy, nut_trt, ppt_trt, background) %>%
  summarise(se_density = sd(density_m2)/sqrt(length(density_m2)),
            density_m2 = mean(density_m2)) %>%
  ggplot(aes(nut_trt, density_m2)) +
  geom_errorbar(aes(ymax = density_m2 + se_density, ymin = density_m2 - se_density), width = 0.2) +
  geom_point() +
  ggtitle("Mean competitor density +- 1SE") +
  facet_grid(background~ppt_trt, scales = "free_y")

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
  ggtitle("Variability in 50x50cm competitor density by subsample size") +
  facet_grid(background~subsample_cm, scales = "free_y")



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
       subtitle = "Average of all blocks (2 and 4) and subplots pooled (n=14 per phyto") +
  coord_flip() +
  facet_grid(ppt_trt~nut_trt)


# does phyto recruitment correlate with density (i.e. if AVBA phytos didn't recruit in XC drought, is AVBA comp density low there too?)
filter(phytos, !is.na(stems)) %>%
  group_by(block, nut_trt, ppt_trt, phyto) %>%
  summarise(se_stems = sd(stems)/sqrt(length(stems)),
            mean_stems = mean(stems)) %>%
  ungroup() %>%
  left_join(dplyr::select(mean_density, block:ppt_trt, background, mean_density_halfm2, mean_density_1m2), by = c("block", "nut_trt", "ppt_trt", "phyto" = "background")) %>%
  ggplot(aes(mean_density_halfm2, mean_stems)) +
  geom_errorbar(aes(ymax = mean_stems + se_stems, ymin = mean_stems - se_stems), width = 0.2) +
  ## > not enough comp surveys yet for std errors
  #geom_errorbarh(aes(xmax = mean_density_halfm2 + se_density, xmin = mean_density_halfm2 + se_density), width = 0.2) +
  geom_point(aes(col = ppt_trt)) +
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
# > AVBA has most -- 33% of subplots have no AVBA phytos
# > 20% subplots have no ERBO *at phytos* but ERBO is often present in background
# others are generally there in background as well

# how many 0s in comp dens?
nrow(subset(mean_density, mean_density_halfm2 == 0)) # no zeros in comp density

# 2) look for high counts in density
density_check <- dplyr::select(mean_density, block:mean_density_halfm2) %>%
  left_join(dplyr::select(seed_lt, -Seed), by = c("background" = "ID")) %>%
  # infill AVBA with AVFA dat
  mutate(max_density_halfm2 = ifelse(background == "AVBA", seed_lt$max_density_halfm2[seed_lt$ID == "AVFA"], max_density_halfm2)) %>%
  filter(mean_density_halfm2 > max_density_halfm2)
# review
View(arrange(density_check, background, block, nut_trt))
# > AVBA weighs less than AVFA so 70s and 80s numbers may be fine.. but also B3 and B4 have a lot of background AVBA. e.g. 132 count is an AVBA-invaded plot
# > ERBO also abundant in background for blocks 1-3
# > TRHI high in background in lower blocks.. and is high in fertilized plot in B4, but not by much..
