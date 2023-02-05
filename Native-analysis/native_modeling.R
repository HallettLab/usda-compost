# native analyses models
# author(s): ctw
# questions: caitlin.t.white@colorado.edu


# -- SETUP -----
rm(list = ls())
# load needed libraries
library(tidyverse)
library(vegan)
library(car)
library(lme4)
library(lmerTest)
library(multcomp)
library(MASS)
library(cowplot)
library(pscl)
library(glmmTMB) # for random effects in ZINF model

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

# specify dropbox pathway (varies by user)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/"
}else{
  ## probs LMH and AS? (fix if not -- rproj is two levels down in github repo)
  datpath <- "~/Dropbox/USDA-compost/Data/"
}

# climvar datpath for native experiment ther
climvar_path <- "~/Dropbox/USDA-climvar/Data/"

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

# climvar natives (not sure if clean)
#climvar_nats <- read.csv("", na.strings = na_vals)


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


# -- ANOVA or GLMMs ----
# plot 33 not seeded because not herbicided (block 4, NW plot)
# maybe try linear MM over ANOVA with block as random
# try ANOVA first

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

# plot to check it looks right
simplecov %>%
  mutate(totcov = ifelse(grepl("Native", fullgrp) & seedtrt == "Unseeded", NA, totcov)) %>%
  ggplot(aes(herbicide, totcov, fill = fullgrp)) + #interaction(fullgrp, seedtrt, herbicide)
  geom_boxplot(aes(group = interaction(fullgrp, herbicide)), varwidth = T, position = position_dodge(width = 0.5), alpha = 0.5) +
  geom_point(aes(col = fullgrp), position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("orchid", "seagreen1", "purple4", "seagreen4")) +
  scale_color_manual(values = c("orchid", "seagreen1", "purple4", "seagreen4")) +
  #facet_grid(seedtrt ~ fulltrt) +
  facet_wrap(~ nut_trt + seedtrt +ppt_trt, ncol = 6) +
  #facet_grid(fulltrt ~ seedtrt) #+
  theme(#axis.text.x = element_text(angle = 90))
        legend.position = "top",
        legend.direction = "horizontal",
        strip.background = element_rect(fill = "transparent"))

# compare background exotics in seeded vs. unseeded plots to see if any diff
lm_exotics <- lmerTest::lmer(totcov ~ fxnlgrp * ppt_trt + nut_trt + herbicide + seedtrt + (1|hillpos), data = subset(simplecov, grepl("Back", fullgrp))) # & seedtrt == "Unseeded"
summary(lm_exotics)
Anova(lm_exotics) # similar if use type III
ls_means(lm_exotics, pairwise = T)
lmerTest::ranova(lm_exotics)
# herbicide makes a difference, forb v. grass, ppt, and interaction between fnxl group x ppt_trt, and fxnlgrp x nut_trt
# nutrient amendment alone does not have main effect, seeding trt doesn't matter, no 3-way interaction between fxnl group, ppt and nut_trt


# compare native seeded in herb vs. non-herb to see if any differences
simplecov$nativity <- with(simplecov, ifelse(grepl("Nati", fullgrp), "Native", "Exotic"))
simplecov$nativity <- factor(simplecov$nativity)
lm_natives <- lmerTest::lmer(totcov ~ herbicide + nativity * fxnlgrp* ppt_trt + nut_trt + (1|hillpos), data = subset(simplecov, seedtrt == "Native seeded"))
summary(lm_natives)
Anova(lm_natives)
ls_means(lm_natives, pairwise = T)
lmerTest::ranova(lm_natives)


# -- SPECIES RESPONSE MODELS ----
# build 5 separate LMMs for each native species seeded
# compare with all species in one model for due diligence

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

# look at distribution
gather(seededspp, spp, cov, BRCA:TRCI) %>%
  ggplot(aes(spp, cov, fill = herbicide, group = paste(herbicide, spp))) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), pch = 21) +
  labs(y = "Abundance (%)", x = "Seeded native species") +
  theme(legend.position = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        legend.title = element_blank())

# mostly just FEMI, BRCA, and maybe NEMA that did ok-ish?
# look at it by treatment
gather(seededspp, spp, cov, BRCA:TRCI) %>%
  mutate(present = ifelse(cov >0, 1, 0)) %>%
  left_join(trtkey) %>%
  group_by(ppt_trt, nut_trt, herbicide, spp) %>%
  summarise(present = sum(present)) %>%
  ungroup() %>%
  mutate(ppt_trt = factor(ppt_trt, levels = c("W", "XC", "D"))) %>%
  ggplot(aes(nut_trt, present, fill = herbicide, group = paste(herbicide, spp))) +
  #geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), pch = 21)+
  labs(y = "# plots present", x = "Nutrient treatment") +
  facet_grid(ppt_trt~spp) +
  theme(legend.position = "none")
# what are counts of non-0?
sapply(seededspp[nats], function(x) sum(x >0))
# pct recruitment:
(sapply(seededspp[nats], function(x) sum(x >0))/nrow(seededspp))*100
# exclude TRCI because only came up in 1 place
sapply(seededspp[nats], function(x) summary(x[x>0])) # NEMA, FEMI, and BRCA definitely did the best. still keep ESCA because present frequently as BRCA or NEMA. (reason for logistic reg?)
# how many possible plots total?
nrow(seededspp)
# did some only recruit in herbicided?
sapply(seededspp[seededspp$herbicide == "Herbicided",nats], function(x) sum(x >0))
sapply(seededspp[seededspp$herbicide != "Herbicided",nats], function(x) sum(x >0)) # all but TRCI recruited regardles of herbicide

# join to neighborhood
neighborhood <- merge(neighborhood, seededspp)
neighborhood_counts <- mutate_at(neighborhood, .vars = nats, function(x) ifelse(x %in% c(1:9), x+1, ifelse(x > 0 & x<1, 1, x)))
neighborhood_PA <- mutate_at(neighborhood, .vars = nats, function(x) ifelse(x>0, 1, 0))
# grams look like they needed herbicide conditions, forbs less affected

# start with lmers on abundance
## FEMI ----
femi_lmer <- lmer(FEMI ~ herbicide + ppt_trt + nut_trt + (1|hillpos), data = neighborhood_counts)
femi_lmer_poisson <- glmer(FEMI ~ herbicide + ppt_trt + nut_trt + (1|hillpos), family = "poisson", data = neighborhood_counts)
femi_lmer_nb <- glmer.nb(FEMI ~ herbicide + ppt_trt + nut_trt + (1|hillpos), data = neighborhood_counts)
femi_lmer_nb_grass <- glmer.nb(FEMI ~ grass + ppt_trt + nut_trt + (1|hillpos), data = neighborhood_counts)
femi_lmer_bin <- glmer(FEMI ~  grass + ppt_trt + nut_trt + (1|hillpos), family = "binomial", data = neighborhood_PA)
femi_lmer_bin2 <- glmer(FEMI ~  grass + herbicide + ppt_trt + nut_trt + (1|hillpos), family = "binomial", data = neighborhood_PA)
femi_lmer_bin3 <- glmer(FEMI ~  herbicide + ppt_trt + nut_trt + (1|hillpos), family = "binomial", data = neighborhood_PA)
summary(femi_lmer)
summary(femi_lmer_poisson)
summary(femi_lmer_nb)
summary(femi_lmer_nb_grass) # throws error, mb because grass is in different units than categorical? warns about scale issues
summary(femi_lmer_bin) # bin will run with both grass and herbicide -- when only include grass not herb, grass and drought are signif, when herb in just herb signif
Anova(femi_lmer_bin)
anova(femi_lmer, femi_lmer_poisson)
anova(femi_lmer, femi_lmer_nb)
anova(femi_lmer_bin, femi_lmer_bin2) # with herb is better model
anova(femi_lmer_bin3, femi_lmer_bin2) # including grass not helpful
anova(femi_lmer_poisson, femi_lmer_nb) # negbin model seems like best choice based on AIC? signif terms agree with gaussian lmer (but dist def not normal)
#another comparison for poisson v negbin
plot(fitted(femi_lmer_nb), resid(femi_lmer_nb))
plot(fitted(femi_lmer_poisson), resid(femi_lmer_poisson)) # resids worse for poisson, negbin better choice
plot(fitted(femi_lmer), resid(femi_lmer)) # bad bad
ls_means(femi_lmer, pairwise = T)
ls_means(femi_lmer_nb, pairwise = T)
lmerTest::ranova(femi_lmer)

# try zinf instead, compare w hurdle
femi_zinf <- zeroinfl(FEMI ~ herbicide + ppt_trt + nut_trt, data = neighborhood_counts, dist = "negbin", link = "logit")
femi_hurdle <- hurdle(FEMI ~ herbicide + ppt_trt + nut_trt, data = neighborhood_counts, dist = "negbin", link = "logit")
summary(femi_zinf)
summary(femi_hurdle)
summary(femi_lmer_nb)
# does this match w what data show?
# look at distribution
gather(seededspp, spp, cov, BRCA:TRCI) %>%
  subset(spp == "FEMI") %>%
  left_join(trtkey) %>%
  ggplot(aes(spp, cov, fill = herbicide, group = paste(herbicide, spp))) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), pch = 21) +
  facet_grid(ppt_trt ~ nut_trt)
# agrees.. does better in drought, and herbicide.. but better in drought patchy (not *so* much better than other ppts) so that's why drought marginal signif

## BRCA ----
brca_lmer_nb <- glmer.nb(BRCA ~ grass + ppt_trt + nut_trt + (1|hillpos), data = neighborhood_counts)
brca_glm_nb <- MASS::glm.nb(BRCA ~ herbicide + ppt_trt + nut_trt, data = neighborhood_counts)
brca_glm_nb2 <- MASS::glm.nb(BRCA ~ grass + ppt_trt + nut_trt, data = neighborhood_counts)
brca_glm <- glm(BRCA ~ herbicide + ppt_trt + nut_trt, data = neighborhood_counts)
brca_lmer_bin <- glmer(BRCA ~ herbicide + ppt_trt + nut_trt + (1|hillpos), family = "binomial", data = neighborhood_PA)
summary(brca_lmer_nb)
summary(brca_glm_nb)
summary(brca_glm_nb2)
anova(brca_lmer_nb1, brca_lmer_nb2)
anova(brca_glm_nb, brca_glm_nb2)
anova(brca_glm, brca_glm_nb)
plot(fitted(brca_glm_nb), resid(brca_glm_nb))
plot(fitted(brca_glm_nb2), resid(brca_glm_nb2)) # residuals seem worse including grass over herbicide
plot(fitted(brca_glm), resid(brca_glm)) # nope
plot(fitted(brca_lmer_nb), resid(brca_lmer_nb))


## NEMA ----
nema_lmer_nb <- glmer.nb(NEMA ~ herbicide + ppt_trt + nut_trt + (1|hillpos), data = neighborhood_counts)
summary(nema_lmer_nb)
nema_lmer_nb2 <- glmer.nb(NEMA ~ neighbors + ppt_trt + nut_trt + (1|hillpos), data = neighborhood_counts)
summary(nema_lmer_nb2)
anova(nema_lmer_nb, nema_lmer_nb2)
nema_lmer_bin <- glmer(NEMA ~ herbicide + ppt_trt + nut_trt + (1|hillpos), family = "binomial", data = neighborhood_PA)
summary(nema_lmer_bin)


## ESCA ----
esca_lmer_nb <- glmer.nb(ESCA ~ herbicide + ppt_trt + nut_trt + (1|hillpos), data = neighborhood_counts)
summary(esca_lmer_nb)
esca_lmer_nb2 <- glmer.nb(ESCA ~ grass + ppt_trt + nut_trt + (1|hillpos), data = neighborhood_counts)
summary(esca_lmer_nb2)
esca_glm_nb2 <- MASS::glm.nb(ESCA ~ neighbors + ppt_trt + nut_trt, data = neighborhood_counts)
summary(esca_glm_nb2)
esca_glm_nb <- MASS::glm.nb(ESCA ~ herbicide + ppt_trt + nut_trt, data = neighborhood_counts)
esca_glm_bin <- glm(ESCA ~ herbicide + ppt_trt + nut_trt, family = "binomial", data = neighborhood_PA)
esca_lmer_bin <- glmer(ESCA ~ grass + ppt_trt + nut_trt + (1|hillpos), family = "binomial", data = neighborhood_PA)
summary(esca_lmer_bin) # singular fit
summary(esca_lmer_nb2)
summary(esca_glm_nb)
summary(esca_glm_nb2)
summary(esca_glm_bin)
anova(esca_glm_nb, esca_glm_nb2)
Anova(esca_glm_nb2)
Anova(esca_glm_bin)
multcomp::glht(esca_glm_bin, linfct =multcomp::mcp(herbicide="Tukey"))
plot(fitted(esca_glm_nb), resid(esca_glm_nb))
plot(fitted(esca_glm_nb2), resid(esca_glm_nb2))
plot(fitted(esca_lmer_nb), resid(esca_lmer_nb))

## All species model ----- 
# try model with all species since having trouble
neighborhood_counts_long <- gather(neighborhood_counts, spp, cov, BRCA:TRCI) %>%
  subset(spp != "TRCI")

all_glm_nb <- MASS::glm.nb(cov ~ spp * herbicide * ppt_trt + nut_trt, data = neighborhood_counts_long)
summary(all_glm_nb)
all_glm_nb2 <- MASS::glm.nb(cov ~ spp * neighbors * ppt_trt + nut_trt, data = neighborhood_counts_long)
summary(all_glm_nb2)
all_glm_nb_pos <- MASS::glm.nb(cov ~ spp * herbicide * ppt_trt + nut_trt + hillpos, data = neighborhood_counts_long)
summary(all_glm_nb_pos)
Anova(all_glm_nb_pos)
Anova(all_glm_nb) # no convergence
Anova(all_glm_nb2)
all_glm_nb_simple <- MASS::glm.nb(cov ~ spp + herbicide + ppt_trt + nut_trt, data = neighborhood_counts_long)
summary(all_glm_nb_simple)
anova(all_glm_nb_simple, all_glm_nb)
all_glm_pos <- glm(cov ~ spp + herbicide + ppt_trt + nut_trt, family = "poisson", data = neighborhood_counts_long)
summary(all_glm_pos)

#compare poisson v nb fits
pchisq(2 * (logLik(all_glm_nb_simple) - logLik(all_glm_pos)), df =1, lower.tail = F) #nb is better fit


# glmers 
# simple
all_glmer_nb <- glmer.nb(cov ~ herbicide + ppt_trt + spp + (1| block), data = neighborhood_counts_long)
summary(all_glmer_nb)
isSingular(all_glmer_nb) # when species is included in model, singular fit
# interactive
all_glmer_nb3x <- MASS::glm.nb(cov ~ spp * herbicide * ppt_trt + nut_trt, data = neighborhood_counts_long)

# -- CLEANED UP NB GLM MODELS ----
# make hillpos factor
neighborhood_counts
# femi
femi_glm_nb <- MASS::glm.nb(FEMI ~ herbicide + ppt_trt + nut_trt, data = neighborhood_counts)
femi_glm_nb_hill <- MASS::glm.nb(FEMI ~ herbicide + ppt_trt + nut_trt + hillpos, data = neighborhood_counts)
femi_glm_nb_hill_grass <- MASS::glm.nb(FEMI ~ grass + ppt_trt + nut_trt + hillpos, data = neighborhood_counts)
femi_glm_nb_hill_x <- MASS::glm.nb(FEMI ~ herbicide * ppt_trt * nut_trt + hillpos, data = neighborhood_counts)
femi_glm_nb_neigh_hill <- MASS::glm.nb(FEMI ~ herbicide + neighbors + ppt_trt + nut_trt + hillpos, data = neighborhood_counts)
summary(femi_glm_nb_hill)
summary(femi_glm_nb_hill_grass)
summary(femi_glm_nb_hill_x)
anova(femi_glm_nb_hill, femi_glm_nb_hill_x) # go with simple
summary(femi_glm_nb_neigh_hill)
anova(femi_glm_nb, femi_glm_nb_hill) # hill more signif :(
anova(femi_glm_nb_hill, femi_glm_nb_neigh_hill) # including neighbor cov does not improve model
Anova(femi_glm_nb_hill) # hillpos and herb most signif
femi_glm_bin <- glm(FEMI ~ herbicide + ppt_trt + nut_trt + hillpos, family = "binomial", data = neighborhood_PA)
summary(femi_glm_nb)
summary(femi_glm_bin)
Anova(femi_glm_bin)


# brca
brca_glm_nb <- MASS::glm.nb(BRCA ~ herbicide + ppt_trt + nut_trt, data = neighborhood_counts)
brca_glm_nb_neigh <- MASS::glm.nb(BRCA ~ herbicide + neighbors + ppt_trt + nut_trt, data = neighborhood_counts)
brca_lmer_nb <- glmer.nb(BRCA ~ herbicide + ppt_trt + nut_trt + (1|hillpos), data = neighborhood_counts) # singular error
brca_glm_nb_x <- MASS::glm.nb(BRCA ~ herbicide * ppt_trt * nut_trt, data = neighborhood_counts)
brca_glm_nb_hill <- MASS::glm.nb(BRCA ~ herbicide + ppt_trt + nut_trt + hillpos, data = neighborhood_counts)
summary(brca_glm_nb)
summary(brca_glm_nb_neigh) # it's really drought that matters
summary(brca_lmer_nb) # says the same thing anyway
summary(brca_glm_nb_x)
summary(brca_glm_nb_hill) # error on this model, either way they say the same thing -- herbicide and drought matter
anova(brca_glm_nb, brca_glm_nb_x) # not strong enough to choose
anova(brca_glm_nb, brca_glm_nb_neigh) # including neighbor not an improvement
Anova(brca_glm_nb_x, type = "III")
Anova(brca_glm_nb) # says nutrient is signif.. F is marginally signif
Anova(brca_glm_nb_neigh, type = "III")
#presence
brca_glm_bin <- glm(BRCA ~ herbicide + ppt_trt + nut_trt + hillpos, family = "binomial", data = neighborhood_PA)
summary(brca_glm_bin)

# nema
nema_glm_nb <- MASS::glm.nb(NEMA ~ herbicide + ppt_trt + nut_trt + hillpos, data = neighborhood_counts)
nema_glm_nb_grass <- MASS::glm.nb(NEMA ~ herbicide + grass + ppt_trt + nut_trt + hillpos, data = neighborhood_counts)
nema_glm_nb_x <- MASS::glm.nb(NEMA ~ herbicide + ppt_trt * nut_trt * hillpos, data = neighborhood_counts)
nema_glm_nb_grass <- MASS::glm.nb(NEMA ~ herbicide + grass + ppt_trt + nut_trt + hillpos, data = neighborhood_counts)
nema_glm_bin <- glm(NEMA ~ herbicide + ppt_trt + nut_trt + hillpos, family = "binomial", data = neighborhood_PA)
summary(nema_glm_nb)
summary(nema_glm_nb_grass) # no difference in what's signif
summary(nema_glm_bin)
summary(nema_glm_nb_x) #no, same things matter, no interactions
summary(nema_glm_nb_grass) # hillpos better
anova(nema_glm_nb, nema_glm_nb_x) # interaction not better
Anova(nema_glm_nb_grass) # just hillpos, ppt trt and somewhat nut trt -- when try with grass only not signif
Anova(nema_glm_nb)
Anova(nema_glm_nb_x) # suggestive of interaction between ppt and hillpos

# esca
esca_glm_nb <- MASS::glm.nb(ESCA ~ herbicide + ppt_trt + nut_trt + hillpos, data = neighborhood_counts)
esca_glm_nb_neighbor <- MASS::glm.nb(ESCA ~ neighbors + ppt_trt + nut_trt + hillpos, data = neighborhood_counts)
esca_glm_nb_x <- MASS::glm.nb(ESCA ~ herbicide * ppt_trt + nut_trt + hillpos, data = neighborhood_counts)
esca_glm_bin <- glm(ESCA ~ herbicide + ppt_trt + nut_trt, family = "binomial", data = neighborhood_PA)
summary(esca_glm_nb)
summary(esca_glm_nb_neighbor) #nope
summary(esca_glm_nb_x)
anova(esca_glm_nb, esca_glm_nb_x) # interaction doesn't improve
summary(esca_glm_bin)
Anova(esca_glm_nb_neighbor)
Anova(esca_glm_nb) # what does it mean when hillpos matters and that's it? could be soil.. or grasses occupying different canopy height allow for esca?

# does esca respond to forb cover (competition?)
esca_glm_nb_grass <- MASS::glm.nb(ESCA ~ grass + ppt_trt + nut_trt + hillpos, data = neighborhood_counts)
summary(esca_glm_nb_grass)
Anova(esca_glm_nb_grass) # grass and hillpos both matter.. but they are correlated

# kns said try functional neighbor (e.g., grams) as proportion of total cover so not correlated w neighbors
neighborhood_counts_long$gram_prop <- (neighborhood_counts_long$grass/neighborhood_counts_long$neighbors)*100
# .. the units will be funky tho.. make rel so same units as cover

# try all spp in one model w zeroinf or hurdle
all_zinf <- zeroinfl(cov ~ spp + herbicide + ppt_trt + nut_trt, data = neighborhood_counts_long, dist = "negbin")
summary(all_zinf)
all_zinf_x <- zeroinfl(cov ~ spp + nut_trt + ppt_trt + neighbors, data = neighborhood_counts_long, dist = "negbin")
summary(all_zinf_x) # won't run w zeroinfl
# w gram prop
all_zinf_full <- zeroinfl(cov ~ spp * neighbors * ppt_trt + nut_trt, data = neighborhood_counts_long, dist = "negbin")
summary(all_zinf_full)

predict_zinf_x <- cbind(all_zinf_x$model, 
                        fitted = all_zinf_x$fitted.values, 
                        pred_count = predict(all_zinf_x, type = "count"),
                        pred_zero = round(predict(all_zinf_x, type = "zero"),4),
                        res = all_zinf_x$residuals)

# hurdle
all_hurdle <- hurdle(cov ~ spp + herbicide + ppt_trt + nut_trt, data = neighborhood_counts_long, dist = "negbin")
summary(all_hurdle)
all_hurdle_x <- hurdle(cov ~ spp * neighbors * nut_trt + ppt_trt , data = neighborhood_counts_long, dist = "negbin")
summary(all_hurdle_x) # will run w hurdle
all_hurdle_x_p <- hurdle(cov ~ spp * neighbors * nut_trt + ppt_trt , data = neighborhood_counts_long, dist = "poisson")
summary(all_hurdle_x_p)

# w gram prop
all_hurdle_full <- hurdle(cov ~ spp * neighbors + herbicide + ppt_trt + nut_trt , data = neighborhood_counts_long, dist = "negbin")
summary(all_hurdle_full) # 

predict(all_hurdle_x)
predict_hurdle_x <- cbind(all_hurdle_x$model, 
                        fitted = all_hurdle_x$fitted.values, 
                        pred_count = predict(all_hurdle_x, type = "count"),
                        pred_zero = predict(all_hurdle_x, type = "zero"),
                        res = all_hurdle_x$residuals)

# test if hurdle or zinf superior to nb glm
vuong(all_glm_nb_simple, all_zinf) # not by BIC corrected
vuong(all_glm_nb_simple, all_hurdle) # not by AIC corrected

glmmTMB()


# -- MULTIVARIATE TO ASSESS NEIGHBORHOOD RESPONSE TO SEEDING ----- 

# prep species abundance dataset with native seeded removed
# need companion environmental dataset with treatments per plot (ppt, nut, herbicided or not, seeded or not, hillslope position)

# assess tot cov before running NMDSes to see how removing natives seeded impacts totcov
# or can just run on relativized and unrelativized to compare if any major diffs
allbutnats <- subset(natlong, !code4 %in% nats, select = c(plot:mon, code4,pct_cover)) %>%
  # choose max cov within season
  group_by(plot, herbicide, seedtrt, code4) %>%
  summarise(pct_cover = max(pct_cover), .groups = "drop_last")

# prep cover -- choose yr 2021 only, and max cov by species (not sure how I can adjust size bc cov still relative)
sppcomp2021 <- subset(coverlong, yr ==2021 & plot != 33, select = c(plot:sample_event, code4, pct_cover)) %>%
  group_by(plot, code4) %>%
  summarise(pct_cover = max(pct_cover), .groups = "drop_last") %>%
  mutate(herbicide = "Non-herbicided", seedtrt = "Unseeded")

# join to all but nats
allbutnats <- rbind(allbutnats, sppcomp2021) %>%
  distinct() %>%
  # remove rare (present in fewer than 5% of plots)
  group_by(code4) %>%
  mutate(freq = length(pct_cover[pct_cover > 0])) %>%
  ungroup() %>%
  mutate(plots = length(unique(paste(plot, herbicide, seedtrt))),
         prop_pres = freq/plots) %>%
  subset(prop_pres >=0.05, select = -c(freq, plots, prop_pres)) %>%
  spread(code4, pct_cover, fill = 0) %>%
  #create unique ID code for each subplot
  mutate(natplotid = 1:nrow(.))

# make enviro dat
allbutnats_env <- dplyr::select(allbutnats, natplotid, plot:seedtrt) %>%
  # join main experiment trtment info
  left_join(trtkey) %>%
  # make hillpos col again
  mutate(hillpos = ifelse(block < 3, "low", "high"),
         hillpos = factor(hillpos, levels = c("low", "high")))

# prep dat for nmds and permanova
summary(rownames(allbutnats) == allbutnats$natplotid) # rowname is already what natplotid is, so don't need to keep col
allbutnats_cov <- dplyr::select(allbutnats, select = -c(plot:seedtrt, natplotid))
# check plotcov
allbutnats_env$totcov <- apply(allbutnats_cov,1, sum)
# add native seeded cover
allbutnats_env2 <- left_join(allbutnats_env, subset(neighborhood, select =c(plot, herbicide, seedtrt, BRCA:NEMA))) %>%
  # fill NAs for seeded spp with 0? or -1 for not there (no chance)
  replace_na(list(FEMI = -1, BRCA = -1, NEMA = -1, ESCA = -1))
# make PA dat too to try to draw hulls
allbutnats_envPA <- mutate_at(allbutnats_env2, .vars = c("ESCA", "FEMI", "BRCA", "NEMA"), function(x) ifelse(x>0, "P","A"))

allbutnats_bc <- vegdist(allbutnats_cov, method = "bray")
allbutnats_nmds <- metaMDS(allbutnats_cov, k = 3, trymax = 100)
# try relativized also
allbutnats_nmds_rel <- metaMDS(decostand(allbutnats_cov, method = "total"), trymax = 100)
plot(allbutnats_nmds)
plot(allbutnats_nmds_rel)
stressplot(allbutnats_nmds)
stressplot(allbutnats_nmds_rel)


ordiplot(allbutnats_nmds, type = "n")
orditorp(allbutnats_nmds, display="species",col="red", air=0.01)
orditorp(allbutnats_nmds, display="sites",col="grey50", air=0.01)
ordisurf(allbutnats_nmds,allbutnats_env2$FEMI,col="forestgreen")
ordisurf(allbutnats_nmds,allbutnats_env2$NEMA,col="chocolate3")
ordisurf(allbutnats_nmds,allbutnats_env2$ESCA,col="orchid")
ordisurf(allbutnats_nmds,allbutnats_env2$BRCA, display = "sites", col="blue") #brca is linear! on axes 1 and 2
ordihull(allbutnats_nmds, paste(allbutnats_env2$seedtrt, allbutnats_env2$herbicide), label = T, col = 1:4, lwd=3)

# axes 1 and 2
ordiplot(allbutnats_nmds, type = "n")
orditorp(allbutnats_nmds, display="species",col="red", air=0.01)
orditorp(allbutnats_nmds, display="sites",col="grey50", air=0.01)
ordiellipse(allbutnats_nmds, paste(allbutnats_env2$seedtrt, allbutnats_env2$herbicide), label = T, col = 1:4, lwd=3)
ordiellipse(allbutnats_nmds, allbutnats_envPA$BRCA, label = T, lwd=3, col = 1:2, choices = c(1,3), kind = "se")

# axes 1 and 3
ordiplot(allbutnats_nmds, type = "n", choices = c(1,3))
orditorp(allbutnats_nmds, display="species",col="red", air=0.01, choices = c(1,3))
orditorp(allbutnats_nmds, display="sites",col="grey50", air=0.01, choices = c(1,3))
ordiellipse(allbutnats_nmds, paste(allbutnats_env2$seedtrt, allbutnats_env2$herbicide), label = T, col = 1:4, lwd=3, choices = c(1,3))
ordisurf(allbutnats_nmds,allbutnats_env2$FEMI,col="forestgreen", choices = c(1,3)) # a little more linear
ordisurf(allbutnats_nmds,allbutnats_env2$ESCA,col="orchid", choices = c(1,3)) # a little more linear
ordisurf(allbutnats_nmds,allbutnats_env2$NEMA,col="chocolate3", choices = c(1,3))
ordiellipse(allbutnats_nmds, allbutnats_envPA$FEMI, label = T, lwd=3, col = 1:2, choices = c(1,3))
ordiellipse(allbutnats_nmds, allbutnats_envPA$ESCA, label = T, lwd=3, col = 5:6, choices = c(1,3), kind = "se")
ordiellipse(allbutnats_nmds, allbutnats_envPA$ESCA, label = T, lwd=3, col = 5:6, choices = c(1,3))
ordiellipse(allbutnats_nmds, allbutnats_envPA$herbicide, label = T, lwd=3, col = 3:4, choices = c(1,3), kind = "se")

# axes 2 and 3
ordiplot(allbutnats_nmds, type = "n", choices = c(2,3))
orditorp(allbutnats_nmds, display="species",col="red", air=0.01, choices = c(2,3))
orditorp(allbutnats_nmds, display="sites",col="grey50", air=0.01, choices = c(2,3))
ordiellipse(allbutnats_nmds, paste(allbutnats_env2$seedtrt, allbutnats_env2$herbicide), label = T, col = 1:4, lwd=3, choices = c(2,3))
ordisurf(allbutnats_nmds,allbutnats_env2$FEMI,col="forestgreen", choices = c(2,3)) # funny (or not) femi and esca capture in same space bc respond similar in models too
ordisurf(allbutnats_nmds,allbutnats_env2$ESCA,col="orchid", choices = c(2,3))
ordisurf(allbutnats_nmds,allbutnats_env2$NEMA,col="chocolate3", choices = c(2,3)) #linear!
ordiellipse(allbutnats_nmds, allbutnats_envPA$NEMA, label = T, lwd=3, choices = c(2,3), kind = "sd")
ordihull(allbutnats_nmds, allbutnats_envPA$seedtrt, label = T, col = 1:2,lwd=3, choices = c(2,3))
allbutnats_envPA$fulltrt <- with(allbutnats_envPA, paste(seedtrt, herbicide, ppt_trt, sep = "_"))
ordiellipse(allbutnats_nmds, allbutnats_envPA$fulltrt, label = T, col = 1:10,lwd=3, choices = c(2,3), kind = "se")

# relativized plot
ordiplot(allbutnats_nmds_rel, type = "n")
orditorp(allbutnats_nmds_rel, display="species",col="red", air=0.01)
orditorp(allbutnats_nmds_rel, display="sites",col="grey50", air=0.01)
ordiellipse(allbutnats_nmds_rel, allbutnats_envPA$NEMA, label = T, lwd=2)
ordiellipse(allbutnats_nmds_rel, allbutnats_envPA$FEMI, label = T, lwd=3, col = 2:3, kind = "se")
ordiellipse(allbutnats_nmds_rel, allbutnats_envPA$ESCA, label = T, lwd=3, col = 4:5, kind = "ehull")
ordisurf(allbutnats_nmds_rel,allbutnats_env2$NEMA,col="chocolate3") #linear!
ordisurf(allbutnats_nmds_rel,allbutnats_env2$FEMI,col="forestgreen")
ordisurf(allbutnats_nmds_rel,allbutnats_env2$ESCA,col="orchid")
ordisurf(allbutnats_nmds_rel,allbutnats_env2$BRCA,col="blue")
ordiellipse(allbutnats_nmds_rel, allbutnats_envPA$fulltrt, label = T, col = 1:10,lwd=1, kind = "se")
ordihull(allbutnats_nmds_rel, allbutnats_envPA$seedtrt, label = T, col = 1:10,lwd=1)

# permanova to see if native seeded, or herbicide explain comm diss at all
adonis2(allbutnats_bc ~ seedtrt + herbicide +ppt_trt + nut_trt, data = allbutnats_envPA, permutations = how(blocks = allbutnats_env2$block))
# for now I'm saying herbicide, ppt_trt, nut_trt, and FEMI explain community diss the most.. seeding treatment does to
adonis2(allbutnats_bc ~ seedtrt + herbicide + FEMI + ppt_trt + nut_trt, data = allbutnats_envPA, permutations = how(blocks = allbutnats_env2$block))
adonis2(allbutnats_bc ~ seedtrt * nut_trt + herbicide + FEMI + ppt_trt, data = allbutnats_envPA, permutations = how(blocks = allbutnats_env2$block))

# one more nmds using fxnl grps




# notes:
latlong_proj <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"
boulder_coords <- sp::SpatialPoints(coords = data.frame(lat = 40.0101190898142, long = -105.24232781534424))
sp::proj4string(boulder_coords) <- latlong_proj
sfrec_coords <- sp::SpatialPoints(coords = data.frame(lat = 39.253229939065406, long = -121.31334268578593))
sp::proj4string(boulder_coords) <- latlong_proj
sp::spDists(boulder_coords, sfrec_coords, longlat = T)


# -- GJAM ----
# try coarse functional group to assess overall response to mgmt
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
natcoarse_mlist <- list(ng=6000, burnin=1000, typeNames = 'CA')
#betaPrior = prior, effort = elist <- these were used in tutorial, but looking at j clark's vignette doesn't look like I have to use them

# run model
natcoarse_gjam <- gjam(formula = natcoarse_model, 
                       xdata = natcoarse_xdata,
                       ydata = natcoarse_ydata,
                       modelList = natcoarse_mlist)
summary(natcoarse_gjam)
gjamPlot(natcoarse_gjam, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, SAVEPLOTS = T, 
                                         outFolder = "/Users/scarlet/Documents/suding_lab/Compost/gjam_analyses/fxnl_group_pptxnut/"))

str(natcoarse_xdata)

# try to plot interaction
xWC <- xWF <- xDC <- xDF <- natcoarse_gjam$inputs$xUnstand[1,]
xWC[c("pptTrt2w", "nutTrt1fert", "nutTrt2comp", "pptTrt2w:nutTrt2comp")]<- c(1,0,1,0)
xWC
xDC[c("pptTrt1d", "nutTrt1fert", "nutTrt2comp", "pptTrt1d:nutTrt2comp")] <- c(1,0,1,0)
testfit <- gjamIIE(natcoarse_gjam, xDC)
gjamIIEplot(testfit,response = "BackgroundExoticForb", 
            effectMu = c("main", "direct", "ind"), effectSd = c("main", "direct", "ind"))
natcoarse_xdata
