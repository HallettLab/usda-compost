# native analyses models
# author(s): ctw
# questions: caitlin.t.white@colorado.edu

# want to know..
# 1) how did seeded natives do? generally and by species seeded
# 2) is there any difference between native cover in herbicide vs. not
# 3) did nutrient addition (esp compost) + wet advantage exotics over natives?
# 4 [maybe]) Did seeding natives matter at all for changing exotic cover (e.g. suppressing non-native cover?)

# notes:
# plot 33 not seeded because not herbicided (block 4, NW plot)

# this is a nested (split-split-plot) experiment design.
# to account for all error terms would need to have error for: block, whole plot (nut), sub-plot (ppt), and maybe herbicide depending on how analyzing data
# lack of 4NW plot makes this an imbalanced design, 4 reps for all treat combos except NW (3 only)
# because of low rep size, cannot nest error as ideally would in split-plot ANOVA (not enough variation in data within treatment groups for it)

# keep block as random effect and be cognizant of bias when interpreting results (models may be overconfident in estimates)



# -- SETUP -----
# load needed libraries
library(tidyverse)
library(cowplot) # for arranging plots
library(car) # for type II and III anova tests
library(lme4)
library(lmerTest)
library(multcomp)
library(MASS)
library(pscl) # for zinf or hurdle
library(glmmTMB) # for random effects in ZINF model

# source prepped data
# > this should set default settings for datpath, plot theme, no factors by default
source("Native-analysis/native_prep_cover.R", )



# -- DID SEED + HERB TRTS WORK GENERALLY? ----
# start with simple cover (seeded natives vs. background exotics, where nfixers lumped in "forb group")
# plot to check it looks right -- warnings will be for "missing" seeded cover in unseeded subplots and that's fine, ignore
coarsecov_all_2021_simple %>%
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

# look at distribution by group
with(coarsecov_all_2021_simple, boxplot(totcov~fullgrp))


# compare background exotics in seeded vs. unseeded plots to see if any diff
# ignore nutrient and ppt trts for now but control for them in error
lm_exotics <- lmerTest::lmer(totcov ~ herbicide * fxnlgrp + seedtrt + (1|wholeplotID) +(1|ppt_trt:nut_trt), 
                             data = subset(coarsecov_all_2021_simple, grepl("Back", fullgrp)))
summary(lm_exotics)
Anova(lm_exotics) # seed addition didn't have much impact on exotic cover relative to itself in unseeded subplots, herbicide mattered more
ls_means(lm_exotics, pairwise = T)
lmerTest::ranova(lm_exotics)
# herbicide, fxnl grp all signif; seeding has no impact on amount of exotic grass and forb cover in or outside of seeded sub-sub-plots (even when include as interaction)
# generally exotic cover greater in non-herb, and grass dominates (by 70% cover on average)
# exotic forb cover is greater in herbicided, exotic grass reduced (tends to have 40% more cover in herbicided)
# model it as negbin with counts
nb_exotics <- glmer.nb(round(totcov,0) ~ herbicide * fullgrp * seedtrt + (1|wholeplotID) +(1|ppt_trt:nut_trt), 
                       data = subset(coarsecov_all_2021_simple, grepl("Back", fullgrp)))
summary(nb_exotics) # similar results
emmeans_allexos <- data.frame(emmeans::emmeans(nb_exotics, spec = c("herbicide", "fullgrp", "seedtrt")))
ggplot(emmeans_allexos, aes(herbicide, exp(emmean), group = paste(herbicide, seedtrt), col = seedtrt)) +
  geom_errorbar(aes(ymin = exp(asymp.LCL), ymax = exp(asymp.UCL)), position = position_dodge(width = 0.5), width = 0) +
  geom_point(position = position_dodge(width = 0.5)) +
  facet_wrap(~fullgrp)


# check out exotic grass cover only -- just want to know if treatments reduced dominant cover
nb_exoticgram <- glmer.nb(round(totcov,0) ~ herbicide + nut_trt + ppt_trt +  (1|block) + (1|seedtrt:subplotID), 
         data = subset(coarsecov_all_2021_simple, grepl("Exotic_G", fullgrp)))
summary(nb_exoticgram) # now yes, as expected: herb and drought reduced grass cover -- but singular fit
Anova(nb_exoticgram) # seeding and nutrient amendments do little overall; responsive to herbicide + ppt trt
emmeans_exogram <- data.frame(emmeans::emmeans(nb_exoticgram, spec = c("herbicide", "nut_trt", "ppt_trt")))
ggplot(emmeans_exogram, aes(herbicide, emmean, group = paste(herbicide, nut_trt, ppt_trt), col = paste(nut_trt, ppt_trt))) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge(width = 0.5), width = 0) +
  geom_point(position = position_dodge(width = 0.5))

# repeat for exotic forb
pois_exoticforb <- glmer(round(totcov,0) ~ herbicide + seedtrt + nut_trt + ppt_trt +  (1|block) + (1|seedtrt:subplotID),
                       family = "poisson",
                          data = subset(coarsecov_all_2021_simple, grepl("Exotic_F", fullgrp)))
summary(pois_exoticforb) # herb advantaged exotic forbs, compost reduced, drought reduced.. interesting .. perhaps bc of native grass response
Anova(pois_exoticforb) # responsive to herbicide + nutrient treat
exoforb <- data.frame(emmeans::emmeans(pois_exoticforb, spec = c("herbicide", "nut_trt", "ppt_trt")))
ggplot(exoforb, aes(herbicide, emmean, group = paste(herbicide, nut_trt, ppt_trt), col = paste(nut_trt, ppt_trt))) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge(width = 0.5), width = 0) +
  geom_point(position = position_dodge(width = 0.5))
  

# compare native seeded in herb vs. non-herb to see if any differences
lm_natives <- lmerTest::lmer(totcov ~ herbicide + fullgrp + ppt_trt + nut_trt + (1|block), 
                             data = subset(coarsecov_all_2021_simple, seedtrt == "Native seeded"))
summary(lm_natives)
Anova(lm_natives)
ls_means(lm_natives, pairwise = T)
lmerTest::ranova(lm_natives)

# round to run negbin but will not count trace amounts of natives that were there
nb_natives <- glmer.nb(round(totcov,0) ~ herbicide + fullgrp * ppt_trt + nut_trt + (1|block), 
                      data = subset(coarsecov_all_2021_simple, grepl("Native", fullgrp) & seedtrt == "Native seeded"))
# singular fit warning for negbin model
pois_natives <- glmer(round(totcov,0) ~ herbicide + nut_trt + ppt_trt + fullgrp + (1|block), 
                    family = "poisson",
                       data = subset(coarsecov_all_2021_simple, grepl("Native", fullgrp) & seedtrt == "Native seeded"))
summary(pois_natives) # herbicide helps, ppt trt also helps
# not being able to run model with random error terms and full interaction is a problem..
plot(pois_natives)
plot(nb_natives) # better
pchisq(2 * (logLik(nb_natives) - logLik(pois_natives)), df =1, lower.tail = F) #nb is better fit
summary(nb_natives)
Anova(nb_natives)
# plot predicted means
emnatives <- data.frame(emmeans::emmeans(nb_natives, spec = c("herbicide", "nut_trt", "ppt_trt", "fullgrp")))
ggplot(mutate(emnatives, ppt_trt = factor(ppt_trt, levels = c("D", "XC", "W"))), 
       aes(herbicide, exp(emmean), group = paste(herbicide, nut_trt, ppt_trt), col = ppt_trt)) +
  #geom_hline(aes(yintercept = 0), lty = 2, col = "grey") +
  geom_errorbar(aes(ymin = exp(asymp.LCL), ymax = exp(asymp.UCL)), position = position_dodge(width = 0.5), width = 0) +
  geom_point(aes(shape = nut_trt), position = position_dodge(width = 0.5)) +
  scale_color_manual(values = RColorBrewer::brewer.pal(4,"Blues")[2:4]) +
  facet_wrap(~fullgrp)
# anything below 0 yields an NaN number (interp: would not occur)


nb_natives_seeded <- glmer.nb(round(totcov) ~ herbicide * fullgrp * ppt_trt + nut_trt +(1|block), 
                                  data = subset(coarsecov_all_2021_simple, seedtrt == "Native seeded"))
summary(nb_natives_seeded)
Anova(nb_natives_seeded)
confint(nb_natives_seeded)
emnatives_all <- data.frame(emmeans::emmeans(nb_natives_seeded, spec = c("herbicide", "nut_trt", "ppt_trt", "fullgrp")))
ggplot(mutate(emnatives_all, ppt_trt = factor(ppt_trt, levels = c("D", "XC", "W"))), aes(herbicide, exp(emmean), group = paste(herbicide, nut_trt, ppt_trt), col = ppt_trt)) +
  #geom_hline(aes(yintercept = 0), lty = 2, col = "grey") +
  geom_errorbar(aes(ymin = exp(asymp.LCL), ymax = exp(asymp.UCL)), position = position_dodge(width = 0.5), width = 0) +
  geom_point(aes(shape = nut_trt, fill = ppt_trt), col = "grey30", position = position_dodge(width = 0.5)) +
  labs(y = "Predicted cover (%)", x = "Pre-seeding treatment") +
  scale_color_manual(values = RColorBrewer::brewer.pal(4,"Blues")[2:4], guide = ) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4,"Blues")[2:4]) +
  scale_shape_manual(values = c(21, 24, 22)) +
  facet_wrap(~fullgrp, scales = "free")

# try to diff groups to see if can model effect of all that way..
natonly_coarsecov_simple <- subset(coarsecov_all_2021_simple, grepl("Native", seedtrt), 
                                   select = c(plot, herbicide, fullgrp, block:nativity, totcov)) %>%
  mutate(herbicide= gsub("Non-h", "NonH", herbicide)) %>%
  spread(herbicide, totcov, fill = 0) %>%
  mutate(nh_less_h = NonHerbicided-Herbicided,
         fullgrp = factor(fullgrp),
         fullgrp =relevel(fullgrp, "Background_Exotic_Grass"))

# look at dist
ggplot(natonly_coarsecov_simple) +
  geom_density(aes(nh_less_h)) +
  facet_wrap(~fullgrp) # exo forbs and nat grams have long left tail..

# run model anyway..
lm_diff <- lmerTest::lmer(nh_less_h ~ fullgrp * nut_trt * ppt_trt + (1|block), 
                             data = natonly_coarsecov_simple)
summary(lm_diff)
Anova(lm_diff) # seed addition didn't have much impact on exotic cover relative to itself in unseeded subplots, herbicide mattered more
ls_means(lm_diff, pairwise = T)
lmerTest::ranova(lm_diff) #no..


# -- SPECIES RESPONSE MODELS ----
# build 5 separate LMMs for each native species seeded
# compare with all species in one model for due diligence

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

# grams look like they needed herbicide conditions, forbs less affected

# start with lmers on abundance
## FEMI ----
femi_lmer <- lmer(FEMI ~ herbicide * ppt_trt + nut_trt  + (1|block), data = neighborhood_counts)
femi_lmer_poisson <- glmer(FEMI ~ herbicide + ppt_trt + nut_trt + (1|block), family = "poisson", data = neighborhood_counts)
femi_lmer_nb <- glmer.nb(FEMI ~ herbicide + ppt_trt + nut_trt + (1|block), data = neighborhood_counts)
femi_lmer_nb_herbonly <- glmer.nb(FEMI ~ ppt_trt + nut_trt + (1|block), data = subset(neighborhood_counts, grepl("^Herb", herbicide)))
femi_lmer_nb_grass <- glmer.nb(FEMI ~ ExoticGrass + ppt_trt + nut_trt + (1|block), data = neighborhood_counts)
femi_lmer_bin <- glmer(FEMI ~  ExoticGrass + ppt_trt + nut_trt + (1|block), family = "binomial", data = neighborhood_PA)
femi_lmer_bin2 <- glmer(FEMI ~  grass + herbicide + ppt_trt + nut_trt + (1|hillpos), family = "binomial", data = neighborhood_PA)
femi_lmer_bin3 <- glmer(FEMI ~  herbicide + ppt_trt + nut_trt + (1|hillpos), family = "binomial", data = neighborhood_PA)
summary(femi_lmer)
summary(femi_lmer_poisson)
summary(femi_lmer_nb)
summary(femi_lmer_nb_herbonly)
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
femi_hurdle <- hurdle(FEMI ~ herbicide * ppt_trt + nut_trt, data = neighborhood_counts, dist = "negbin", link = "logit")
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
summary(all_glmer_nb3x)

# -- CLEANED UP NB GLM MODELS ----
# make hillpos factor
neighborhood_counts
# femi
femi_glm_nb <- MASS::glm.nb(FEMI ~ herbicide + ppt_trt + nut_trt + factor(block), data = neighborhood_counts)
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

# compare 3 model to model with all in it
allspp_hurdle_herb_only <- hurdle(cov ~ spp + ppt_trt * nut_trt, data = subset(neighborhood_counts_long, herbicide == "Herbicided"), dist = "negbin")
allspp_hurdle_nonherb_only <- hurdle(cov ~ spp + ppt_trt * nut_trt, data = subset(neighborhood_counts_long, herbicide != "Herbicided"), dist = "negbin")
summary(allspp_hurdle_herb_only)
summary(allspp_hurdle_nonherb_only)

# diff herb from non herb
nat_difference <- subset(neighborhood_counts_long, grepl("Herbi", herbicide), select = c(plot, seedtrt:hillpos, spp, cov)) %>%
  rename(herb_cov = cov) %>%
  left_join(subset(neighborhood_counts_long, grepl("Non", herbicide), select = c(plot, spp, cov))) %>%
  rename(nonherb_cov = cov) %>%
  mutate(herb_less_nonherb = herb_cov - nonherb_cov,
         rr_herbNH = log(nonherb_cov/herb_cov))
hist(nat_difference$rr_herbNH) # I guess it's normalish?..central mode is overepped and has slight left tail [did better in herbcided]
ggplot(nat_difference, aes(spp, rr_herbNH)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = 0), col = "red") +
  geom_point()

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

