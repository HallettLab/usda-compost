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
library(fitdistrplus) # for assessing distributions
library(car) # for type II and III anova tests
library(lme4)
library(lmerTest)
library(multcomp)
library(MASS)
library(glmmTMB) # for random effects in ZINF model
library(bbmle) # for AIC tab to compare quasi log-likelihoods
library(DHARMa)

# source prepped cover and trtkey data
# > sets default settings for datpath, plot theme, no factors by default
source("Native-analysis/native_prep_cover.R", print.eval = F)


# -- ASSESS DISTRIBUTIONS -----
# to know which type of distribution to specify in models
# > note: fxnl groups in tidy are only each fxnl group (e.g., forbs and n-fixers not summed)
plotdist(tidyfxnl0$totcov_pct)
plotdist(tidyfxnl0_seeded$totcov_pct) # about half are 0
plotdist(tidysp0$pct_cover)
plotdist(tidysp0_seeded$pct_cover) # over 80% are 0
# look at prop 0 (rule of thumb: anything above 60% 0 candidate for zero inf model)
table(tidyfxnl0$totcov_pct == 0)/nrow(tidyfxnl0) # 50 % are 0, 50% >0
table(tidyfxnl0_seeded$totcov_pct == 0)/nrow(tidyfxnl0_seeded) # 42% are 0, 58% not 0 -- probably don't need zero inf
table(tidysp0_seeded$pct_cover == 0)/nrow(tidysp0_seeded) # 76% are 0, 24% not 0 -- possible use zero inf

# by fxnl groups: using pct cover and count
descdist(tidyfxnl0$totcov_pct, boot = 10000) # beta
descdist(tidyfxnl0$totcov_count, boot = 10000, discrete = T) # poisson is closest
# check native seeded only
descdist(tidyfxnl0$totcov_pct[tidyfxnl0$seedtrt != "Unseeded"], boot = 10000) #beta
descdist(tidyfxnl0$totcov_count[tidyfxnl0$seedtrt != "Unseeded"], boot = 10000, discrete = T) # poisson

# species level
descdist(tidysp0$pct_cover, boot = 10000) # beta dist
descdist(tidysp0$count_cover, boot = 10000, discrete = T) # negative binomial

# look at native seeded spp
descdist(widesp_count$FEMI, boot = 10000, discrete = T) # negbin
descdist(widesp_count$FEMI[widesp_count$seedtrt == "Native seeded"], boot = 10000, discrete = T) # close to poisson in seeded only plots, but still negbin
descdist(widesp_count$ESCA, boot = 10000, discrete = T) # also close to poisson, but can in negbin space
descdist(widesp_count$ESCA[widesp_count$seedtrt == "Native seeded"], boot = 10000, discrete = T) 
descdist(widesp_count$BRCA, boot = 10000, discrete = T) # negbin
descdist(widesp_count$BRCA[widesp_count$seedtrt == "Native seeded"], boot = 10000, discrete = T) # negbin
descdist(widesp_count$NEMA, boot = 10000, discrete = T) # negbin
descdist(widesp_count$NEMA[widesp_count$seedtrt == "Native seeded"], boot = 10000, discrete = T) # negbin

# check abundant spp
## AVBA
descdist(round(widesp_count$AVBA), boot = 10000, discrete = T) # poisson?
hist(widesp_count$AVBA, breaks = 100) # not normal
descdist(round(widesp_pct$AVBA), boot = 10000, discrete = F) # beta
## TACA
descdist(round(widesp_count$TACA), boot = 10000, discrete = T) # poisson? like AVBA
hist(widesp_count$TACA, breaks = 100) # not normal
## ERBO
descdist(round(widesp_count$ERBO), boot = 10000, discrete = T) # negbin
descdist(seedtrt_ydata$FEMI, boot = 10000)
fit_fxnlp <- fitdist(tidyfxnl0$totcov_count, method = "mle", discrete = T, distr = "pois")
summary(fit_fxnlp)
fit_fxnlnb <- fitdist(tidyfxnl0$totcov_count, method = "mle", discrete = T, distr = "nbinom")
summary(fit_fxnlnb)


# check mean vs variance
mean(tidyfxnl0$totcov_count) # 13.32545
var(tidyfxnl0$totcov_count) # 911.3632


# -- PREP PCT COVER FOR BETA FAM ----
# needs to be between 0 and 1
tidyfxnl0$totcov_prop <- (tidyfxnl0$totcov_pct)/100
tidyfxnl0_seeded$totcov_prop <- tidyfxnl0_seeded$totcov_pct/100

# beta dist needs values between 0 and 1. some fxnl groups sum to over 100% over season
# add small amount to 0 and subtract small amount from max value
tidyfxnl0$totcov_scale01 <- with(tidyfxnl0, ifelse(totcov_pct == 0, 0, totcov_pct/max(totcov_pct)))
# check values before adjust min and max
head(sort(unique(tidyfxnl0$totcov_scale01))) # smallest value is at e-5 scale; largest is at 0.922 before 1
tidyfxnl0$totcov_scale01beta <- tidyfxnl0$totcov_scale01
tidyfxnl0$totcov_scale01beta[tidyfxnl0$totcov_scale01==0] <- 1e-7
tidyfxnl0$totcov_scale01beta[tidyfxnl0$totcov_scale01==1] <- 1-0.01
range(tidyfxnl0$totcov_scale01); range(tidyfxnl0$totcov_scale01beta)
#tidyfxnl0$totcov_scale01t <- with(tidyfxnl0, ifelse(totcov_pct == 0, 0.0000000001, totcov_scale01))
plot(tidyfxnl0$totcov_pct ~ tidyfxnl0$totcov_scale01)

fit_fxnlbeta <- fitdist(tidyfxnl0$totcov_scale01beta, distr = "beta")
summary(fit_fxnlbeta)

# plot diagnostic fits
par(mfrow = c(2,2))
denscomp(fit_fxnlbeta) # not sure how to interp this
cdfcomp(fit_fxnlbeta) # fit looks good for beta
ppcomp(fit_fxnlbeta) # more 0s than expected, reason to use negbin
qqcomp(fit_fxnlbeta)

# negbin fit on 'count' cover
denscomp(fit_fxnlnb) # nb seems like a better fit than poisson
cdfcomp(fit_fxnlnb) # not a perfect fit but better than poisson 
ppcomp(fit_fxnlnb)
qqcomp(fit_fxnlnb)

# poisson
denscomp(fit_fxnlp)
cdfcomp(fit_fxnlp) # fit looks good for beta
ppcomp(fit_fxnlbeta) # more 0s than expected, reason to use negbin
qqcomp(fit_fxnlp)


# beta is a an okay distribution if want to use pct cover; otherwise use negbin


# -- DID SEED + HERB TRTS WORK GENERALLY? ----
# screen fxnl groups for low occurrence (to drop from analysis)
sort(colSums(widefxnl_pa[grepl("Exotic|Native", names(widefxnl_pa))])/nrow(widefxnl_pa))
# can drop native n-fixers (seeded or native)

# start with simple cover (seeded natives vs. background exotics, where nfixers lumped in "forb group")
# plot to check it looks right -- warnings will be for "missing" seeded cover in unseeded subplots and that's fine, ignore
tidyfxnl0 %>%
  #mutate(totcov = ifelse(grepl("Native", fullgrp) & seedtrt == "Unseeded", NA, totcov)) %>%
  ggplot(aes(herbicide, totcov_pct, fill = coarse_fxnl)) + #interaction(fullgrp, seedtrt, herbicide)
  geom_boxplot(aes(group = interaction(coarse_fxnl, herbicide)), varwidth = T, position = position_dodge(width = 0.5), alpha = 0.5) +
  geom_point(aes(col = coarse_fxnl), position = position_dodge(width = 0.5)) +
  # scale_fill_manual(values = c("orchid", "seagreen1", "purple4", "seagreen4")) +
  # scale_color_manual(values = c("orchid", "seagreen1", "purple4", "seagreen4")) +
  #facet_grid(seedtrt ~ coarse_fxnl) +
  facet_wrap(~ nut_trt + seedtrt +ppt_trt, ncol = 6) +
  #facet_grid(coarse_fxnl ~ seedtrt) #+
  theme(#axis.text.x = element_text(angle = 90))
        legend.position = "top",
        legend.direction = "horizontal",
        strip.background = element_rect(fill = "transparent"))

# look at distribution by group
with(tidyfxnl0_seeded, boxplot(totcov_pct~coarse_fxnl))

#tidyfxnl0$plot <- factor(tidyfxnl0$plot)
# relevel coarse_fxnl so background exotic grams are reference (otherwise herbicide has a positive effect because exo forbs are the reference)

tidyfxnl0$coarse_fxnl <- factor(tidyfxnl0$coarse_fxnl) 
tidyfxnl0$coarse_fxnl <- relevel(tidyfxnl0$coarse_fxnl, ref = "Background Exotic Grass") 
levels(tidyfxnl0$coarse_fxnl)
# compare background exotics in seeded vs. unseeded plots to see if any diff
# ignore nutrient and ppt trts for now but control for them in error
glm_beta_exotics <- glmmTMB(totcov_scale01 ~ seedtrt + herbicide * coarse_fxnl + (1|block/wholeplotID/subplotID), 
                      family = ordbeta,
                      data = subset(tidyfxnl0, grepl("Background Exo", coarse_fxnl))) # select = c(plot, block, coarse_fxnl, totcov_scale01, herbicide)
summary(glm_beta_exotics)

glm_nb1_exotics <- glmmTMB(totcov_count ~ seedtrt + herbicide * coarse_fxnl +(1|block/wholeplotID/subplotID), 
                      family = nbinom1(),
                      data = subset(tidyfxnl0, grepl("Background Exo", coarse_fxnl)))
summary(glm_nb1_exotics)

glm_nb2_exotics <- glmmTMB(totcov_count ~ seedtrt + herbicide * coarse_fxnl +(1|block/wholeplotID/subplotID), 
                           family = nbinom2,
                           data = subset(tidyfxnl0, grepl("Background Exo", coarse_fxnl)))
summary(glm_nb2_exotics) # theta is better for nbinom2

# run on all fxnl groups
# using pct rescaled to 0<y<1
glm_beta_allfxnl <- glmmTMB(totcov_scale01 ~ seedtrt + herbicide * coarse_fxnl + (1|block/wholeplotID/subplotID), 
                            family = ordbeta, # beta_family
                            # drop seeded N fixer because it barely came up
                            data = subset(tidyfxnl0, !grepl("Native N-fix", coarse_fxnl))) # select = c(plot, block, coarse_fxnl, totcov_scale01, herbicide)
summary(glm_beta_allfxnl)

# using count with poisson
glm_nb2_allfxnl <- glmmTMB(totcov_count ~ seedtrt + herbicide * coarse_fxnl + (1|block/wholeplotID/subplotID), 
                            family = nbinom2,
                            # drop seeded N fixer because it barely came up
                            data = subset(tidyfxnl0, !grepl("Native N-fix", coarse_fxnl)))
summary(glm_nb2_allfxnl)

# update to hurdle model to see if zero inflation model is better
# > any 0 in dataset is true (don't think ashley and carmen did missed species or were there are the wrong time, no false zeroes)
glm_nb2_allfxnl_hurdle <- glmmTMB(totcov_count ~ seedtrt + herbicide * coarse_fxnl + (1|block/wholeplotID/subplotID), 
                                  family = truncated_nbinom2,
                                  # seedtrt influences whether natives there or not, herbicide should also influence finding natives
                                  zi = ~ seedtrt + herbicide,
                                  # drop seeded N fixer because it barely came up
                                  data = subset(tidyfxnl0, !grepl("Native N-fix", coarse_fxnl)))
summary(glm_nb2_allfxnl_hurdle)
glm_nb2_allfxnl_zinb <- glmmTMB(totcov_count ~ seedtrt + herbicide * coarse_fxnl + (1|block/wholeplotID/subplotID), 
                                  family = nbinom2,
                                  # seedtrt influences whether natives there or not, herbicide should also influence finding natives
                                  zi = ~ seedtrt + herbicide,
                                  # drop seeded N fixer because it barely came up
                                  data = subset(tidyfxnl0, !grepl("Native N-fix", coarse_fxnl)))
summary(glm_nb2_allfxnl_zinb)

# test poisson with hurdle
glm_pois_allfxnl_hurdle <- glmmTMB(totcov_count ~ seedtrt + herbicide * coarse_fxnl + (1|block/wholeplotID/subplotID), 
                                family = truncated_poisson(),
                                # seedtrt influences whether natives there or not, herbicide should also influence finding natives
                                zi = ~ seedtrt + herbicide,
                                # drop seeded N fixer because it barely came up
                                data = subset(tidyfxnl0, !grepl("Native N-fix", coarse_fxnl)))
summary(glm_pois_allfxnl_hurdle)

# compare hurdle to non-hurdle negbin
AICtab(glm_nb2_allfxnl, glm_nb2_allfxnl_hurdle,glm_pois_allfxnl_hurdle, glm_nb2_allfxnl_zinb, base = T, logLik = T) # poisson definitely not the right fit
# hurdle is worst.. zinb is best (but model doesn't fit data in concept; only slightly better than non-zero inflated model)

# more model eval with dharma, compare beta v. nb
fxnl_nb2_hurdle_simres <- simulateResiduals(glm_nb2_allfxnl_hurdle)
plot(fxnl_nb2_hurdle_simres)
testOutliers(fxnl_nb2_hurdle_simres, type = "bootstrap")

# zinb model
fxnl_nb2_zinb_simres <- simulateResiduals(glm_nb2_allfxnl_zinb)
plot(fxnl_nb2_zinb_simres)
testOutliers(fxnl_nb2_zinb_simres, type = "bootstrap")
plotResiduals(fxnl_nb2_zinb_simres)
testZeroInflation(fxnl_nb2_zinb_simres)

# negative binomial model
fxnl_nb2_simres <- simulateResiduals(glm_nb2_allfxnl)
plot(fxnl_nb2_simres)
plotResiduals(fxnl_nb2_simres)
testOutliers(fxnl_nb2_simres, type = "bootstrap")
testZeroInflation(fxnl_nb2_simres)



fxnl_beta <- simulateResiduals(glm_beta_allfxnl) # ordered beta distribution [y can be 0 or 1] does not look quite as bad as beta_family
plot(fxnl_beta)



# -- TEST EFFECTS ON SEEDED SPP -----
widesp_fxnlneighbors_count <- left_join(widefxnl_count, widesp_count[c("plot", "seedtrt", "herbicide", nats)])
widesp_fxnlneighbors_pct <- left_join(widefxnl_pct, widesp_pct[c("plot", "seedtrt", "herbicide", nats)])
# make proportional for ordered beta dist
widesp_fxnlneighbors_pct$FEMI_prop <- widesp_fxnlneighbors_pct$FEMI/100
widesp_fxnlneighbors_pct$BRCA_prop <- widesp_fxnlneighbors_pct$BRCA/100
widesp_fxnlneighbors_pct$NEMA_prop <- widesp_fxnlneighbors_pct$NEMA/100
widesp_fxnlneighbors_pct$ESCA_prop <- widesp_fxnlneighbors_pct$ESCA/100

# use neighbor grams and neighbor forbs as predictors
# just within seeded plots
# look at effects of env on species first
femi_nb2_herbicide <- glmmTMB(FEMI ~ herbicide + nut_trt * ppt_trt, family = nbinom2() ,data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
femi_nb2_allneighbors <- glmmTMB(FEMI ~ herbicide + AllNeighbors + nut_trt * ppt_trt, 
                                 family = nbinom2(),
                                 data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
femi_nb2_exoneighbors <- glmmTMB(FEMI ~ herbicide + ExoticGrass + ExoForbNfix + nut_trt * ppt_trt, 
                                 family = nbinom2(),
                                 data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
femi_nb2_reblock <- glmmTMB(FEMI ~ herbicide + nut_trt * ppt_trt + (1|block), 
                            family = nbinom2(),
                            data = subset(widesp_count, seedtrt != "Unseeded"))
femi_nb2_reblockwp <- glmmTMB(FEMI ~ herbicide + nut_trt * ppt_trt + (1|block/wholeplotID), 
                            family = nbinom2(),
                            data = subset(widesp_count, seedtrt != "Unseeded"))
femi_nb2_refullnest <- glmmTMB(FEMI ~ herbicide + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                              family = nbinom2(),
                              data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
femi_nb2_exograss_refullnest <- glmmTMB(FEMI ~ herbicide + ExoticGrass + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                               family = nbinom2(),
                               data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
femi_nb2_exoneighbors_refullnest <- glmmTMB(FEMI ~ herbicide + ExoticGrass + ExoForbNfix + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                                        family = nbinom2(),
                                        data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
summary(femi_nb2_herbicide)
summary(femi_nb2_allneighbors)
summary(femi_nb2_exoneighbors)
summary(femi_nb2_reblock)
summary(femi_nb2_reblockwp)
summary(femi_nb2_refullnest)
summary(femi_nb2_exograss_refullnest)
summary(femi_nb2_exoneighbors_refullnest)

AICtab(femi_nb2_herbicide, femi_nb2_allneighbors, femi_nb2_exoneighbors,
       femi_nb2_reblock, femi_nb2_reblockwp, femi_nb2_refullnest, 
       femi_nb2_exograss_refullnest, femi_nb2_exoneighbors_refullnest,
       base =T, logLik = T)


# test using pct cover
femiprop_beta_refullnest <- glmmTMB(FEMI_prop ~ herbicide + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                               family = ordbeta,
                               data = subset(widesp_fxnlneighbors_pct, seedtrt != "Unseeded"))
femiprop_beta_reblock <- glmmTMB(FEMI_prop ~ herbicide + nut_trt * ppt_trt + (1|block), 
                                    family = ordbeta,
                                    data = subset(widesp_fxnlneighbors_pct, seedtrt != "Unseeded"))
femiprop_beta_reblockwp <- glmmTMB(FEMI_prop ~ herbicide + nut_trt * ppt_trt + (1|block/wholeplotID), 
                                 family = ordbeta,
                                 data = subset(widesp_fxnlneighbors_pct, seedtrt != "Unseeded"))
summary(femiprop_beta_refullnest) # can't use full nesting
summary(femiprop_beta_reblock)
summary(femiprop_beta_reblockwp) # effects are about the same
# FEMI doesn't care about neighbors. it just needs a clean slate to start (if had litter pre-growth that might be helpful).

AICtab(femiprop_beta_reblock, femiprop_beta_reblockwp, base = T, logLik = T)

femi_beta_simres <- simulateResiduals(femiprop_beta_reblockwp) # ordered beta distribution [y can be 0 or 1] does not look quite as bad as beta_family
plot(femi_beta_simres) # no issues
femi_nb2_refullnest_simres <- simulateResiduals(femiprop_beta_reblockwp) # ordered beta distribution [y can be 0 or 1] does not look quite as bad as beta_family
plot(femi_nb2_refullnest_simres) # no issues here either


# -- compare with BRCA: negbin and beta -----  
brca_nb2_refullnest <- glmmTMB(BRCA ~ herbicide + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                               family = nbinom2(),
                               data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
brca_nb2_exogram_refullnest <- glmmTMB(BRCA ~ herbicide + ExoticGrass + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                               family = nbinom2(),
                               data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
brca_nb2_exogf_refullnest <- glmmTMB(BRCA ~ herbicide + ExoticGrass + ExoticForb + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                                       family = nbinom2(),
                                       data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
brca_nb2_exogfn_refullnest <- glmmTMB(BRCA ~ herbicide + ExoticGrass + ExoForbNfix + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                                     family = nbinom2(),
                                     data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
summary(brca_nb2_refullnest)
summary(brca_nb2_exogram_refullnest) # if use all neighbors does not produce
summary(brca_nb2_exogf_refullnest)
summary(brca_nb2_exogfn_refullnest)
AICtab(brca_nb2_refullnest, brca_nb2_exogram_refullnest, brca_nb2_exogf_refullnest, brca_nb2_exogfn_refullnest, base = T)

brca_beta_refullnest <- glmmTMB(BRCA_prop ~ herbicide + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                               family = ordbeta,
                               data = subset(widesp_fxnlneighbors_pct, seedtrt != "Unseeded"))
brca_beta_exofn_refullnest <- glmmTMB(BRCA_prop ~ herbicide + ExoForbNfix + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                                family = ordbeta,
                                data = subset(widesp_fxnlneighbors_pct, seedtrt != "Unseeded"))
brca_beta_exogfn_refullnest <- glmmTMB(BRCA_prop ~ herbicide + ExoticGrass + ExoForbNfix + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                                      family = ordbeta,
                                      data = subset(widesp_fxnlneighbors_pct, seedtrt != "Unseeded"))
summary(brca_beta_refullnest)
summary(brca_beta_exofn_refullnest)
summary(brca_beta_exogfn_refullnest)

AICtab(brca_beta_refullnest, brca_beta_exofn_refullnest, brca_beta_exogfn_refullnest, base = T, logLik = T)
ppcomp(fitdist(widesp_fxnlneighbors_count$BRCA, distr = "nbinom"))
cdfcomp(fitdist(widesp_fxnlneighbors_pct$BRCA_prop, method = "mme", distr = "beta"), plotstyle = "ggplot") # beta looks right but there are many 0s
cdfcomp(fitdist(widesp_fxnlneighbors_pct$FEMI_prop, method = "mme", distr = "beta"), plotstyle = "ggplot")
cdfcomp(fitdist(widesp_fxnlneighbors_count$FEMI, distr = "negbin"), plotstyle = "ggplot")

# run brca with hurdle
brca_nb2_refullnest_hurdle <-  glmmTMB(BRCA ~ herbicide + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                                       ziformula = ~ herbicide,
                                       family = truncated_nbinom2(),
                                       data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
summary(brca_nb2_refullnest_hurdle)

brca_nb2_exogfn_refullnest_hurdle <- glmmTMB(BRCA ~ herbicide + ExoticGrass + ExoForbNfix + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                                       ziformula = ~ herbicide,
                                       family = truncated_nbinom2(),
                                       data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
summary(brca_nb2_exogfn_refullnest_hurdle)

brca_nb2_refullnest_zinb <- glmmTMB(BRCA ~ herbicide + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                                             ziformula = ~ herbicide,
                                             family = nbinom2(),
                                             data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
summary(brca_nb2_refullnest_zinb) # doesn't run with exotic grams and forbs

AICtab(brca_nb2_refullnest, brca_nb2_exogfn_refullnest, 
       brca_nb2_refullnest_hurdle, brca_nb2_refullnest_zinb, brca_nb2_exogfn_refullnest_hurdle,
       base = T, logLik = T) # hurdle model is worse.. but not zinb

# ESCA ----
esca_nb2_refullnest <- glmmTMB(ESCA ~ herbicide + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                               family = nbinom2(),
                               data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
esca_nb2_exogfn_refullnest_zinb <- glmmTMB(ESCA ~ herbicide + ExoticForb + nut_trt * ppt_trt + (1|block/wholeplotID), 
                                      ziformula = ~ herbicide,
                                      family = nbinom2(),
                                      data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
esca_pois_exogfn_refullnest_zinb <- glmmTMB(ESCA ~ herbicide + ExoticGrass + ExoForbNfix + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                                           ziformula = ~ herbicide,
                                           family = poisson,
                                           data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
esca_beta_refullnest <- glmmTMB(ESCA_prop ~ herbicide + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                               family = ordbeta,
                               data = subset(widesp_fxnlneighbors_pct, seedtrt != "Unseeded"))
summary(esca_nb2_refullnest)
summary(esca_nb2_exogfn_refullnest_zinb)
summary(esca_pois_exogfn_refullnest_zinb)
summary(esca_beta_refullnest)

nema_nb2_refullnest <- glmmTMB(NEMA ~ herbicide + nut_trt * ppt_trt + (1|block/wholeplotID), # won't run with subplot included
                               family = nbinom2(),
                               data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
nema_nb2_refullnest_zinb <- glmmTMB(NEMA ~ herbicide + nut_trt * ppt_trt + block + (1|wholeplotID/subplotID),
                               ziformula = ~ herbicide,
                               family = nbinom2(),
                               data = subset(widesp_fxnlneighbors_count, seedtrt != "Unseeded"))
nema_beta_refullnest <- glmmTMB(NEMA_prop ~ herbicide + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), # won't run with subplot included
                               family = ordbeta,
                               data = subset(widesp_fxnlneighbors_pct, seedtrt != "Unseeded"))
summary(nema_beta_refullnest) # runs, but won't calculate st errors
summary(nema_nb2_refullnest)
summary(nema_nb2_refullnest_zinb) # runs, won't calculate st errors


# -- TRY ALL SPECIES MODEL ----
# don't anticipate this will work with full nested error structure, but see what can be modeled
tidynats0_seeded <- subset(tidysp0_seeded, code4 %in% nats[nats != "TRCI"]) %>%
  left_join(subset(widefxnl_pct, select = c(plot:ppt_trt, seedtrt:ExoticNfixer, ExoForbNfix:AllNeighbors))) %>%
  mutate(code4 = factor(code4, levels = c("TRCI", "ESCA", "NEMA", "BRCA", "FEMI")))
                    
# use tidysp0, subsetting to just nats seeded. join neighbors of interest as columns. start with no neighbors in model to begin

allsp_nb <- glmmTMB(formula = count_cover ~ herbicide + code4 + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                    family = nbinom2,
                    # subset to just native spp seeded, order levels so low abundance forbs are reference [otherwise it's BRCA]
                    data = subset(mutate(tidysp0_seeded, code4 = factor(code4, levels = c("TRCI", "ESCA", "NEMA", "BRCA", "FEMI"))),
                                         code4 %in% nats[nats != "TRCI"]))

# ran with code4 as additive term; TRCI is reference (which is basically 0)
summary(allsp_nb)

# three-way interaction between sp code and enviro trts
allspx_nb <- glmmTMB(formula = count_cover ~ herbicide + code4 * nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                    family = nbinom2,
                    data = subset(mutate(tidysp0_seeded, code4 = factor(code4, levels = c("TRCI", "ESCA", "NEMA", "BRCA", "FEMI"))),
                                  code4 %in% nats[nats != "TRCI"]))

# interaction model runs.. not much is signif
summary(allspx_nb)

# three-way interaction between sp code and enviro trts
allspx_nb2_allneigh <- glmmTMB(formula = count_cover ~ herbicide + AllNeighbors + code4 * nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                       family = nbinom2,
                       data = tidynats0_seeded)

summary(allspx_nb2_allneigh)

# three way interaction with grams specifically
allspx_nb2_exog <- glmmTMB(formula = count_cover ~ herbicide + ExoticGrass + code4 * nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                               family = nbinom2,
                               data = tidynats0_seeded)

summary(allspx_nb2_exog)

# grams + forb neighbors
# note: too many terms in model with,
allspx_nb2_exofn <- glmmTMB(formula = count_cover ~ herbicide + ExoForbNfix + code4 * nut_trt * ppt_trt + (1|block/wholeplotID/subplotID), 
                           family = nbinom2,
                           data = tidynats0_seeded)

summary(allspx_nb2_exofn) 

# try model without herbicide, but with neighbors (in this case, herbcide should be part of nested error, but probably won't fit)
allspx_nb2_fullnest <- glmmTMB(formula = count_cover ~ code4 * nut_trt * ppt_trt + (1|block/wholeplotID/subplotID/herbicide), 
                            family = nbinom2,
                            data = tidynats0_seeded)

summary(allspx_nb2_fullnest) 

allspx_nb2_fullnest_exoneigh <- update(allspx_nb2_fullnest, update.formula(count_cover ~ code4 * nut_trt * ppt_trt + (1|block/wholeplotID/subplotID/herbicide),
                                                                           . ~ ExoticNeighbors + .))
summary(allspx_nb2_fullnest_exoneigh)

allspx_nb2_fullnest_exogfn <- update(allspx_nb2_fullnest, update.formula(count_cover ~ code4 * nut_trt * ppt_trt + (1|block/wholeplotID/subplotID/herbicide),
                                                                           . ~ ExoticGrass + ExoForbNfix + code4 * nut_trt * ppt_trt + (1|block/wholeplotID/subplotID/herbicide)))

summary(allspx_nb2_fullnest_exogfn) 

AICtab(allsp_nb, allspx_nb, allspx_nb2_allneigh, allspx_nb2_exog, allspx_nb2_exofn, 
       allspx_nb2_fullnest, allspx_nb2_fullnest_exoneigh, allspx_nb2_fullnest_exogfn,
       base = T, logLik = T)
# like others above, best model is the one without any neighbors, but next best are any of the models with any neighbor (agnostic to who)
# > what i'm taking from this is neighbor ID really does not matter that much other than herbicide knocks back dominant spp competition (in this case, non-native grass)
# disperion is close to 1 for all, so less need for zero-inflated model

# one more: interact native spp with grams instead
allspx_nb2_fullnest_codexneighbor <- update(allspx_nb2_fullnest, update.formula(count_cover ~ code4 * nut_trt * ppt_trt + (1|block/wholeplotID/subplotID/herbicide),
                                                                         . ~ ExoticNeighbors * code4 + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID/herbicide)))

summary(allspx_nb2_fullnest_codexneighbor) 

allspx_nb2_fullnest_codexgram <- update(allspx_nb2_fullnest, update.formula(count_cover ~ code4 * nut_trt * ppt_trt + (1|block/wholeplotID/subplotID/herbicide),
                                                                                . ~ ExoticGrass * code4 + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID/herbicide)))

summary(allspx_nb2_fullnest_codexgram) 

AICtab(allsp_nb, allspx_nb, allspx_nb2_allneigh, allspx_nb2_exog, allspx_nb2_exofn, 
       allspx_nb2_fullnest, allspx_nb2_fullnest_codexneighbor, allspx_nb2_fullnest_exoneigh,
       allspx_nb2_fullnest_codexgram, allspx_nb2_fullnest_exogfn,
       base = T, logLik = T)

allspx_nb2_simres <- simulateResiduals(allspx_nb)
testDispersion(allspx_nb2_simres) 
testResiduals(allspx_nb2_simres) # simulated residual values are not more extreme than expected (dispersion okay for model used)
plotResiduals(allspx_nb2_simres) # surprisingly looks ok!

# test model with code interacting with neighbors
allspx_nb2_codexneigh_simres <- simulateResiduals(allspx_nb2_fullnest_codexneighbor)
testResiduals(allspx_nb2_codexneigh_simres)
plotResiduals(allspx_nb2_codexneigh_simres) # slightly better predictions, but has outlier (but is also bootstrapped)


# -- TEST LITTER ----
# if neighbor ID does not matter, is it just litter (whether there is opporunity to emerge and grow?)
# compile april visit mean litter depth, just in seeded; also if any gopher holes disturbed plots

natlitter <- group_by(natlong, plot, seedtrt, herbicide) %>%
  mutate(gopher = any(grepl("gopher", notes, ignore.case = T))) %>%
  ungroup() %>%
  subset(mon == 4, select = c(plot, seedtrt, herbicide, gopher, litter_depth_cm)) %>%
  distinct()

# test for influence of litter and/or gopher disturbance (only impacted 11 plots)
tidynats0_seeded <- left_join(tidynats0_seeded, natlitter)

# leave out neighbors
allspx_nb2_littergoph <- update(allspx_nb, 
                                update.formula(count_cover ~ herbicide + code4 * nut_trt *ppt_trt + (1 | block/wholeplotID/subplotID),
                                               . ~ litter_depth_cm + gopher + .), data = tidynats0_seeded)
summary(allspx_nb2_littergoph)  
  
allspx_nb2_litter_codexneighbors <- update(allspx_nb2_fullnest_codexneighbor, 
                                update.formula(allspx_nb2_fullnest_codexneighbor$modelInfo$allForm$formula,
                                               . ~ litter_depth_cm + .), data = tidynats0_seeded)
summary(allspx_nb2_litter_codexneighbors)
AICtab(allspx_nb, allspx_nb2_fullnest_codexneighbor, allspx_nb2_littergoph, allspx_nb2_litter_codexneighbors, base = T, logLik = T)

# allspx_nb is still most parsimonious model, with code4 * exoneighbors is slightly better


# re-run with beta dist to confirm results, totcov_pct/max(totcov_pct)
tidynats0_seeded$cover_scale01 <- with(tidynats0_seeded, (pct_cover - min(pct_cover))/(max(pct_cover)- min(pct_cover)))
allspx_ordbeta <- glmmTMB(formula = cover_scale01 ~ herbicide + code4 * nut_trt * ppt_trt + (1|block/wholeplotID/subplotID),
                         family = ordbeta,
                         data = tidynats0_seeded)
summary(allspx_ordbeta) # same things are signif

# test with code 4 * exo neighbors interaction
allspx_ordbeta_codexneighbors <- glmmTMB(formula = cover_scale01 ~ ExoticNeighbors * code4 + nut_trt * ppt_trt + (1|block/wholeplotID/subplotID/herbicide),
                          family = ordbeta,
                          data = tidynats0_seeded) # doesn't like it
summary(allspx_ordbeta_codeneighbors)


# ## -- OLD CODE -----
Anova(glm_exotics) # seed addition didn't have much impact on exotic cover relative to itself in unseeded subplots, herbicide mattered more
lmerTest::ls_means(glm_poisson_exotics, pairwise = T)
lmerTest::ranova(glm_exotics)
# herbicide, fxnl grp all signif; seeding has no impact on amount of exotic grass and forb cover in or outside of seeded sub-sub-plots (even when include as interaction)
# generally exotic cover greater in non-herb, and grass dominates (by 70% cover on average)
# exotic forb cover is greater in herbicided, exotic grass reduced (tends to have 40% more cover in herbicided)
# model it as negbin with counts
nb_exotics <- glmer.nb(totcov_count ~ herbicide + coarse_fxnl + seedtrt + (1|block), 
                       data = subset(tidyfxnl0, grepl("Background Ex", coarse_fxnl)))
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

