# native multivariate analyses
# author(s): ctw
# questions: caitlin.t.white@colorado.edu


# notes from gavin simpson tutorials (most are verbatim):

# 'lc' scores are like the regression coefficients in a constrained ordination (closest to regression/fitted scores)
# > lc reflects purely fitted values of the environment, whereas sites scores reflect variation between sites given the species composition
# correlation = T helps standardize scores if you did not scale response variables prior to running ordination model
# hills number gives correlations for cca, correlation gives correlation-based scores for rda
# hills number also tells you effective sample size from renyi entropy (is exponentiated renyi)
# if you don't scale response variable prior to model, it reflect covariances more, and so largest abundance spp will dominate/drive the results

# partial constrained ordinations remove the effect of variables entered in z argument (or Condition(vars)),
# then fits constrained ordination on the unconstrained variation leftover
# 'partial' out nuisance variables [like block]
# vegan doesn't consider partialed 'bits' in calculating eigenvalues because it first removes their influence before fitting a cca on whatever vars entered as constraints
# VIF = how much variance of some factors get inflated by presence of correlated variables
# ^things can get large effects and look important because they are just trying to offset variables they are correlated with 

# with permutation, we're curious if F is large compared to null expectation F
# with anova testing signif of axes, partials out variation explained by axis 1, then tests signif of next axis, partials it out, repeat
# ^ helpful to show which axes are good to report/show in results
# by = margin gives partial effect (how much extra info is explained by addition of a term)
# ^ what is the additional effect given all the other terms present

# other notes: 
# on hellinger transformation (online def from 'hellinger {labdsv}' R fxn [I'm using the fxn in vegan tho, is same):
# Calculates a sample total standardization (all values in a row are divided by the row sum), and then takes the square root of the values.
# Hellinger standardization is a convex standardization that simultaneously helps minimize effects of vastly different sample total abundances.



# -- SETUP -----
# load needed libraries
library(tidyverse)
library(vegan)
library(FD)

# source prepped data
# > this should set default settings for datpath, plot theme, no factors by default
#source("Native-analysis/native_prep_cover.R", )
source("Native-analysis/native_prep_traits.R", print.eval = F, echo = F) # this script loads veg data needed
# ^note: warning message about NA coercian is okay. just converting notes about no sample available to NA, rest of numeric vals become numberic

# remove things not needed for simple enviro
rm(list = ls()[grepl("^asla|aslatr|aslta", ls())])


# -- PREP FXNL DATA FOR ORDINATIONS -----
# separate plant dat from plot dat
wide_env <- subset(widefxnl_count, select = c(plot:yr)) %>%
  # add hillpos
  mutate(hillpos = ifelse(block %in% c(1,2), "downhill", "uphill"),
         hillpos = factor(hillpos)) %>%
  # make rowids: plot + seedtrt + herbicide
  mutate(rowid = paste(plot, seedtrt, herbicide, sep = "_"),
         rowid = gsub("herbi.*", "Herb", rowid, ignore.case = T),
         rowid = gsub("Non-|Un", "No", rowid, ignore.case = T),
         rowid = gsub("seeded", "Seed", rowid, ignore.case = T),
         rowid = gsub("Native ", "Nat", rowid, ignore.case = T))
rownames(wide_env) <- wide_env$rowid
# make ypos and xpos factors
wide_env$xpos <- factor(wide_env$xpos)
wide_env$ypos <- factor(wide_env$ypos)

# cca by default weights by rowsums and colsums -- for this reason, should probably include bare ground for herbicided plots. take may value.
# determine pct bare (relevant for droughted or herbicided plots)
amb_groundcov21 <- subset(coverlong, yr == 2021 & sample_event == 2, select = c(plot, pct_bare)) %>%
  distinct() %>%
  mutate(herbicide = "Non-herbicided", seedtrt = "Unseeded")
nat_groundcov21 <- subset(natlong, select = c(plot, mon, seedtrt, herbicide, pct_bare))  %>%
  distinct() %>%
  group_by(plot, seedtrt, herbicide) %>%
  mutate(pct_bare = ifelse(length(mon) == 1, pct_bare, pct_bare[mon == 5])) %>%
  # drop mon
  select(-mon) %>%
  distinct()
# row-bind ambient bare cover with seeded/herbicided subplots
may_pctbare <- rbind(nat_groundcov21, amb_groundcov21[names(nat_groundcov21)]) %>%
  # join enviro data so can assign rowid
  left_join(wide_env[c("plot", "seedtrt", "herbicide", "rowid")])
  
# keep individual fxnl groups
# check frequency of nfixers (if should lump with forbs)
colSums(subset(widefxnl_pa, select = -c(plot:yr)))
# split by hebicide and seedtrt
sapply(split(subset(widefxnl_pa, select = -c(plot:yr)), paste(widefxnl_pa$seedtrt, widefxnl_pa$herbicide)), colSums)
# can drop native nfixer from all (or combine w native forb), only keep native grass for unseeded plots

# keep seeded grps so easy to use for ordinations when when to keep (drop in command for ambient analyses)
widefxnl_pct_ordgrps <- left_join(widefxnl_pct, may_pctbare)
# assign rowid
rownames(widefxnl_pct_ordgrps) <- widefxnl_pct_ordgrps$rowid
# drop cols not needed -- use combined forb-nfix groups
widefxnl_pct_ordgrps <- subset(widefxnl_pct_ordgrps, select = c(ExoForbNfix, ExoticGrass, NatForbNfix, NativeGrass, SeededNativeForb, SeededNativeGrass, pct_bare))
# be sure widefxnl dataset in same order as env
widefxnl_pct_ordgrps <- widefxnl_pct_ordgrps[rownames(wide_env),]

# for convenience, make dfs for unseeded and seeded sub-subplots only
# drop background native grass from seeded plots -- only occurs in 1 native seeded subsubplot
widefxnl_pct_ordgrps_seeded <- widefxnl_pct_ordgrps[grepl("NatSeed", rownames(widefxnl_pct_ordgrps)),names(widefxnl_pct_ordgrps) != "NativeGrass"]
# drop seeded cols in unseeded plots df
widefxnl_pct_ordgrps_unseed <- widefxnl_pct_ordgrps[grepl("NoSeed", rownames(widefxnl_pct_ordgrps)),!grepl("^Seeded", names(widefxnl_pct_ordgrps))]
# subset enviro dat to seeded and unseeded
wide_env_seeded <- wide_env[grepl("NatSeed", rownames(wide_env)),]
wide_env_unseed <- wide_env[grepl("NoSeed", rownames(wide_env)),]



# -- PREP SPP DATA FOR ORDINATIONS -----
# CA and CCA are more approp than pca and rda because can accomodate non-linear responses, is robust to assumptions not being met, and long gradients
# PCA is better for short gradients (on latent scale)

ordsp_pct_all <- left_join(wide_env[c("plot", "seedtrt", "herbicide", "rowid")], widesp_pct)
rownames(ordsp_pct_all) <- ordsp_pct_all$rowid
ordsp_pct_all <- subset(ordsp_pct_all, select = -c(plot:yr))
# df for unseeded only (ambient), drop anything not present
ordsp_pct_unseed <- ordsp_pct_all[grepl("NoSeed", rownames(ordsp_pct_all)),]
ordsp_pct_unseed <- ordsp_pct_unseed[,colSums(ordsp_pct_unseed) > 0]
# repeat for seeded
ordsp_pct_seeded <- ordsp_pct_all[grepl("NatSeed", rownames(ordsp_pct_all)),]
ordsp_pct_seeded <- ordsp_pct_seeded[,colSums(ordsp_pct_seeded) > 0]

# set pct plots sp should be present
plotpres <- 0.05


# -- aggregate anything present but rare in spp level data ----
freq_spcheck_all <- sapply(ordsp_pct_all, function(x) sum(x>0) > (plotpres*nrow(ordsp_pct_all)))
lumpsp_all <- freq_spcheck_all[!freq_spcheck_all | names(freq_spcheck_all) %in% c("AST3", "UNBU")]
# what are they?
sort(names(lumpsp_all))
# aggregate
allplots_lowsp <- ordsp_pct_all[names(lumpsp_all)] %>%
  rownames_to_column("rowid") %>%
  gather(code4, pct, names(.)[names(.) != "rowid"]) %>%
  left_join(distinct(natambsp_cov[c("code4", "fxnl_grp", "nativity", "coarse_fxnl")])) %>%
  group_by(rowid, coarse_fxnl) %>%
  reframe(pct_cov = sum(pct)) %>%
  distinct() %>%
  # shorten group names
  mutate(coarse_fxnl = gsub("Background |-| ", "", coarse_fxnl),
         coarse_fxnl = gsub("Native", "Nat", coarse_fxnl),
         coarse_fxnl = gsub("ixer", "ix", coarse_fxnl)) %>%
  spread(coarse_fxnl, pct_cov) %>%
  # combine nfix with forb because few obs
  mutate(ExoForbNfix = ExoticForb + ExoticNfix,
         NatForbNfix = NatForb + NatNfix) %>%
  # join pct bare
  left_join(may_pctbare[c("rowid", "pct_bare")]) %>%
  data.frame()
sapply(allplots_lowsp[,-1], function(x) sum(x > 0))
sapply(allplots_lowsp[,-1], function(x) sort(unique(x[x>0]))) # 5% of all plots is 7.05. exotic grass is in 7, but lumped cover seems negligible.
# ^keep just exo and nat forbnfix group + pct_bare
# winnow, and order correctly
rownames(allplots_lowsp) <- allplots_lowsp$rowid
allplots_lowsp <- allplots_lowsp[rownames(ordsp_pct_all),]

# drop low freq and abundance spp from all subsubplots, cbind aggregated groups
ordsp_pct_all <- ordsp_pct_all[!names(ordsp_pct_all) %in% names(lumpsp_all)] %>%
  cbind(allplots_lowsp[c("ExoForbNfix", "NatForbNfix", "pct_bare")])



# -- prep unseeded subsubplots ----
# check frequency, lump things that don't appear often
freq_spcheck_unseed <- sapply(ordsp_pct_unseed, function(x) sum(x>0) > (plotpres*nrow(ordsp_pct_unseed)))
lumpsp_unseed <- freq_spcheck_unseed[!freq_spcheck_unseed | names(freq_spcheck_unseed) %in% c("AST3", "UNBU")]
unseed_lowsp <- ordsp_pct_unseed[names(lumpsp_unseed)] %>%
  rownames_to_column("rowid") %>%
  gather(code4, pct, names(.)[names(.) != "rowid"]) %>%
  left_join(distinct(natambsp_cov[c("code4", "fxnl_grp", "nativity", "coarse_fxnl")])) %>%
  group_by(rowid, coarse_fxnl) %>%
  reframe(pct_cov = sum(pct)) %>%
  distinct() %>%
  # shorten group names
  mutate(coarse_fxnl = gsub("Background |-| ", "", coarse_fxnl),
         coarse_fxnl = gsub("Native", "Nat", coarse_fxnl),
         coarse_fxnl = gsub("ixer", "ix", coarse_fxnl)) %>%
  spread(coarse_fxnl, pct_cov) %>%
  # combine nfix with forb because few obs
  mutate(ExoForbNfix = ExoticForb + ExoticNfix,
         NatForbNfix = NatForb + NatNfix) %>%
  # join pct bare
  left_join(may_pctbare[c("rowid", "pct_bare")]) %>%
  data.frame()
sapply(unseed_lowsp[,-1], function(x) sum(x > 0)) 
sapply(unseed_lowsp[,-1], function(x) sort(unique(x))) # it's really just the forbnfix cols + pct bare to keep. grass cover is negligible.
# winnow, and order correctly
rownames(unseed_lowsp) <- unseed_lowsp$rowid
unseed_lowsp <-unseed_lowsp[rownames(ordsp_pct_unseed),]

# drop low freq and abundance spp from unseed sp pct cov, cbind aggregated groups
ordsp_pct_unseed <- ordsp_pct_unseed[!names(ordsp_pct_unseed) %in% names(lumpsp_unseed)] %>%
  cbind(unseed_lowsp[c("ExoForbNfix", "NatForbNfix", "pct_bare")])

# check that rownames in same order as env data
summary(rownames(ordsp_pct_unseed) == rownames(wide_env)[!grepl("NatSe", rownames(wide_env))])



# -- prep seeded subsubplots ----
# add ast3 and unbu to low abd since unknowns
freq_spcheck_seeded <- sapply(ordsp_pct_seeded, function(x) sum(x>0) > (plotpres*nrow(ordsp_pct_seeded)))
lumpsp_seeded <- freq_spcheck_seeded[!freq_spcheck_seeded | names(freq_spcheck_seeded) %in% c("AST3", "UNBU")]

# aggregate low freq-abd species in seeded sub-subplots
seeded_lowsp <- ordsp_pct_seeded[names(lumpsp_seeded)] %>%
  rownames_to_column("rowid") %>%
  gather(code4, pct, names(.)[names(.) != "rowid"]) %>%
  left_join(distinct(natambsp_cov[c("code4", "fxnl_grp", "nativity", "coarse_fxnl")])) %>%
  group_by(rowid, coarse_fxnl) %>%
  reframe(pct_cov = sum(pct)) %>%
  distinct() %>%
  # shorten group names
  mutate(coarse_fxnl = gsub("Background |-| ", "", coarse_fxnl),
         coarse_fxnl = gsub("Native", "Nat", coarse_fxnl),
         coarse_fxnl = gsub("ixer", "ix", coarse_fxnl)) %>%
  spread(coarse_fxnl, pct_cov) %>%
  # combine nfix with forb for consistency
  mutate(ExoForbNfix = ExoticForb + ExoticNfix,
         NatForbNfix = NatForb + NatNfix) %>%
  # join pct bare
  left_join(may_pctbare[c("rowid", "pct_bare")]) %>%
  data.frame()
sapply(seeded_lowsp[,-1], function(x) sum(x > 0)) # it's really just the forbnfix cols + pct bare to keep. grass cover is negligible (manual check)
sapply(seeded_lowsp[,-1], function(x) sort(unique(x))) # 30% seems like a lot for nat forbs
# what's the max cov for the spp pulled out?
sapply(ordsp_pct_seeded[names(lumpsp_seeded)], max) # it's popcorn flower
# winnow, and order correctly
rownames(seeded_lowsp) <- seeded_lowsp$rowid
seeded_lowsp <-seeded_lowsp[rownames(ordsp_pct_seeded),]

# drop low freq and abundance spp from seeded sp pct cov, cbind aggregated groups
ordsp_pct_seeded <- ordsp_pct_seeded[!names(ordsp_pct_seeded) %in% names(lumpsp_seeded)] %>%
  cbind(seeded_lowsp[c("ExoForbNfix", "NatForbNfix", "pct_bare")])


# check that rownames in same order as env data
summary(rownames(ordsp_pct_seeded) == rownames(wide_env)[!grepl("NoSee", rownames(wide_env))])





# -- ASSESS VEG DISTRIBUTIONS -----
# also make distance-based matrices for convenience here if need

# hellinger transformation on all subsubplots
widefxnl_hel_ordgrps <- decostand(widefxnl_pct_ordgrps, method = "hellinger")
# calculate relative for comparison (curiosity)
widefxnl_rel_ordgrps <- decostand(widefxnl_pct_ordgrps, method = "total")
rbind(cbind(covtype = "raw", widefxnl_pct_ordgrps),
      cbind(covtype = "relative", widefxnl_rel_ordgrps),
      cbind(covtype = "hellinger", widefxnl_hel_ordgrps)) %>%
  rownames_to_column("rowid") %>%
  gather(spp, cover, names(.)[!names(.) %in% c("rowid", "covtype")]) %>%
  # clean up rowid (1 appended when row-binded)
  mutate(rowid = gsub("Herb[0-9]$", "Herb", rowid)) %>%
  spread(covtype, cover) %>%
  ggplot() +
  geom_point(aes(raw, hellinger), pch = 1) +
  geom_point(aes(raw, relative), col = "orchid", pch = 1) +
  theme_minimal() +
  labs(y = "transformed cover", x = "raw pct cover",
       subtitle = "compare relative (orchid) v. hellinger (black) cover with raw pct cover") +
  facet_wrap(~spp, nrow = 2)
# ^ seems like it gives a little more weight (higher values) to things that are less abundant
# notes did say hellinger gives less importance to highly abundant spp
# compare rowsums because hellinger sums > 1
plot(rowSums(widefxnl_hel_ordgrps) ~ rowSums(widefxnl_pct_ordgrps))

# apply hellinger transformation to seeded and unseeded plots
widefxnl_hel_ordgrps_seeded <- decostand(widefxnl_pct_ordgrps_seeded, method = "hellinger")
widefxnl_hel_ordgrps_unseed<- decostand(widefxnl_pct_ordgrps_unseed, method = "hellinger")

# repeat for species level veg data
ordsp_hel_all <- decostand(ordsp_pct_all, method = "hellinger")
ordsp_rel_all <- decostand(ordsp_pct_all, method = "total")
ordsp_hel_seeded <- decostand(ordsp_pct_seeded, method = "hellinger")
ordsp_hel_unseed <- decostand(ordsp_pct_unseed, method = "hellinger")

# consider transformation vs. raw pct at spp level
rbind(cbind(covtype = "raw", ordsp_pct_all),
      cbind(covtype = "relative", ordsp_rel_all),
      cbind(covtype = "hellinger", ordsp_hel_all)) %>%
  rownames_to_column("rowid") %>%
  gather(spp, cover, names(.)[!names(.) %in% c("rowid", "covtype")]) %>%
  # clean up rowid (1 appended when row-binded)
  mutate(rowid = gsub("Herb[0-9]$", "Herb", rowid)) %>%
  spread(covtype, cover) %>%
  ggplot() +
  geom_point(aes(raw, hellinger), pch = 1) +
  geom_point(aes(raw, relative), col = "orchid", pch = 1) +
  theme_minimal() +
  labs(y = "transformed cover", x = "raw pct cover",
       subtitle = "compare sp-level relative (orchid) v. hellinger (black) cover with raw pct cover") +
  facet_wrap(~spp) # many spp panels if uncomment


 
# -- FXNL GROUP CA AND CCA ------
# ca can take abundance data; if use hellinger transformation and euclidean then can use pca

# -- unconstrained -----
# hellinger transformed
fxnl_hel_ca_unseed <- cca(widefxnl_hel_ordgrps_unseed ~ 1, data = wide_env_unseed)
plot(fxnl_hel_ca_unseed)
# constrast with raw pct
fxnl_pct_ca_unseed <- cca(widefxnl_pct_ordgrps_unseed ~ 1, data = wide_env_unseed)
plot(fxnl_pct_ca_unseed) # seems different
# take euclidean distance in pca comparison -- scaling compares correlations of sites/species instead of covariance
fxnl_hel_pca_unseed <- rda(widefxnl_hel_ordgrps_unseed, scale = T)

# constrast with ca
biplot(fxnl_hel_pca_unseed, scaling = "symmetric")
biplot(fxnl_hel_pca_unseed, display = "sites", scaling = "sites")
biplot(fxnl_hel_pca_unseed, display = "species", scaling = "species")

# compare unconstrained axes and variation loaded
summary(fxnl_pct_ca_unseed) # 4 axes; 83% variation in data captured on first two axes
summary(fxnl_hel_ca_unseed) # 4 axes, less variance with hellinger transformations; 70% on first two axes
summary(fxnl_hel_pca_unseed) # 5 axes, first two capture ~ 69% data variation

scores(fxnl_pct_ca_unseed, display = "species")
scores(fxnl_hel_ca_unseed, display = "species")
scores(fxnl_hel_pca_unseed, display = "species")

plot(fxnl_pct_ca_unseed, display = "species")
plot(fxnl_hel_ca_unseed, display = "species") # difference is where exotic forbs and native grams ordinate
plot(fxnl_hel_pca_unseed, display = "species") # agree on native forb and nfixers plotting closest to bare ground


# -- constrained ordinations ---- 
# nest design structure: ppt plot nested in nut plot nested in block
fxnl_pct_unseed_cca_fullnesterror <- cca(widefxnl_pct_ordgrps_unseed ~ Condition(block/xpos/ypos), data = wide_env_unseed)
fxnl_pct_unseed_cca_fullnesterror # nested structure explains about ~60% data variance; 40% unexplained
plot(fxnl_pct_unseed_cca_fullnesterror, scaling = "symmetric") # CA 2 looks like a herbicide axis

# compare with hellinger transformed
fxnl_hel_unseed_cca_fullnesterror <- cca(widefxnl_hel_ordgrps_unseed ~ Condition(block/xpos/ypos), data = wide_env_unseed)
summary(fxnl_hel_unseed_cca_fullnesterror) # conditional structure explains 67%; 33% unexplained 
plot(fxnl_hel_unseed_cca_fullnesterror, scaling = "symmetric") # strange diamond pattern in middle

# pca to compare
fxnl_hel_unseed_pca_fullnesterror <- rda(widefxnl_hel_ordgrps_unseed ~ Condition(block/xpos/ypos), data = wide_env_unseed)
summary(fxnl_hel_unseed_pca_fullnesterror) # conditional structure explains 62%

# look at env influence only -- no conditioning
fxnl_pct_unseed_cca_envx <- cca(widefxnl_pct_ordgrps_unseed ~ herbicide * nut_trt* ppt_trt, data = wide_env_unseed)
summary(fxnl_pct_unseed_cca_envx) # 38% of data, 61% unexplained when have herbicide as interaction; 28% explained by constraints when herbicide does not interact
fxnl_pct_unseed_cca_env <- cca(widefxnl_pct_ordgrps_unseed ~ herbicide + nut_trt* ppt_trt, data = wide_env_unseed)
fxnl_pct_unseed_cca_envadd <- cca(widefxnl_pct_ordgrps_unseed ~ herbicide + nut_trt + ppt_trt, data = wide_env_unseed)
# permutest signif of models
anova(fxnl_pct_ca_unseed, fxnl_pct_unseed_cca_envadd, fxnl_pct_unseed_cca_env, fxnl_pct_unseed_cca_envx)
# ^ suggests additive best model, permuting all
fxnlenv_eval <- ordistep(fxnl_pct_ca_unseed,
                         scope = list(lower = formula(fxnl_pct_ca_unseed), upper = formula(fxnl_pct_unseed_cca_envx)),
                         direction = "both",
                         # allow whole plots to shuffle within blocks -- cannot permute because design is unbalanced (one of the block 4 plots was not herbicided)
                         #permutations = how(within = Within("none"), plots = Plots(wide_env_unseed$xpos, type = "free"), blocks = wide_env_unseed$block),
                          trace = F)
fxnlenv_eval # suggests only need herbicide + ppt..
fxnlenv_eval$anova

drop1(fxnl_pct_unseed_cca_envx) # drop interaction
drop1(fxnl_pct_unseed_cca_envadd) # also drop nut_trt
anova(fxnl_pct_unseed_cca_envadd, by = "margin", strata = wide_env_unseed$block) # nutrient addition does not matter

plot(fxnl_pct_unseed_cca_env, scaling = "symmetric") # looks like ppt is axis 2 (wet pos, drought neg); nut trt is cca 1
plot(fxnl_pct_unseed_cca_envadd, display = c("sp", "bp"), scaling = "species")


#  -- full model with conditioning and environmental effects ---
fxnl_unseed_cca_envxcond <- cca(widefxnl_pct_ordgrps_unseed ~ herbicide + nut_trt * ppt_trt + Condition(block), #Condition(block + block:xpos + block:xpos:ypos),
                                  data = wide_env_unseed)
# ^if condition with full nested error structure, will not consider nut_trt and ppt_trt as factors, only herbicide
# ^of conditioned vars, first two axes capture the most variation (seems like herbicide on axis 1, maybe also ppt; nut on axis 2)
summary(fxnl_unseed_cca_envxcond)
# unconstrained accounts for 60% of variation; herb + env 26% and block 11%
plot(fxnl_unseed_cca_envxcond, display = c("sp", "bp"), main = "herbicide + nut x ppt (condtioned on block)",  scaling = "symmetric")
plot(fxnl_unseed_cca_envxcond, display = "wa", main = "sites scores", scaling = "symmetric")
# check out model
vif.cca(fxnl_unseed_cca_envxcond) # no inflated terms (wasn't expecting there to be)
RsquareAdj(fxnl_unseed_cca_envxcond) # adj r sq is 17%
anova(fxnl_unseed_cca_envxcond, by = "axis",
      permutations = how(plots = Plots(wide_env_unseed$xpos, type = "none"), blocks = wide_env_unseed$block)) # only cca 1 and 2 are signif
anova(fxnl_unseed_cca_envxcond, by = "term", 
      permutations = how(blocks = wide_env_unseed$block)) # herb and ppt signif; marginal interaction;
# ^ if reorder so ppt comes before nut, then nut is also marginally signif
anova(fxnl_unseed_cca_envxcond, by = "margin", 
      permutations = how(blocks = wide_env_unseed$block)) # interaction has marginal contribution

# test additive model to see if nutrient amendment is important
fxnl_unseed_cca_envcond <- cca(widefxnl_pct_ordgrps_unseed ~ herbicide + ppt_trt + nut_trt + Condition(block), #Condition(block + block:xpos + block:xpos:ypos),
                                data = wide_env_unseed)
anova(fxnl_unseed_cca_envcond, by = "axis",
      permutations = how(blocks = wide_env_unseed$block)) # first two axes important still
anova(fxnl_unseed_cca_envcond, by = "term", 
      permutations = how(blocks = wide_env_unseed$block)) # herb and ppt signif; nut_trt marginal
anova(fxnl_unseed_cca_envcond, by = "margin", 
      permutations = how(blocks = wide_env_unseed$block)) # nutrient addition adds nothing discernable

# -- hellinger full env-conditioned models -----
fxnl_hel_unseed_cca_envxcond <- cca(widefxnl_hel_ordgrps_unseed ~ herbicide + nut_trt * ppt_trt + Condition(block), #Condition(block + block:xpos + block:xpos:ypos),
                               data = wide_env_unseed)
plot(fxnl_hel_unseed_cca_envxcond, display = c("bp", "sp"), scaling = "symmetric")
summary(fxnl_hel_unseed_cca_envxcond)
anova(fxnl_hel_unseed_cca_envxcond, by = "axis",
      permutations = how(blocks = wide_env_unseed$block)) # first two axes important still
anova(fxnl_hel_unseed_cca_envxcond, by = "term", 
      permutations = how(plots = Plots(strata = wide_env_unseed$xpos, type = "none"), blocks = wide_env_unseed$block)) # herb and ppt signif; nut_trt marginal
anova(fxnl_hel_unseed_cca_envxcond, by = "margin", 
      permutations = how(plots = Plots(strata = wide_env_unseed$xpos, type = "none"), blocks = wide_env_unseed$block)) # nutrient addition adds nothing discernable
# marginal signif for env interaction at functional group level
goodness(fxnl_hel_unseed_cca_envxcond)

# additive model
fxnl_hel_unseed_cca_envcond <- cca(widefxnl_hel_ordgrps_unseed ~ herbicide + nut_trt + ppt_trt + Condition(block), #Condition(block + block:xpos + block:xpos:ypos),
                                    data = wide_env_unseed)
anova(fxnl_hel_unseed_cca_envcond, by = "margin", 
      permutations = how(plots = Plots(strata = wide_env_unseed$xpos, type = "none"), blocks = wide_env_unseed$block)) # nutrient addition marginally signif

# contrast with rda (pca)
fxnl_hel_unseed_rda_envxcond <- rda(widefxnl_hel_ordgrps_unseed ~ herbicide + nut_trt * ppt_trt + Condition(block), #Condition(block + block:xpos + block:xpos:ypos),
                                    data = wide_env_unseed)
summary(fxnl_hel_unseed_rda_envxcond) # proportion similar to cca
plot(fxnl_hel_unseed_rda_envxcond, display = c("bp", "sp"), scaling = "species")


# -- pull results to make nicer plots -----
fxnl_hel_unseed_envxcond_cca_df <- scores(fxnl_hel_unseed_cca_envxcond, tidy = T) %>%
  left_join(wide_env, by = c("label" = "rowid"))

ggplot(subset(fxnl_hel_unseed_envxcond_cca_df, score == "sites"), aes(CCA1, CCA2)) +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey50") +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey50") +
  geom_point(aes(fill = paste(ppt_trt, nut_trt), shape = herbicide, group = paste(herbicide)), size = 3, alpha = 0.75) +
  scale_shape_manual(values = c(21, 23)) +
  scale_fill_manual(name = "PPT x NUT", values = c(paste0("goldenrod", c(4,3,1)), 
                                                   paste0("royalblue", c(3,1)), "powderblue", paste0("palegreen", c(4,3,1)))) +
  facet_grid(ppt_trt ~block) +
  theme(strip.background = element_rect(fill = "transparent")) +
  labs(subtitle = "Constrained and conditioned CA on functional groups, paneled by block x ppt trt") +
  guides(fill = guide_legend(override.aes = list(shape = 21)))


ggplot(fxnl_hel_unseed_envxcond_cca_df, aes(CCA1, CCA2)) +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey50") +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey50") +
  geom_segment(data = subset(fxnl_hel_unseed_envxcond_cca_df, score %in% c("biplot", "factorbiplot")),
               aes(x = 0, xend = CCA1, y = 0, yend = CCA2),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"),
               color = "grey20") + 
  ggrepel::geom_text_repel(data = subset(fxnl_hel_unseed_envxcond_cca_df, score %in% c("biplot", "factorbiplot")), 
                           box.padding = unit(10, "pt"),min.segment.length = 10,
                           aes(label = gsub("herbicide", "", label)), 
                           size = 3) +
  geom_text(data = subset(fxnl_hel_unseed_envxcond_cca_df, score %in% c("species")), aes(label = label), size = 5, col = "forestgreen") +
  # geom_text(data = subset(fxnl_hel_unseed_envxcond_cca_df, score %in% c("regression") & grepl("block", label) & !is.na(CCA1)), aes(label = gsub("wholeplotID", "", label)), 
  #           size = 3, col = "purple", position = position_jitter(width = 0.02, height = 0.02)) +
  scale_x_continuous(expand = c(0.05, 0.05))


# -- all plots ----
fxnl_hel_all_cca_envxcond <- cca(widefxnl_hel_ordgrps ~ seedtrt + herbicide + nut_trt * ppt_trt + Condition(block), #Condition(block + block:xpos + block:xpos:ypos),
                                    data = wide_env)
summary(fxnl_hel_all_cca_envxcond)
plot(fxnl_hel_all_cca_envxcond, display = c("bp", "sp"), scaling = "symmetric")
plot(fxnl_hel_all_cca_envxcond, display = c("wa"), scaling = "sites")

anova(fxnl_hel_all_cca_envxcond, by = "axis",
      permutations = how(plots = Plots(strata = wide_env$xpos, type = "none"), blocks = wide_env$block)) # first three signif
anova(fxnl_hel_all_cca_envxcond, by = "margin",
      permutations = how(plots = Plots(strata = wide_env$xpos, type = "none"), blocks = wide_env$block)) # all signif
anova(fxnl_hel_all_cca_envxcond, by = "term",
      permutations = how(plots = Plots(strata = wide_env$xpos, type = "none"), blocks = wide_env$block)) # all

fxnl_hel_all_cca_envcond <- cca(widefxnl_hel_ordgrps ~ seedtrt + herbicide + ppt_trt + nut_trt + Condition(block), #Condition(block + block:xpos + block:xpos:ypos),
                                 data = wide_env)
summary(fxnl_hel_all_cca_envcond) # not much change in variance explained
plot(fxnl_hel_all_cca_envcond, display = c("bp", "sp"), scaling = "symmetric")

# pull for nicer plot
fxnl_hel_all_envxcond_cca_df <- scores(fxnl_hel_all_cca_envxcond, tidy = T) %>%
  left_join(wide_env, by = c("label" = "rowid"))

# site loadings
ggplot(subset(fxnl_hel_all_envxcond_cca_df, score == "sites"), aes(CCA1, CCA2)) +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey50") +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey50") +
  geom_point(aes(fill = paste(ppt_trt, nut_trt), shape = paste(herbicide, seedtrt), group = paste(herbicide)), size = 3, alpha = 0.75) +
  scale_shape_manual(name = "Herb + seeding", values = c(21, 22, 24, 25)) +
  scale_fill_manual(name = "PPT x NUT", values = c(paste0("goldenrod", c(4,3,1)), 
                                                   paste0("royalblue", c(3,1)), "powderblue", paste0("palegreen", c(4,3,1)))) +
  facet_grid(ppt_trt ~block) +
  theme(strip.background = element_rect(fill = "transparent")) +
  labs(subtitle = "Constrained and conditioned CA on functional groups, paneled by block x ppt trt, all sub-sub-subplots") +
  guides(fill = guide_legend(override.aes = list(shape = 21)))

# fxnl group loadings and treatment effects
ggplot(fxnl_hel_all_envxcond_cca_df, aes(CCA1, CCA2)) +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey50") +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey50") +
  geom_segment(data = subset(fxnl_hel_all_envxcond_cca_df, score %in% c("biplot", "factorbiplot")),
               aes(x = 0, xend = CCA1, y = 0, yend = CCA2),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"),
               color = "grey20") + 
  # yneg labels
  ggrepel::geom_text_repel(data = subset(fxnl_hel_all_envxcond_cca_df, score %in% c("biplot", "factorbiplot") & CCA2 <0), 
                           box.padding = unit(10, "pt"),min.segment.length = 10,nudge_y = -0.05,
                           aes(label = gsub("herbicide|seedtrt", "", label)), 
                           size = 3) +
  # ypos labels
  ggrepel::geom_text_repel(data = subset(fxnl_hel_all_envxcond_cca_df, score %in% c("biplot", "factorbiplot") & CCA2 >= 0), 
                           box.padding = unit(10, "pt"),min.segment.length = 10, nudge_y = 0.05,
                           aes(label = gsub("herbicide|seedtrt", "", label)), 
                           size = 3) +
  geom_text(data = subset(fxnl_hel_all_envxcond_cca_df, score %in% c("species")), aes(label = label), 
            size = 4, col = "forestgreen", hjust = 0, vjust = 0) +
  # geom_text(data = subset(fxnl_hel_unseed_envxcond_cca_df, score %in% c("regression") & grepl("block", label) & !is.na(CCA1)), aes(label = gsub("wholeplotID", "", label)), 
  #           size = 3, col = "purple", position = position_jitter(width = 0.02, height = 0.02)) +
  scale_x_continuous(expand = c(0.1, 0.1))


# -- seeded plots only ----
fxnl_hel_seeded_cca_envxcond <- cca(widefxnl_hel_ordgrps_seeded ~ herbicide + ppt_trt * nut_trt + Condition(block), #Condition(block + block:xpos + block:xpos:ypos),
                                 data = wide_env_seeded)
summary(fxnl_hel_seeded_cca_envxcond) # similar variance breakdown
anova(fxnl_hel_seeded_cca_envxcond, by = "axis",
      permutations = how(blocks = wide_env_seeded$block)) # first two signif
anova(fxnl_hel_seeded_cca_envxcond, by = "term",
      permutations = how(blocks = wide_env_seeded$block)) # interaction marginal
anova(fxnl_hel_seeded_cca_envxcond, by = "margin",
      permutations = how(blocks = wide_env_seeded$block)) # interaction does add something

eigenvals(fxnl_hel_seeded_cca_envxcond) # first two unconditioned axes explain more than first two (signif) constrained axes
plot(fxnl_hel_seeded_cca_envxcond, display = c("bp", "sp"), scaling = "species")



# -- SPP LEVEL UNCONSTRAINED AND CONSTRAINED CCA -----
sp_pct_unseed_ca <- cca(ordsp_pct_unseed ~ 1, data = wide_env_unseed)
summary(sp_pct_unseed_ca) # first axis captures 15% of variation in data, lots of unexplained; 37 axes
plot(sp_pct_unseed_ca, scaling = "symmetric") # site clustering on pos half of CA1

# test on hellinger transformed
sp_hel_unseed_ca <- cca(ordsp_hel_unseed ~ 1, data = wide_env_unseed)
summary(sp_hel_unseed_ca) # PC1 captures 38%.. ; 35 axes
plot(sp_hel_unseed_ca, scaling = "symmetric") # sites more spread out

# test pca
sp_hel_unseed_rda <- rda(ordsp_hel_unseed ~ 1, data = wide_env_unseed)
summary(sp_hel_unseed_rda) # PC1 captures 33%.. ; 36 axes
biplot(sp_hel_unseed_rda) # more like pct-based cca; PC1 seems to split hill pos, and maybe PC2 splits block 3 and block 4


# -- test enviro-conditioned model -----
# stick with cca
sp_hel_unseed_cca_envxcond <- cca(formula = ordsp_hel_unseed ~ herbicide + ppt_trt * nut_trt + Condition(block), 
                              data = wide_env_unseed)
plot(sp_hel_unseed_cca_envxcond, display = c("bp", "sp"), scaling = "symmetric")
summary(sp_hel_unseed_cca_envxcond) # 35 unconstrained axes..

# permutest signif
anova(sp_hel_unseed_cca_envxcond, by = "axis", 
      permutations = how(block = wide_env_unseed$block)) # first 2 CCA axes signif
anova(sp_hel_unseed_cca_envxcond, by = "margin", 
      permutations = how(block = wide_env_unseed$block)) # interaction not signif

# additive only 
sp_hel_unseed_cca_envcond <- cca(formula = ordsp_hel_unseed ~ herbicide + ppt_trt + nut_trt + Condition(block), 
                              data = wide_env_unseed)
plot(sp_hel_unseed_cca_envcond, display = c("bp", "sp"), scaling = "symmetric")
summary(sp_hel_unseed_cca_envcond) # block captures much more variation for spp compared to functional groups
anova(sp_hel_unseed_cca_envcond, by = "margin", 
      permutations = how(block = wide_env_unseed$block)) # all add something
eval_sp_unseed <- ordistep(sp_hel_unseed_ca, scope = list(upper = sp_hel_unseed_cca_envxcond), 
         permutations = how(block = wide_env_unseed$block), trace = F)
eval_sp_unseed # additive model is best, without conditioning
eval_sp_unseed$anova
# reduced model is better
ordiR2step(sp_hel_unseed_cca_envcond, sp_hel_unseed_cca_envxcond, trace = F)


# -- all and seeded plots -----
# base model 
sp_hel_all_ca <- cca(formula = ordsp_hel_all ~ 1, data = wide_env)

# additive only 
sp_hel_all_cca_envcond <- cca(formula = ordsp_hel_all ~ seedtrt + herbicide + ppt_trt + nut_trt + Condition(block), 
                                 data = wide_env)
# with ppt x nut
sp_hel_all_cca_envxcond <- cca(formula = ordsp_hel_all ~ seedtrt + herbicide + ppt_trt * nut_trt + Condition(block), 
                              data = wide_env)
eval_sp_all <- ordistep(sp_hel_all_ca, scope = list(upper = sp_hel_all_cca_envxcond),direction = "both", 
                        permutations = how(blocks = wide_env$block), trace = F)
eval_sp_all
eval_sp_all$anova # full model best
anova(sp_hel_all_cca_envxcond, by = "axis", permutations = how(blocks = wide_env$block)) # first 4 cc axes signif
anova(sp_hel_all_cca_envxcond, by = "margin", permutations = how(blocks = wide_env$block)) # all add something

# plot
plot(sp_hel_all_cca_envxcond, display = c("bp", "sp"), scaling = "species")
plot(sp_hel_all_cca_envxcond, display = c("wa"), scaling = "species")


# -- seeded plots -----
sp_hel_seeded_cca_envcond <- cca(formula = ordsp_hel_seeded ~ herbicide + ppt_trt + nut_trt + Condition(block), 
                              data = wide_env_seeded)
# with ppt x nut
sp_hel_seeded_cca_envxcond <- cca(formula = ordsp_hel_seeded ~ herbicide + nut_trt * ppt_trt + Condition(block), 
                               data = wide_env_seeded)
summary(sp_hel_seeded_cca_envxcond)
plot(sp_hel_seeded_cca_envxcond, display = c("bp", "sp"), scaling = "symmetric",
     main = "CCA: seeded plot spp ~ herbicide + nut_trt * ppt_trt + Cond(block)\n(symmetric scaling)")
plot.cca(sp_hel_seeded_cca_envxcond, axis.bp = T, scaling = "symmetric")
anova(sp_hel_seeded_cca_envcond, sp_hel_seeded_cca_envxcond) # interactive model better
anova(sp_hel_seeded_cca_envxcond, by = "axis", permutations = how(blocks = wide_env_seeded$block)) # only first two axes
anova(sp_hel_seeded_cca_envxcond, by = "margin", permutations = how(plots = Plots(wide_env_seeded$xpos, type = "none"), 
                                                             blocks = wide_env_seeded$block)) 
summary(sp_hel_seeded_cca_envxcond)
# never permute samples between blocks, only within
# permute can permute whole and/or split plot; whole plots have to have same number of obs [must be balaned design]
# can keep whole plot fixed and shuffle split plots within to test signif of effects within whole plots
# min p is 1/(num of perms, including the observed)
# exact test will eval all possible permutations (good when have few perm options via constant = T)

how()
# within = lowest level, plots = interemediate, blocks = outer level
how(within= Within(), plots = Plots(), blocks = Blocks())

adonis2(ordsp_hel_unseed ~ herbicide + nut_trt + ppt_trt,
        strata = wide_env_unseed$block,
        by = "margin",
       data = wide_env_unseed)

permutest(betadisper(vegdist(ordsp_hel_unseed), wide_env_unseed$block), pairwise = T) # block 4 different from others but block 3
# different at block level, but not at wholeplot level
plot(betadisper(vegdist(ordsp_hel_unseed), wide_env_unseed$block)) 

boxplot(betadisper(vegdist(ordsp_hel_unseed), wide_env_unseed$block))

permutest(betadisper(vegdist(ordsp_hel_all),wide_env$block),
          pairwise = T) # block 4 still different from other blocks at all spp level
boxplot(betadisper(vegdist(ordsp_hel_all), wide_env$block))

permutest(betadisper(vegdist(ordsp_hel_all),wide_env$wholeplotID),pairwise = T)
# but variation within wholeplot is not too different ITO dispersion
boxplot(betadisper(vegdist(ordsp_hel_all), wide_env$wholeplotID))



# -- PCA ON TRAITS ------
traitmat <- subset(natamb_avgtraits_wide, select = c(code4, Coarse.root.diameter.mm:Total.biomass.g, fxnl_grp, nativity, duration.abbr)) %>%
  # join alternate code4 and rename to code 4
  rbind(rename(subset(natamb_avgtraits_wide, !is.na(code4_alternate), select = c(code4_alternate, Coarse.root.diameter.mm:Total.biomass.g, fxnl_grp, nativity, duration.abbr)), code4 = code4_alternate)) %>%
  # make fxnl grp and duration factor
  mutate(fxnl_grp = factor(fxnl_grp, levels = c("Grass", "Forb", "N-fixer")),
         # alpha order is correct factor level (short lived to long lived)
         duration.abbr = factor(duration.abbr),
         nativity = factor(nativity)) %>%
  # be sure arranged alphabetically
  arrange(code4) %>%
  data.frame() 
# set rownames
rownames(traitmat) <- traitmat$code4

# for convenience:
# trait matrix for full dataset (all 4 sub-sub plots) in the order of ydat
traitmat_full <- traitmat[rownames(traitmat) %in% names(ordsp_pct_all),]
traitmat_seeded <- traitmat[rownames(traitmat) %in% names(ordsp_pct_seeded),]

# note cols in vegdat that don't have traits
notrait_spp <- unique(c(names(ordsp_pct_all)[!names(ordsp_pct_all) %in% traitmat_full$code4],
                        names(ordsp_pct_seeded)[!names(ordsp_pct_seeded) %in% rownames(traitmat_seeded)]))


# clean up enviro before proceed
#rm(list = ls()[grepl("asla|asltat|notraits_prop|seedsum|checkcolS|plots_sum", ls())])

# check that distribution of trait dat is roughly normal
traitmat_full[,!grepl("code4|fxnl_grp|nativity|durat", names(traitmat_full))] %>%
  gather(met, val) %>%
  ggplot() +
  geom_density(aes(val)) +
  theme_minimal() +
  facet_wrap(~gsub("[.]", " ", met), scales = "free", labeller = label_wrap_gen(width= 15))
# LDMC is okay, RMF is approx ok
# most everything else is skewed right except Propoprtion fine roots (skewed right)
# specific root length data look more unusual .. large value is crassula, which should have a super small value 
# others are two native flowers: popcorn flower and tritelia laxa
histogram(1/(traitmat_full$Coarse.root.specific.length.cm.g)) # inverse transform makes it a bit better
histogram(1/(traitmat_full$Fine.root.specific.length.cm.g))

# see if logging helps others
traitmat_full[,!grepl("code4|fxnl_grp|nativity|durat", names(traitmat_full))] %>%
  gather(met, val) %>%
  group_by(met) %>%
  mutate(val = 1/(1-val)) %>% #sqrt and log, subtract from constant for left skewed
  ggplot() +
  geom_density(aes(val)) +
  facet_wrap(~met, scales = "free")

# ln works for sla, height
# sqrt works for dry leaf mass, fine root length, length, leaf area, root dry biomass
# ^ also sqrt Fresh leaf mass, Root.dry.bio, Total.bio (better than log but not ideal)
# log 10 works for coarse root length, shoot dry bio
# inv works for root volume, root spec length, coarse root diam, spec root length

invtran <- names(traitmat_full)[grepl("specific.len|volume|oarse.root.di", names(traitmat_full))]
log10tran <- names(traitmat_full)[grepl("oarse.root.le", names(traitmat_full))]
lntran <- names(traitmat_full)[grepl("SLA|Heigh|Shoot.dry", names(traitmat_full))]
sqtran <- names(traitmat_full)[grepl("Dry.le|Fine.root.l|Length.m|Leaf.Are|Root.dry|Total.biomass|Fresh.le", names(traitmat_full))]
neginvtran <- names(traitmat_full)[grepl("Propor", names(traitmat_full))]

traitmat_full_norm <- traitmat %>%
  mutate_at(invtran, function(x) 1/x) %>%
  mutate_at(neginvtran, function(x) 1/(1-x)) %>%
  mutate_at(log10tran, function(x) log(x, 10)) %>%
  mutate_at(lntran, function(x) log(x)) %>%
  mutate_at(sqtran, function(x) sqrt(x))

traitmat_full_norm[,!grepl("code4|fxnl_grp|nativity|durat", names(traitmat_full))] %>%
  gather(met, val) %>%
  ggplot() +
  geom_density(aes(val)) +
  theme_minimal() +
  facet_wrap(~gsub("[.]", " ", met), scales = "free", labeller = label_wrap_gen(width= 15))

# repull seeded plot traits
traitmat_seeded_norm <- traitmat_full_norm[rownames(traitmat_full_norm) %in% names(ordsp_pct_seeded),]

boxplot(decostand(traitmat_full_norm[!grepl("code4|fxnl|nativ|durat", names(traitmat_full_norm))], method = "standardize")) # it's about the same
# see what's highly correlated
corrplot::corrplot(
  cor(decostand(traitmat_full_norm[!grepl("code4|fxnl|nativ|durat", names(traitmat_full_norm))], 
                method = "standardize")),
  method = "square", addrect = 6, tl.cex = 0.9, 
  order = "hclust"
  )


fulltrait_pca <- rda(traitmat_full_norm[!grepl("code4|fxnl|nativ|durat", names(traitmat_full_norm))], scale = T) # standardizes without decostand needed
summary(fulltrait_pca)
# compare with base R function
fulltrait_prcomp <- prcomp(traitmat_full_norm[!grepl("code4|fxnl|nativ|durat", names(traitmat_full_norm))], center = T, scale. = T)
fulltrait_prcomp
# these should be about the same, unless base R is looking more at covar
biplot(fulltrait_pca)
biplot(fulltrait_prcomp) # scales are a little different, but otherwise the same
summary(fulltrait_pca)
# what is on the 3rd axis?
screeplot(fulltrait_pca, main = "PC1 captures 53% variation in data, PC2 13%")
plot(fulltrait_prcomp) # this is the same


biplot(fulltrait_pca, choices = c(2,3)) # it's sort of acquisitive to conservative, but pc 1 and 2 show as much
biplot(fulltrait_pca, choices = c(2,3), scaling = "none") 

#fulltrait_pca_reduced <- rda(traitmat_full_norm[grepl("LDMC|SLA|Height|SLA|root.diam|root.specif|Prop|Root.dens", names(traitmat_full_norm))], scale = T) # standardizes without decostand needed
reducedtrait_pca <- rda(traitmat_full_norm[grepl("SLA|Fresh|RMF|LDM|Height|SLA|root.diam|Fine.root.specif|volu", names(traitmat_full_norm))], scale = T) # standardizes without decostand needed
biplot(reducedtrait_pca, display = "species")
summary(reducedtrait_pca) # PC1 captures 27%, PC2 cumulative up to 50%
biplot(reducedtrait_pca, display = "sites")
biplot(reducedtrait_pca, scaling = "symmetric")


traitpc_df <- scores(reducedtrait_pca,choices = c(1:3), tidy = T) %>% # scores(fulltrait_pca,choices = c(1:3), tidy = T) %>%
  left_join(traitmat_full[c("code4", "fxnl_grp", "nativity", "duration.abbr")], by = c("label" = "code4")) %>%
  left_join(distinct(natambsp_cov[c("code4", "native_seeded")]), by = c("label" = "code4")) %>%
  mutate(label = ifelse(score == "species", gsub("cm3|cm2|cm.g|cm|[.]g|.mm|[.]$", "", label), label))
  # add min, max for pc1 and pc2 for plotting traits as arrows from 0
  # mutate(pc1_min = ifelse(PC1 > 0, 0, PC1),
  #        pc1_max = ifelse(PC1 > 0, PC1, 0))

exoforb4 <- RColorBrewer::brewer.pal(4, "Reds")
natforb3 <- RColorBrewer::brewer.pal(3, "Purples")#[2:8]
forb7 <- RColorBrewer::brewer.pal(8, "PuRd")[c(1:4, 6:8)]
grass5 <- RColorBrewer::brewer.pal(5, "Greens")#[2:6]
nfix <- RColorBrewer::brewer.pal(3, "Oranges")#[2:4]
spp <- ggplot(subset(traitpc_df, !is.na(fxnl_grp)), aes(PC1, PC2)) +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey50") +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey50") +
  #ggrepel::geom_text_repel(aes(label = label), size = 3, box.padding = unit(10, "pt")) +
  #geom_point(aes(fill = paste(fxnl_grp, nativity, duration.abbr), size = nativity),col = "grey50", pch = 21) +
  geom_label(aes(fill = paste(fxnl_grp, nativity, duration.abbr), col = paste(fxnl_grp, nativity, duration.abbr), label = label, size = native_seeded)) +
  labs(y = NULL) +
  scale_fill_manual(name = "fxnl grp", values = c(forb7, grass5, nfix)) +
  scale_color_manual(name = "fxnl grp", values = c(rep("grey20", 4), rep("white",3),
                                                   rep("grey20", 2), rep("white", 3),
                                                   rep("grey20", 2), "white")) +
  scale_size_manual( values = c(3,5)) +
  scale_x_continuous(limits = c(min(traitpc_df$PC1)-0.5, max(traitpc_df$PC1)+1)) +
  scale_y_continuous(limits = c(min(traitpc_df$PC2), max(traitpc_df$PC2)))


trts <- ggplot(subset(traitpc_df, score == "species")) +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey50") +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey50") +
  #geom_line(aes(PC1, PC2, group = label)) +
  geom_segment(aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"),
               color = "grey20") + 
  #geom_point(data = subset(traitpc_df, score == "species"), aes(PC1, PC2)) +
  ggrepel::geom_text_repel(data = subset(traitpc_df, score == "species" & PC1 > 0 & PC2 <0), aes(PC1, PC2, label = label),
                           box.padding = unit(5, "pt"),vjust = 0,nudge_x = 0.05,size = 3, min.segment.length = 10) +
  # upper right quad
  ggrepel::geom_text_repel(data = subset(traitpc_df, score == "species" & PC1 > 0 & PC2 >0), aes(PC1, PC2, label = label),
                           hjust = 0, vjust = 0, 
                           box.padding = unit(1, "pt"),
                           point.padding = unit(3, "pt"),
                           #nudge_y = 0.05, 
                           nudge_x = 0.1, 
                           size = 3, min.segment.length = 10) +
  # lower left
  ggrepel::geom_text_repel(data = subset(traitpc_df, score == "species" & PC1 < 0 & PC2 <0), aes(PC1, PC2, label = label),
                           nudge_y = -0.05, nudge_x = -0.05, size = 3, min.segment.length = 10) +
  # upper left
  ggrepel::geom_text_repel(data = subset(traitpc_df, score == "species" & PC1 < 0 & PC2 >0), aes(PC1, PC2, label = label),
                           nudge_y = 0.05, nudge_x = -0.05, size = 3, min.segment.length = 10) +
  labs(x = "PC1", y = "PC2") +
  scale_x_continuous(limits = c(min(traitpc_df$PC1)-0.5, max(traitpc_df$PC1)+1)) +
  scale_y_continuous(limits = c(min(traitpc_df$PC2), max(traitpc_df$PC2)))


traits_pcafig <- cowplot::plot_grid(trts, spp, nrow = 1,axis = "l",
                   rel_widths = c(1, 1.3),
                   align = "v")
ggsave(plot = traits_pcafig, filename = "allspp21_reducedtraitsfig.png", 
       path = "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/chapter/figs/prelim/",
       width = 8, height = 4, units = "in", scale = 1.5)

# RMF, Height, SLA, LDMC, Coarse root diameter, Coarse and fine root spec length, prop fine roots, root density

# check correlation on reduced
corrplot::corrplot(cor(traitmat_full_norm[grepl("SLA|Fresh|RMF|LDM|Height|SLA|root.diam|Fine.root.specif|volu", names(traitmat_full_norm))]),
                  p.mat = corrplot::cor.mtest(traitmat_full_norm[grepl("SLA|Fresh|RMF|LDM|Height|SLA|root.diam|Fine.root.specif|volu", names(traitmat_full_norm))])$p,
                   method = "square", order = "hclust")


# try varimax rotate
mydata <- traitmat_full_norm[grepl("SLA|Fresh|RMF|LDM|Height|SLA|root.diam|Fine.root.specif|volu", names(traitmat_full_norm))]
#mydata <- traitmat_full_norm[!grepl("code4|fxnl|nativ|durat", names(traitmat_full_norm))]
rloading <- varimax(scores(reducedtrait_pca, display = "species"))$loadings #fulltrait_pca
iloading = t(pracma::pinv(rloading))
scores = scale(mydata) %*% iloading
r = range(c(rloading, scores))
plot(scores, xlim = r, ylim= r, xlab= "PC1 ", ylab= "PC2 ")
arrows(0,0, rloading[,1], rloading[,2], length=0.1, col=2)
text(rloading[,1], rloading[,2], labels = colnames(mydata), pos=3, col=2)
text(scores[,1], scores[,2], labels = rownames(mydata), pos = 3)
abline(h=0, lty=3)
abline(v=0, lty=3)

# recreate fig with rotated pca
rottraitpc_df <- cbind(data.frame(scores), score = "sites") %>%
  rename(PC1 = X1, PC2 = X2) %>%
  rbind(cbind(data.frame(rloading[1:(length(rloading)/2),]), score = "species")) %>%
  rownames_to_column("label") %>%
  left_join(traitmat_full_norm[c("code4", "fxnl_grp", "nativity", "duration.abbr")], by = c("label" = "code4")) %>%
  left_join(distinct(natambsp_cov[c("code4", "native_seeded")]), by = c("label" = "code4"))

rotspp <- ggplot(subset(rottraitpc_df, !is.na(fxnl_grp)), aes(PC1, PC2)) +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey50") +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey50") +
  #ggrepel::geom_text_repel(aes(label = label), size = 3, box.padding = unit(10, "pt")) +
  geom_point(aes(fill = paste(fxnl_grp, nativity, duration.abbr), size = nativity),col = "grey50", pch = 21) +
  geom_label(data = subset(rottraitpc_df, native_seeded== "Seeded"), aes(fill = paste(fxnl_grp, nativity, duration.abbr), label = label), color = "white", size = 4, #
  #            #alpha = 0.75
              show.legend = F) +
  labs(y = NULL) +
  scale_fill_manual(name = "fxnl grp", values = c(forb7, grass5, nfix)) +
  scale_color_manual(name = "fxnl grp", values = c(rep("grey20", 4), rep("white",3),
                                                   rep("grey20", 2), rep("white", 3),
                                                   rep("grey20", 2), "white")) +
  guides(color = "none", text = "none", fill = guide_legend(override.aes = list(size = 3))) +
  scale_size_manual( values = c(3,5)) +
  scale_x_continuous(limits = c(min(rottraitpc_df$PC1)-0.5, max(rottraitpc_df$PC1)+1)) +
  scale_y_continuous(limits = c(min(rottraitpc_df$PC2), max(rottraitpc_df$PC2)))

rottrts <- ggplot(subset(rottraitpc_df, score == "species")) +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey50") +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey50") +
  #geom_line(aes(PC1, PC2, group = label)) +
  geom_segment(aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"),
               color = "grey20") + 
  #geom_point(data = subset(rottraitpc_df, score == "species"), aes(PC1, PC2)) +
  ggrepel::geom_text_repel(data = subset(rottraitpc_df, score == "species" & PC1 > 0 & PC2 <0), aes(PC1, PC2, label = label),
                           box.padding = unit(5, "pt"),vjust = 0,nudge_x = 0.05,size = 3, min.segment.length = 10) +
  # upper right quad
  ggrepel::geom_text_repel(data = subset(rottraitpc_df, score == "species" & PC1 > 0 & PC2 >0), aes(PC1, PC2, label = label),
                           hjust = 0, vjust = 0, 
                           box.padding = unit(1, "pt"),
                           point.padding = unit(3, "pt"),
                           #nudge_y = 0.05, 
                           nudge_x = 0.1, 
                           size = 3, min.segment.length = 10) +
  # lower left
  ggrepel::geom_text_repel(data = subset(rottraitpc_df, score == "species" & PC1 < 0 & PC2 <0), aes(PC1, PC2, label = label),
                           nudge_y = -0.05, nudge_x = -0.05, size = 3, min.segment.length = 10) +
  # upper left
  ggrepel::geom_text_repel(data = subset(rottraitpc_df, score == "species" & PC1 < 0 & PC2 >0), aes(PC1, PC2, label = label),
                           nudge_y = 0.05, nudge_x = -0.05, size = 3, min.segment.length = 10) +
  labs(x = "PC1", y = "PC2") +
  scale_x_continuous(limits = c(min(rottraitpc_df$PC1)-0.5, max(rottraitpc_df$PC1)+1)) +
  scale_y_continuous(limits = c(min(rottraitpc_df$PC2), max(rottraitpc_df$PC2)))

rottraits_pcafig <- cowplot::plot_grid(rottrts, rotspp, nrow = 1,axis = "l",
                                    rel_widths = c(1, 1.3),
                                    align = "v")
ggsave(plot = rottraits_pcafig, filename = "allspp21_reducedtraitsfig_varimaxrot.png", 
       path = "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/chapter/figs/prelim/",
       width = 8, height = 4, units = "in", scale = 1.5)


# -- TEST TRAIT INFLUENCE ON SPP ORDINATIONS -----
# using RDA and normalized trait values

ordsp_hel_seeded_traits <- ordsp_hel_seeded[names(ordsp_hel_seeded) %in% rownames(traitmat_seeded_norm)]
# check that order of trait matrix is same as order of vegdat
summary(rownames(traitmat_seeded_norm) == names(ordsp_hel_seeded_traits)) # yes


# calculate cwms
seeded_cwms <- functcomp(x = traitmat_seeded_norm, a = as.matrix(ordsp_hel_seeded_traits), CWM.type = "dom")
# pull out only traits
seeded_cwms_traits <- seeded_cwms[!names(seeded_cwms) %in% c("code4", "fxnl_grp", "nativity", "duration.abbr")]
# check correlations of CWMs
corrplot::corrplot(cor(seeded_cwms_traits), 
                   method = "square", 
                   order = "hclust")
# check distribution of CWMs (from transformed)
rownames_to_column(seeded_cwms_traits, "rowid") %>%
  gather(trait, val, names(.)[!grepl("rowid", names(.))]) %>%
  ggplot() +
  geom_density(aes(val)) +
  facet_wrap(~trait, scales = "free") # distributions looks like the transformed traits.. good candidate for rda

# to be sure, consider cwms from non-normalized trait (raw trait) data
seeded_cwms_rawtraits <- functcomp(x = traitmat_seeded, a = as.matrix(ordsp_hel_seeded_traits), CWM.type = "dom")
rownames_to_column(seeded_cwms_rawtraits, "rowid") %>%
  # drop character/factor traits
  dplyr::select(-c("code4", "fxnl_grp", "nativity", "duration.abbr")) %>%
  gather(trait, val, names(.)[!grepl("rowid", names(.))]) %>%
  # standardize so on comparable scales
  group_by(trait) %>%
  mutate(val = scale(val)) %>%
  ggplot() +
  geom_density(aes(val)) +
  facet_wrap(~trait, scales = "free") # i suppose it looks closer to the raw trait distributions, but not too terrible for running with cca
# ^but! reason to use normalized trait cwms is so scales aren't very off 

# run cca and rda
seeded_cwms_cca_envxcond <- cca(seeded_cwms_traits ~ herbicide + ppt_trt * nut_trt + Condition(block), data = wide_env_seeded)
summary(seeded_cwms_cca_envxcond)
seeded_cwms_rda_envxcond <- rda(seeded_cwms_traits ~ herbicide + ppt_trt * nut_trt + Condition(block),
                                scale = T,  # scale 'species' [trait] cols
                                data = wide_env_seeded)
summary(seeded_cwms_rda_envxcond)
plot(seeded_cwms_rda_envxcond, display = c("bp", "sp"), scaling = "species")
plot(seeded_cwms_cca_envxcond, display = c("sp"), scaling = "species")
plot(seeded_cwms_cca_envxcond, display = c("bp", "sp"), scaling = "symmetric")

# compare with cwms from raw trait dat
seeded_cwmsraw_rda_envxcond <- rda(seeded_cwms_rawtraits[!grepl("code4|fxnl|abbr|nati", names(seeded_cwms_rawtraits))] ~ herbicide + ppt_trt * nut_trt + Condition(block),
                                scale = T, data = wide_env_seeded)

summary(seeded_cwmsraw_rda_envxcond) # block explains 35% of data variance, environmental effects 17% -- which lie on 9 axes
# unconstrained axes fall on 18
# first two rda axes explain the most variation of the variation constraints explain
# BUT PC1 captures more variance than RDA1; PC2 and 3 explain more than RDA2
plot(seeded_cwmsraw_rda_envxcond, display = c("bp", "sp"), scaling = "species")

# run on reduced trait mod: SLA|Fresh|RMF|LDM|Height|SLA|root.diam|Fine.root.specif|volu
seeded_selectcwms_rda_envxcond <- rda(seeded_cwms_traits[grepl("SLA|Fresh|RMF|LDM|Height|SLA|root.diam|Fine.root.specific|volu", names(seeded_cwms_traits))] ~ herbicide + nut_trt * ppt_trt + Condition(block),
                                scale = T,  # scale 'species' [trait] cols
                                data = wide_env_seeded)
summary(seeded_selectcwms_rda_envxcond)
seeded_selectcwms_rda_envxcond
plot(seeded_selectcwms_rda_envxcond, display = c("bp", "sp"), scaling = "species", 
     #cex = 0.8,
     main = "Select trait CWMs ~ herbicide + nut_trt * ppt + Cond(block)\n(symmetric scaling)")
plot(seeded_cwms_rda_envxcond, display = c("bp", "sp"), scaling = "sites")
anova(seeded_selectcwms_rda_envxcond, by = "axis", permutations = how(blocks = wide_env_seeded$block)) # just rda 1 matters?
anova(seeded_selectcwms_rda_envxcond, by = "margin", permutations = how(blocks = wide_env_seeded$block)) # interaction matters
anova(seeded_selectcwms_rda_envxcond, by = "term", permutations = how(blocks = wide_env_seeded$block))
# ^unlike with species ID, ppt_trt is not signif (alone, interaction yes)
# so traits matter for how species respond to joint effects, but for species ID more responsive to direct effects


# -- OLD CODE -----
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

# plot all but nats in ggplot so can see plot trtments better
nmds_allbutnats_gg <- cbind(allbutnats_env2, allbutnats_nmds$points)

ggplot(nmds_allbutnats_gg, aes(MDS1, MDS2, col = herbicide)) +
  geom_label(aes(label = block)) +
  facet_grid(ppt_trt~nut_trt)

# one more nmds using fxnl grps -- compare experimental subplots only 
mds_nats_df <- subset(natlong, seedtrt != "Unseeded", select = c(plot:mon, code4,pct_cover)) %>%
  # choose max cov within season
  group_by(plot, herbicide, seedtrt, code4) %>%
  summarise(pct_cover = max(pct_cover), .groups = "drop_last") %>%
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

mds_nats_nmds <- metaMDS(subset(mds_nats_df, select = -c(plot:seedtrt, natplotid)), k = 2, trymax = 100)
mds_nats_nmds
plot(mds_nats_nmds)

# create bc matrix of all spp, only seededplots
allspp_bc <- vegdist(subset(mds_nats_df, select = -c(plot:seedtrt, natplotid)), method = "bray")
seeded_env <- subset(mds_nats_df, select = c(plot:seedtrt, natplotid)) %>%
  left_join(distinct(natlong[natlong$seedtrt != "Unseeded",], plot, block, ppt_trt, nut_trt, herbicide))
# check signif of treatments
adonis2(allspp_bc ~ herbicide * ppt_trt *nut_trt, data = seeded_env, permutations = how(blocks = seeded_env$block))

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Blocks:  seeded_env$block 
# Permutation: free
# Number of permutations: 199

adonis2(formula = allspp_bc ~ herbicide * nut_trt * ppt_trt, data = seeded_env, permutations = how(blocks = seeded_env$block))
adonis2(formula = allspp_bc ~nut_trt * ppt_trt *herbicide, data = seeded_env, permutations = how(blocks = seeded_env$block))
# > doesn't matter if change order of terms
# Df SumOfSqs      R2      F Pr(>F)   
# herbicide                  1   0.8831 0.04454 3.2606  0.005 **
# nut_trt                    2   1.4821 0.07475 2.7361  0.005 **
# ppt_trt                    2   0.6851 0.03455 1.2648  0.035 * 
# herbicide:nut_trt          2   0.6313 0.03184 1.1654  0.040 * 
# herbicide:ppt_trt          2   0.5236 0.02641 0.9667  0.125   
# nut_trt:ppt_trt            4   0.9785 0.04935 0.9032  0.170   
# herbicide:nut_trt:ppt_trt  4   0.5601 0.02825 0.5170  0.870   
# Residual                  52  14.0833 0.71031                 
# Total                     69  19.8270 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis2(formula = allspp_bc ~ herbicide + nut_trt * ppt_trt, data = seeded_env, permutations = how(blocks = seeded_env$block))

# Df SumOfSqs      R2      F Pr(>F)   
# herbicide          1   0.8831 0.04454 3.3910  0.005 **
# nut_trt            2   1.4821 0.07475 2.8456  0.005 **
# ppt_trt            2   0.6851 0.03455 1.3154  0.025 * 
# herbicide:nut_trt  2   0.6313 0.03184 1.2121  0.035 * 
# Residual          62  16.1455 0.81432                 
# Total             69  19.8270 1.0000

mds_nats_results <- left_join(mds_nats_df[c("plot", "herbicide", "seedtrt")], trtkey) %>%
  cbind(mds_nats_nmds$points)

mds_nats_spp <- data.frame(mds_nats_nmds$species) %>%
  mutate(code4 = rownames(.)) %>%
  left_join(distinct(subset(spplist, select = -c(species, compost_synonym)))) %>%
  mutate(nativity = gsub("Unknown", "Exotic", nativity))

ggplot(mds_nats_results, aes(MDS1, MDS2)) +
  geom_label(aes(label = block)) + # nmds separates on axis 1 more by hillpos than anything else
  facet_grid(ppt_trt~nut_trt)

ggplot(mds_nats_results, aes(MDS1, MDS2)) +
  geom_label(aes(label = block)) + # nmds separates on axis 1 more by hillpos than anything else
  facet_grid(ppt_trt~nut_trt)

ggplot(data = mds_nats_results, aes(MDS1, MDS2)) +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey") +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey") +
  geom_point(aes(shape = paste(block >2, herbicide)), alpha = 0.9, size = 3) +
  geom_text(data = subset(mds_nats_spp, !code4 %in% nats), aes(MDS1, MDS2, color = paste(nativity,fxnl_grp), label = code4), size = 2.5) +
  geom_text(data = subset(mds_nats_spp, code4 %in% nats), aes(MDS1, MDS2, color = paste(nativity,fxnl_grp), label = code4), size = 3) +
  scale_shape_manual(name = "hillpos + herbicide", values = c(16,1, 17, 2), labels = c("lo, herbicide", "lo, control", "hi, herbicide", "hi, control")) +
  scale_color_manual(name = NULL, values = c('Exotic Forb' = "red", 'Native Forb' = "purple", 'Exotic Grass' = "forestgreen", 'Native Grass' = "dodgerblue", "Exotic N-fixer" = "orange")) +
  facet_grid(ppt_trt~nut_trt)
ggsave("/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/compost_natdiv_spp_mds.pdf")

# make hulls (quickplot.. spider will be easier on eyes to show maybe)
hulldat <- data.frame()
for(h in unique(mds_nats_results$herbicide)){
  for(f in unique(mds_nats_results$fulltrt)){
    grphull <- mds_nats_results[mds_nats_results$fulltrt == f & mds_nats_results$herbicide == h, ][chull(mds_nats_results[mds_nats_results$fulltrt == f & mds_nats_results$herbicide == h, 
                                                                                                                          c("MDS1", "MDS2")]), ]
    # add to df
    hulldat <- rbind(hulldat, grphull)  
  }
}

# kns says show species in one panel, treatments in other

# make pretty names for spp nmds
mds_nats_spp <- mutate(mds_nats_spp, 
                       pretty_name = ifelse(code4 %in% nats, paste(paste0(substr(genus,1,1), "."), epithet),paste(nativity, fxnl_grp)),
                       pretty_name = gsub("Exotic", "Non-native", pretty_name),
                       pretty_name = gsub("Grass", "grasses", pretty_name),
                       pretty_name = gsub("Forb", "forbs", pretty_name),
                       pretty_name = gsub("N-fixer", "N-fixers", pretty_name),
                       pretty_name = gsub("Native fo", "Background native fo", pretty_name),
                       pretty_name = factor(pretty_name, 
                                            levels = c("Non-native forbs", "Non-native grasses", "Non-native N-fixers",
                                                       "Background native forbs", "E. californica", "N. maculata", 
                                                       "B. carinatus", "F. microstachys"))
                       )
unique(mds_nats_spp$pretty_name)
#mds_nats_spp$pretty_name <- with(mds_)

mds_spp <- ggplot(data = mutate(mds_nats_results, 
                                ppt_trt = factor(ppt_trt, levels = c("D", "XC", "W"))
                                #ppt_trt = relevel(ppt_trt, "D")
                                ), aes(MDS1, MDS2)) +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey") +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey") +
  #geom_polygon(data=hulldat, aes(fill=ppt_trt,group=ppt_trt, col = ppt_trt), alpha=0.35) +
  #geom_point(aes(shape = block<3), col = "grey30", alpha = 0.7, size = 3) +
  #geom_point(data = subset(mds_nats_spp, !code4 %in% nats), aes(MDS1, MDS2, color = pretty_name)) +
    geom_point(data = mds_nats_spp, aes(MDS1, MDS2, fill = pretty_name), pch = 21, col = "grey30", alpha = 0.85, size = 3) +
  #geom_text(data = subset(mds_nats_spp, !code4 %in% nats), aes(MDS1, MDS2, color = paste(nativity,fxnl_grp), label = code4), size = 3) +
  ggrepel::geom_text_repel(data = subset(mds_nats_spp, code4 %in% nats), 
                           aes(MDS1, MDS2, color = pretty_name, label = pretty_name), show.legend = FALSE, size = 4) +
  scale_color_manual(name = NULL, values = plant_cols2) +
    scale_fill_manual(name = NULL, values = plant_cols2) +
  #theme( legend.text = element_text(size = 10))
  #scale_color_manual(name = NULL, values = c("orchid", "seagreen2", "red", "purple4", "seagreen4")) +
  theme(legend.position = c(0.01,0.99),
        legend.text = element_text(size = 10),
        legend.justification = c("left", "top"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.1, "pt"))
                       #c('Exotic Forb' = "red", 'Native Forb' = "purple", 'Exotic Grass' = "forestgreen", 'Native Grass' = "dodgerblue", "Exotic N-fixer" = "orange"))

# plot hulls
mds_plots <- ggplot(data = mutate(mds_nats_results, 
                                  herbicide = factor(herbicide, levels = c("Non-herbicided", "Herbicided")),
              ppt_trt = factor(ppt_trt, levels = c("D", "XC", "W")),
              #ppt_trt = relevel(ppt_trt, "D"),
              nut_trt = factor(nut_trt, levels = c("C", "F", "N"), labels = c("Compost", "Fertilizer", "No amendment"))),
       aes(MDS1, MDS2, col = ppt_trt, fill = ppt_trt)) +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey") +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey") +
  geom_polygon(data=mutate(hulldat, 
                           herbicide = factor(herbicide, levels = c("Non-herbicided", "Herbicided")),
                           ppt_trt = factor(ppt_trt, levels = c("D", "XC", "W")),
                           nut_trt = factor(nut_trt, levels = c("C", "F", "N"), labels = c("Compost", "Fertilizer", "No amendment"))
                           ), 
               aes(fill=ppt_trt,group=ppt_trt, col = ppt_trt), alpha=0.35) +
  geom_point(aes(shape = block<3), col = "grey30", alpha = 0.7, size = 3) + #col = "grey30",
  geom_point(aes(shape = block<3), alpha = 0.7, size = 3) + #col = "grey30",
  #geom_point(alpha = 0.7, size = 3) +
  scale_shape_manual(name = "Hillslope", values = c(24, 25), labels = c("Up", "Down"), guides(color = "none")) +
  #scale_color_manual(name = NULL, values = c('Exotic Forb' = "red", 'Native Forb' = "purple", 'Exotic Grass' = "forestgreen", 'Native Grass' = "dodgerblue", "Exotic N-fixer" = "orange")) +
  scale_color_manual(name = "Precipitation\ntreatment", 
                     labels = c("D" = "Drought", "XC" = "Ambient", "W" = "Wet"),
                     values = RColorBrewer::brewer.pal(4,"Blues")[2:4]) +
  scale_fill_manual(name = "Precipitation\ntreatment", 
                    labels = c("D" = "Drought", "XC" = "Ambient", "W" = "Wet"), 
                    values = RColorBrewer::brewer.pal(4,"Blues")[2:4]) +
  #guides( = "none") +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 10)) +
  facet_grid(herbicide~nut_trt)

cowplot::plot_grid(mds_plots, mds_spp,
                   nrow = 1,
                   rel_widths = c(1.1,0.9))

ggsave(paste0(figpath,"natfxnl_seededspp_nmds2.pdf"), units = "in", width = 10.2, height = 5)


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




