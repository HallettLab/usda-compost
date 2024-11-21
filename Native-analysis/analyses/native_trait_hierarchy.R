# move trait modeling from native_ordinations to its own script
# but source native ordinations so have all the prepped compostition and trait data frame



# -- SETUP -----
# load needed libraries
library(tidyverse)
library(vegan)
library(FD)
library(glmmTMB)
library(lme4)
library(emmeans)
library(effects)
library(effectsize)
library(performance)

source("Native-analysis/analyses/native_ordinations.R")

# clean up environment
# > anything functions, any ordination analyses
rm(list = ls()[grepl("^fxnl|_rda|_ca|_cca", ls())])



# -- REVIEW DISTRIBUTIONS BEFORE MODELING ----
# >> see native_mixedmods.R for full workup on distributions
# > KNS advises model first without herbicide and environmental drivers; few control (to allow max variation across obs)
# > then add in manipulated drivers to test whether traits still important 

# how much of focals are 0%?
with(subset(tidysp0_seeded, code4 %in% nats), sapply(split(pct_cover, code4), function(x) table(x==0)))
with(subset(tidysp0_seeded, code4 %in% nats), aggregate(pct_cover ~ paste(code4, herbicide), FUN = function(x) table(x==0)))     
# ^ first col under pct_cover is present, second absent


# =======================================================
# -- HOW SIMILAR ARE TARGET SPP TO NEIGHBORS? -----
# compare target abundance to Fdis by plot (if targets similar to neighbors, expect negative relationship)
# can also try trait hierarchical distance as in Hao 2023 Sci of Total Env

# (1) Functional metrics on SPECIES TRAITS (not relative to target species) ------------ 
# calculate Fdis for plots

# community trait calculations based on presence absence
traitdb_seeded <-dbFD(x = traitmat_seeded_norm[,-1],
                      a = ordsp_pct_seeded[names(ordsp_hel_seeded_traits)],
                      #a = ordsp_hel_seeded_traits,
                      w.abun = F,  # not weighted by abundance
                      stand.x = T)

# abundance weighted
traitdb_seeded_abdwgt <-dbFD(x = traitmat_seeded_norm, a = ordsp_hel_seeded_traits,
                             w.abun = T, 
                             stand.x = T)

seeded_fdis <- data.frame(FDis = traitdb_seeded$FDis,
                          FDiv = traitdb_seeded$FDiv,
                          CWM = traitdb_seeded$CWM) %>% # can't use abundance weighted or else focal abundance *should* go up with fdis
  rownames_to_column("rowid") %>% 
  #data.frame() %>%
  merge(mutate(rowid= rownames(ordsp_hel_seeded_traits), data.frame(ordsp_hel_seeded_traits[nats[nats != "TRCI"]]))) %>%
  left_join(wide_env_seeded) %>%
  # join count dat
  left_join(rename_at(widesp_count[c("subplotID", "seedtrt", "herbicide", nats[nats != "TRCI"])], .vars = nats[nats != "TRCI"], function(x) paste0(x,"_ct")))

# how related are the different traits (ITO correlations?)
tidy_commfd_pres <- subset(seeded_fdis,select = c(rowid:CWM.Total.biomass.g, plot:herbicide)) %>%
  gather(traitmet, val, grep("^FDi|^CWM", names(.)))

ggplot(tidy_commfd_pres, aes(fulltrt, val, col = herbicide, group = herbicide)) +
  stat_summary() +
  coord_flip() +
  facet_wrap(~traitmet, scales = "free_x")

ggplot(seeded_fdis, aes(FDis, FDiv)) +
  geom_point(aes(col = nut_trt)) +
  geom_smooth(method = "lm") + # no difference in slopes, it is positively related
  facet_grid(herbicide ~ ppt_trt) # similar relationship


# -- femi -----
# does femi grow with things that are similar or dissimilar to it?
femi_fullnb_fdislme <- lme4::glmer.nb(FEMI_ct ~ FDis + ppt_trt + nut_trt + herbicide + (1|block), data = seeded_fdis)
femi_fullnb_fdis <- glmmTMB(FEMI_ct ~ FDis + ppt_trt + nut_trt + herbicide + (1|block), family = nbinom2(), data = seeded_fdis)
femi_rednb_fdis <- glmmTMB(FEMI_ct ~ FDis + herbicide + (1|block/xpos), family = nbinom2(), data = seeded_fdis) # can't include ypos or else doesn't converge
summary(femi_fullnb_fdislme); summary(femi_fullnb_fdis); summary(femi_rednb_fdis)
anova(femi_fullnb_fdis, femi_rednb_fdis) # full is no better; reduced penalized for nestedness

# visualize
ggplot(seeded_fdis, aes(FDis, FEMI_ct)) +
  geom_point(aes(shape = herbicide, col = ppt_trt), size = 2) +
  facet_wrap(~ ppt_trt) +
  #geom_smooth(method = "lm",linetype= 2, col = "grey30") +
  geom_smooth(aes(group = paste(herbicide, ppt_trt), linetype = herbicide), method = "lm",formula = 'y ~ poly(x, 2)') +
  scale_shape_manual(values = c(16, 1)) +
  ggtitle("FEMI: Fdis not weighted by spp abundances")


# -- brca ----
brca_fullnb_fdis <- glmmTMB(BRCA_ct ~ FDis * ppt_trt * nut_trt + (1|block), family = nbinom2(), data = seeded_fdis)
brca_rednb_fdis <- glmmTMB(BRCA_ct ~ FDis + ppt_trt + herbicide + (1|block), family = nbinom2(), data = seeded_fdis) # can't include xpos or else NaN for std error
summary(brca_fullnb_fdis); summary(brca_rednb_fdis)
anova(brca_fullnb_fdis, brca_rednb_fdis) # full better.. but 

ggplot(seeded_fdis, aes(FDis, BRCA_ct)) +
  geom_point(aes(shape = herbicide, col = nut_trt), size = 2) +
  facet_grid(. ~ ppt_trt, scales = "free") +
  geom_smooth(aes(col = nut_trt, group = nut_trt), method = "lm") + #formula = 'y ~ poly(x,2)'
  scale_shape_manual(values = c(1, 16)) +
  ggtitle("BRCA: Fdis not weighted by spp abundances")


# -- esca ----
esca_fullnb_fdis <- glmmTMB(ESCA_ct ~ FDis + ppt_trt + nut_trt + herbicide + (1|block), family = nbinom2(), data = seeded_fdis)
esca_rednb_fdis <- glmmTMB(ESCA_ct ~ FDis + nut_trt + herbicide + (1|block), family = nbinom2(), data = seeded_fdis) # can't include xpos or else NaN for std error
esca_subnb_fdis <- glmmTMB(ESCA_ct ~ FDis + nut_trt + herbicide + (1|block), family = nbinom2(), data = subset(seeded_fdis, FDis > 0.16)) # can't include xpos or else NaN for std error
summary(esca_fullnb_fdis); 
summary(esca_rednb_fdis)
summary(esca_subnb_fdis)
anova(esca_fullnb_fdis, esca_rednb_fdis) # full better.. but 

ggplot(seeded_fdis, aes(FDis, ESCA_ct)) +
  geom_point(aes(shape = herbicide), size = 2) +
  #facet_grid(. ~ ppt_trt, scales = "free") +
  geom_smooth(col = "grey50", method = "lm", linetype = 2) +
  geom_smooth(data = subset(seeded_fdis, FDis > 0.16), col = "grey30", method = "lm") +
  scale_shape_manual(values = c(1, 16)) +
  ggtitle("ESCA: Fdis not weighted by spp abundances")
# ^ if don't include outlier in herbicide drought, then esca does increase as traits in plot are more dispersed


# -- nema ----
nema_fullnb_fdis <- glmmTMB(NEMA_ct ~ FDis + ppt_trt + nut_trt + herbicide + (1|block), family = nbinom2(), data = seeded_fdis)
nema_rednb_fdis <- glmmTMB(NEMA_ct ~ FDis * nut_trt + herbicide + (1|block), family = nbinom2(), data = seeded_fdis) # can't include xpos or else NaN for std error
nema_subnb_fdis <- glmmTMB(NEMA_ct ~ FDis * nut_trt + ppt_trt + herbicide + (1|block), family = nbinom2(), data = seeded_fdis) # can't include xpos or else NaN for std error
summary(nema_fullnb_fdis); 
summary(nema_rednb_fdis)
summary(nema_subnb_fdis)
anova(nema_fullnb_fdis, nema_subnb_fdis) # with interaction maginally better 

ggplot(seeded_fdis, aes(FDis, NEMA_ct)) +
  geom_point(aes(shape = herbicide, col = ppt_trt), size = 2) +
  facet_grid(. ~ nut_trt, scales = "free") +
  geom_smooth(method = "lm",linetype= 2, col = "grey30") +
  scale_shape_manual(values = c(1, 16)) +
  ggtitle("NEMA: Fdis not weighted by spp abundances")



# -- multi-trait PC score ----
spp_pcscores <- traitpc_df$PC1[traitpc_df$score == "sites" & traitpc_df$label %in% traitmat_seeded_norm$code4]
names(spp_pcscores) <- traitpc_df$label[traitpc_df$label %in% traitmat_seeded_norm$code4]
pcdisp <- dbFD(x = spp_pcscores, ordsp_hel_seeded_traits) # don't standardize, don't weight by abundance

seeded_pcfdis <- data.frame(FDis = pcdisp$FDis, # FRic is the range of convex hulls that captures all spp when there is only 1 trait
                            CWM.PC = pcdisp$CWM$Trait) %>%
  rownames_to_column("rowid") %>% 
  #data.frame() %>%
  merge(mutate(rowid= rownames(ordsp_hel_seeded_traits), data.frame(ordsp_hel_seeded_traits[nats[nats!= "TRCI"]]))) %>%
  left_join(wide_env_seeded) %>%
  # join count dat
  left_join(rename_at(widesp_count[c("subplotID", "seedtrt", "herbicide", nats[nats != "TRCI"])], .vars = nats[nats != "TRCI"], function(x) paste0(x,"_ct")))


ggplot(seeded_pcfdis, aes(FDis, FEMI_ct)) +
  geom_point(aes(shape = herbicide), size = 2) +
  #facet_grid(herbicide ~ .) +
  geom_smooth(method = "lm",linetype= 2, col = "grey30") +
  scale_shape_manual(values = c(1, 16)) +
  labs("FEMI: PC1 Fdis not weighted by spp abundances or standardized")

ggplot(seeded_pcfdis, aes(CWM.PC, FEMI_ct)) +
  geom_point(aes(shape = herbicide, col = FDis), size = 2) +
  facet_grid(herbicide ~ .) +
  geom_smooth(method = "lm",linetype= 2, col = "grey30") +
  scale_shape_manual(values = c(1, 16)) +
  labs("FEMI: PC1 Fdis not weighted by spp abundances or standardized")

ggplot(seeded_pcfdis, aes(CWM.PC, FDis)) +
  geom_point(aes(col = log(FEMI_ct+0.001)), size = 3) +
  facet_grid(herbicide ~ ppt_trt) +
  geom_smooth(method = "lm",linetype= 2, col = "grey30") +
  labs("FEMI: PC1 Fdis not weighted by spp abundances or standardized")


summary(glmmTMB(FEMI_ct ~ CWM.PC + herbicide + (1|block/xpos/ypos), family = nbinom2(), data = seeded_pcfdis))

summary(glmmTMB(BRCA_ct ~ FDis + CWM.PC + ppt_trt + nut_trt + (1|block/herbicide), family = nbinom2(), data = seeded_pcfdis))
cor(seeded_pcfdis$FDis, seeded_pcfdis$CWM.PC) # not correlated

ggplot(seeded_pcfdis, aes(FDis, BRCA_ct)) +
  geom_point(aes(shape = herbicide), size = 2) +
  facet_grid(. ~ ppt_trt, scales = "free") +
  geom_smooth(method = "lm",linetype= 1, col = "grey30") +
  scale_shape_manual(values = c(1, 16)) +
  labs("BRCA: PC1 Fdis not weighted by spp abundances or standardized")

ggplot(seeded_pcfdis, aes(CWM.PC, BRCA_ct)) +
  geom_point(aes(col = ppt_trt), size = 2) +
  #facet_grid(. ~ ppt_trt) +
  geom_smooth(method = "lm",linetype= 1, col = "grey30") +
  scale_shape_manual(values = c(1, 16)) +
  labs("BRCA: PC1 Fdis not weighted by spp abundances or standardized")


summary(glmmTMB(ESCA_ct ~ CWM.PC + herbicide + (1|block/xpos/ypos), family = nbinom2(), data = seeded_pcfdis))
ggplot(seeded_pcfdis, aes(CWM.PC, ESCA_ct)) +
  geom_point(aes(shape = herbicide, col = ppt_trt), size = 2) +
  facet_grid(. ~ ppt_trt) +
  geom_smooth(method = "lm",linetype= 2, col = "grey30") +
  scale_shape_manual(values = c(1, 16)) +
  labs("FEMI: PC1 Fdis not weighted by spp abundances or standardized")


summary(glmmTMB(NEMA_ct ~ FDis * ppt_trt + (1|block/xpos), family = nbinom2(), data = seeded_pcfdis))
summary(glmmTMB(NEMA_ct ~ CWM.PC * ppt_trt + herbicide + (1|block/xpos), family = nbinom2(), data = seeded_pcfdis))


ggplot(seeded_pcfdis, aes(FDis, NEMA_ct)) +
  geom_point(aes(shape = herbicide), size = 2) +
  facet_grid(. ~ ppt_trt) +
  geom_smooth(method = "lm",linetype= 2, col = "grey30") +
  scale_shape_manual(values = c(1, 16)) +
  labs("FEMI: PC1 Fdis not weighted by spp abundances or standardized")

ggplot(seeded_pcfdis, aes(CWM.PC, NEMA_ct)) +
  geom_point(aes(shape = herbicide), size = 2) +
  facet_grid(. ~ ppt_trt) +
  geom_smooth(method = "lm",linetype= 2, col = "grey30") +
  scale_shape_manual(values = c(1, 16)) +
  labs("FEMI: PC1 Fdis not weighted by spp abundances or standardized")




# -- (2) CALCULATE TRAIT HIERARCHIES ----

# absolute trait hiearchy = [ti - tc] (where i = target, c = competitor, and t is the trait value)
# trait hierarchy = non-absolute value of that
# Hao et al 2023: Trait hierarchical distance is trait hiearchy weighted by abundance (except they did tc-ti)
# ^CWM on trait hierarchy


# try for FEMI:
# difference spp' traits from femi value, then scale for modeling (per Yin et al 2021)
# > can only include numeric traits here
traits_normscale <- subset(traitmat_seeded_norm, select = c(code4:Total.biomass.g)) %>%
  # add in multi-trait pc score
  left_join(data.frame(code4 = names(spp_pcscores), PC1 = spp_pcscores)) %>%
  # scale traits and PC scoer (mean center, divided by 1 sd)
  mutate_at(.vars = names(.)[!names(.) == "code4"], function(x) as.numeric(scale(x))) %>%
  gather(trait, value, 2:ncol(.)) %>%
  spread(code4, value)

# long form all but femi to diff
femi_traitdiffs <- traits_normscale %>%
  gather(code4, value, !grep("trait|FEMI", names(.))) %>%
  mutate(femi_diff = FEMI - value,
         femi_absdiff = abs(FEMI - value))

summary(femi_traitdiffs)
ggplot(femi_traitdiffs) +
  geom_vline(aes(xintercept = 0), col = "red") +
  geom_density(aes(x = femi_diff)) + # on using scale: these traits are already normalized, hard to interpret if scale?
  facet_wrap(~trait, scales = "free")
# interp for hiearchy:
# relative to femi, values in positive direction (to the right of 0) should be bigger, taller, longer, etc.
# value to the left are smaller, shorter, thinner, lighter
# the abs diff: bigger number = more different from femi trait (more dissimilar)


# function to compile functional metrics
funchier <- function(sp = "FEMI", diss = F){
  df1 <- traits_normscale
  # rename target as focal for generic functions
  names(df1)[names(df1) == sp] <- "focal"
  # take diffs
  focal_traitdiffs <- df1 %>%
    gather(code4, value, !grep("trait|focal", names(.))) %>%
    mutate(focal_diff = focal - value) # target minus neighbor value (everything is relative to target)
  # if dissimilarity is TRUE then use absolute
  if(diss){
    focal_traitdiffs$focal_diff <- abs(focal_traitdiffs$focal_diff)
  }
  
  # use actual difference for hierarchy, choose traits used in PCA
  focal_hier <- subset(focal_traitdiffs, 
                       grepl(paste0("^PC|", str_flatten(traitpc_df$label[traitpc_df$score == "species"], collapse = "|")), trait), 
                       select = c(trait, code4, focal_diff)) %>%
    spread(trait, focal_diff) %>%
    data.frame()
  # set code4 to rownames
  rownames(focal_hier) <- focal_hier$code4
  
  # calculate functional db without PC
  hier_funcdb <- dbFD(x = focal_hier[!grepl("code|^PC", names(focal_hier))], # drop code 4 in first col
                      a = ordsp_hel_seeded_traits[!names(ordsp_hel_seeded_traits) == sp],
                      w.abun = T, # weight centroids and distances by abundance
                      stand.x = F # don't standardize? (because already standardized)
  )
  
  # repeat for PC score
  pc_funcdb <- dbFD(x = focal_hier["PC1"], # drop code 4 in first col
                    a = ordsp_hel_seeded_traits[!names(ordsp_hel_seeded_traits) == sp],
                    w.abun = T, # weight centroids and distances by abundance
                    stand.x = F # don't standardize? (because already standardized)
  )
  
  # compile functional metrics with abundance of focal sp
  hier_df <- data.frame(FDis = hier_funcdb$FDis, 
                        FDiv = hier_funcdb$FDiv,
                        FRic = hier_funcdb$FRic, # convex hull volume, can't be weighted
                        CWM = hier_funcdb$CWM,
                        CWM.PC = pc_funcdb$CWM,
                        PC.FDis = pc_funcdb$FDis,
                        PC.FRic = pc_funcdb$FRic) %>%
    rename(PC.CWM = PC1) %>%
    rownames_to_column("rowid") %>% 
    merge(mutate(rowid= rownames(ordsp_hel_seeded_traits), data.frame(ordsp_pct_seeded[sp]))) %>%
    left_join(wide_env_seeded) %>%
    data.frame()
  
  #return
  return(hier_df)
  
}

# calculate dissimilarity and hierarchy values for each species
femi_hier <- funchier()
femi_diss <- funchier(diss = T)
# brca
brca_hier <- funchier(sp = "BRCA")
brca_diss <- funchier(sp = "BRCA", diss = T)
# esca
esca_hier <- funchier(sp = "ESCA")
esca_diss <- funchier(sp = "ESCA", diss = T)
# nema
nema_hier <- funchier(sp = "NEMA")
nema_diss <- funchier(sp = "NEMA", diss = T)


# rbind all then join plot info
all_hier <- rbind(cbind(spp = "FEMI", rename(femi_hier, cover = FEMI)),
                  cbind(spp = "BRCA", rename(brca_hier, cover = BRCA)),
                  cbind(spp = "ESCA", rename(esca_hier, cover = ESCA)),
                  cbind(spp = "NEMA", rename(nema_hier, cover = NEMA))
) %>%
  left_join(subset(target_long, select = -c(xpos, ypos)))

all_diss <- rbind(cbind(spp = "FEMI", rename(femi_diss, cover = FEMI)),
                  cbind(spp = "BRCA", rename(brca_diss, cover = BRCA)),
                  cbind(spp = "ESCA", rename(esca_diss, cover = ESCA)),
                  cbind(spp = "NEMA", rename(nema_diss, cover = NEMA))
) %>%
  left_join(subset(target_long, select = -c(xpos, ypos)))


# pull trait cols to test
test_traits <- names(all_diss)[grepl("CWM|FDis", names(all_diss))]
# remove fresh leaf mass
test_traits <- test_traits[!grepl("Fresh.leaf|PC.FD", test_traits)]


# make trait dats long for storing glmms in list format
all_hier_long <- cbind(type = "hierarchy", all_hier) %>%
  data.frame() %>%
  subset(select = -c(grep("FRic|FDiv|Fresh|PC.FD", names(.)))) %>%
  gather(met, val, test_traits)

all_diss_long <- cbind(type = "dissimilarity", all_diss) %>%
  data.frame() %>%
  subset(select = -c(grep("FRic|FDiv|Fresh|PC.FD", names(.)))) %>%
  gather(met, val, test_traits)



# spp trait dissimilarities -----
form1 <- as.formula(count ~ val * (herbicide + nut_trt + ppt_trt) + (1|block/wholeplotID))
form2 <- as.formula(count ~ val * (herbicide + nut_trt + ppt_trt) + (1|block))
inflform <- as.formula(~1)
get_glmms <- function(df = all_hier_long, s = "FEMI", form = form1, fam = "nbinom2", zform = NA){
  tempdat <- subset(df, spp == s)
  if(length(zform)>1){
  glmm_list <- lapply(split(tempdat[c("count", "herbicide", "nut_trt", "ppt_trt", "met", "val", "block", "wholeplotID")], tempdat$met), 
                                 function(x) m <- glmmTMB(formula = form, 
                                                          ziformula = zform, 
                                                          data = x, family = fam))
  }else{
    glmm_list <- lapply(split(tempdat[c("count", "herbicide", "nut_trt", "ppt_trt", "met", "val", "block", "wholeplotID")], tempdat$met), 
                        function(x) m <- glmmTMB(formula = form, 
                                                 data = x, family = fam))
  }
  names(glmm_list) <- unique(tempdat$met)
  return(glmm_list)
  
}

femi_hier_glmms <- get_glmms()
femi_diss_glmms <- get_glmms(df = all_diss_long)
femi_diss_glmms_nb1 <- get_glmms(df = all_diss_long, fam = "nbinom1")
femi_diss_glmms_pois <- get_glmms(df = all_diss_long, fam = "poisson")
femi_diss_glmms_zip <- get_glmms(df = all_diss_long, zform = inflform, fam = "poisson")
femi_diss_glmms_zinb1 <- get_glmms(df = all_diss_long, zform = inflform, fam = "nbinom1") # convergence issue
femi_diss_glmms_zinb2 <- get_glmms(df = all_diss_long, zform = inflform, fam = "nbinom2")
femi_diss_glmms_hurd1 <- get_glmms(df = all_diss_long, zform = inflform, fam = "truncated_nbinom1") # hurdle
femi_diss_glmms_noenv <- get_glmms(df = all_diss_long, form = as.formula(count ~ val + herbicide + (1|block/nut_trt/ppt_trt)))
femi_diss_glmms_addval <- get_glmms(df = all_diss_long, form = as.formula(count ~ herbicide * (nut_trt + ppt_trt) +val + (1|block)))
femi_diss_glmms_addvalp <- get_glmms(df = all_diss_long, form = as.formula(count ~ herbicide * (nut_trt + ppt_trt) +val + (1|block)), fam = "poisson")
femi_diss_glmms_addvalnb1 <- get_glmms(df = all_diss_long, form = as.formula(count ~ herbicide * (nut_trt + ppt_trt) +val + (1|block)), fam = "nbinom1")
femi_hier_glmms_addval <- get_glmms(form = as.formula(count ~ herbicide * (nut_trt + ppt_trt) +val + (1|block/wholeplotID)))



DHARMa::plotResiduals(femi_hier_glmms$CWM.Coarse.root.diameter.mm)

lapply(femi_hier_glmms, function(x)plot(simulateResiduals(x))) # LDMC major issues, Height has issues for 25th q, PC.CWM has a poly shape      

# -- EVAL FITS -----
# just for femi dissimilarity, then move on to others
lapply(femi_diss_glmms, summary)
lapply(femi_diss_glmms, testResiduals) # quantile resid issues in many. poly shape


# repeat for nbinom1
lapply(femi_diss_glmms_nb1, summary)
lapply(femi_diss_glmms_nb1, testResiduals)
par(mfrow = c(3,3))
# compare residuals
lapply(femi_diss_glmms, plotResiduals) # nbinom2 issues for coarse root diam, root volum, sla, pc.cwm
lapply(femi_diss_glmms_nb1, plotResiduals) # issues for CWM.PC and LDMC, a little for height
lapply(femi_diss_glmms_pois, plotResiduals) # poisson fit is definitely worse
lapply(femi_diss_glmms_hurd1, plotResiduals) # hurdle seems worse
lapply(femi_diss_glmms_zinb1, plotResiduals) # convergence error is on PC.CWM
lapply(femi_diss_glmms_zinb2, plotResiduals)
lapply(femi_diss_glmms_zip, plotResiduals) # this is not too bad
lapply(femi_diss_glmms_noenv, plotResiduals) # no. accounting for env important
lapply(femi_diss_glmms_addval, plotResiduals) # this one is good! no interaction between trait diss and treatments
lapply(femi_diss_glmms_addvalp, plotResiduals) # no. negbin is better.
lapply(femi_diss_glmms_addvalnb1, plotResiduals) # nbinom2 is slightly better than nbinom1
lapply(femi_hier_glmms_addval, plotResiduals) # slightly off for 2nd to last

lapply(femi_diss_glmms_zip, summary)
lapply(femi_diss_glmms_zip, testDispersion)
lapply(femi_diss_glmms_zinb1, testDispersion)
lapply(femi_diss_glmms_nb1, testDispersion) # seems like better dispersion using nbinom1 than zip..
lapply(femi_diss_glmms_addval, testDispersion) 
lapply(femi_diss_glmms_addvalnb1, testDispersion) 

lapply(femi_diss_glmms_addval, summary) # SLA dissimilarity helps FEMI
lapply(femi_hier_glmms_addval, summary) # hierarchy: SLA**, PC, RMF (negative), Coarse root diam


# test hierarchy for femi


plot(emmeans(femi_diss_glmms_zip$CWM.Root.volume.cm3, ~ nut_trt + ppt_trt | herbicide + val))

ggplot(subset(femi_diss), aes(CWM.SLA.cm2.g, FEMI)) +
  geom_point(aes(col = herbicide)) +
  geom_smooth(method = "lm") +
  facet_wrap(~herbicide)

ggplot(subset(femi_hier), aes(CWM.Height.cm, FEMI)) +
    geom_point(aes(col = herbicide)) +
    geom_smooth(method = "lm") +
  facet_wrap(~nut_trt)
  



# -- CONTINUE W OTHER SPP ------
form3 <- as.formula(count ~ herbicide * (ppt_trt + nut_trt) + val + (1|block/wholeplotID))
form4 <- as.formula(count ~ herbicide * (ppt_trt + nut_trt + val) + (1|block/wholeplotID))

# -- brca -----
brca_hier_glmms1 <- get_glmms(s = "BRCA") # best
brca_hier_glmms <- get_glmms(s = "BRCA", form = form3)
brca_hier_glmms4 <- get_glmms(s = "BRCA", form = form4) 
brca_diss_glmms <- get_glmms(df = all_diss_long, s = "BRCA", form = form3) # two model convergence problems
brca_diss_glmms1 <- get_glmms(df = all_diss_long, s = "BRCA") 

lapply(brca_hier_glmms1, plotResiduals)
lapply(brca_hier_glmms, plotResiduals)
lapply(brca_hier_glmms4, plotResiduals) 
lapply(brca_diss_glmms, plotResiduals)
lapply(brca_diss_glmms1, plotResiduals)

AIC(brca_hier_glmms1$FDis, brca_hier_glmms$FDis, brca_hier_glmms4$FDis) # AIC lower for simpler model

lapply(brca_hier_glmms1, summary) # rmf matters (negatively) except pos in fert, comp, herb
# > PC matters (pos) in herbicide, RMF, fine root, coarse root,
lapply(brca_hier_glmms, summary)
lapply(brca_hier_glmms4, summary)
lapply(brca_diss_glmms, summary) # no trait advantage brca
lapply(brca_diss_glmms1, summary) # coarse root diameter (comp) seems to matter for brca; SLA; RMf F and C; LDMC compost


# -- esca -----
form5 <- as.formula(count ~ herbicide + nut_trt + ppt_trt + val + (1|block/wholeplotID))

esca_hier_glmms <- get_glmms(s = "ESCA") # 6 model convergence problems
esca_hier_glmmsnb1 <- get_glmms(s = "ESCA", fam = "nbinom1") # lots of issue w convergence
esca_hier_glmmsp <- get_glmms(s = "ESCA", fam = "poisson")
esca_hier_glmms3 <- get_glmms(s = "ESCA", form = form3)
esca_hier_glmms4<- get_glmms(s = "ESCA", form = form4)
esca_hier_glmms5<- get_glmms(s = "ESCA", form = form5)
esca_hier_glmms_zinb2 <- get_glmms(s = "ESCA", zform = inflform) # convergence issues

esca_diss_glmms <- get_glmms(df = all_diss_long, s = "ESCA")
esca_diss_glmms3 <- get_glmms(df = all_diss_long, s = "ESCA", form = form3)
esca_diss_glmms4 <- get_glmms(df = all_diss_long, s = "ESCA", form = form4)

lapply(esca_hier_glmms, plotResiduals)
lapply(esca_hier_glmmsnb1, plotResiduals) # no
lapply(esca_hier_glmmsp, plotResiduals) # nb2 fit still better
lapply(esca_hier_glmms3, plotResiduals) # no
lapply(esca_hier_glmms4, plotResiduals) # no
lapply(esca_hier_glmms5, plotResiduals) # no
lapply(esca_hier_glmms_zinb2, plotResiduals) # normal nb2 still least problematic

lapply(esca_diss_glmms, plotResiduals) # problems with 2nd trait only
lapply(esca_diss_glmms3, plotResiduals) # probems with several others but less severe
lapply(esca_diss_glmms4, plotResiduals) # definite no

lapply(esca_hier_glmms, summary) # neg in drought only (interactions): fine root sl, ldmc (weak), sla, trait pc
lapply(esca_diss_glmms, summary) # LDMC (positive), fine root (pos)

lapply(esca_hier_glmms5, summary)

lapply(esca_diss_glmms3, summary)


nema_hier_glmms <- get_glmms(s = "NEMA") # 4 model convergence issues
nema_diss_glmms <- get_glmms(df = all_diss_long, s = "NEMA")
lapply(nema_diss_glmms, summary)
lapply(nema_diss_glmms, plotResiduals)


# (2.1) FEMI trait hierarchy --------
# troubleshoot FEMI dissimilarity
lapply(femi_diss_glmms, summary)
lapply(femi_diss_glmms, testResiduals) # quantile resid issues in many. poly shape
lapply(femi_diss_glmms, plotQQunif) # passes dispersion test
femi_diss_glmms_form1nb1 <- get_glmms(df = all_diss_long, fam = "nbinom1")
femi_diss_glmms_form1nb1 <- get_glmms(df = all_diss_long, fam = "nbinom1") 
femi_diss_glmms_form1genp <- get_glmms(df = all_diss_long, fam = "genpois") 
femi_diss_glmms_form1zip <- get_glmms(df = all_diss_long, fam = "poisson") 
femi_diss_glmms_form1zinb1 <- get_glmms(df = all_diss_long, fam = "nbinom1") 

testzip <- lapply()
lapply(femi_hier_glmms, plotResiduals) # Coarse root has issues, Height S curve issues, LDMC better but still issue, RMF has issues (poly shape), PC.CWM okay     
lapply(femi_diss_glmms_form1zip, plotResiduals) # generally best looking resids
lapply(femi_diss_glmms_form1zip, plotQQunif)
lapply(femi_diss_glmms_form1zinb1, plotResiduals) # good except for PC.CWM
lapply(femi_diss_glmms_form1nb1, plotResiduals) # LDMC major issues, Height has issues for 25th q, PC.CWM has a poly shape      
lapply(femi_diss_glmms_form1genp, plotResiduals) # still a few issue but better than either nb fits
lapply(femi_diss_glmms_form1genp, DHARMa::testDispersion)
lapply(femi_diss_glmms_form1nb1, DHARMa::testDispersion)
lapply(femi_diss_glmms, DHARMa::testDispersion)
lapply(femi_diss_glmms, DHARMa::plotSimulatedResiduals)
plot(simulateResiduals(femi_diss_glmms$CWM.LDMC))

# compare for height
AIC(femi_diss_glmms$CWM.Height.cm, femi_diss_glmms_form1genp$CWM.Height.cm, femi_diss_glmms_form1nb1$CWM.Height.cm,
    femi_diss_glmms_form1zip$CWM.Height.cm, femi_diss_glmms_form1zinb1$CWM.Height.cm)
# general poisson has lowest AIC


femi_diss_traitxtrts <- lapply(split(all_diss_long[c("FEMI_ct", "herbicide", "nut_trt", "ppt_trt", "met", "val", "block", "wholeplotID")], femidiss_long$met), 
                               function(x) m <- glmmTMB(formula = FEMI_ct ~ val * (herbicide + nut_trt + ppt_trt) + (1|block/wholeplotID), 
                                                   data = x, family = "nbinom2"))
names(femi_diss_traitxtrts) <- test_taits
femi_hier_traitxtrts <- lapply(split(femihier_long[c("FEMI_ct", "herbicide", "nut_trt", "ppt_trt", "met", "val", "block", "wholeplotID")], femihier_long$met), 
                               function(x) m <-(glmmTMB(formula = FEMI_ct ~ val * (herbicide + nut_trt + ppt_trt) + (1|block/wholeplotID), 
                                                   data = x, family = "nbinom2")))
names(femi_hier_traitxtrts) <- test_traits

lapply(femi_hier_traitxtrts, function(x) broom.mixed::tidy(x, conf.int = T))
lapply(femi_diss_traitxtrts, function(x) broom.mixed::tidy(x, conf.int = T))
lapply(femi_hier_traitxtrts, coef)

plot(emmeans::emmeans(femi_diss_traitxtrts$CWM.LDMC, ~ herbicide + ppt_trt + nut_trt | val))
plot(emmeans::emmeans(femi_hier_traitxtrts$CWM.LDMC, ~ ppt_trt + nut_trt | val +herbicide, ))


# test dissimilarity first
# > does Fdisp for traits and Fdis PCA give you the same info?
plot(femi_diss$FDis ~ femi_diss$PC.FDis) # not quite

cor.test(~femi_diss$FDis + femi_diss$PC.FDis) # rho = .3445938 95% CI 0.1192778 0.5361551


summary(glmmTMB(FEMI_ct ~ FDis * (nut_trt + ppt_trt + herbicide) + (1|block), data = femi_diss, family = "nbinom2"))



# are neighbor species morfdisp()# are neighbor species more similar to or different from FEMI?
ggplot(femi_hier, aes(FDis, FEMI, col = herbicide)) + #
  geom_point(aes(col = herbicide)) +
  geom_smooth(method = "lm")
  facet_grid(nut_trt ~ ppt_trt)

ggplot(femi_hier, aes(CWM.SLA.cm2.g, FEMI, col = herbicide)) +
  geom_point(aes(col = herbicide)) +
  facet_grid(ppt_trt ~ .)

# plot all CWMs, ignore treatments
gather(femi_hier, traitmet, val, grep("^F[A-Z][a-z]|^CWM|^PC", names(femi_hier))) %>%
  subset(grepl("Herb", herbicide)) %>%
  #subset(FEMI > 0) %>%
  ggplot(aes(val, FEMI_ct)) +
  #ggplot(aes(val, log(FEMI + 0.0001))) +
  #ggplot(aes(val, scale(FEMI))) +
  #geom_vline(aes(xintercept = 0), col = "black", lwd = 1.5) +
  geom_point(aes(col = ppt_trt), size = 3, alpha = .75) + #col = 
  geom_smooth(aes(col = ppt_trt, fill = ppt_trt), method = "lm") + # formula = "y ~ poly(x,2)" 
  stat_smooth(aes(group = herbicide, lty = herbicide),col = "grey30", method = "lm", formula = "y ~ poly(x,2)") +
  #stat_smooth(aes(),col = "grey30", method = "lm") +
  scale_color_manual(values = c("seagreen2", "chocolate", "blue")) +
  scale_fill_manual(values = c("seagreen2", "chocolate", "blue")) +
  ggtitle("FEMI abundance ~ functional trait hierarchies for select traits and multi-trait PC1") + 
  facet_wrap(~traitmet, scales = "free_x") +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(face = "bold"))


# assess functional dissimilarity (does femi grow when it's different from it's neighbors?)
gather(femi_diss, traitmet, val, grep("^F[A-Z][a-z]|^CWM|^PC", names(femi_diss))) %>%
  subset(grepl("Herb", herbicide)) %>%
  #subset(FEMI > 0) %>%
  ggplot(aes(val, FEMI_ct)) +
  #ggplot(aes(val, log(FEMI+0.0001))) +
  #ggplot(aes(val, scale(FEMI))) +
  #geom_vline(aes(xintercept = 0), col = "black", lwd = 1.5) +
  geom_point(aes(col = ppt_trt), size = 3, alpha = .75) + #col = 
  geom_smooth(aes(col = ppt_trt, fill = ppt_trt, lty = herbicide), se = T, method = "lm") + #formula = "y ~ poly(x,2)"
  stat_smooth(aes(),col = "grey30", method = "lm", formula = "y ~ poly(x,2)") +
  #stat_smooth(aes(),col = "grey30", method = "lm") +
  scale_color_manual(values = c("seagreen2", "chocolate", "blue")) +
  scale_fill_manual(values = c("seagreen2", "chocolate", "blue")) +
  ggtitle("FEMI abundance ~ functional trait dissimilarity for select traits and multi-trait PC1") + 
  facet_wrap(~traitmet, scales = "free_x") +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(face = "bold"))



summary(lme4::glmer.nb(FEMI_ct ~ poly(PC.CWM,2) + (1|block:nut_trt:ppt_trt), data = femi_hier)) #+ herbicide
summary(glmmTMB::glmmTMB((FEMI/100) ~ poly(PC.CWM,2) + (1|block/nut_trt/ppt_trt), family = ordbeta(), data = femi_hier))
summary(glmmTMB::glmmTMB((FEMI/100) ~ poly(PC.CWM,2) + (1|block), family = ordbeta(), data = femi_diss))
summary(glmmTMB::glmmTMB(formula = (FEMI/100) ~ poly(PC.CWM,2) + herbicide + (1|block), 
                         family = ordbeta(), data = femi_diss))
summary(glmmTMB::glmmTMB((FEMI/100) ~ poly(PC.CWM,2) + herbicide + (1|block), family = ordbeta(), data = femi_diss))
summary(glmmTMB::glmmTMB(FEMI_ct ~ poly(PC.CWM,2) + (1|block), family = nbinom2(), data = femi_hier))
summary(glmmTMB::glmmTMB(FEMI_ct ~ poly(PC.CWM,2) + herbicide + (1|block), family = nbinom2(), data = femi_hier))
summary(glmmTMB::glmmTMB(FEMI_ct ~ PC.CWM + herbicide + (1|block/xpos), family = nbinom2(), data = femi_hier))
summary(glmmTMB::glmmTMB(FEMI_ct ~ FDis + ppt_trt + herbicide + (1|block/xpos), family = nbinom2(), data = femi_hier))
summary(glmmTMB::glmmTMB(FEMI_ct ~ poly(CWM.Root.volume.cm3,2) + ppt_trt + herbicide + (1|block/xpos), family = nbinom2(), data = femi_hier))


summary(glmmTMB(FEMI_ct ~ CWM.Coarse.root.diameter.mm + CWM.Fine.root.specific.length.cm.g + 
                  CWM.Height.cm + CWM.LDMC +CWM.SLA.cm2.g + CWM.RMF + CWM.Root.volume.cm3 + (1|block),
                family = nbinom2, data = femi_hier))

summary(glmmTMB(FEMI_ct ~  (CWM.Height.cm + CWM.SLA.cm2.g) * herbicide + (1|block/xpos),
                family = nbinom2, data = femi_hier))

summary(glmmTMB(FEMI_ct ~  (CWM.Height.cm + CWM.SLA.cm2.g) * herbicide + nut_trt + (1|block),
                family = nbinom2, data = femi_hier))

summary(glmmTMB(FEMI_ct ~  (CWM.Height.cm) * herbicide + nut_trt + (1|block), #+ CWM.SLA.cm2.g
                family = nbinom2, data = femi_diss))

summary(glmer.nb(FEMI_ct ~  CWM.Height.cm * herbicide + nut_trt + (1|block), #+ CWM.SLA.cm2.g
                data = femi_hier))


ggplot(femi_hier, aes(CWM.Height.cm, FEMI_ct)) +
  geom_point(aes(col = herbicide), size = 3, alpha = .75) + #col = 
  geom_smooth(aes(col = herbicide, fill = herbicide), method = "lm") + #  , formula = "y ~ poly(x,2)"
  scale_color_manual(values = c("seagreen2", "chocolate")) +
  scale_fill_manual(values = c("seagreen2", "chocolate")) +
  ggtitle("FEMI abundance ~ functional trait hierarchies for select traits and multi-trait PC1")

summary(glmmTMB(FEMI_ct ~  FDis * herbicide * nut_trt + (1|block),
                family = nbinom2, data = femi_diss))


# (2.2) BRCA trait hierarchy ----
brca_hier <- funchier(sp = "BRCA")
brca_diss <- funchier(sp = "BRCA", diss = T)
# stick with count to be consistent with gllvm community modeling
brca_hier <- left_join(brca_hier, rename(widesp_count[c("subplotID", "seedtrt", "herbicide", "BRCA")], BRCA_ct = BRCA))
brca_diss <- left_join(brca_diss, rename(widesp_count[c("subplotID", "seedtrt", "herbicide", "BRCA")], BRCA_ct = BRCA))


# quick check nb vs zinb
brca_fullmod <- glmmTMB(BRCA_ct ~  FDis + herbicide + ppt_trt + nut_trt + (1|block),family = nbinom2, data = brca_hier)
brca_fullmod_zinb <- glmmTMB(BRCA_ct ~  FDis + herbicide + ppt_trt + nut_trt + (1|block),
                             zi = ~ ppt_trt + nut_trt + herbicide,
                             family = nbinom2, data = brca_hier)

brca_fullmod <- glmer.nb(BRCA_ct ~  FDis + herbicide + ppt_trt + nut_trt + (1|block), data = brca_hier)
brca_fullmod_zinb <- glmmTMB(BRCA_ct ~  FDis + herbicide + ppt_trt + nut_trt + (1|block),
                             zi = ~ ppt_trt + nut_trt + herbicide,
                             family = nbinom2, data = brca_hier)

nonnest2::vuongtest(brca_fullmod, brca_fullmod_zinb)


brca_fullmodlme <- glmer.nb(BRCA_ct ~  FDis + herbicide + ppt_trt + nut_trt + (1|block),family = nbinom2, data = brca_hier)
brca_redmod <- glmmTMB(BRCA_ct ~  FDis + herbicide + ppt_trt + (1|block/xpos),family = nbinom2, data = brca_hier)
brca_redmod_diss <- glmmTMB(BRCA_ct ~  FDis + herbicide + ppt_trt + (1|block/xpos),family = nbinom2, data = brca_diss)
summary(brca_fullmod); summary(brca_fullmodlme)
summary(brca_redmod); summary(brca_redmod_diss)
anova(brca_fullmod, brca_redmod) # same


brca_fullmod_pc <- glmmTMB(BRCA_ct ~  PC.CWM + herbicide + ppt_trt + nut_trt + (1|block),family = nbinom2, data = brca_hier)
brca_redmod_pc <- glmmTMB(BRCA_ct ~  PC.CWM + herbicide + ppt_trt + (1|block/xpos),family = nbinom2, data = brca_hier)
brca_redmod_pcdiss <- glmmTMB(BRCA_ct ~  poly(PC.CWM,2) + herbicide + ppt_trt + (1|block/xpos),family = nbinom2, data = brca_diss)
summary(brca_fullmod_pc); summary(brca_redmod_pc)
summary(brca_redmod_pcdiss)

brca_redmod_ldmc <- glmmTMB(BRCA_ct ~  CWM.LDMC + herbicide + ppt_trt + (1|block/xpos),family = nbinom2, data = brca_hier)
brca_redmod_rootvol <- glmmTMB(BRCA_ct ~  CWM.Root.volume.cm3 + herbicide + (1|block/xpos/ypos),family = nbinom2, data = brca_hier)
brca_redmod_ldmc <- glmmTMB(BRCA_ct ~  poly(CWM.LDMC,2) + herbicide + (1|block/xpos/ypos),family = nbinom2, data = brca_hier)
brca_redmod_ldmc_full <- glmmTMB(BRCA_ct ~  poly(CWM.LDMC,2) + ppt_trt + herbicide + (1|block/xpos),family = nbinom2, data = brca_hier)
brca_redmod_ldmc_fulldiss <- glmmTMB(BRCA_ct ~  poly(CWM.LDMC,2) + ppt_trt + herbicide + (1|block/xpos),family = nbinom2, data = brca_diss)
summary(brca_redmod_ldmc); summary(brca_redmod_rootvol); 
summary(brca_redmod_ldmc); summary(brca_redmod_ldmc_full)
car::Anova(brca_redmod_ldmc_fulldiss, type = 2)

# plot all CWMs, ignore treatments
# trait hierachies
gather(brca_hier, traitmet, val, grep("^F[A-Z][a-z]|^CWM|^PC", names(brca_hier))) %>%
  #subset(grepl("Herb", herbicide)) %>%
  #subset(BRCA > 0) %>%
  ggplot(aes(val, BRCA_ct)) +
  #ggplot(aes(val, log(BRCA+0.0001))) +
  #ggplot(aes(val, scale(BRCA))) +
  geom_vline(aes(xintercept = 0), col = "black", lwd = 1.5) +
  geom_point(aes(col = ppt_trt), size = 3, alpha = .75) + 
  #geom_point(aes(col = ppt_trt, shape = nut_trt), size = 3, alpha = .75) + #col = 
  geom_smooth(aes(col = ppt_trt, fill = ppt_trt),lty = 2, se = F, method = "lm") +
  geom_smooth(aes(col = ppt_trt, fill = ppt_trt), method = "lm", formula = "y ~ poly(x,2)") +
  #stat_smooth(aes(),col = "grey30", method = "lm", formula = "y ~ poly(x,2)") + #group = herbicide, lty = herbicide
  #stat_smooth(aes(),col = "grey30", method = "lm") +
  scale_color_manual(values = c("seagreen2", "chocolate", "blue")) +
  scale_fill_manual(values = c("seagreen2", "chocolate", "blue")) +
  facet_wrap(~traitmet, scales = "free_x")


# assess functional dissimilarity (does brca grow when it's different from it's neighbors?)
gather(brca_diss, traitmet, val, grep("^F[A-Z][a-z]|^CWM|^PC", names(brca_hier))) %>%
  #subset(grepl("Herb", herbicide)) %>%
  #subset(BRCA > 0) %>%
  ggplot(aes(val, BRCA_ct)) +
  #ggplot(aes(val, log(BRCA+0.0001))) +
  #ggplot(aes(val, scale(BRCA))) +
  geom_vline(aes(xintercept = 0), col = "black", lwd = 1.5) +
  geom_point(aes(col = ppt_trt), size = 3, alpha = .75) + #col = 
  #geom_point(aes(col = ppt_trt, shape = nut_trt), size = 3, alpha = .75) + #col = 
  #geom_smooth(aes(col = ppt_trt, fill = ppt_trt), method = "lm") +
  stat_smooth(aes(col = ppt_trt, fill = ppt_trt), method = "lm", formula = "y ~ poly(x,2)") +
  stat_smooth(aes(),col = "grey30", method = "lm", formula = "y ~ poly(x,2)") + #group = herbicide, lty = herbicide
  #stat_smooth(aes(),col = "grey30", method = "lm") +
  scale_color_manual(values = c("seagreen2", "chocolate", "blue")) +
  scale_fill_manual(values = c("seagreen2", "chocolate", "blue")) +
  facet_wrap(~traitmet, scales = "free_x")


# (2.3) ESCA trait hierarchy -----

esca_hier <- funchier(sp = "ESCA")
esca_diss <- funchier(sp = "ESCA", diss = T)
# stick with count to be consistent with gllvm community modeling
esca_hier <- left_join(esca_hier, rename(widesp_count[c("subplotID", "seedtrt", "herbicide", "ESCA")], ESCA_ct = ESCA))
esca_diss <- left_join(esca_diss, rename(widesp_count[c("subplotID", "seedtrt", "herbicide", "ESCA")], ESCA_ct = ESCA))


# test trait hierarchy on root volume, coarse root diameter, cwmpc
esca_hier_pcmod <- glmmTMB(ESCA_ct ~ PC.CWM + (1|block/xpos/herbicide), family = nbinom2(), data = esca_hier) # can't run w ypos
esca_hier_pcmod_full <- glmmTMB(ESCA_ct ~ PC.CWM + herbicide + ppt_trt + nut_trt + (1|block), family = nbinom2(), data = esca_hier)
anova(esca_hier_pcmod_full, esca_hier_pcmod)
car::Anova(esca_hier_pcmod_full, type = 2)
esca_hier_pcmod_red <- glmmTMB(ESCA_ct ~ PC.CWM + ppt_trt + (1|block/nut_trt/herbicide), family = nbinom2(), data = esca_hier)
summary(esca_hier_pcmod_red) # as neighbors are more positive along pc axis (more aquisitive), esca declines
# fdispersion based on pc score difference
esca_hier_pcfdismod <- glmmTMB(ESCA_ct ~ poly(PC.FDis,2) + (1|block/ypos/herbicide), family = nbinom2(), data = esca_hier)
summary(esca_hier_pcfdismod) # esca has lowest abundances in intermediate PC dispersion convex hull ranges, better when range smaller or larger on trait hierarchy
esca_diss_pcfdismod <- glmmTMB(ESCA_ct ~ poly(PC.FDis,2) + (1|block/ypos/herbicide), family = nbinom2(), data = esca_diss)
summary(esca_diss_pcfdismod) # not as strong with absolute dissimilarity. trait hierarchy matters

esca_diss_pcmod <- glmmTMB(ESCA_ct ~ PC.CWM + (1|block/xpos/herbicide), family = nbinom2(), data = esca_diss) # can't run w ypos
esca_diss_pcmod_full <- glmmTMB(ESCA_ct ~ PC.CWM + herbicide + ppt_trt + nut_trt + (1|block), family = nbinom2(), data = esca_diss)
anova(esca_diss_pcmod_full, esca_diss_pcmod) # doesn't matter
car::Anova(esca_diss_pcmod_full, type = 2) # hierarchy captures more of the variance than dissimilarity
esca_diss_pcmod_red <- glmmTMB(ESCA_ct ~ PC.CWM + ppt_trt + (1|block/nut_trt/herbicide), family = nbinom2(), data = esca_diss)
summary(esca_diss_pcmod_red) # as neighbors are increasingly dissimilar from esca, esca declines (but esca is one of the more conservative species)

# Func dissimilarity relative to ESCA
esca_diss_fdismod <- glmmTMB(ESCA_ct ~ FDis * ppt_trt + (1|block/xpos), family = nbinom2(), data = esca_diss) # can't run w herbicide
summary(esca_diss_fdismod) # the more dissimilar from ESCA, the more ESCA declines?
esca_hier_fdismod <- glmmTMB(ESCA_ct ~ FDis * ppt_trt + (1|block/herbicide), family = nbinom2(), data = esca_hier) # can't run w xpos
summary(esca_hier_fdismod) 
car::Anova(esca_hier_fdismod, type= 3)
car::Anova(esca_diss_fdismod, type= 3) # hierarchy captures more of the variance than dissimilarity

# proceed with hierarchy
esca_hier_fdismod_full <- glmmTMB(ESCA_ct ~ FDis + ppt_trt + (1|block/xpos/herbicide), family = nbinom2(), data = esca_hier) # can't run w ypos
summary(esca_hier_fdismod_full)
car::Anova(esca_hier_fdismod_full, type = 2)
esca_hier_fdismod_full <- glmmTMB(ESCA_ct ~ FDis + herbicide + (1|block/xpos), family = nbinom2(), data = esca_hier) # can't run w ypos
summary(esca_hier_fdismod_full)
car::Anova(esca_hier_fdismod_full, type = 2)
anova(esca_diss_fdismod_full, esca_diss_fdismod) # simpler model is better

# > trait hierarchy is stronger for esca
# look at root and aboveground traits with trait hierarchy
esca_hier_indivtraits <- glmmTMB(ESCA_ct ~ CWM.RMF + CWM.SLA.cm2.g + CWM.Height.cm + CWM.Coarse.root.diameter.mm +
                                   CWM.Fine.root.specific.length.cm.g + 
                                   herbicide + nut_trt + ppt_trt + (1|block), family = nbinom2(), data = esca_hier) 
summary(esca_hier_indivtraits)
car::Anova(esca_hier_indivtraits, type = 2)
esca_hier_indivtraits_red <- glmmTMB(ESCA_ct ~ (CWM.Height.cm + CWM.Fine.root.specific.length.cm.g ) + ppt_trt + (1|block/nut_trt/herbicide), 
                                 family = nbinom2(), data = esca_hier) 
summary(esca_hier_indivtraits_red)
# simple root model
esca_hier_indivtraits_simple <- glmmTMB(ESCA_ct ~ CWM.Fine.root.specific.length.cm.g + ppt_trt + (1|block/xpos/herbicide), 
                                     family = nbinom2(), data = esca_hier) 
summary(esca_hier_indivtraits_simple)
car::Anova(esca_hier_indivtraits_simple, type = 2)
esca_diss_indivtraits_simple <- glmmTMB(ESCA_ct ~ CWM.Fine.root.specific.length.cm.g + ppt_trt + (1|block/xpos/herbicide), # can't run with xpos
                                        family = nbinom2(), data = esca_diss) 
summary(esca_diss_indivtraits_simple)
car::Anova(esca_diss_indivtraits_simple, type = 2) # hierarchy captures more variance

# plot all CWMs, ignore treatments
# trait hierachies
gather(esca_hier, traitmet, val, grep("^F[A-Z][a-z]|^CWM|^PC", names(esca_hier))) %>%
  subset(grepl("Herb", herbicide)) %>%
  #subset(ESCA > 0) %>%
  ggplot(aes(val, ESCA_ct)) +
  #ggplot(aes(val, log(ESCA+0.0001))) +
  #ggplot(aes(val, scale(ESCA))) +
  geom_vline(aes(xintercept = 0), col = "black", lwd = 1.5) +
  geom_point(aes(col = ppt_trt, shape = herbicide), size = 3, alpha = .75) + #col = 
  #geom_smooth(aes(col = ppt_trt, fill = ppt_trt, lty = herbicide), se = F, method = "lm") +
  stat_smooth(aes(col = ppt_trt, fill = ppt_trt), method = "lm", formula = "y ~ poly(x,2)") +
  #stat_smooth(aes(),col = "grey30", method = "lm", formula = "y ~ poly(x,2)") + #group = herbicide, lty = herbicide
  #stat_smooth(aes(),col = "grey30", method = "lm") +
  scale_color_manual(values = c("seagreen2", "chocolate", "blue")) +
  scale_fill_manual(values = c("seagreen2", "chocolate", "blue")) +
  scale_shape_manual(values = c("Non-herbicided" = 1, "Herbicided" = 19)) +
  facet_wrap(~traitmet, scales = "free_x")


# just signif things
subset(esca_hier, select = grep("rowid|PC.CW|FDis|CWM.Fine|block|herbici|ppt_|nut_tr|ESCA", names(esca_hier))) %>%
  gather(traitmet, val, grep("^F[A-Z][a-z]|^CWM|^PC", names(.))) %>%
  #subset(ESCA > 0) %>%
  ggplot(aes(val, ESCA_ct)) +
  #ggplot(aes(val, log(ESCA+0.0001))) +
  #ggplot(aes(val, scale(ESCA))) +
  geom_vline(aes(xintercept = 0), col = "black", lwd = 1.5) +
  geom_point(aes(col = ppt_trt, shape = herbicide), size = 3, alpha = .75) + #col = 
  #geom_smooth(aes(col = ppt_trt, fill = ppt_trt, lty = herbicide), se = F, method = "lm") +
  stat_smooth(aes(col = ppt_trt, fill = ppt_trt), method = "lm", formula = "y ~ poly(x,2)") +
  #stat_smooth(aes(),col = "grey30", method = "lm", formula = "y ~ poly(x,2)") + #group = herbicide, lty = herbicide
  #stat_smooth(aes(),col = "grey30", method = "lm") +
  scale_color_manual(values = c("seagreen2", "chocolate", "blue")) +
  scale_fill_manual(values = c("seagreen2", "chocolate", "blue")) +
  scale_shape_manual(values = c("Non-herbicided" = 1, "Herbicided" = 19)) +
  ggtitle("ESCA: trait hierarchies that drive abundance") +
  facet_wrap(~traitmet, scales = "free_x")

subset(esca_diss, select = grep("rowid|PC.CW|FDis|CWM.Fine|block|herbici|ppt_|nut_tr|ESCA", names(esca_hier))) %>%
  gather(traitmet, val, grep("^F[A-Z][a-z]|^CWM|^PC", names(.))) %>%
  #subset(ESCA > 0) %>%
  ggplot(aes(val, log(ESCA_ct+.01))) +
  #ggplot(aes(val, log(ESCA+0.0001))) +
  #ggplot(aes(val, scale(ESCA))) +
  geom_vline(aes(xintercept = 0), col = "black", lwd = 1.5) +
  geom_point(aes(col = ppt_trt, shape = herbicide), size = 3, alpha = .75) + #col = 
  geom_smooth(col = "black", se = T, method = "lm") + #lty = herbicide
  geom_smooth(aes(col = ppt_trt, fill = ppt_trt), se = T, method = "lm") + #lty = herbicide
  #stat_smooth(aes(col = ppt_trt, fill = ppt_trt), method = "lm", formula = "y ~ poly(x,2)") +
  #stat_smooth(aes(),col = "grey30", method = "lm", formula = "y ~ poly(x,2)") + #group = herbicide, lty = herbicide
  #stat_smooth(aes(),col = "grey30", method = "lm") +
  scale_color_manual(values = c("seagreen2", "chocolate", "blue")) +
  scale_fill_manual(values = c("seagreen2", "chocolate", "blue")) +
  scale_shape_manual(values = c("Non-herbicided" = 1, "Herbicided" = 19)) +
  ggtitle("ESCA: trait dissimilarities that drive abundance") +
  facet_wrap(~traitmet, scales = "free_x")


# assess functional dissimilarity (does esca grow when it's different from it's neighbors?)
gather(esca_diss, traitmet, val, grep("^F[A-Z][a-z]|^CWM|^PC", names(esca_hier))) %>%
  subset(grepl("Herb", herbicide)) %>%
  #subset(ESCA > 0) %>%
  ggplot(aes(val, ESCA)) +
  #ggplot(aes(val, scale(ESCA))) +
  geom_vline(aes(xintercept = 0), col = "black", lwd = 1.5) +
  geom_point(aes(col = ppt_trt, shape = nut_trt), size = 3, alpha = .75) + #col = 
  #geom_smooth(aes(col = ppt_trt, fill = ppt_trt), method = "lm") +
  #geom_smooth(aes(col = ppt_trt, fill = ppt_trt), method = "lm", formula = "y ~ poly(x,2)") +
  stat_smooth(aes(),col = "grey30", method = "lm", formula = "y ~ poly(x,2)") + #group = herbicide, lty = herbicide
  #stat_smooth(aes(),col = "grey30", method = "lm") +
  scale_color_manual(values = c("seagreen2", "chocolate", "blue")) +
  scale_fill_manual(values = c("seagreen2", "chocolate", "blue")) +
  facet_wrap(~traitmet, scales = "free_x")



# (2.3) NEMA trait hierarchy -----
nema_hier <- funchier(sp = "NEMA")
nema_diss <- funchier(sp = "NEMA", diss = T)
# stick with count to be consistent with gllvm community modeling
nema_hier <- left_join(nema_hier, rename(widesp_count[c("subplotID", "seedtrt", "herbicide", "NEMA")], NEMA_ct = NEMA))
nema_diss <- left_join(nema_diss, rename(widesp_count[c("subplotID", "seedtrt", "herbicide", "NEMA")], NEMA_ct = NEMA))

# -- pc1 differences ----
nema_hier_pcmod <- glmmTMB(NEMA_ct ~ PC.CWM + (1|block/xpos/ypos), family = nbinom2(), data = nema_hier) # can't run w ypos
nema_hier_pcmod_full <- glmmTMB(NEMA_ct ~ PC.CWM + herbicide + ppt_trt + nut_trt + (1|block), family = nbinom2(), data = nema_hier)
anova(nema_hier_pcmod_full, nema_hier_pcmod) # full is much better?? 
summary(nema_hier_pcmod_full)
car::Anova(nema_hier_pcmod_full, type = 2) # ppt trt is the most important
nema_hier_pcmod_red <- glmmTMB(NEMA_ct ~ PC.CWM * ppt_trt + nut_trt + herbicide+ (1|block), family = nbinom2(), data = nema_hier)
summary(nema_hier_pcmod_red)
car::Anova(nema_hier_pcmod_red, type = 2)
# compare with dissimilarity
nema_diss_pcmod_red <- glmmTMB(NEMA_ct ~ PC.CWM * ppt_trt + nut_trt + herbicide+ (1|block), family = nbinom2(), data = nema_diss)
summary(nema_diss_pcmod_red) # trait hierarchy more important for PC CWM

# -- FDis ----
nema_hier_fdismod <- glmmTMB(NEMA_ct ~ FDis + (1|block/ypos), family = nbinom2(), data = nema_hier) # can't run w ypos
nema_hier_fdismodx_full <- glmmTMB(NEMA_ct ~ FDis * ppt_trt + nut_trt + (1|block/herbicide), family = nbinom2(), data = nema_hier)
summary(nema_hier_fdismodx_full)
car::Anova(nema_hier_fdismodx_full) # FDis interaction marginally signif
nema_hier_fdismod_full <- glmmTMB(NEMA_ct ~ FDis + ppt_trt + nut_trt + (1|block/herbicide), family = nbinom2(), data = nema_hier)
summary(nema_hier_fdismod_full)
anova(nema_hier_fdismodx_full, nema_hier_fdismod_full) # interaction marginally better

# look at dissimilarity rather than hierarchy
nema_diss_fdismodx_full <- glmmTMB(NEMA_ct ~ FDis * ppt_trt + nut_trt + (1|block/herbicide), family = nbinom2(), data = nema_diss)
summary(nema_diss_fdismodx_full) # no
nema_diss_fdismod_full <- glmmTMB(NEMA_ct ~ FDis + ppt_trt + nut_trt + (1|block/herbicide), family = nbinom2(), data = nema_diss)
summary(nema_diss_fdismod_full)
car::Anova(nema_diss_fdismod_full) # absolute dissimilarity not meaningful for dissimilarity; trait hierarchy more meaningful

# -- FDiv ----
nema_hier_fdivmod_full <- glmmTMB(NEMA_ct ~ FDiv + ppt_trt + nut_trt + (1|block), family = nbinom2(), data = nema_hier) # can't run with interaction
summary(nema_hier_fdivmod_full) # not meaningful for trait hierarchies

nema_diss_fdivmodx_full <- glmmTMB(NEMA_ct ~ FDiv * ppt_trt + nut_trt + (1|block), family = nbinom2(), data = nema_diss)
summary(nema_diss_fdivmodx_full)
car::Anova(nema_diss_fdivmodx_full) # FDiv not meaninful as main effect, but moderates ppt effect
# trait dissimilarity stronger for functional divergence (more neighbors are on the boundary edge of convex) 


# -- individual traits ----
nema_hier_individual_full <- glmmTMB(NEMA_ct ~ CWM.SLA.cm2.g + CWM.Height.cm + CWM.Coarse.root.diameter.mm + CWM.Root.volume.cm3 +
                                       CWM.Fine.root.specific.length.cm.g + ppt_trt + nut_trt + herbicide + (1|block), family = nbinom2(), 
                                     data = nema_hier)
summary(nema_hier_individual_full)
car::Anova(nema_hier_individual_full, type = 2)
nema_hier_individual_red <- glmmTMB(NEMA_ct ~ CWM.SLA.cm2.g + CWM.Height.cm + CWM.Coarse.root.diameter.mm + CWM.Root.volume.cm3 +
                                       CWM.Fine.root.specific.length.cm.g + ppt_trt + nut_trt + herbicide + (1|block), family = nbinom2(), 
                                     data = nema_hier)
summary(nema_hier_individual_red)
nema_hier_individual_red <- glmmTMB(NEMA_ct ~ (CWM.Coarse.root.diameter.mm + CWM.Root.volume.cm3 +
                                      CWM.Fine.root.specific.length.cm.g) * ppt_trt + nut_trt + herbicide + (1|block), family = nbinom2(), 
                                    data = nema_hier)
summary(nema_hier_individual_red)
car::Anova(nema_hier_individual_red, type = 2)
nema_hier_individual_red <- glmmTMB(NEMA_ct ~ (CWM.Coarse.root.diameter.mm + CWM.Fine.root.specific.length.cm.g) * ppt_trt + nut_trt + herbicide + (1|block), family = nbinom2(), 
                                    data = nema_hier)
summary(nema_hier_individual_red)
car::Anova(nema_hier_individual_red, type = 2)

nema_hier_individualx_red <- glmmTMB(NEMA_ct ~ CWM.Fine.root.specific.length.cm.g * ppt_trt + nut_trt + herbicide + (1|block), family = nbinom2(), 
                                    data = nema_hier)
summary(nema_hier_individualx_red)
car::Anova(nema_hier_individualx_red, type = 2)
nema_hier_individual_red <- glmmTMB(NEMA_ct ~ CWM.Fine.root.specific.length.cm.g + ppt_trt + nut_trt + herbicide + (1|block), family = nbinom2(), 
                                    data = nema_hier)
summary(nema_hier_individual_red)
car::Anova(nema_hier_individual_red, type = 2)
anova(nema_hier_individualx_red, nema_hier_individual_red) # with interaction is better

# compare with trait dissimilarity
nema_diss_individualx_red <- glmmTMB(NEMA_ct ~ CWM.Fine.root.specific.length.cm.g * ppt_trt + nut_trt + herbicide + (1|block), family = nbinom2(), 
                                     data = nema_diss)
summary(nema_diss_individualx_red)
nema_diss_individual_red <- glmmTMB(NEMA_ct ~ CWM.Fine.root.specific.length.cm.g + ppt_trt + nut_trt + herbicide + (1|block), family = nbinom2(), 
                                     data = nema_diss)
summary(nema_diss_individual_red) # fine root is only meaningful using hierarchy
nema_diss_ldmc_red <- glmmTMB(NEMA_ct ~ poly(CWM.LDMC,2) * ppt_trt + nut_trt + herbicide + (1|block), family = nbinom2(), 
                                    data = nema_diss)
summary(nema_diss_ldmc_red)
car::Anova(nema_diss_ldmc_red, type = 2)
nema_diss_ldmc_red <- glmmTMB(NEMA_ct ~ poly(CWM.LDMC,2) * ppt_trt + nut_trt + (1|block/herbicide), family = nbinom2(), 
                              data = nema_diss)
summary(nema_diss_ldmc_red)
car::Anova(nema_diss_ldmc_red, type = 2) # this is the only individual trait model meaningful for dissimilarity
nema_diss_ldmc_redlme <- glmer.nb(NEMA_ct ~ poly(CWM.LDMC,2) * ppt_trt + nut_trt + (1|block/herbicide), data = nema_diss)
summary(nema_diss_ldmc_redlme)
nema_diss_ldmc_simple <- glmmTMB(NEMA_ct ~ poly(CWM.LDMC,2) + ppt_trt + nut_trt + (1|block), family = nbinom2(), 
                              data = nema_diss)
summary(nema_diss_ldmc_simple) # LDMC has a nonlinear relationship (not signif if just linear), specifically in wet plots
anova(nema_diss_ldmc_red, nema_diss_ldmc_simple) # reduced (intereacting with ppt trt is better)

emmeans(nema_hier_fdismod_full, ~ FDis | ppt_trt, component = "response")
plot(predictorEffects(nema_hier_pcmod_red, partial.residuals=TRUE))
plot(Effect("CWM.LDMC", nema_diss_ldmc_redlme, partial.residuals=TRUE))
plot(Effect("ppt_trt", nema_diss_ldmc_redlme, partial.residuals=TRUE,
     partial.residual=list(pch="+", col="#FF00FF80")))
standardize_parameters(nema_hier_fdismod_full, two_sd = T, include_response = T)
standardize_parameters(nema_hier_fdismod_full, two_sd = T, method = "basic")
effectsize(nema_diss_ldmc_simple)
standardize_parameters(nema_diss_ldmc_simple, two_sd = T, method = "basic")
pairs(emmeans(nema_diss_ldmc_red, ~ ppt_trt | CWM.LDMC))
effectsize(nema_diss_ldmc_red)
# ? i am not seeing any difference in marginal means from what you can eyeball? standardized effect sizes look the same
insight::standardize_names(parameters::model_parameters(nema_diss_ldmc_red, exponentiate = T))
# plot all CWMs, ignore treatments
# trait hierachies
gather(nema_hier, traitmet, val, grep("^F[A-Z][a-z]|^CWM|^PC", names(nema_hier))) %>%
  #subset(grepl("Herb", herbicide)) %>%
  #subset(NEMA > 0) %>%
  ggplot(aes(val, NEMA_ct)) +
  #ggplot(aes(val, log(NEMA+0.0001))) +
  #ggplot(aes(val, scale(NEMA))) +
  geom_vline(aes(xintercept = 0), col = "black", lwd = 1.5) +
  geom_point(aes(col = ppt_trt), size = 3, alpha = .75) + #col =  shape = herbicide
  geom_smooth(aes(col = ppt_trt, fill = ppt_trt), se = F, method = "lm") +
  #geom_smooth(aes(col = ppt_trt, fill = ppt_trt), method = "lm", formula = "y ~ poly(x,2)") +
  #stat_smooth(aes(),col = "grey30", method = "lm", formula = "y ~ poly(x,2)") + #group = herbicide, lty = herbicide
  stat_smooth(aes(),col = "grey30", method = "lm") +
  scale_color_manual(values = c("seagreen2", "chocolate", "blue")) +
  scale_fill_manual(values = c("seagreen2", "chocolate", "blue")) +
  scale_shape_manual(values = c(1,19)) +
  ggtitle("NEMA: trait hierarhcy (focal - neighbor trait)") +
  facet_wrap(~traitmet, scales = "free_x")


# assess functional dissimilarity (does nema grow when it's different from it's neighbors?)
gather(nema_diss, traitmet, val, grep("^F[A-Z][a-z]|^CWM|^PC", names(nema_hier))) %>%
  #subset(grepl("Herb", herbicide)) %>%
  #subset(NEMA > 0) %>%
  ggplot(aes(val, NEMA_ct)) +
  #ggplot(aes(val, scale(NEMA))) +
  #geom_vline(aes(xintercept = 0), col = "black", lwd = 1.5) +
  geom_point(aes(col = ppt_trt), size = 3, alpha = .75) + #col = 
  #geom_smooth(aes(col = ppt_trt, fill = ppt_trt), method = "lm") +
  geom_smooth(aes(col = ppt_trt, fill = ppt_trt), method = "lm", formula = "y ~ poly(x,2)") +
  #stat_smooth(aes(),col = "grey30", method = "lm", formula = "y ~ poly(x,2)") + #group = herbicide, lty = herbicide
  #stat_smooth(aes(),col = "grey30", method = "lm") +
  scale_color_manual(values = c("seagreen2", "chocolate", "blue")) +
  scale_fill_manual(values = c("seagreen2", "chocolate", "blue")) +
  ggtitle("NEMA: trait dissimilarity (absolute trait difference)") +
  facet_wrap(~traitmet, scales = "free_x")



# (3) MODEL TRAIT HIERARCHY -----
VarCorr(nema_hier_fdismod_full)

# adjusted pvalues for multiple hypothesis testing

# determine marginal importance


# (4) CHAPTER MODELS -----
# just run mods for trait dissimilarity, separate mods for each spp, and per KNS one w  only trait, then trait 

femi_diss_nb <- glmmTMB(FEMI_ct ~ FDis + herbicide + (1|block/xpos/ypos), family = nbinom2(), data = femi_diss)
femi_diss_env_nb <- glmmTMB(FEMI_ct ~ FDis + ppt_trt + nut_trt + herbicide + (1|block), family = nbinom2(), data = femi_diss)
femi_diss_envx_nb <- glmmTMB(FEMI_ct ~ FDis + (ppt_trt + nut_trt) * herbicide + (1|block), family = nbinom2(), data = femi_diss)
summary(femi_diss_nb)
summary(femi_diss_env_nb)
summary(femi_diss_envx_nb)

femi_hier_nb <- glmmTMB(FEMI_ct ~ FDis + herbicide + (1|block/xpos/ypos), family = nbinom2(), data = femi_hier)
femi_hier_env_nb <- glmmTMB(FEMI_ct ~ FDis + ppt_trt + nut_trt + herbicide + (1|block), family = nbinom2(), data = femi_hier)
femi_hier_envx_nb <- glmmTMB(FEMI_ct ~ FDis + (ppt_trt * nut_trt) + herbicide + (1|block), family = nbinom2(), data = femi_hier)
femi_hier_envxonly_nb <- glmmTMB(FEMI_ct ~ (ppt_trt + nut_trt) * herbicide + (1|block), family = nbinom2(), data = femi_hier)
summary(femi_hier_nb)
summary(femi_hier_env_nb)
summary(femi_hier_envx_nb)
summary(femi_hier_envxonly_nb)

anova(femi_hier_envx_nb, femi_hier_envxonly_nb, femi_hier_nb)


# integrated trait pca: disimilarity v trait hierarchy
# BRCA

# FEMI
femi_diss_pccwm.envx_nb <- glmmTMB(FEMI_ct ~ PC.CWM + (ppt_trt + nut_trt) * herbicide + (1|block), family = nbinom2(), data = femi_diss)
femi_hier_pccwm.envx_nb <- glmmTMB(FEMI_ct ~ PC.CWM + (ppt_trt + nut_trt) * herbicide + (1|block), family = nbinom2(), data = femi_hier)

# NEMA

# ESCA
esca_diss_pccwm.envx_nb <- glmmTMB(ESCA_ct ~ PC.CWM + (ppt_trt + nut_trt) * herbicide + (1|block), family = nbinom2(), data = esca_diss)
esca_hier_pccwm.envx_nb <- glmmTMB(ESCA_ct ~ PC.CWM + (ppt_trt + nut_trt) * herbicide + (1|block), family = nbinom2(), data = esca_hier)
summary(esca_diss_pccwm.envx_nb)
summary(esca_hier_pccwm.envx_nb)

nema_diss_nb
esca_diss_nb
brca_diss_nb