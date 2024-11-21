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

# major changes to file:
# feb 2023:
# >> analysis and figs in support of esa 2023 abstract submission
# aug 2023: 
# >> add code for figs in esa 2023 presentation, including multivariate
# jan 2024: 
# >> revise_/cleanup analysis code for chapter prep (ctw has better understanding of gjam now, chain convergence checks to do)
# >> move multivar code back to ordinations script


# -- SETUP -----
# load needed libraries
library(tidyverse)
library(cowplot)
library(gjam)
library(coda)

# get cover data
source("Native-analysis/native_prep_cover.R", print.eval = F) # default options prefs set, theme is test

# specify path for writing out models and figs
outpath <- "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/chapter/gjam/"



# -- DATA PREP ----




# -- VIZ NATIVE RECRUITS ----
# check abundance of each nat per treatment
left_join(seededspp, trtkey) %>%
  gather(code4, cover, nats) %>%
  group_by(herbicide, block > 2) %>%
  ggplot(aes(code4, cover, group = paste(herbicide, block >2), col = block > 2)) +
  geom_jitter(aes(shape = herbicide), width = 0.2,height = 0, alpha = 0.5, size = 2) +
  stat_summary(aes(shape = herbicide), alpha = 0.7) +
  facet_grid(ppt_trt~nut_trt) +
  coord_flip()

# plot by species to see cover change across treatments relative for each spp
left_join(seededspp, trtkey) %>%
  gather(code4, cover, nats) %>%
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
  reframe(max_cov = max(pct_cover)) %>%
  distinct() %>%
  mutate(expgroup = paste(plot, fulltrt, herbicide, seedtrt, sep = "_")) %>%
  subset(seedtrt != "Unseeded")

# plot overall rank abundance
global_ranks <- group_by(sppcover, code4) %>%
  reframe(totcov = sum(max_cov),
            nplots = length(unique(expgroup)),
            prop_plots = nplots/length(unique(sppcover$expgroup))) %>%
  distinct() %>%
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
  reframe(totcov = sum(max_cov),
            nplots = length(unique(expgroup)),
            prop_plots = nplots/allplots) %>% # 105 = total # unique experiment subplots
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
  reframe(totcov = sum(max_cov),
            nplots = length(unique(expgroup)),
            prop_plots = nplots/allplots) %>% # 105 = total # unique experiment subplots
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

# -- MODEL EVALUATION FXNS -----
# eval functions

# function to compile geweke diagnostics as data frame
ctw_gewdf <- function(gew){
  data.frame(z = gew$z) %>%
    rownames_to_column() %>%
    separate(rowname, into = c("spp", "term"), sep = "_") %>%
    # drop other if present
    subset(!grepl("other", spp)) %>%
    return()
}

# function to create custom geweke diagnostic plots
ctw_gewfig <- function(gewdf, cap, termfac, dropvars = "keepall", lpos = c(0.9,0.1)){
  p <- mutate(gewdf, term = factor(term, termfac)) %>%
    subset(!grepl(dropvars, term)) %>%
    ggplot(aes(term, z)) +
    geom_hline(aes(yintercept = -2), lty = 2) +
    geom_hline(aes(yintercept = 2), lty = 2) +
    geom_point(aes(shape = abs(z) > 2, color = abs(z) > 2)) +
    scale_shape_manual(values = c(1,19)) +
    scale_color_manual(values = c("black","red")) +
    facet_wrap(~spp) +
    labs(caption = cap, y = "Geweke diagnostic z-score", x = "Model covariate") +
    coord_flip() +
    theme(legend.position = lpos)
  # return plot
  return(p)
}
# function to compile gelman diagnostics as data frame
ctw_geldf <- function(gel){
  data.frame(gel$psrf) %>%
    rownames_to_column() %>%
    separate(rowname, into = c("spp", "term"), sep = "_") %>%
    # drop other if present
    subset(!grepl("other", spp)) %>%
    return()
}
# function to plot gelman diagnostics
ctw_gelfig <- function(geldf, cap, termfac = rev(varorder), dropvars = "keepall"){
  mutate(geldf, term = factor(term, levels = termfac)) %>%
    subset(!grepl(dropvars, term)) %>%
    ggplot() +
    geom_pointrange(aes(x = term, y  = Point.est., ymin = Point.est., ymax = `Upper.C.I.`, shape = Upper.C.I. >1.2, col = Upper.C.I. >1.2)) +
    scale_shape_manual(values = c(1,19)) +
    scale_color_manual(values = c("black","red")) +
    facet_wrap(~spp) +
    coord_flip() +
    labs(caption = cap) +
    theme(legend.position = c(0.9, 0.1))
}
# function to extract an mcmc chain
get_mcmc <- function(gmod, chain = "bgibbs"){
  tmp <- mcmc(gmod$chains[[chain]][(gmod$modelList$burnin + 1):gmod$modelList$ng,])
  return(tmp)
}


# compile tidy chains by all spp by for variable of interest
chaindf <- function(mlist, v = "herbicide"){ 
  # extract colnames
  mcnames <- varnames(mlist[1])
  # pull variable and spp names
  tmplist <- mlist[,grep(v, mcnames)]
  tmplist <- lapply(tmplist, data.frame)
  tmplist <- lapply(tmplist, function(x) gather(x,"spp", "val", 1:ncol(x)))
  tmpdf <- data.frame()
  for(i in 1:length(tmplist)){
    tmpdf <- rbind(tmpdf, cbind(chain = i, tmplist[[i]]))
  }
  tmpdf$term <- gsub("^[A-Z]+_", "", tmpdf$spp)
  tmpdf$spp <- gsub("_.*", "", tmpdf$spp)
  tmpdf$chain <- as.character(tmpdf$chain)
  return(tmpdf)
}
ctw_densplot <- function(df, subt = NULL){
  p <- ggplot(df) +
    geom_density(aes(val, fill = chain), alpha = 0.5) +
    geom_vline(xintercept = 0, lwd = 0.75) +
    labs(x = "Beta estimate", subtitle = subt) +
    theme(legend.position = c(0.97, 0.01),
          legend.direction = "horizontal",
          legend.justification = c("right","bottom")) +
    facet_wrap(~spp, scales = "free_x")
  return(p)
} 



# -- 1. Coarse fxnl group response to treatments ----
# > note seeded natives are split from background natives here

# 4 types of cover groups per plot
natcoarse_ydata <- subset(natcoarse_all, select = c(plot, seedtrt, herbicide, coarse_fxnl, totcov)) %>%
  # add ambient cover 2021, unseeded
  rbind(subset(covcoarse, select = c(plot, seedtrt, herbicide, coarse_fxnl, totcov))) %>%
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

# only seeded spp will have seeded prefix, else assume ambient community
names(natcoarse_ydata) <- gsub("Background", "", names(natcoarse_ydata))
names(natcoarse_ydata) <- gsub("Exotic", "Intro", names(natcoarse_ydata))
names(natcoarse_ydata) <- gsub("Native", "Nat", names(natcoarse_ydata))

natcoarse_xdata <- subset(natcoarse_ydata, select = c(rowid, plot, seedtrt, herbicide)) %>%
  left_join(trtkey) %>%
  # convert env vars to factors
  mutate(nut_trt = factor(nut_trt, levels = c("N", "F", "C"), labels = c("0noNut", "1fert", "2comp")),
         ppt_trt = factor(ppt_trt, levels = c("XC", "D", "W"), labels = c("0xc", "1d", "2w")),
         #seedtrt = factor(seedtrt, levels = c("Unseeded", "Native seeded"), labels = c("0unseed", "1seed")),
         #herbicide = gsub("Non-herb", "NonHerb", herbicide),
         #herbicide = factor(herbicide,ordered = T,  levels = c("NonHerbicided", "Herbicided"), labels = c("0noHerb", "1herb")),
         # make seedtrt and herbicide dummy instead
         seedtrt = as.numeric(grepl("Nat", seedtrt)),
         herbicide = as.numeric(grepl("^Herb", herbicide)),
         # create hill position -- upslope is grassier than bottom two blocks
         uphill = ifelse(block <3, 0, 1), 
         # make block and plot factors
         block = factor(block),
         plot = factor(plot),
         # note area of each plot -- ambient 1x1m, otherwise 0.5x0.5m
         area = ifelse(herbicide == "0noHerb" & seedtrt == "0unseed", 1, 0.5*0.5)
         )
rownames(natcoarse_xdata) <- natcoarse_xdata$rowid
# clean up colnames
colnames(natcoarse_xdata) <- gsub("_trt", "Trt", colnames(natcoarse_xdata))

# assign rownames and check colsums on fxnl grps
rownames(natcoarse_ydata) <- natcoarse_ydata$rowid 
natcoarse_ydata <- subset(natcoarse_ydata, select = -c(plot, seedtrt, herbicide, rowid))  
sumcheck <- sapply(natcoarse_ydata, function(x) sum(x > 0))
sort(sumcheck)
# move anything not in at least 10% of plots to other
natcoarse_ydata <- data.frame(gjamTrimY(natcoarse_ydata, minObs = round(nrow(natcoarse_ydata)*0.05))$y)
# check ranges on each 
summary(natcoarse_ydata)
# relativize cover data to make it fractional composition
natcoarse_ydata_rel <- vegan::decostand(natcoarse_ydata, method = "total")
# make pa dataset
natcoarse_ydata_pa <- mutate_all(natcoarse_ydata, function(x) ifelse(x > 0, 1, 0))
# plot continuous cover and relativized side by side
par(mfrow = c(1,2))
boxplot(natcoarse_ydata, main = "CA data")
boxplot(natcoarse_ydata_rel, main = "FC data")
par(mfrow = c(1,1))
# check frequencies
sort(colSums(natcoarse_ydata_pa)/nrow(natcoarse_ydata_pa))

# run continuous cover model -- run 2 chains
allcoarse_gjam_r3 <- gjam(formula = as.formula(~ seedtrt + herbicide + pptTrt * nutTrt), 
                       xdata = natcoarse_xdata,
                       ydata = natcoarse_ydata,
                       modelList = list(ng=20000, burnin=10000, typeNames = 'CA', random = "block"))
summary(allcoarse_gjam_r1)
# write out plots
gjamPlot(allcoarse_gjam_r1, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, SAVEPLOTS = T, 
                                         outFolder = paste0(outpath,"figs/allcoarse_blockrand_20kx10kr1_resummed/")))
# save models to avoid re-run
saveRDS(allcoarse_gjam_r3, paste0(outpath, "mods/allcoarse_gjam_ca20x10k_r3.rdata"))

# test if chain converged at 15k x 5k runs
allcoarse_gjam_r1$fit[1:3] # 20kx10k: RMSPE is 15.75752, and the DIC is 23262; 15kx5k: RMSPE 15.77179 DIC 24324.14
allcoarse_gjam_r2$fit[1:3] # 20kx10k: RMSPE is 15.76134, and the DIC is 22053; 15kx5k: RMSPE 15.77781 DIC 80289.68
allcoarse_gjam_r3$fit[1:3] # 20kx10k: RMSPE is 15.8, and the DIC is 23747
allcoarse_gjam_nogram_r1$fit[1:3] # RMSPE is 16.9, and the DIC is 20322
allcoarse_gjam_nogram_r2$fit[1:3] # RMSPE is 16.9, and the DIC is 20276
allcoarse_gjam_nogram_r3$fit[1:3] # RMSPE is 16.9, and the DIC is 20798

# pull standardized bgibbs chains
allcoarse1_bgibbs <- get_mcmc(allcoarse_gjam_r1)
allcoarse1_nogram_bgibbs <- get_mcmc(allcoarse_gjam_nogram_r1)
allcoarse2_bgibbs <- get_mcmc(allcoarse_gjam_r2)
allcoarse2_nogram_bgibbs <- get_mcmc(allcoarse_gjam_nogram_r2)
allcoarse3_bgibbs <- get_mcmc(allcoarse_gjam_r3)
allcoarse3_nogram_bgibbs <- get_mcmc(allcoarse_gjam_nogram_r3)
# geweke diagnostics 
allcoarse1_gewdf <- ctw_gewdf(geweke.diag(allcoarse1_bgibbs))
allcoarse1_nogram_gewdf <- ctw_gewdf(geweke.diag(allcoarse1_nogram_bgibbs))
allcoarse2_gewdf <- ctw_gewdf(geweke.diag(allcoarse2_bgibbs))
allcoarse2_nogram_gewdf <- ctw_gewdf(geweke.diag(allcoarse2_nogram_bgibbs))
allcoarse3_gewdf <- ctw_gewdf(geweke.diag(allcoarse3_bgibbs))
# specify order for variables -- as compiled
varorder <- allcoarse1_gewdf$term[allcoarse1_gewdf$spp == "NatGrass"]
allcoarse1_gewfig <- ctw_gewfig(allcoarse1_gewdf, termfac = rev(varorder), lpos = "right", dropvars = "intercept",
                                cap = "all plots coarse functional, 20kx10k, chain 1")
allcoarse1_nogram_gewfig <- ctw_gewfig(allcoarse1_nogram_gewdf, termfac = rev(varorder), lpos = "right", cap = "all plots coarse functional, no nat gram, 20kx10k, chain 1")
allcoarse2_gewfig <- ctw_gewfig(allcoarse2_gewdf, termfac = rev(varorder), lpos = "right", cap = "all plots coarse functional, 20kx10k, chain 2")
allcoarse2_nogram_gewfig <- ctw_gewfig(allcoarse1_nogram_gewdf, termfac = rev(varorder), lpos = "right", cap = "all plots coarse functional, no nat gram, 20x10k, chain 2")
allcoarse3_gewfig <- ctw_gewfig(allcoarse3_gewdf, termfac = rev(varorder), lpos = "right", cap = "all plots coarse functional, 20kx10k, chain 3")
# gelman diagnostic
allcoarse_mlist <- mcmc.list(allcoarse1_bgibbs, allcoarse2_bgibbs, allcoarse3_bgibbs)
allcoarse_geldf <- ctw_geldf(gelman.diag(allcoarse_mlist, multivariate = F)) # throws error on cholesky decomp [suggest params correlated.. which, yes bc factorial exp], stackoverflow says set multivar estimate = F
allcoarse_gelfig <- ctw_gelfig(allcoarse_geldf, termfac = varorder,dropvars = "intercept", cap = "all plots coarse functional, 20kx10k")

allcoarse_nogram_mlist <- mcmc.list(allcoarse1_nogram_bgibbs, allcoarse2_nogram_bgibbs, allcoarse3_nogram_bgibbs)
allcoarse_nogram_geldf <- ctw_geldf(gelman.diag(allcoarse_nogram_mlist, multivariate = F)) # throws error on cholesky decomp [suggest params correlated.. which, yes bc factorial exp], stackoverflow says set multivar estimate = F
allcoarse_nogram_gelfig <- ctw_gelfig(allcoarse_nogram_geldf, termfac = varorder, cap = "all plots coarse functional, no nat gram, 20kx10k")



# run as FC
allcoarse_gjam_fc <- gjam(formula = as.formula(~ seedtrt + herbicide + pptTrt * nutTrt), 
                                 xdata = natcoarse_xdata,
                                 ydata = natcoarse_ydata_rel,
                                 modelList = list(ng=20000, burnin=10000, typeNames = 'FC', random = "block"))
summary(allcoarse_gjam_fc) # The RMSPE is 0.128, and the DIC is 161423
gjamPlot(allcoarse_gjam_fc, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, SAVEPLOTS = T, 
                                            outFolder = paste0(outpath, "figs/allcoarse_blockrand_fc_20kx10kr1/")))
# > results from fractional cover are about the same as using continuous cover; stick with continuous

# run with hillpos to see if that changes effects of treatments (upslope is grassier)
# run continuous cover model -- run 2 chains
allcoarse_hillpos_gjam_r1 <- gjam(formula = as.formula(~ seedtrt + herbicide + uphill + pptTrt * nutTrt), 
                          xdata = natcoarse_xdata,
                          ydata = natcoarse_ydata,
                          modelList = list(ng=20000, burnin=10000, typeNames = 'CA', random = "block"))
summary(allcoarse_hillpos_gjam_r1) # about the same prediction error, DIC higher: RMSPE is 15.7, and the DIC is 31837
gjamPlot(allcoarse_hillpos_gjam_r1, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, SAVEPLOTS = T, 
                                            outFolder = paste0(outpath, "figs/allcoarse_hillpos_blockrand_20kx10kr1/")))
# including hill position slightly changes the clustering in correlation, but otherwise does nothing
# > does not help with prediction, beta estimates for other effects do not change. ignore.

# since can't nest error, try block as factor to capture differences in community comp by block, use plotid as random
allcoarse_blockeffect_gjam_r1 <- gjam(formula = as.formula(~ block + seedtrt + herbicide + pptTrt * nutTrt), 
                          xdata = natcoarse_xdata,
                          ydata = natcoarse_ydata,
                          modelList = list(ng=20000, burnin=10000, typeNames = 'CA', random = "plot"))
summary(allcoarse_blockeffect_gjam_r1)
gjamPlot(allcoarse_blockeffect_gjam_r1, plotPars = list(PLOTALLY = T, GRIDPLOTS = T))



# -- 2. Individual native seeded with coarse fxnl group neighbors -----

# keep natcoarse_ydata, but drop seeded cols and cbind specific natives
# also drop other -- few obs that are > 0
natspp_fxnlneighbors_ydata <- natcoarse_ydata[!grepl("Seeded|other", names(natcoarse_ydata))] %>%
  cbind(natcoarse_xdata[c("plot", "seedtrt", "herbicide")]) %>%
  # add individual seeded spp
  left_join(subset(
    mutate(neighborhood_allplots, 
           seedtrt = as.numeric(grepl("Nati", seedtrt)),
           herbicide = as.numeric(grepl("^Herb", herbicide)),
           plot = factor(plot)), 
    select = c(plot, seedtrt, herbicide, BRCA:NEMA))) %>%
  # drop keys used to join
  subset(select = -c(plot, seedtrt, herbicide)) %>%
  # infill missing nat spp in plot 33 with 0
  mutate_at(nats[nats != "TRCI"], function(x) ifelse(is.na(x), 0, x)) 
# specify rownames again
rownames(natspp_fxnlneighbors_ydata) <- rownames(natcoarse_ydata)

# check colsums
sort(sapply(natspp_fxnlneighbors_ydata, function(x) sum(x > 0))) # ok

# run model
natspp_fxnlneighbors_gjam_run3 <- gjam(formula = as.formula(~ seedtrt + herbicide + pptTrt * nutTrt), 
                       xdata = natcoarse_xdata,
                       ydata = natspp_fxnlneighbors_ydata,
                       modelList = list(ng=40000, burnin=20000, typeNames = 'CA', random = "block"))
summary(natspp_fxnlneighbors_gjam_run1)
gjamPlot(natspp_fxnlneighbors_gjam_run1, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, SAVEPLOTS = T, 
                                         outFolder = paste0(outpath, "figs/natspp_fxnlneighbors_40kx20kr1/")))
# save runs
saveRDS(natspp_fxnlneighbors_gjam_run3, paste0(outpath, "mods/natspp_fxnlneighbors_gjam_ca40x20k_run3.rdata"))

# check rmspe, dic, and variance contributions
natspp_fxnlneighbors_gjam_run1$fit[1:3] # 40x20: 13.40883, DIC: 47725.95; 30kx15k: RMSPE is 14.9, and the DIC is 591,875,271; # 20kx10k: RMSPE is 14.9, and the DIC is 711011167
natspp_fxnlneighbors_gjam_run2$fit[1:3] # 40x20: 13.40703, DIC: 28459.63; 30kx15k: RMSPE is 14.9, and the DIC is 279,763,364; 20x10 RMSPE is 14.9 and DIC 513,403,925
natspp_fxnlneighbors_gjam_run3$fit[1:3] # 40k x 20k: 13.40788, DIC: 57440.2

# check chain convergence -- start with rhat diag
natspp_fxnl1_bgibbs <- get_mcmc(natspp_fxnlneighbors_gjam_run1)
natspp_fxnl2_bgibbs <- get_mcmc(natspp_fxnlneighbors_gjam_run2)
natspp_fxnl3_bgibbs <- get_mcmc(natspp_fxnlneighbors_gjam_run3)
natspp_fxnl_geldf <- ctw_geldf(gelman.diag(mcmc.list(natspp_fxnl1_bgibbs, natspp_fxnl2_bgibbs, natspp_fxnl3_bgibbs)))
natspp_fxnl_gelfig <- ctw_gelfig(natspp_fxnl_geldf, cap = "CA model: natspp + fxnl neighbors, all subplots, 40k x 20k burnin, 3 chains")

# geweke test
natspp_fxnl1_gewdf <- ctw_gewdf(geweke.diag(natspp_fxnl1_bgibbs))
natspp_fxnl2_gewdf <- ctw_gewdf(geweke.diag(natspp_fxnl2_bgibbs))
natspp_fxnl3_gewdf <- ctw_gewdf(geweke.diag(natspp_fxnl2_bgibbs))
natspp_fxnl1_gewfig <- ctw_gewfig(natspp_fxnl1_gewdf, termfac = rev(varorder), cap = "all subplots, 40k x 20k burnin, chain 1", lpos = "top")
natspp_fxnl2_gewfig <- ctw_gewfig(natspp_fxnl2_gewdf, termfac = rev(varorder), cap = "all subplots, 40k x 20k burnin, chain 2", lpos = "top")
natspp_fxnl3_gewfig <- ctw_gewfig(natspp_fxnl3_gewdf, termfac = rev(varorder), cap = "all subplots, 40k x 20k burnin, chain 3", lpos = "top")


natspp_fxnlneighbors_gjam_run1$parameters$betaStandXTable %>%
  data.frame() %>%
  rownames_to_column() %>%
  separate(rowname, c("spp", "term"), sep = "_") %>%
  View()

allcoarse_gjam$parameters$betaStandXTable %>%
  data.frame() %>%
  rownames_to_column() %>%
  separate(rowname, c("spp", "term"), sep = "_") %>%
  View()


# plot distribution of prediction to observed values as density plots
natspp_fxnl1_preds <- data.frame(natspp_fxnlneighbors_gjam_run1$prediction$yPresentMu) %>%
  rownames_to_column("rowid") %>%
  mutate(herbicide = ifelse(grepl("NoHerb", rowid), 0, 1), 
         seedtrt = ifelse(grepl("Unseed", rowid), 0, 1),
         plot = factor(parse_number(rowid))) %>%
  # drop rowid then gather spp 
  dplyr::select(-rowid) %>%
  gather(spp, preds, names(natspp_fxnlneighbors_ydata))

# tidy ydata
natspp_fxnlneighbors_tidy <- natspp_fxnlneighbors_ydata %>%
  rownames_to_column("rowid") %>%
  mutate(herbicide = ifelse(grepl("NoHerb", rowid), 0, 1), 
         seedtrt = ifelse(grepl("Unseed", rowid), 0, 1),
         plot = factor(parse_number(rowid))) %>%
  # drop rowid then gather spp 
  dplyr::select(-rowid) %>%
  gather(spp, obs, names(natspp_fxnlneighbors_ydata))

natspp_fxnlneighbors_predsobs <- merge(natspp_fxnl1_preds, natspp_fxnlneighbors_tidy) %>%
  merge(natcoarse_xdata) %>%
  mutate(herbicide = factor(herbicide, labels = c("NoHerb", "Herb")),
         seedtrt = factor(seedtrt, labels = c("Unseed", "Seeded")))

# as density plot
ggplot(natspp_fxnlneighbors_predsobs) +
  geom_density(aes(preds-obs, fill = fulltrt), alpha = 0.5) +
  facet_wrap(~paste(spp, herbicide), scales = "free")

# just one species at a time
subset(natspp_fxnlneighbors_predsobs, spp == "FEMI" & seedtrt == "Seeded") %>%
  gather(type, cover, preds:obs) %>%
  ggplot() +
  geom_density(aes(cover, fill = paste(herbicide, seedtrt, type)), alpha = 0.5) +
  facet_wrap(~fulltrt)

# means all together
subset(natspp_fxnlneighbors_predsobs, seedtrt == "Seeded") %>%
  #gather(type, cover, preds:obs) %>%
  ggplot() +
  stat_summary(aes(fulltrt, preds-obs, color = paste(herbicide, seedtrt), group = paste(herbicide, seedtrt)), position = position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = 0)) +
  facet_wrap(~spp) +
  theme_minimal() +
  coord_flip()
# these are generally not great.. maybe due to different scales in cover (e.g., single digits vs nearly 100)


# run with block as effect and plot as random
natspp_fxnlneighbors_blockeffect_gjam_run1 <- gjam(formula = as.formula(~ block + seedtrt + herbicide + pptTrt * nutTrt), 
                                       xdata = natcoarse_xdata,
                                       ydata = natspp_fxnlneighbors_ydata,
                                       modelList = list(ng=40000, burnin=20000, typeNames = 'CA', random = "block"))

summary(natspp_fxnlneighbors_blockeffect_gjam_run1)
gjamPlot(natspp_fxnlneighbors_blockeffect_gjam_run1, plotPars = list(PLOTALLY = T, GRIDPLOTS = T))

# run with just the seeded plots
natspp_fxnlneighbors_blockeffect_seedsonly_gjam_run1 <- gjam(formula = as.formula(~ block + herbicide + pptTrt * nutTrt), 
                                                   xdata = natcoarse_xdata[grep("NatSeed", rownames(natcoarse_xdata)), ],
                                                   ydata = natspp_fxnlneighbors_ydata[grep("NatSeed", rownames(natcoarse_xdata)),],
                                                   modelList = list(ng=40000, burnin=20000, typeNames = 'CA', random = "block"))
summary(natspp_fxnlneighbors_blockeffect_seedsonly_gjam_run1)
gjamPlot(natspp_fxnlneighbors_blockeffect_seedsonly_gjam_run1, plotPars = list(PLOTALLY = T, GRIDPLOTS = T))


natspp_fxnlneighbors_blockeffect_seedsonlydimred_gjam_run1 <- gjam(formula = as.formula(~ block + herbicide + pptTrt * nutTrt), 
                                                             xdata = mutate(natcoarse_xdata[grep("NatSeed", rownames(natcoarse_xdata)), ], herbicide = factor(herbicide)),
                                                             ydata = natspp_fxnlneighbors_ydata[grep("NatSeed", rownames(natcoarse_xdata)),],
                                                             modelList = list(ng=50000, burnin=25000, 
                                                                              reductList = list(N= 6, r = 4),
                                                                              typeNames = 'CA', random = "block"))
summary(natspp_fxnlneighbors_blockeffect_seedsonlydimred_gjam_run1)
gjamPlot(natspp_fxnlneighbors_blockeffect_seedsonlydimred_gjam_run1, plotPars = list(PLOTALLY = T, GRIDPLOTS = T))

# -- 3. Individual native seeded with most abundant/frequent neighbors -----
# > ideal if neighbors have trait corresponding trait data to pair (see what's available)

# > prep same way as for fxnl cover
# prep just seeded-trt specdata since ambient comp plot will have more species in it (not fair comparison)
# choose max pct_cov per species per plot per herbicide and seeding trt, unless grass, then sum fxnl groups
seedtrt_ydata_tidy <- subset(natlong, grepl("Native", seedtrt)) %>%
  # drop TRCI -- only came up in 1 plot
  subset(code4 != "TRCI") %>%
  dplyr::select(plot:mon, code4:pct_cover) %>%
  spread(code4, pct_cover, fill = 0) %>%
  # gather any spp codes [codenames won't have lowcase character]
  gather(code4, pct_cover, names(.)[!grepl("[a-z]{2,}", names(.))]) %>%
  left_join(distinct(natlong, code4, fxnl_grp, nativity, native_seeded)) %>%
  mutate(nativity = ifelse(nativity == "Unknown", "Exotic", nativity),
         native_seeded = ifelse(native_seeded == "No", "Background", "Seeded"),
         coarse_fxnl = paste(native_seeded, nativity, fxnl_grp,sep = " ")) %>%
  # by not including month as grouping factor, take max cover val of species over season
  group_by(plot, seedtrt, herbicide, fulltrt, block, ppt_trt, nut_trt, fxnl_grp, nativity, native_seeded, coarse_fxnl, code4) %>%
  reframe(nobs = length(mon),
          cover = ifelse(nobs == 2 & coarse_fxnl == "Background Exotic Grass", 
                         mean(pct_cover, na.rm = T), 
                         #max(pct_cover, na.rm = T), 
                         max(pct_cover, na.rm = T))) %>%
  distinct() %>%
  # getting cholesky error below. try converting anything >0 and <0.05 to 0.25. Some values got smaller than trace value (0.01) because of averaging across months? 
  mutate(cover = ifelse(cover >0 & cover<=0.25, 0.25, 
                            # otherwise round to nearest 1 decimal to put on original scale of recording)
                            round(cover,1))) %>%
  # change nut_trt control to XC
  mutate(nut_trt = dplyr::recode(nut_trt, N = "XC"),
         fulltrt = gsub("N", "XC", fulltrt))

# ydata wide
seedtrt_ydata <- subset(seedtrt_ydata_tidy, grepl("Native", seedtrt), select = c(plot, seedtrt, herbicide, code4, cover)) %>%
  # make cover count data (trace amounts --> 1)
  # decimals = str_extract(cover, "(?<=[.])[0-9]+"),
  # tenths = as.numeric(substr(decimals, 1, 1)),
  # hundredths = as.numeric(substr(decimals, 2, 2)),
  # addct = ifelse(tenths == 5 & !is.na(tenths), 1+ hundredths, hundredths))
  spread(code4, cover, fill = 0) %>%
  mutate(abbrv_herb = ifelse(grepl("^Herbi", herbicide), "Herbi", "NoHerbi"),
         abbrv_seed2 = ifelse(grepl("Native", seedtrt), "NatSeed", "Unseed")) %>%
  unite(rowid, abbrv_herb, abbrv_seed2, sep = "", remove = T) %>%
  mutate(rowid = paste(rowid, plot, sep = "_")) %>%
  data.frame()
# set rownames
rownames(seedtrt_ydata) <- seedtrt_ydata$rowid
# pull xdata
seedtrt_xdata <- natcoarse_xdata[rownames(natcoarse_xdata) %in% rownames(seedtrt_ydata),]

# aggregate low abundance, low frequency to other by functional group
# check total abundance over exp by sp
seedtrt_abds <- colSums(seedtrt_ydata[, grep("[A-Z]+", names(seedtrt_ydata))])
sort(seedtrt_abds)
# check frequency by sp
seedtrt_freqs <- colSums(seedtrt_ydata[, grep("[A-Z]+", names(seedtrt_ydata))] >0) 
sort(seedtrt_freqs) # set cutoff at 10 plots
lowfrq <- seedtrt_freqs[seedtrt_freqs < 10]
# check that remaining spp have tot abundances generally over 10%
seedtrt_abds[!names(seedtrt_abds) %in% names(lowfrq)]
# check median cover
seedtrt_medians <- sapply(seedtrt_ydata[, grep("[A-Z]+", names(seedtrt_ydata))], median)
sort(seedtrt_medians) # not many have a median value > 0
seedtrt_means <- sapply(seedtrt_ydata[, grep("[A-Z]+", names(seedtrt_ydata))], mean)
sort(round(seedtrt_means,4))
seedtrt_max <- sapply(seedtrt_ydata[, grep("[A-Z]+", names(seedtrt_ydata))], max)
par(mfrow = c(1,2))
plot(seedtrt_means ~ seedtrt_freqs)
plot(seedtrt_max ~ seedtrt_freqs)

# see if there are certain species that only occur in types of treatments before use 10 as cutoff
seedtrt_sppfreq_bytrt <- group_by(seedtrt_ydata_tidy, code4) %>%
  # scan for things that have more than minmal cover
  mutate(nglobal = sum(cover > 1),
         nhits = sum(cover > 0),
         maxcov = max(cover),
         other = nglobal <5) %>%
  group_by(ppt_trt, nut_trt, code4, nglobal, nhits, maxcov, other) %>%
  reframe(npres = sum(cover > 1),)
# screen out plots in fewer than 10 total plots with a max cov 5 or less (ESCA has >1% cov in 10 plots; is present in more than that)
# seedtrt_sppfreq_bytrt$other <- with(seedtrt_sppfreq_bytrt, (nglobal < 10 & maxcov <5) | nglobal < 5)
# get cholesky error using ^, try nglobal must be 15 or more
seedtrt_sppfreq_bytrt$other <- with(seedtrt_sppfreq_bytrt, nglobal < 15 & !code4 %in% nats) # if don't specify this esca gets dropped
grpspp <- subset(spplist, code4 %in% with(seedtrt_sppfreq_bytrt, unique(code4[other])))
# try changing all trace to 0.025 to see if that helps with sigma (trace value is arbitrary)

grpnatforbs <- seedtrt_ydata[,with(grpspp, code4[grepl("Native", nativity) & !grepl("Gras", fxnl_grp)])]
grpnatgrams <- data.frame(seedtrt_ydata[with(grpspp, code4[grepl("Native", nativity) & grepl("Gras", fxnl_grp)])]) # stpu is just present in one plot. ignore for species level
grpexoforbs <- seedtrt_ydata[,with(grpspp, code4[!grepl("Native", nativity) & grepl("Forb", fxnl_grp)])]
grpexonfix <- seedtrt_ydata[,with(grpspp, code4[!grepl("Native", nativity) & grepl("^N", fxnl_grp)])]
grpexograms <- seedtrt_ydata[,with(grpspp, code4[!grepl("Native", nativity) & grepl("Gras", fxnl_grp)])]

# sum across
grpnatforbs$NatForbs <- rowSums(grpnatforbs)
grpexoforbs$IntroForbs <- rowSums(grpexoforbs)
grpexonfix$IntroNFix <- rowSums(grpexonfix)
grpexograms$IntroGrass <- rowSums(grpexograms)

seedtrt_ydata_common <- seedtrt_ydata[with(seedtrt_sppfreq_bytrt, unique(code4[!other]))] %>%
  cbind(grpnatforbs[c("NatForbs")], grpexonfix[c("IntroNFix")], grpexoforbs[c("IntroForbs")], grpexograms[c("IntroGrass")])

# check rownames in ydata match xdata before run
summary(rownames(seedtrt_ydata_common) == rownames(seedtrt_xdata)) # good to go

# run individual species gjam model

# run model -- start with as many runs as nat spp x fxnl needed
# > this will take longer to run, don't predict x to speed up
natseed_spp_gjam_run5<- gjam(formula = as.formula(~ herbicide + pptTrt * nutTrt), 
                                       xdata = mutate(seedtrt_xdata, herbicide = factor(herbicide)), 
                                       # test on most common spp to see if can run
                                       ydata = seedtrt_ydata_common,# nats[-grep("TRCI", nats)]))],
                                       modelList = list(ng=20000, burnin=10000, 
                                                        PREDICTX = F,
                                                        # dim red
                                                        reductList = list(N = 15, r = 10), ## dimred error is because herbicide was not factor
                                                        typeNames = 'CA', random = "block"
                                                        ))
summary(natseed_spp_gjam_run5)
data.frame(natseed_sppfc_gjam_run1$parameters$betaStandXTable) %>%
  rownames_to_column() %>%
  separate(rowname, into = c("spp", "term"), sep = "_") %>%
  View()
# cannot run at 20kx10k without dim red, get cholesky error: Error in chol.default(SS) : 
# the leading minor of order 25 is not positive definite .. which means cannot determine true covariance matrix (not enough variation in some spp.. values close to 0 may be causing problem?)

gjamPlot(natseed_spp_gjam_run1, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, SAVEPLOTS = T, 
                                                         outFolder = "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/chapter/natseed_spp_gjam_r1/"))
gjamPlot(natseed_spp_gjam_run2, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, SAVEPLOTS = T, 
                                                outFolder = "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/chapter/natseed_spp_gjam_r2/"))
gjamPlot(natseed_spp_gjam_run3, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, SAVEPLOTS = T, 
                                                outFolder = "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/chapter/natseed_spp_gjam_r3/"))
gjamPlot(natseed_spp_gjam_run4, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, SAVEPLOTS = T, 
                                                outFolder = "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/chapter/natseed_spp_gjam_r4/"))
gjamPlot(natseed_spp_gjam_run5, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, SAVEPLOTS = T, 
                                                outFolder = "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/chapter/natseed_spp_gjam_r5/"))

# save mods to avoid rerun
saveRDS(natseed_spp_gjam_run1, paste0(outpath, "mods/natseed_spp_gjam_ca30kx15k_run1.rdata"))

# check chains
natseed_spp1_bgibbs <- get_mcmc(natseed_spp_gjam_run1)
natseed_spp2_bgibbs <- get_mcmc(natseed_spp_gjam_run2)

# geweke diagnostics
natseed_spp1_gewdf <- ctw_gewdf(geweke.diag(natseed_spp1_bgibbs))
natseed_spp2_gewdf <- ctw_gewdf(geweke.diag(natseed_spp2_bgibbs))
natseed_spp1_gewfig <- ctw_gewfig(natseed_spp1_gewdf, termfac = rev(varorder), dropvars = "inter", cap = "seeded trts only, chain 1", lpos = "top")
natseed_spp2_gewfig <- ctw_gewfig(natseed_spp2_gewdf, termfac = rev(varorder), dropvars = "inter", cap = "seeded trts only, chain 2", lpos = "top")

# gelman diagnostics
natseed_gelman <- gelman.diag(mcmc.list(natseed_spp1_bgibbs, natseed_spp2_bgibbs))
natseed_geldf <- ctw_geldf(natseed_gelman) 
natseed_gelfig <- ctw_gelfig(natseed_geldf, cap = "seeded trts only, 2 chains", dropvars = "inter") + theme(legend.position = "top")

# try as fractional cover to see if that helps model converge
seedtrt_ydata_commonfc <- vegan::decostand(seedtrt_ydata_common, method = "total")
# fc gjam
natseed_sppfc_gjam_run2 <- gjam(formula = as.formula(~ herbicide + pptTrt * nutTrt), 
                              xdata = mutate(seedtrt_xdata, herbicide = factor(herbicide)),
                              # test on most common spp to see if can run
                              ydata = seedtrt_ydata_commonfc,#[unique(c(with(seedtrt_sppfreq_bytrt,code4[nglobal >= 20]), nats[-grep("TRCI", nats)]))],
                              modelList = list(ng=20000, burnin=10000, 
                                               PREDICTX = F,
                                               # dim red
                                               reductList = list(N = 15, r = 10),
                                               typeNames = 'FC', random = "block"))

summary(natseed_sppfc_gjam_run2)
gjamPlot(natseed_sppfc_gjam_run1, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, SAVEPLOTS = T, 
                                                outFolder = "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/chapter/natseed_sppfc_gjam_r1/"))
gjamPlot(natseed_sppfc_gjam_run2, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, SAVEPLOTS = T, 
                                                  outFolder = "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/chapter/natseed_sppfc_gjam_r2_dimred/"))




## ========== ESA 2023 CODE BELOW ===============
# -- ESA 2023 model: Don't split hillslope, block as random effect -----
## this is the model to use
# native spp seeded (less TRCI) + 4 types of cover groups per plot
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
  mutate(nut_trt = factor(nut_trt, levels = c("N", "F", "C")),
         ppt_trt = factor(ppt_trt, levels = c("XC", "D", "W")),
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
  ggplot(aes(predictor, val)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))

beta_df <- data.frame(natcoarse_gjam$parameters$betaMu) %>%
  mutate(predictor = rownames(.)) %>%
  gather(response, beta, 1:(ncol(.)-1))
betaSE_df <- data.frame(natcoarse_gjam$parameters$betaSe) %>%
  mutate(predictor = rownames(.)) %>%
  gather(response, betaSE, 1:(ncol(.)-1))
beta_df <- left_join(beta_df, betaSE_df)
# clean up treatments
beta_df <- mutate(beta_df, ppt_trt = ifelse(grepl("1d", predictor), "D", 
                                          ifelse(grepl("2w", predictor), "W", "XC")),
                nut_trt = ifelse(grepl("2comp", predictor), "C",
                                 ifelse(grepl("1f", predictor), "F", "XC")),
                herbicide = ifelse(grepl("1herb", predictor), "Herbicided", "Non-Herbicided"))
# NA anything that is intercept
beta_df[beta_df$predictor == "intercept", c("nut_trt", "ppt_trt", "herbicide")] <- NA


ggplot(subset(beta_df, predictor != "intercept"), aes(predictor, beta)) +
  geom_hline(aes(yintercept = 0), lty = 2, col = "red") +
  geom_errorbar(aes(ymax = betaSE + beta, ymin = beta-betaSE)) +
  geom_point() +
  #theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(response~ .) +
  coord_flip()

natcoarse_gjam$parameters$betaStandXWmu
gjamOrdination(natcoarse_gjam, method = "NMDS")



# -- FIGS FOR ESA ------
# path to figs
figpath <- "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/esa2023_presfigs/" 

# plot coefficients for each species ()
beta_df <- data.frame(natcoarse_gjam$parameters$betaMu) %>%
  mutate(predictor = rownames(.)) %>%
  gather(response, beta, 1:(ncol(.)-1))
betaSE_df <- data.frame(natcoarse_gjam$parameters$betaSe) %>%
  mutate(predictor = rownames(.)) %>%
  gather(response, betaSE, 1:(ncol(.)-1))
beta_df <- left_join(beta_df, betaSE_df)

betacoeff_df <- data.frame(natcoarse_gjam$parameters$betaStandXWTable)
betacoeff_df$rowid <- rownames(betacoeff_df)
betacoeff_df$response <- str_extract(betacoeff_df$rowid, "^[:alpha:]+_")
betacoeff_df$response <- gsub("_", "", betacoeff_df$response)
betacoeff_df$predictor <- gsub("^.*_", "", betacoeff_df$rowid)
betacoeff_df$group <- ifelse(grepl("Background", betacoeff_df$response), "Background", "Native Seeded")
betacoeff_df$lifeform <- ifelse(grepl("BRCA|FEMI|Gra", betacoeff_df$response), "Grass", "Forb")


# -- plot results -----
# make plant colors
plantnames <- unique(betacoeff_df$response)
natgrasses <- plantnames[grepl("FEM|BR", plantnames)]
natgrass_cols <- RColorBrewer::brewer.pal((length(natgrasses)+2), "Greens")[3:(length(natgrasses)+2)]
exogr <- "BackgroundExoticGrass"
exogr_cols <- RColorBrewer::brewer.pal(6, "YlGnBu")[4]
natforbs <- plantnames[grepl("NativeFor|NEMA|ESCA", plantnames)]
natforbs_cols <- RColorBrewer::brewer.pal((length(natforbs)+1), "Purples")[2:(length(natforbs)+1)]
exoforb <- "BackgroundExoticForb"
exoforb_cols <- RColorBrewer::brewer.pal(6, "YlOrBr")[5]
nfix <- plantnames[grepl("Nfix", plantnames)]
nfix_cols <- RColorBrewer::brewer.pal(6, "YlOrBr")[3]

plant_cols <- c(natgrass_cols, exogr_cols, natforbs_cols, exoforb_cols, nfix_cols)
plants1 <- c(natgrasses, exogr, natforbs, exoforb, nfix)
names(plant_cols) <- plants1

plant_cols2 <- plant_cols
plants2 <- c("B. carinatus", "F. microstachys", "Non-native grasses", "Background native forbs", "E. californica",
                          "N. maculata", "Non-native forbs", "Non-native N-fixers")
names(plant_cols2) <- plants2

values = c('Exotic Forb' = "red", 
           'Native Forb' = "purple", 
           'Exotic Grass' = "forestgreen", 
           'Native Grass' = "dodgerblue", 
           "Exotic N-fixer" = "orange")

signif_betas <- subset(betacoeff_df, grepl("[*]", sig95)) %>%
  # drop no herbicide because it's flip of herbicide
  subset(!grepl("0noHerb", predictor)) %>%
  # note whether competitive effect or environmental
  mutate(effect = ifelse(grepl("herb", predictor), "Competitive", "Environmental"))

ggplot(signif_betas, aes(predictor,Estimate)) +
  geom_point() +
  facet_wrap(~effect) +
  coord_flip()

subset(signif_betas, effect == "Competitive") %>%
  mutate(response = factor(response, 
                           levels = c("BackgroundExoticGrass", "BackgroundExoticForb", "FEMI", "BRCA"),
         labels = c("Non-native grasses", "Non-native forbs", "F. microstachys", "B. carinatus")
         )) %>%
ggplot(aes(response,Estimate)) +
  geom_hline(aes(yintercept = 0), alpha = 0.5) +
  labs(y = "Effect of herbicide (with 95% credible interval)",
       x = NULL) +
  geom_errorbar(aes(ymin = CI_025, ymax = CI_975), width = 0.1) +
  geom_point(aes(fill = response), pch = 21, size = 4) +
  scale_fill_manual(values = c(exogr_cols, exoforb_cols, natgrass_cols)) +
  theme(axis.text = element_text(size = 11),
        #axis.title = element_text(size = ),
        legend.position = "none") +
  coord_flip()
ggsave(paste(figpath, "gjam_herbsig.pdf"), units = "in", height = 4, width = 5)

# also show raw?



specNames <- colnames(natcoarse_ydata)
specColor_natdiv <- unname(plant_cols[specNames])

gjamPlot(natcoarse_gjam, plotPars = list(GRIDPLOTS=T, specColor = specColor_natdiv, #cex = 0.75,
                                         #PLOTALLY = T, #SAVEPLOTS = T,
                                         outfolder = "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/esa2023_presfigs/"))
# save gjammodel
saveRDS(natcoarse_gjam, "/Users/scarlet/Documents/suding_lab/SFREC/Compost/natdiv/gjam_model_esa.rdata")



# -- make ordination using coarse fxnl groups and natives ---
library(vegan)

natcoarse_nmds <- metaMDS(natcoarse_ydata, k = 2, trymax = 100)
plot(natcoarse_nmds)
adonis2(natcoarse_ydata ~ herbicide + nutTrt * pptTrt, data = natcoarse_xdata, permutations = how(blocks = natcoarse_xdata$block))

mds_natcoarse_results <- rename(natcoarse_xdata[c("plot", "herbicide", "seedtrt", "hillpos", "fulltrt", "nutTrt", "pptTrt")]) %>%
  left_join(subset(trtkey)) %>%
  cbind(natcoarse_nmds$points) %>%
  mutate(herbicide_pretty = ifelse(herbicide == "0noHerb", "Non-herbicided", "Herbicided"),
         herbicide_pretty = factor(herbicide_pretty, levels = c("Non-herbicided", "Herbicided")))

mds_natcoarse_spp <- data.frame(natcoarse_nmds$species) %>%
  mutate(spp = rownames(.)) %>%
  mutate(pretty_name = gsub("Background", "", spp),
         pretty_name = gsub("Grass", "grasses", pretty_name),
         pretty_name = gsub("Forb", "forbs", pretty_name),
         pretty_name = gsub("Nfixer", "N-fixers", pretty_name),
         pretty_name = gsub("Exotic", "Non-native ", pretty_name))
# clean up names
mds_natcoarse_spp$pretty_name[mds_natcoarse_spp$spp == "FEMI"] <- "F. microstachys"
mds_natcoarse_spp$pretty_name[mds_natcoarse_spp$spp == "BRCA"] <- "B. carinatus"
mds_natcoarse_spp$pretty_name[mds_natcoarse_spp$spp == "NEMA"] <- "N. maculata"
mds_natcoarse_spp$pretty_name[mds_natcoarse_spp$spp == "ESCA"] <- "E. californica"
mds_natcoarse_spp$pretty_name[mds_natcoarse_spp$pretty_name == "Nativeforbs"] <- "Background native forbs"

# make hulls (quickplot.. spider will be easier on eyes to show maybe)
natcoarse_hulldat <- data.frame()
for(h in unique(mds_natcoarse_results$herbicide)){
  for(f in unique(mds_natcoarse_results$fulltrt)){
    grphull <- mds_natcoarse_results[mds_natcoarse_results$fulltrt == f & mds_natcoarse_results$herbicide == h, ][chull(mds_natcoarse_results[mds_natcoarse_results$fulltrt == f & mds_natcoarse_results$herbicide == h, 
                                                                                                                          c("MDS1", "MDS2")]), ]
    # add to df
    natcoarse_hulldat <- rbind(natcoarse_hulldat, grphull)  
  }
}
# add herbicide pretty
natcoarse_hulldat <- left_join(natcoarse_hulldat, distinct(mds_natcoarse_results, herbicide, herbicide_pretty))

# make treatments factors
mds_natcoarse_results <- mutate(mds_natcoarse_results, 
       ppt_trt = factor(ppt_trt, levels = c("D", "XC", "W")),
       #ppt_trt = relevel(ppt_trt, "D"),
       nut_trt = factor(nut_trt, levels = c("C", "F", "N"), labels = c("Compost", "Fertilizer", "No amendment")))
natcoarse_hulldat <- mutate(natcoarse_hulldat, 
                            ppt_trt = factor(ppt_trt, levels = c("D", "XC", "W")),
                            nut_trt = factor(nut_trt, levels = c("C", "F", "N"), labels = c("Compost", "Fertilizer", "No amendment")))
# plot treatment hulls
mds_natcoarse_trthulls <- ggplot(mds_natcoarse_results,
       aes(MDS1, MDS2, col = ppt_trt, fill = ppt_trt)) +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey") +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey") +
  geom_polygon(data=natcoarse_hulldat,  aes(fill=ppt_trt,group=ppt_trt, col = ppt_trt), alpha=0.35) +
  geom_point(aes(shape = block<3), col = "grey30", alpha = 0.5, size = 3) + #col = "grey30",
  geom_point(aes(shape = block<3), alpha = 0.5, size = 3) + #col = "grey30",
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
        legend.text = element_text(size = 10)) +
  facet_grid(herbicide_pretty~nut_trt)

# plot species
mds_natcoarse_sppfig <- ggplot(mds_natcoarse_results, aes(MDS1, MDS2)) +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey") +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey") +
  geom_polygon(data=natcoarse_hulldat,  fill = "transparent", col = "transparent", alpha=0.05) +
  geom_point(data = mds_natcoarse_spp, aes(MDS1, MDS2, color = spp), size = 3) +
  #geom_text(data = subset(mds_natcoarse_spp, grepl("Backg", spp)), aes(MDS1, MDS2, color = spp, label = spp), size = 3) +
  #geom_label(data = subset(mds_natcoarse_spp, !grepl("Backg", spp)), aes(MDS1, MDS2, color = spp, label = spp), size = 3) +
  ggrepel::geom_text_repel(data = mds_natcoarse_spp, aes(MDS1, MDS2, col = spp, label = pretty_name, fontface = nchar(spp))) +
  #ggrepel::geom_text_repel(data = subset(mds_natcoarse_spp, nchar(spp) > 4), aes(MDS1, MDS2, col = spp, label = pretty_name)) +
  #ggrepel::geom_text_repel(data = subset(mds_natcoarse_spp, nchar(spp) == 4), aes(MDS1, MDS2, col = spp, label = pretty_name), fontface = "italic") +
  scale_color_manual(name = NULL, values = plant_cols) +
  scale_fill_manual(name = NULL, values = plant_cols) +
  theme(legend.position = "none")
        # legend.justification = c("left", "top"),
        # legend.key.size = unit(0.1, "pt"))

cowplot::plot_grid(mds_natcoarse_trthulls,
                   mds_natcoarse_sppfig, 
                   nrow = 1,
                   rel_widths = c(1.1,0.8))

ggsave(paste0(figpath,"fxnlgrp_natspp_nmds.pdf"), units = "in", width = 10, height = 5)


# -- plot cover ----
# functional groups only first, then fxnl groups and native species IDs
natcoarse_seedtrts <- subset(natcoarse_all, seedtrt != "Unseeded") %>%
  # drop bg native nfix and bg natgram because rare
  subset(!grepl("Background Native G|Native N-", coarse_fxnl)) %>%
  mutate(ppt_trt = factor(ppt_trt, levels = c("D","XC", "W")))
# rbind natives seeded

seededspp_long <- left_join(seededspp, distinct(subset(natcoarse_seedtrts, select =c(plot:nut_trt)))) %>%
  gather(coarse_fxnl, totcov, BRCA:TRCI) %>%
  subset(!grepl("TRCI", coarse_fxnl)) %>%
  #dplyr::select(names(natcoarse_seedtrts)) %>%
  rbind(subset(natcoarse_seedtrts, !grepl("Seeded", coarse_fxnl), select = names(.))) %>%
  mutate(pretty_name = trimws(gsub("Background", "", coarse_fxnl)),
         pretty_name = gsub("Grass", "grasses", pretty_name),
         pretty_name = gsub("Forb", "forbs", pretty_name),
         pretty_name = gsub("N-fixer", "N-fixers", pretty_name),
         pretty_name = gsub("Exotic", "Non-native", pretty_name))
# clean up names
seededspp_long$pretty_name[seededspp_long$coarse_fxnl == "FEMI"] <- "F. microstachys"
seededspp_long$pretty_name[seededspp_long$coarse_fxnl == "BRCA"] <- "B. carinatus"
seededspp_long$pretty_name[seededspp_long$coarse_fxnl == "NEMA"] <- "N. maculata"
seededspp_long$pretty_name[seededspp_long$coarse_fxnl == "ESCA"] <- "E. californica"
seededspp_long$pretty_name[seededspp_long$pretty_name == "Native forbs"] <- "Background native forbs"


ggplot(subset(natcoarse_seedtrts, grepl("Back", coarse_fxnl)), aes(ppt_trt, totcov, fill = coarse_fxnl, group = paste(coarse_fxnl))) +
  stat_summary(geom = "bar", position = position_dodge()) +
  stat_summary(geom = "errorbar", position = position_dodge()) +
  # scale_color_manual(name = "Precipitation\ntreatment", 
  #                    labels = c("D" = "Drought", "XC" = "Ambient", "W" = "Wet"),
  #                    values = RColorBrewer::brewer.pal(4,"Blues")[2:4]) +
  facet_grid(herbicide~nut_trt)

# don't show herbicide
ggplot(subset(natcoarse_seedtrts, grepl("Back", coarse_fxnl) & !grepl("Non", herbicide)), aes(nut_trt, totcov, fill = coarse_fxnl, group = paste(coarse_fxnl, nut_trt))) +
  #geom_boxplot() +
  stat_summary(geom = "bar", position = position_dodge()) +
  # scale_color_manual(name = "Precipitation\ntreatment", 
  #                    labels = c("D" = "Drought", "XC" = "Ambient", "W" = "Wet"),
  #                    values = RColorBrewer::brewer.pal(4,"Blues")[2:4]) +
  facet_grid(~ppt_trt)

ggplot(subset(natcoarse_seedtrts, !grepl("Back", coarse_fxnl) & !grepl("Non", herbicide)), aes(nut_trt, totcov, fill = coarse_fxnl, group = paste(coarse_fxnl, nut_trt))) +
  #geom_boxplot() +
  stat_summary(geom = "bar", position = position_dodge()) +
  # scale_color_manual(name = "Precipitation\ntreatment", 
  #                    labels = c("D" = "Drought", "XC" = "Ambient", "W" = "Wet"),
  #                    values = RColorBrewer::brewer.pal(4,"Blues")[2:4]) +
  facet_grid(~ppt_trt)

ggplot(subset(natcoarse_seedtrts, !grepl("Back", coarse_fxnl) & grepl("Non", herbicide)), aes(nut_trt, totcov, fill = coarse_fxnl, group = paste(coarse_fxnl, nut_trt))) +
  #geom_boxplot() +
  stat_summary(geom = "bar", position = position_dodge()) +
  # scale_color_manual(name = "Precipitation\ntreatment", 
  #                    labels = c("D" = "Drought", "XC" = "Ambient", "W" = "Wet"),
  #                    values = RColorBrewer::brewer.pal(4,"Blues")[2:4]) +
  facet_grid(~ppt_trt)






# -- anova for not seeded ----
summary(aov(totcov ~ coarse_fxnl * ppt_trt * herbicide * nut_trt, data = natcoarse_seedtrts))

summary(aov(totcov ~ herbicide * nut_trt * ppt_trt, data = subset(natcoarse_seedtrts, coarse_fxnl == "Background Exotic Grass"))) # only precip and herb matter for exograms (main effects only)
summary(aov(totcov ~ herbicide * nut_trt * ppt_trt, 
            data = subset(natcoarse_seedtrts, coarse_fxnl == "Background Exotic Forb"))) # herbicide and nutrient treatment matter for exo forbs (including their interaction)

summary(aov(totcov ~ herbicide * nut_trt * ppt_trt, 
            data = subset(natcoarse_seedtrts, coarse_fxnl == "Background Exotic N-fixer"))) #ppt marginally matters

summary(aov(totcov ~ herbicide * nut_trt * ppt_trt, 
            data = subset(natcoarse_seedtrts, coarse_fxnl == "Seeded Native Forb"))) # just ppt matters


summary(aov(totcov ~ herbicide * nut_trt * ppt_trt, 
            data = subset(natcoarse_seedtrts, coarse_fxnl == "Seeded Native Grass"))) # herbicide, ppt and their interaction

summary(aov(FEMI ~ herbicide * nut_trt * ppt_trt, 
            data = subset(neighborhood_counts))) # just herb and interaction with herb:ppt

summary(aov(BRCA ~ herbicide * nut_trt * ppt_trt, 
            data = subset(neighborhood_counts))) # herb, nut, ppt, nut:ppt

summary(aov(NEMA ~ herbicide * nut_trt * ppt_trt, 
            data = subset(neighborhood_counts))) # ppt strong, marginal nut

summary(aov(ESCA ~ herbicide * nut_trt * ppt_trt, 
            data = subset(neighborhood_counts))) # marginal 3-way. native seeded forbs just happy to be there

summary(aov(totcov ~ herbicide * nut_trt * ppt_trt, 
            data = subset(natcoarse_seedtrts, coarse_fxnl == "Background Native Forb"))) # bg nat forbs don't care


# just show herbicided. show exogr and exoforb + seeded nats

natcoarse_seedtrts_short <- subset(natcoarse_seedtrts, grepl("Exotic [G|F]|Seeded", coarse_fxnl)) %>%
 mutate(pretty_names = trimws(gsub("Background", "", coarse_fxnl)),
        pretty_names = gsub("Forb", "forbs", pretty_names),
        pretty_names = gsub("Grass", "grasses", pretty_names),
        pretty_names = gsub("Exotic", "Non-native", pretty_names),
        pretty_names = gsub("Native", "native", pretty_names)
        ) %>%
  mutate(nut_trt = factor(nut_trt, levels = c("C", "F", "XC"), labels = c("Compost", "Fertilizer", "No amendment")),
         herbicide = factor(herbicide, levels = c("Non-herbicided", "Herbicided")))
simple_plantcols <- plant_cols2[grepl("B. |N. |Non-native [g|f]", names(plant_cols2))]
names(simple_plantcols)[grep("B. ", names(simple_plantcols))] <- "Seeded native grasses"
names(simple_plantcols)[grep("N. ", names(simple_plantcols))] <- "Seeded native forbs"
#simple_plantcols #<-

ggplot(natcoarse_seedtrts_short, aes(ppt_trt, totcov, col = pretty_names, group = coarse_fxnl)) +
  stat_summary(position = position_dodge(width = 0.5)) +
  scale_color_manual(name = "Functional group", values = simple_plantcols) +
  facet_grid(herbicide~nut_trt, scales = "free_y") +
  labs(y = "Mean cover (%)", x = "Precipitation treatment") +
  theme(strip.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 11))

ggsave(paste0(figpath, "simple_cover_summaries.pdf"), unit = "in", width = 6.5, height = 4)

# all groups (it's a lot)
seededspp_long %>%
  mutate(pretty_name = factor(pretty_name,
                              levels = c("Non-native forbs", "Non-native grasses", "Non-native N-fixers",
                                         "Background native forbs",
                                         "E. californica", "N. maculata",
                                         "B. carinatus", "F. microstachys")),
         nut_trt = factor(nut_trt, levels = c("C", "F", "XC"), labels = c("Compost", "Fertilizer", "No amendment"))) %>%
  subset(herbicide == "Herbicided") %>%
ggplot(aes(as.numeric(ppt_trt), totcov, col = pretty_name, group = pretty_name)) +
  geom_vline(aes(xintercept = 1.5), linetype = 3, col = "grey50") +
  geom_vline(aes(xintercept = 2.5), linetype = 3, col = "grey50") +
  stat_summary(position = position_dodge(width = 0.75)) +
  scale_x_continuous(breaks =c(1,2,3), labels = c("D", "XC", "W")) +
  scale_color_manual(name = NULL, values = plant_cols2) +
  facet_grid(herbicide~nut_trt, scales = "free_y") +
  labs(y = "Mean cover (%)", x = "Precipitation treatment") +
  theme(strip.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 11))

ggsave(paste0(figpath, "spp_cover_summaries.pdf"), unit = "in", width = 10, height = 3)


# try to cover natives
testherb <- seededspp_long %>%
  mutate(pretty_name = factor(pretty_name,
                              levels = c("Non-native forbs", "Non-native grasses", "Non-native N-fixers",
                                         "Background native forbs",
                                         "E. californica", "N. maculata",
                                         "B. carinatus", "F. microstachys")),
         nut_trt = factor(nut_trt, levels = c("C", "F", "XC"), labels = c("Compost", "Fertilizer", "No amendment"))) %>%
  subset(herbicide == "Herbicided") #%>%

ggplot(testherb, aes(as.numeric(ppt_trt), totcov, col = pretty_name, group = pretty_name)) +
  geom_vline(aes(xintercept = 1.5), linetype = 3, col = "grey50") +
  geom_vline(aes(xintercept = 2.5), linetype = 3, col = "grey50") +
  stat_summary(aes(alpha = grepl("Non-nat|Back", pretty_name)), position = position_dodge(width = 0.75)) +
  #stat_summary(data = subset(testherb, grepl("Non-nat|Back", pretty_name)), position = position_dodge(width = 0.75)) +
  scale_x_continuous(breaks =c(1,2,3), labels = c("D", "XC", "W")) +
  scale_color_manual(name = NULL, values = plant_cols2) +
  scale_alpha_manual(values = c(`FALSE`=0,`TRUE` = 1), guide = "none") + 
  facet_grid(herbicide~nut_trt, scales = "free_y") +
  labs(y = "Mean cover (%)", x = "Precipitation treatment") +
  theme(strip.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 11))

ggsave(paste0(figpath, "spp_cover_summaries_ghostnats.pdf"), unit = "in", width = 10, height = 3)


# try to cover native grams
ggplot(testherb, aes(as.numeric(ppt_trt), totcov, col = pretty_name, group = pretty_name)) +
  geom_vline(aes(xintercept = 1.5), linetype = 3, col = "grey50") +
  geom_vline(aes(xintercept = 2.5), linetype = 3, col = "grey50") +
  stat_summary(aes(alpha = !grepl("B. |F. ", pretty_name)), position = position_dodge(width = 0.75)) +
  #stat_summary(data = subset(testherb, grepl("Non-nat|Back", pretty_name)), position = position_dodge(width = 0.75)) +
  scale_x_continuous(breaks =c(1,2,3), labels = c("D", "XC", "W")) +
  scale_color_manual(name = NULL, values = plant_cols2) +
  scale_alpha_manual(values = c(`FALSE`=0,`TRUE` = 1), guide = "none") + 
  facet_grid(herbicide~nut_trt, scales = "free_y") +
  labs(y = "Mean cover (%)", x = "Precipitation treatment") +
  theme(strip.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 11))

ggsave(paste0(figpath, "spp_cover_summaries_ghostnatgrams.pdf"), unit = "in", width = 10, height = 3)
