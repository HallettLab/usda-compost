# native recruitment composition modeling with gllvm


# -- SETUP -----
#devtools::install_github("https://github.com/BertvanderVeen/gllvm/tree/roweff") #use this version to use nested random effects
library(gllvm)
library(vegan) # need to call gllvm::ordiplot


# source prepped cover and trtkey data
# > sets default settings for datpath, plot theme, no factors by default, loads tidyverse
#source("Native-analysis/native_prep_cover.R", print.eval = F, echo = F) # trait script loads all datasets needed
# source prepped trait data
source("Native-analysis/native_prep_traits.R", print.eval = F, echo = F)
# ^note: warning message about NA coercian is okay. just converting notes about no sample available to NA, rest of numeric vals become numberic


# -- PREP Y, X, AND TR DATA  ----
# ydat = species (row = site, cols = sp abundances)
# xdat = environment (row = site, cols = env vars)
# TR = trait data (sp as rownames, cols = traits)

# first make rowid key
rowkey <- subset(widesp_count, select = c(plot:herbicide)) %>%
  # make rowids: plot + seedtrt + herbicide
  mutate(rowid = paste(plot, seedtrt, herbicide, sep = "_"),
         rowid = gsub("herbi.*", "Herb", rowid, ignore.case = T),
         rowid = gsub("Non-|Un", "No", rowid, ignore.case = T),
         rowid = gsub("seeded", "Seed", rowid, ignore.case = T),
         rowid = gsub("Native ", "Nat", rowid, ignore.case = T))

# make xdat -- needs to have env predictors and any random effects to use
xdat <- rowkey
rownames(xdat) <- rowkey$rowid
# drop anything not using in models
xdat <- xdat[!names(xdat) %in% c("rowid")] #"xpos", "ypos" <-- keep these if want to use as randos

# species in count form
ydat_count <- merge(widesp_count, rowkey) # merge to be sure

# assign rownames and drop all else but species
rownames(ydat_count) <- ydat_count$rowid
ydat_count <- ydat_count[names(ydat_count)[grepl("^[A-Z]+", names(ydat_count))]]
# check colsums
checkcolSums <- sapply(ydat_count, function(x) sum(x>0))
sort(checkcolSums)/nrow(ydat_count)

# for first test, drop anything not in at least 5% of plots
nrow(ydat_count)*.05 # at least 7 plots

# pull rare (<5% plots to group by fxnl group)
# > to not over-inflate count on rare, sum on pct cover then transform to count
ydat_raresp <- names(ydat_count)[checkcolSums < nrow(ydat_count)*.05]
# keep pct cover df separate to recycle for beta distribution
ydat_rare <- subset(tidysp0, code4 %in% ydat_raresp) %>%
  # join spp functional info
  left_join(distinct(natambsp_cov, code4, coarse_fxnl)) %>%
  group_by(plot, seedtrt, herbicide, coarse_fxnl) %>%
  reframe(pct_cover = sum(pct_cover)) %>%
  distinct() %>%
  # shorten fxln group names
  mutate(coarse_fxnl = gsub("Background | |-|er", "", coarse_fxnl),
         # trace gets 1, all else gets rounded; add 0.1 sound 0.5 rounds (R doesn't round it always)
         count_cover = ifelse(pct_cover == 0, 0, 
                              ifelse(pct_cover == 0.01, 1, round(pct_cover + 1.1)))) %>%
  # join rowid
  merge(rowkey)

# convert to count as in cover prep script
ydat_count_rare <- subset(ydat_rare, select = -pct_cover) %>%
  spread(coarse_fxnl, count_cover)
# set rownames
rownames(ydat_count_rare) <- ydat_count_rare$rowid
# drop anything not plant cov
ydat_count_rare <- ydat_count_rare[grepl("Nat|Exo", names(ydat_count_rare))]
# be sure roworder in rare is same as in sp-cover
summary(rownames(ydat_count) == rownames(ydat_count_rare)) # yes

# cbind then recrunch colsums to drop anything in less than 5% of plots
ydat_count <- cbind(ydat_count, ydat_count_rare)
# check colsums
checkcolSums <- sapply(ydat_count, function(x) sum(x>0))
sort(checkcolSums)/nrow(ydat_count)
# combine background native nfix and native forb, and exo forb and exo nfix (note each nfixers only present in 2 plots)
ydat_count <- ydat_count[checkcolSums > nrow(ydat_count)*.05]
ydat_count$NativeForbNfix <- ydat_count$NativeForb + ydat_count_rare$NativeNfix
ydat_count$ExoticForbNfix <- ydat_count$ExoticForb + ydat_count_rare$ExoticNfix
# drop NativeForb and ExoticForb cats
ydat_count <- ydat_count[!names(ydat_count) %in% c("ExoticForb", "NativeForb")]
# review
summary(ydat_count) # ok
# make sure no rare sp in final
summary(ydat_raresp %in% names(ydat_count)) # none

#  sp as beta distribution [0,1]
ydat_beta <- merge(widesp_pct, rowkey)
# since already have ID'd rare can combine aggregageted rare with spp first, then drop all nonspp cols
# create rare groups for beta dist
ydat_beta_rare <- subset(ydat_rare, select = -count_cover) %>%
  spread(coarse_fxnl, pct_cover) %>%
  # sum natforb and nat nfix, exoforb and exonfix
  mutate(NativeForbNfix = NativeForb + NativeNfix,
         ExoticForbNfix = ExoticForb + ExoticNfix)
ydat_beta <- cbind(ydat_beta, ydat_beta_rare)
# set rownames
rownames(ydat_beta) <- ydat_beta$rowid
# drop anything not plant cov
ydat_beta <- ydat_beta[names(ydat_beta)[grepl("^[A-Z]+|Nat|Exo", names(ydat_beta))]]

# check colsums
checkcolSums <- sapply(ydat_beta, function(x) sum(x>0))
sort(checkcolSums)/nrow(ydat_beta)
# combine background native nfix and native forb, and exo forb and exo nfix (note each nfixers only present in 2 plots)
ydat_beta <- ydat_beta[checkcolSums > nrow(ydat_beta)*.05]

# drop NativeForb and ExoticForb cats
ydat_beta <- ydat_beta[!names(ydat_beta) %in% c("ExoticForb", "NativeForb")]

# make all proportional, so between 0 and 1
ydat_beta <- mutate_all(ydat_beta, function(x) x/100)
# make sure all between 0,1
summary(ydat_beta) # all btwn 0 and 1

# specify design structure
sDesign_full <- xdat[c("plot", "block", "xpos", "ypos", "wholeplotID", "subplotID")]


# check that ydat, sDesign, and xdat are in same order
summary(rownames(ydat_count) == rownames(xdat))
head(rownames(ydat_count)) 
head(rownames(ydat_beta))
head(xdat) # rows ordered by factor levels rather than alpha order
ydat_count <- ydat_count[rownames(xdat),]
ydat_beta <- ydat_beta[rownames(xdat),]
# doublecheck
summary(rownames(ydat_beta) == rownames(xdat))
summary(rownames(ydat_beta) == rownames(sDesign_full)) # okay


# -- prep seeded plot y data -----
# > can re-use xdat and traitdat, just need to subset to relevant spp
# need to split out seeded so no columns with all 0s -- this also means re-aggregating fxnl cover for rare spp
# what is 5% of plots?
sum(xdat$seedtrt != "Unseeded") * 0.05 # at least 4 plots

# build count first; build beta only if need (will probably use count)
ydat_count_seeded <- merge(subset(widesp_count, seedtrt != "Unseeded"), rowkey) %>%
  # order in same way as xdat (not sure why got out of order)
  arrange(plot, seedtrt, herbicide)
  
# assign rownames and drop all else but species
rownames(ydat_count_seeded) <- ydat_count_seeded$rowid
ydat_count_seeded <- ydat_count_seeded[names(ydat_count_seeded)[grepl("^[A-Z]+", names(ydat_count_seeded))]]
# check colsums
checkcolSums <- sapply(ydat_count_seeded, function(x) sum(x>0))
sort(checkcolSums)/nrow(ydat_count_seeded)

# pull rare (<5% plots to group by fxnl group)
# > to not over-inflate count on rare, sum on pct cover then transform to count
ydat_raresp_seeded <- names(ydat_count_seeded)[checkcolSums < nrow(ydat_count_seeded)*.05]
# keep pct cover df separate to recycle for beta distribution
ydat_rare_seeded <- subset(tidysp0, seedtrt != "Unseeded" & code4 %in% ydat_raresp_seeded) %>%
  # join spp functional info
  left_join(distinct(natambsp_cov, code4, coarse_fxnl)) %>%
  group_by(plot, seedtrt, herbicide, coarse_fxnl) %>%
  reframe(pct_cover = sum(pct_cover)) %>%
  distinct() %>%
  # shorten fxln group names
  mutate(coarse_fxnl = gsub("Background | |-|er", "", coarse_fxnl),
         # trace gets 1, all else gets rounded; add 0.1 sound 0.5 rounds (R doesn't round it always)
         count_cover = ifelse(pct_cover == 0, 0, 
                              ifelse(pct_cover == 0.01, 1, round(pct_cover + 1.1)))) %>%
  # join rowid
  merge(rowkey)

# convert to count as in cover prep script
ydat_count_rare_seeded <- subset(ydat_rare_seeded, select = -pct_cover) %>%
  spread(coarse_fxnl, count_cover)

# join based on rowid then keep only plant cols
ydat_count_seeded <- mutate(ydat_count_seeded, rowid = rownames(ydat_count_seeded)) %>%
  merge(ydat_count_rare_seeded) %>%
  arrange(plot, seedtrt, herbicide)
# set rownames
rownames(ydat_count_seeded) <- ydat_count_seeded$rowid
# drop anything not plant
ydat_count_seeded <- ydat_count_seeded[names(ydat_count_seeded)[grepl("^[A-Z]+|Nat|Exo", names(ydat_count_seeded))]]

# check colsums
checkcolSums <- sapply(ydat_count_seeded, function(x) sum(x>0))
sort(checkcolSums)/nrow(ydat_count_seeded) 
checkcolSums # sum nfix and forbs for ambient native plants only

# combine background native nfix and native forb, and exo forb and exo nfix (note each nfixers only present in 2 plots)
ydat_count_seeded <- ydat_count_seeded[checkcolSums > nrow(ydat_count_seeded)*.05]
ydat_count_seeded$NativeForbNfix <- ydat_count_seeded$NativeForb + ydat_count_rare_seeded$NativeNfix
# drop NativeForb and ExoticForb cats
ydat_count_seeded <- ydat_count_seeded[!names(ydat_count_seeded) %in% c("NativeForb")]
# review
summary(ydat_count_seeded) # ok
# make sure no rare sp in final
summary(ydat_raresp_seeded %in% names(ydat_count_seeded)) # none

# create xdat and study design for seeded for convenience, and to be sure roworder same
xdat_seeded <- subset(xdat, seedtrt != "Unseeded")
sDesign_seeded <- sDesign_full[grepl("NatSeed", rownames(sDesign_full)),]
# check rownames the same order
summary(rownames(ydat_count_seeded) == rownames(xdat_seeded))
summary(rownames(ydat_count_seeded) == rownames(sDesign_seeded)) # yes



# -- prep trait data ----
# prep trait data so rows in same order as cols in y dat, trait vals will be in cols
# > categorical traits need to be factors; spp = rownames; remember to scale
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
traitmat_full <- traitmat[rownames(traitmat) %in% names(ydat_count),]
traitmat_seeded <- traitmat[rownames(traitmat) %in% names(ydat_count_seeded),]
# check order agrees (should be alphabetical)
rownames(traitmat_full);names(ydat_count) # yes, alphabetical
nrow(traitmat_full) # 34 rows but 37 spp in ydat_count.. may need to drop spp not in trait dat if want 4th corner model
names(ydat_count)[!names(ydat_count) %in% traitmat_full$code4] # unk aster, 2 native plnts, aggregates
names(ydat_count_seeded)[!names(ydat_count_seeded) %in% rownames(traitmat_seeded)] # ast3, napu, and aggregate groups

# note cols in ydat to drop
notrait_spp <- unique(c(names(ydat_count)[!names(ydat_count) %in% traitmat_full$code4],
                 names(ydat_count_seeded)[!names(ydat_count_seeded) %in% rownames(traitmat_seeded)]))
                 

# clean up enviro before proceed
rm(list = ls()[grepl("asla|asltat|notraits_prop|seedsum|checkcolS|plots_sum", ls())])

# normalize traits (math trans)
# > this was worked out in native_ordinations script. just copy-pasting and adapting code.

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

# check all looks good
traitmat_full_norm[,!grepl("code4|fxnl_grp|nativity|durat", names(traitmat_full_norm))] %>%
  gather(met, val) %>%
  ggplot() +
  geom_density(aes(val)) +
  theme_minimal() +
  facet_wrap(~gsub("[.]", " ", met), scales = "free", labeller = label_wrap_gen(width= 15))
# good enough! 

# specify reduced trait set to pull (also determined in native_ordinations.R based on correlations/axes) 
reduced_traits <- names(traitmat_full)[grepl("SLA|Fresh|RMF|LDM|Height|SLA|root.diam|Fine.root.specif|volu", names(traitmat_full))]
reduced_traits # be sure have the ones desired

 
# -- PcoA and db-RDA to set expectations -----
# this may produce neg eigenvalues, but to get a general idea of constrained ordination
# use hellinger transformation -- drop native seeded spp because won't be in unseeded plots
allplots_hel <- decostand(ydat_count[!names(ydat_count) %in% c("FEMI", "BRCA", "ESCA", "NEMA")], method = "hellinger")
allplots_hel_pcoa <- cmdscale(dist(allplots_hel), k=(nrow(allplots_hel)-1), eig=TRUE) # warning message is fine because only interested in first several
allplots_hel_pcoa.ape <- ape::pcoa(dist(allplots_hel))
which(allplots_hel_pcoa$eig <0) # not until 90th
View(allplots_hel_pcoa.ape$values) # first axis explains ~30% of variation in data, 2 and 3 both explain 11%
biplot(allplots_hel_pcoa.ape, Y = scale(allplots_hel)) # axes 1,2
biplot(allplots_hel_pcoa.ape, plot.axes = c(2,3)) # axes 2,3
biplot(allplots_hel_pcoa.ape, plot.axes = c(1,3))

vegan::ordiplot(allplots_hel_pcoa)
ordiellipse(allplots_hel_pcoa, xdat$fulltrt, kind = "sd", label = T)
ordiplot(allplots_hel_pcoa)
ordiellipse(allplots_hel_pcoa, xdat$nut_trt, col = c("chocolate2", "mediumpurple", "seagreen"), kind = "sd", lwd = 3, label = T)
ordiplot(allplots_hel_pcoa)
ordiellipse(allplots_hel_pcoa, xdat$ppt_trt, kind = "sd", col = c("lightblue1", "dodgerblue2", "steelblue4"), lwd = 3, label = T) # precip doesn't do much
ordiplot(allplots_hel_pcoa)
ordiellipse(allplots_hel_pcoa, xdat$fulltrt, kind = "se", col = c(paste0("mediumpurple", 1:3),
                                                                  paste0("sienna", 1:3),
                                                                  paste0("seagreen", 1:3)), lwd = 3, label = T)


ordiplot(allplots_hel_pcoa)
ordiellipse(allplots_hel_pcoa, xdat$block, kind = "sd", 
            col = c("aquamarine", "cadetblue3", "orchid", "mediumpurple4"),
            label = T) # axis 1 splits hillslope (as I've seen in nmds)
ordiplot(allplots_hel_pcoa)
ordihull(allplots_hel_pcoa, xdat$block, label = T, 
         col = c("aquamarine", "cadetblue3", "orchid", "mediumpurple4"),lwd = 2) # still nearly cleanly split capturing all composition in block
title(main = "PCoA: all sub-sub-subplots species abundances (seeded spp removed)\n(euclidean distance on hellinger transformed composition)")
         # try axes 2,3
ordiplot(allplots_hel_pcoa, choices = c(2,3))
ordiellipse(allplots_hel_pcoa, xdat$block, kind = "sd",choices = c(2,3), label = T) # this is something besides hill position.. sort of
ordiplot(allplots_hel_pcoa, choices = c(1,3))
ordiellipse(allplots_hel_pcoa, xdat$block, kind = "sd",choices = c(1,3), label = T)
ordiplot(allplots_hel_pcoa, choices = c(1,2))
ordiellipse(allplots_hel_pcoa, paste(xdat$block, xdat$seedtrt), col = c(paste0("skyblue", c(1,3)),
                                                                          paste0("orchid", c(2,4)),
                                                                          paste0("seagreen", c(1,3)),
                                                                           paste0("chocolate", c(1,3))), 
            lwd = 3, kind = "sd",choices = c(1,2), label = T)
ordiplot(allplots_hel_pcoa, choices = c(2,3))
ordiellipse(allplots_hel_pcoa, paste(xdat$block, xdat$herbicide), col = c(paste0("skyblue", c(1,3)),
                                                                        paste0("orchid", c(2,4)),
                                                                        paste0("seagreen", c(1,3)),
                                                                        paste0("chocolate", c(1,3))), 
            lwd = 3, kind = "sd",choices = c(2,3), label = T)
# dims 2 and sort of split block 3 from 4, and dim3 block 1 from 2?

# see how much wholeplotID makes a difference -- this would ignore water treatment
# plot one block at a time
par(mfrow = c(2,2))
ordiplot(allplots_hel_pcoa$points[xdat$block == 1,])
ordihull(allplots_hel_pcoa$points[xdat$block == 1,], xdat$wholeplotID[xdat$block == 1], label = T) # axis 1 splits hillslope (as I've seen in nmds)
ordiplot(allplots_hel_pcoa$points[xdat$block == 2,])
ordihull(allplots_hel_pcoa$points[xdat$block == 2,], xdat$wholeplotID[xdat$block == 2], label = T)
ordiplot(allplots_hel_pcoa$points[xdat$block == 3,])
ordihull(allplots_hel_pcoa$points[xdat$block == 3,], xdat$wholeplotID[xdat$block == 3], label = T)
ordiplot(allplots_hel_pcoa$points[xdat$block == 4,])
ordihull(allplots_hel_pcoa$points[xdat$block == 4,], xdat$wholeplotID[xdat$block == 4], label = T)

# subplot ID = block x nut trt x ppt trt
ordiplot(allplots_hel_pcoa$points[xdat$block == 1,], choices = c(2,3))
ordiellipse(allplots_hel_pcoa$points[xdat$block == 1,], choices = c(2,3), xdat$subplotID[xdat$block == 1], kind = "se", label = T) # axis 1 splits hillslope (as I've seen in nmds)
ordiplot(allplots_hel_pcoa$points[xdat$block == 2,], choices = c(2,3))
ordiellipse(allplots_hel_pcoa$points[xdat$block == 2,], choices = c(2,3), xdat$subplotID[xdat$block == 2], kind = "se", label = T)
ordiplot(allplots_hel_pcoa$points[xdat$block == 3,], choices = c(2,3))
ordiellipse(allplots_hel_pcoa$points[xdat$block == 3,], choices = c(2,3), xdat$subplotID[xdat$block == 3], kind = "se", label = T)
ordiplot(allplots_hel_pcoa$points[xdat$block == 4,], choices = c(2,3))
ordiellipse(allplots_hel_pcoa$points[xdat$block == 4,],choices = c(2,3), xdat$subplotID[xdat$block == 4], kind = "se", label = T)
# ^looking at relative position of different full trtments across blocks, it seems like compost is having slightly different effect upslope than down
# downslope (probably wetter?), cd and cw overlap, separate from others; upslope (probably drier?), cxc and cw overlap, separate from others 
# everything else is pretty jumbled

# dims 2 and 3 seems to be separating environmental effects by hillslope pos. different things going on. blocks 1 and 2 are slightly more similar than 3 and 4 to one another. 

# constrain by block
blockRDA <- rda(allplots_hel ~ block, data=xdat)
summary(blockRDA)
par(mfrow = c(1,1))
screeplot(blockRDA)
ordiplot(blockRDA, display = "all")

# constrain by env
envRDA <- rda(allplots_hel ~ fulltrt, data=xdat)
summary(envRDA)
screeplot(envRDA)
ordiplot(envRDA, display = "all")

# constrain by local site (block), herbicide and full enviro trt
fullRDA <- rda(allplots_hel ~ block + herbicide + fulltrt, data=xdat)
summary(fullRDA)
screeplot(fullRDA) # local composition [block] has much bigger influence than environment
ordiplot(fullRDA, display = "all", scaling = 2, choices = c(1,2))
ordiplot(fullRDA, display = "all", scaling = 2, choices = c(2,3))
ordiplot(fullRDA, display = "all", scaling = 2, choices = c(2,3))
ordiellipse(fullRDA, xdat$ppt_trt, kind = "se", choices = c(2,3), label = T)
ordiellipse(fullRDA, xdat$nut_trt, kind = "se", col = "orange", lwd = 3, choices = c(2,3), label = T)
spscores <- scores(fullRDA)$species
ordiplot(fullRDA, display = "all", scaling = 2, choices = c(1,2))
ordiellipse(fullRDA, xdat$ppt_trt, kind = "se", col = "dodgerblue", choices = c(2,3), label = T)
ordiellipse(fullRDA, xdat$nut_trt, kind = "se", col = "orange", lwd = 3, choices = c(1,2), label = T)
text(spscores, rownames(spscores))


# quick test on seeded plots
seedplots_hel <- decostand(ydat_count_seeded, method = "hellinger")
seedplots_hel_pcoa <- cmdscale(dist(seedplots_hel), k=(nrow(seedplots_hel)-1), eig=TRUE) # warning message is fine because only interested in first several
seedplots_hel_pcoa.ape <- ape::pcoa(dist(seedplots_hel))

biplot(seedplots_hel_pcoa.ape, Y = ydat_count_seeded) # femi came up in lower non grassy, herbicided
ordiplot(seedplots_hel_pcoa)
ordiellipse(seedplots_hel_pcoa, xdat_seeded$wholeplotID, kind = "sd", label = T)
ordiellipse(seedplots_hel_pcoa, xdat_seeded$herbicide, kind = "sd", col = c("pink1", "pink4"), label = T)
ordiellipse(seedplots_hel_pcoa, paste0(xdat_seeded$wholeplotID, xdat_seeded$herbicide), kind = "ehull", col = c("pink1", "pink4"), label = T)


# -- PRELIM GLLVM MODELS -----
# to determine distribution to use (negbin or beta)
# specify seed for reproducibility
myseed <- 1371
myseed2 <- 81615

# unconstrained ordination model, without random error structure so can plot resids (doesn't work with )
null_nb <- gllvm(y = ydat_count, family = "negative.binomial", seed = myseed)
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(null_nb) # looks pretty good visually
summary(null_nb, dispersion = T, Lvcoefs = T) # AIC:  14377.27 AICc:  14385.86 BIC:  15341.55 LL:  -7042 df:  147 

# still unconstrained model, but with nested design specified in random errom
null_nb_nested <- gllvm(y = ydat_count, family = "negative.binomial", 
                        studyDesign = sDesign_full,
                        row.eff = ~(1|block/wholeplotID/subplotID),
                        seed = myseed
)
summary(null_nb_nested, row.intercepts = T, dispersion = T, spp.intercepts = T, Lvcoefs = T) # AIC:  16237.27 AICc:  16246.22 BIC:  17221.22 LL:  -7969 df:  150 
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(null_nb_nested) # cannot plot when nested error specified (maybe a bug in development version, works in CRAN)

# null model on proportional cover
null_beta <- gllvm(y = ydat_beta, family = "beta",
                   seed = myseed)
plot(null_beta) # bad fits due to high number 0s [maybe because of seeded native spp], stick with negbin

# try with seeded plots only to see if they are better fits
null_beta_seeded <- gllvm(y = ydat_beta_seeded, family = "beta",
                          seed = myseed)
plot(null_beta_seeded) # still not good fits because of high number of 0s, stick with negbin


# > move forward with negative binomial gllvm

# -- screen best number latent variables ----
# above pcoa and rda suggest maybe 3, but 1 axis (diving hillslope) captures the most variation in data 
bestlv <- function(y = ydat_count, X = NULL, numseq = 1:4, yform = as.formula(~ herbicide + seedtrt), smat = sDesign_full, rowform = as.formula(~ (1|block/wholeplotID)), seed = myseed){
  criteria <- NULL
  for(i in numseq){
    print(i) # print to know where you are in loop while waiting
    fiti <- gllvm(y, X, family = "negative.binomial", num.lv = i, sd.errors = FALSE,
                  formula = yform, studyDesign = smat, row.eff = rowform, seed = myseed)
    criteria[i] <- summary(fiti)$AICc
    names(criteria)[i] = i
  }
  return(criteria)
}

bestnull_nb <- bestlv(yform = NULL)
bestnull_nb # 2 latent variables are good! then 1
#      1        2        3        4 
# 14704.41 14693.05 14716.80 14771.02 
# 14823.09 14809.80 14824.02 14878.27  <-- different numbers when set seed(?), still same conclusion: 2 then 1/3 are comp

# test on seeding and herbicide treats (ignore env for now): are 2 lv's still best?
bestseedherb_nb <- bestlv(ydat_count, X = xdat)
bestseedherb_nb 
#       1        2        3        4 
# 13960.32 13988.62 14050.48 14126.30 # 1 LV is best

# run with a different seed to see if 1 LV still best
bestseedherb_nb <- bestlv(ydat_count, X = xdat, rowform = as.formula(~ 1|block/wholeplotID/subplotID), seed = myseed)
bestseedherb_nb 
# 1 still best.. maybe not as much randomness involved
# 1        2        3        4 
# 13967.41 19001.61 14060.54 29141.78  # 1 still best

# run with enviro terms added
bestseedherbenv_nb <- bestlv(ydat_count, X = xdat, yform = as.formula(~ herbicide + seedtrt + nut_trt*ppt_trt))
bestseedherbenv_nb # 1 LV is also best
#     1        2        3        4 
# 13888.62 13959.78 14046.78 14282.83 


# -- GLLVM NEGBIN MODELS -----
# add hillpos to xdat
xdat <- mutate(xdat, hillpos = ifelse(block %in% c(1,2), "down", "up"),
               # alpha order is correct factor order (downslope = ref, just as block 1 is ref)
               hillpos = factor(hillpos))

# without environment, just nested random structure with herbicide and seedtrt
seedherb_nb <- gllvm(y = ydat_count[!names(ydat_count) %in% nats], X = xdat, 
                     family = "negative.binomial", 
                     #num.lv = 3, 
                     formula = as.formula(~ herbicide + seedtrt), # if have seeedtrt and 3 latent vars, have zero in covar matrix
                     # nest subplot if can, but block/wholeplot captures most of the nested similarity
                     studyDesign = sDesign_full[c("block", "wholeplotID", "subplotID")], 
                     row.eff = ~ (1|block/wholeplotID), # with subplot throws singular fit error, also throws error with wholeplotID nested [if have 3 lvs]
                     seed = myseed)
summary(seedherb_nb)
gllvm::coefplot(seedherb_nb, order = F)
gllvm::ordiplot(seedherb_nb) # split by hillpos

# plot residual spp correlations
corrplot::corrplot(getResidualCor(seedherb_nb), method = "square", 
                   main = "spp residual correlation: herbcide + seedtrt, 1|block/wp",
                   mar = c(0,0,2,0),
                   order = "hclust", addrect = 4
)
# presidual covariance
corrplot::corrplot(getResidualCov(seedherb_nb)$cov, is.corr = F, method = "square", 
                   main = "spp residual covariance: herbcide + seedtrt, 1|block/wp",
                   mar = c(0,0,2,0),
                   order = "hclust", addrect = 4
)



# add in env vars as predictors
seedherbenv_nb <- gllvm(y = ydat_count[!names(ydat_count) %in% nats], X = xdat, family = "negative.binomial",
                        formula = ~ herbicide + seedtrt + nut_trt * ppt_trt, 
                        # lv.formula = ~ ppt_trt*nut_trt,
                        # num.lv.c = 2, 
                        studyDesign = sDesign_full[c("block", "wholeplotID", "subplotID")],
                        row.eff = ~ (1|block), #(1|block/wholeplotID),
                        seed = myseed)

summary(seedherbenv_nb)
coefplot(seedherbenv_nb,order = F, mfrow = c(2,5), cex.ylab = 1) #xlim.list = rep(list(c(-200,200), c(-200,200)),10))
par(mfrow = c(1,1))
gllvm::ordiplot(seedherbenv_nb, s.cex = 0.5, biplot = T, arrow.scale = 1, symbols = T)

# add env as constraints
seedherb_hillenvls_nb <- gllvm(y = ydat_count[!names(ydat_count) %in% nats], X = xdat, family = "negative.binomial",
                        # put block since captures so much variation
                        formula = ~ herbicide + seedtrt , 
                        lv.formula = ~ hillpos * ppt_trt*nut_trt,
                        num.lv.c = 2, 
                        studyDesign = sDesign_full[c("block", "wholeplotID", "subplotID")],
                        row.eff = ~ (1|block), #(1|block/wholeplotID),
                        seed = myseed)

summary(seedherb_hillenvls_nb)
coefplot(seedherb_hillenvls_nb, mfrow = c(4,5), order = F, cex.ylab = 0.5, which.Xcoef = c("ppt_trtD:nut_trtC"))
coefplot(seedherb_hillenvls_nb, mfrow = c(1,1), order = F, cex.ylab = 0.5, which.Xcoef = c("ppt_trtD:nut_trtC")) #xlim.list = rep(list(c(-200,200), c(-200,200)),10))
par(mfrow = c(1,1))
gllvm::ordiplot.gllvm(seedherb_hillenvls_nb, s.cex = 0.5,  biplot = T, alpha = 1, arrow.scale = 1, symbols = T, spp.arrows = T)
gllvm::ordiplot(seedherb_hillenvls_nb, s.cex = 0.5, arrow.scale = 1, symbols = F, type = "residual")

# add native spp seeded in
seedherb_hillenvls_nb_wnats <- gllvm(y = ydat_count, X = xdat, family = "negative.binomial",
                               # put block since captures so much variation
                               formula = ~ herbicide + seedtrt , 
                               lv.formula = ~ hillpos * ppt_trt*nut_trt,
                               num.lv.c = 2, 
                               studyDesign = sDesign_full[c("block", "wholeplotID", "subplotID")],
                               row.eff = ~ (1|block), #(1|block/wholeplotID),
                               seed = myseed)

summary(seedherb_hillenvls_nb_wnats)
coefplot(seedherb_hillenvls_nb_wnats, mfrow = c(4,5), order = F, cex.ylab = 0.5) #xlim.list = rep(list(c(-200,200), c(-200,200)),10))
coefplot(seedherb_hillenvls_nb_wnats, mfrow = c(1,1), order = F, cex.ylab = 0.5, which.Xcoef = c("hillposup:ppt_trtD:nut_trtC"))
par(mfrow = c(1,1))
gllvm::ordiplot.gllvm(seedherb_hillenvls_nb, s.cex = 0.5,  biplot = T, alpha = 1, arrow.scale = 1, symbols = T, spp.arrows = T)
gllvm::ordiplot(seedherb_hillenvls_nb, s.cex = 0.5, arrow.scale = 1, symbols = F, type = "residual")


seedherbenv_nb_latblock <- gllvm(y = ydat_count, X = xdat, family = "negative.binomial",
                        formula = ~ herbicide + seedtrt + nut_trt, 
                        lv.formula = ~ ppt_trt,
                        num.lv.c = 2, 
                        studyDesign = sDesign_full[c("block", "wholeplotID", "subplotID")],
                        row.eff = ~ (1|block), #(1|block/wholeplotID),
                        seed = myseed)
ordiplot(seedherbenv_nb_latblock, s.cex = 0.75, jitter = T, biplot = T)
# ^ when use more latent variables the se errors are narrower (e.g., 1 vs. 3)
par(mfrow = c(1,1))
corrplot::corrplot(getResidualCor(seedherbenv_nb), method = "square", order = "hclust", addrect = 4)
corrplot::corrplot(getResidualCov(seedherbenv_nb)$cov,is.corr = F, method = "square", order = "hclust", addrect = 4)
ordiplot(seedherbenv_nb, s.cex = 0.75, biplot = T) # LV 1 is capturing hill position, whether use 1, 2 or 3 LVs

# guidance from ben bolker's glmm FAQ:
# "Treating factors with small numbers of levels as random will in the best case lead to very small and/or imprecise estimates of random effects;" 
# "in the worst case it will lead to various numerical difficulties such as lack of convergence, zero variance estimates, etc.."

seedherbenv_nb_blockfixed <- gllvm(y = ydat_count, X = xdat, num.lv = 2, family = "negative.binomial",
                                   formula = ~ herbicide + seedtrt + ppt_trt * nut_trt, 
                                   # keep wholeplotID as random. has 12 levels  
                                   studyDesign = sDesign[c("wholeplotID", "subplotID")],
                                   row.eff = ~(1|block/wholeplotID),
                                   seed = myseed)

summary(seedherbenv_nb_blockfixed)
plot(seedherbenv_nb_blockfixed)
rownames(sDesign_full)

# split models by hillpos to see if will run, env effect estimats are different
lowhill_plots <- unique(xdat$plot[xdat$block %in% c(1,2)])
lowhill_plots_noseed <- rownames(xdat)[xdat$block %in% c(1,2) & xdat$seedtrt == "Unseeded"]
uphill_plots <- unique(xdat$plot[xdat$block %in% c(3,4)])
parse_number(rownames(ydat_count)) %in% lowhill_plots
View(ydat_count[parse_number(rownames(ydat_count)) %in% lowhill_plots,])

# need to remove ydat with all 0s in low hill pos
lowhill_plotpres <- sapply(ydat_count[parse_number(rownames(ydat_count)) %in% lowhill_plots,], function(x) sum(x > 0) > 4)
lowhill_noseed_plotpres <- sapply(ydat_count[rownames(ydat_count) %in% lowhill_plots_noseed,], function(x) sum(x > 0) > 4)
lowhill_plotpres
lowhill_noseed_plotpres

uphill_plotpres <- sapply(ydat_count[parse_number(rownames(ydat_count)) %in% uphill_plots,], function(x) sum(x > 0) > 4)
uphill_plotpres

seedherbenv_nb_lowhill <- gllvm(y = ydat_count[parse_number(rownames(ydat_count)) %in% lowhill_plots,lowhill_plotpres], 
                                X = xdat[parse_number(rownames(ydat_count)) %in% lowhill_plots,], #num.lv = 2, 
                                family = "negative.binomial",
                                formula = ~ herbicide + seedtrt, 
                                lv.formula = ~ ppt_trt*nut_trt,
                                num.lv.c = 2,
                                # keep wholeplotID as random. has 12 levels  
                                #studyDesign = sDesign_full[parse_number(rownames(ydat_count)) %in% lowhill_plots, c("block", "wholeplotID", "subplotID")],
                                #row.eff = ~(1|block), #~(1|block/wholeplotID),
                                seed = myseed)

coefplot(seedherbenv_nb_lowhill, mfrow = c(4,3), cex.ylab = 1)
ordiplot(seedherbenv_nb_lowhill, type = "conditional", s.cex = 0.75, jitter = T, jitter.amount = 2, rotate = T)
# unseeded plots only
seedherbenv_nb_lowhill_noseed <- gllvm(y = ydat_count[rownames(ydat_count) %in% lowhill_plots_noseed,lowhill_noseed_plotpres], 
                                X = xdat[rownames(ydat_count) %in% lowhill_plots_noseed,], #num.lv = 2, 
                                family = "negative.binomial",
                                formula = ~ herbicide,
                                lv.formula = ~ ppt_trt * nut_trt,
                                num.lv. = 2,
                                # keep wholeplotID as random. has 12 levels  
                                #studyDesign = sDesign_full[parse_number(rownames(ydat_count)) %in% lowhill_plots, c("block", "wholeplotID", "subplotID")],
                                #row.eff = ~(1|block), #~(1|block/wholeplotID),
                                seed = myseed)

ordiplot(seedherbenv_nb_lowhill_noseed, 
         jitter = T, jitter.amount = 0.5, 
         s.cex = 0.75, biplot = T, type = "conditional", rotate = T) # lol


seedherbenv_nb_uphill <- gllvm(y = ydat_count[parse_number(rownames(ydat_count)) %in% uphill_plots,uphill_plotpres], 
                                X = xdat[parse_number(rownames(ydat_count)) %in% uphill_plots,], #num.lv = 2, 
                                family = "negative.binomial",
                                formula = ~ herbicide + seedtrt + ppt_trt * nut_trt, 
                                # keep wholeplotID as random. has 12 levels  
                                #studyDesign = sDesign_full[parse_number(rownames(ydat_count)) %in% lowhill_plots, c("block", "wholeplotID", "subplotID")],
                                #row.eff = ~(1|block), #~(1|block/wholeplotID),
                                seed = myseed)
summary(seedherbenv_nb_uphill)
shortnames <- dimnames(seedherbenv_nb_lowhill$X.design)[[2]]
coefplot(seedherbenv_nb_lowhill, mfrow = c(3,2), cex.ylab = 0.75, which.Xcoef = c(shortnames[!grepl(":", shortnames)]))
coefplot(seedherbenv_nb_uphill, mfrow = c(3,2), cex.ylab = 0.75, which.Xcoef = c(shortnames[!grepl(":", shortnames)]))
seedherbenv_nb_lowhill$params$sigma.lv

ordiplot(seedherbenv_nb_lowhill, biplot = T, spp.arrows = T, s.cex = 0.75)

sppeffects_lowsd <- seedherbenv_nb_lowhill$sd$Xcoef

# -- test role of traits ----
# scale numeric cols in traitmat
traitmat_scaled <- mutate_at(traitmat_full_norm, .vars = names(traitmat_full_norm)[!grepl("code4|durat|fxnl_|nativ", names(traitmat_full_norm))], function(x) as.numeric(scale(x)))

# just env vars, no random effects
seedherb_nb_traits <- gllvm(y = ydat_count[,names(ydat_count) %in% rownames(traitmat)], X = xdat[c("seedtrt", "herbicide", "nut_trt", "ppt_trt")],
                            TR = traitmat_scaled[rownames(traitmat_scaled) %in% names(ydat_count),reduced_traits],#,c("Fine.root.specific.length.cm.g", "Proportion.fine.roots", "RMF", "LDMC", "SLA.cm2.g", "Height.cm")], # drop code4 
                            family = "negative.binomial", #num.lv = 1, 
                            #formula = y ~ (herbicide + seedtrt + nut_trt * ppt_trt) + (Coarse.root.diameter.mm + Fine.root.specific.length.cm.g +Height.cm + LDMC),
                            #studyDesign = sDesign_full[c("block", "wholeplotID", "subplotID")], #row.eff = ~ (1|block/wholeplotID/subplotID),
                            seed = myseed)
summary(seedherb_nb_traits)
seed
coefplot(seedherb_nb_traits)
gllvm::ordiplot(seedherb_nb_traits, biplot = T)

# with traits, no random effect
seedherb_nb_scaletraits <- gllvm(y = ydat_count[,names(ydat_count) %in% rownames(traitmat)], 
                                 X = xdat[c("seedtrt", "herbicide", "nut_trt", "ppt_trt")],
                                 TR = traitmat_scaled[rownames(traitmat_scaled) %in% names(ydat_count),!grepl("code4|biomass.g", names(traitmat_scaled))],#,c("Fine.root.specific.length.cm.g", "Proportion.fine.roots", "RMF", "LDMC", "SLA.cm2.g", "Height.cm")], # drop code4 
                                 family = "negative.binomial", #num.lv = 1, 
                                 formula = y ~ (herbicide + seedtrt + nut_trt * ppt_trt):(Coarse.root.diameter.mm + Coarse.root.length.mm  + Fine.root.specific.length.cm.g + Height.cm + LDMC + SLA.cm2.g),
                                 #studyDesign = sDesign_full[c("block", "wholeplotID", "subplotID")], row.eff = ~ (1|block/wholeplotID),
                                 seed = myseed)

seedherb_nb_scaletraits_re <- gllvm(y = ydat_count[,names(ydat_count) %in% rownames(traitmat)], X = xdat[c("seedtrt", "herbicide", "nut_trt", "ppt_trt")],
                                 TR = traitmat_scaled[rownames(traitmat_scaled) %in% names(ydat_count),!grepl("code4|biomass.g", names(traitmat_scaled))],#,c("Fine.root.specific.length.cm.g", "Proportion.fine.roots", "RMF", "LDMC", "SLA.cm2.g", "Height.cm")], # drop code4 
                                 family = "negative.binomial", #num.lv = 1, 
                                 #formula = y ~ (herbicide + seedtrt + nut_trt * ppt_trt)*(Coarse.root.diameter.mm + Coarse.root.length.mm  + Fine.root.specific.length.cm.g + Height.cm + LDMC + SLA.cm2.g),
                                 studyDesign = sDesign_full[c("block", "wholeplotID", "subplotID")], row.eff = ~ (1|block/wholeplotID),
                                 seed = myseed)


# with unique wholeplot ID nested within block, all variables (env and traits) with all interactions
seedherb_nb_allscaletraits_re <- gllvm(y = ydat_count[,names(ydat_count) %in% rownames(traitmat)], 
                                    X = xdat[c("seedtrt", "herbicide", "nut_trt", "ppt_trt", "fulltrt")],
                                    TR = traitmat_scaled[rownames(traitmat_scaled) %in% names(ydat_count),!grepl("code4|biomass.g", names(traitmat_scaled))],#,c("Fine.root.specific.length.cm.g", "Proportion.fine.roots", "RMF", "LDMC", "SLA.cm2.g", "Height.cm")], # drop code4 
                                    family = "negative.binomial", #num.lv = 1, 
                                    #formula = y ~ (herbicide + seedtrt + nut_trt * ppt_trt) + (Coarse.root.diameter.mm + Fine.root.specific.length.cm.g +Height.cm + LDMC),
                                    #studyDesign = sDesign_full[c("block", "wholeplotID", "subplotID")], 
                                    #row.eff = ~ (1|block/wholeplotID),
                                    seed = myseed)

seedherb_nb_scaletraits_re$fourth.corner

coeffs_envxtrait_re <- names(data.frame(seedherb_nb_scaletraits_re$X.design))
coefplot(seedherb_nb_scaletraits,mar = c(4,10,1,1), order = T, cex.ylab = 0.75)
coefplot.gllvm(seedherb_nb_scaletraits_re, mar = c(4,17,4,1), order = T, cex.ylab = 0.75,
         cex.lab = 1, cex.main = 1, cex.axis = 1,
         main = "All env, trait, and interactions, with random effect (1|block/wholeplotID)")

lattice::levelplot(seedherb_nb_scaletraits$fourth.corner)
lattice::levelplot(seedherb_nb_scaletraits_re$fourth.corner)
dev.off()
corrplot::corrplot(seedherb_nb_scaletraits$fourth.corner, is.corr = F, method = "square")
corrplot::corrplot(seedherb_nb_scaletraits_re$fourth.corner, is.corr = F, method = "square")
ordiplot(seedherb_nb_scaletraits_re, biplot = T, symbols = T)
ordiplot(seedherb_nb_scaletraits_re,s.cex = 0.75)
ordiplot(seedherb_nb_scaletraits_re, symbols = T, biplot = T, spp.arrows = T, cex.spp = 0.75)



# -- test seeded plots only ----
seedsums <- colSums(ydat_count_seeded)
ydat_count_seeded2 <- ydat_count_seeded[seedsums > 10]
seedherb_nb_seeded_scaletraits <- gllvm(y = ydat_count_seeded2[,names(ydat_count_seeded2) %in% rownames(traitmat_scaled)], 
                                        X = xdat[grepl("NatSeed", rownames(xdat)), c("herbicide", "nut_trt", "ppt_trt")],
                                        TR = traitmat_scaled[rownames(traitmat_scaled) %in% names(ydat_count_seeded2),c(reduced_traits)],#,c("Fine.root.specific.length.cm.g", "Proportion.fine.roots", "RMF", "LDMC", "SLA.cm2.g", "Height.cm")], # drop code4 
                                        family = "negative.binomial", num.lv = 2, 
                                        formula = y ~ herbicide + ((nut_trt * ppt_trt)*(Coarse.root.diameter.mm + Fine.root.specific.length.cm.g + Root.volume.cm3)),
                                        #lv.formula = ~ block,
                                        studyDesign = data.frame(block = sDesign_full$block[grepl("NatSeed", rownames(sDesign_full))]), 
                                        corWithin = T,
                                        #randomX = ~ herbicide,
                                        row.eff = ~ (1|block),
                                        gradient.check = T,
                                        seed = myseed)
plot(seedherb_nb_seeded_scaletraits)
corrplot::corrplot(seedherb_nb_seeded_scaletraits$fourth.corner, is.corr = T, method = "square")
coefplot(seedherb_nb_seeded_scaletraits, mar = c(4,10,1,1))
confint_df <- data.frame(gllvm::confint.gllvm(seedherb_nb_seeded_scaletraits)) %>% rownames_to_column("term")
gllvm::AICc(seedherb_nb_seeded_scaletraits)
coeff_df <- data.frame(beta  = coefficients(seedherb_nb_seeded_scaletraits)$B) %>% rownames_to_column("term")
alltog_df <- left_join(coeff_df, confint_df)
alltog_df$keep <- with(alltog_df, (X2.5.. <0 & X97.5.. < 0) | (X2.5.. > 0 & X97.5.. > 0))
data.frame(seedherb_nb_seeded_scaletraits$params$theta) %>%
  rownames_to_column("spp") %>%
  ggplot(aes(LV1, LV2)) +
  geom_text(aes(label = spp))

data.frame(seedherb_nb_seeded_scaletraits$lvs) %>%
  rownames_to_column("site") %>%
  left_join(rowkey, by = c("site" = "rowid")) %>%
  #mutate(site = gsub("_NatSeed", "", site)) %>%
  ggplot(aes(LV1, LV2)) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_vline(aes(xintercept = 0), lty = 2) +
  #geom_text(aes(label = site))
  geom_point(aes(fill = ppt_trt, shape = herbicide), size = 3) +
  ggrepel::geom_text_repel(aes(label = block)) +
  scale_fill_brewer(name = "ppt_trt", palette = "Blues") +
  scale_shape_manual(values = c(21,22)) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  facet_wrap(~ paste(nut_trt, herbicide))


ggplot(alltog_df, aes(term, beta)) +
  geom_errorbar(aes(ymin = X2.5.., ymax = X97.5..), width = 0) +
  geom_point() +
  coord_flip()

gllvm::AICc(seedherb_envxnest_nb_seeded)
seedherb_envxnest_nb_seeded <- gllvm(y = ydat_count_seeded2[,names(ydat_count_seeded2) %in% rownames(traitmat_scaled)], 
                                     X = xdat[grepl("NatSeed", rownames(xdat)), ], family = "negative.binomial",
                                     formula = ~ herbicide + ppt_trt * nut_trt, 
                                     studyDesign = sDesign_full[grepl("NatSeed", rownames(sDesign_full)), ],
                                     #lv.formula = ~ block,
                                     corWithin = T,
                                     row.eff = ~(1|block),
                                     seed = myseed)
summary(seedherb_envxnest_nb_seeded)
gllvm::ordiplot(seedherb_envxnest_nb_seeded, biplot = T, s.cex = 0.75)
corrplot::corrplot(getResidualCor(seedherb_envxnest_nb_seeded), method = "square", 
                   main = c("Residual species correlations"),
                   mar = c(0,0, 5,0),
                   order = "hclust", addrect = 4
)
corrplot::corrplot(getResidualCov(seedherb_envxnest_nb_seeded)$cov, is.corr = F, method = "square", 
                   main = c("Resid spp covar; model = nested random error: block/wholeplotID"),
                   mar = c(0,0, 5,0),
                   order = "hclust", addrect = 4
)
coefplot(seedherb_envxnest_nb_seeded, cex.ylab = 1, mfrow = c(2,5))

# try running full nested model since removed rare-in-seeded spp
# seedherb_envxfullnest_nb_seeded <- gllvm(y = ydat_count_seeded[seedsums>10], X = xdat[grepl("NatSeed", rownames(xdat)), ], family = "negative.binomial",
#                                      formula = ~ herbicide + ppt_trt * nut_trt, studyDesign = sDesign[grepl("NatSeed", rownames(sDesign)), ],
#                                      row.eff = ~(1|block/wholeplotID/subplotID),
#                                      seed = myseed) # singular fit

# if you include block as fixed effect, how does that affect resid covar?
seedherb_envxblock_nb_seeded <- gllvm(y = ydat_count_seeded[seedsums>10], X = xdat[grepl("NatSeed", rownames(xdat)), ], family = "negative.binomial",
                                      formula = ~ block + herbicide + ppt_trt * nut_trt, #studyDesign = data.frame(wholeplotID = sDesign[grepl("NatSeed", rownames(sDesign)),"wholeplotID"]),
                                      #row.eff = ~(1|wholeplotID),
                                      seed = myseed)
plot(seedherb_envxblock_nb_seeded)
summary(seedherb_envxblock_nb_seeded)
coefplot(seedherb_envxblock_nb_seeded, cex.ylab = 1.5, cex.xlab = 1.5, mfrow = c(2,6))
coefplot(seedherb_envxblock_nb_seeded, cex.ylab = 1.5, cex.xlab = 1.5,
         which.Xcoef = c("herbicideHerbicided", "ppt_trtD", "ppt_trtW", "nut_trtF", "nut_trtC"))
coefplot(seedherb_envxblock_nb_seeded, cex.ylab = 1.5, cex.xlab = 1.5,
         which.Xcoef = c( "ppt_trtD:nut_trtF","ppt_trtD:nut_trtC","ppt_trtW:nut_trtF", "ppt_trtW:nut_trtC"))
corrplot::corrplot(getResidualCor(seedherb_envxblock_nb_seeded), is.corr = F, 
                   method = "square", 
                   main = c("Resid spp cor; model = block as fixed, no rand error"),
                   mar = c(0,0, 5,0)
                   #order = "hclust", addrect = 4
)
corrplot::corrplot(getResidualCov(seedherb_envxblock_nb_seeded)$cov, is.corr = F, 
                   method = "square", 
                   main = c("Resid spp covar; model = block as fixed, no rand error"),
                   mar = c(0,0, 5,0)
                   #order = "hclust", addrect = 4
)


ordiplot.gllvm(seedherb_envxblock_nb_seeded)


# -- try constrained ordination model ----
# use herbicide as latent
# -- package example for traits ----
data(antTraits)
View(antTraits$traits)
View(antTraits$abund)

# -- dev package example for nested random effects ----
data(spider)
str(spider$abund)
View(spider$trait)

# nested
model <- gllvm(spider$abund, num.lv = 2, family="poisson", #
               row.eff= ~(1|group/plot), studyDesign = data.frame(group=rep(1:2,14), plot = rep(1:4, each = 7)))
summary(model)
gllvm::getResidualCor(model)
