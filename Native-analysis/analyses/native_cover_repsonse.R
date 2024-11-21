# chapter veg cover analyses
# author: ctw

# 'gllvm' package has been majorly updated since i ran analysis for it in jan 2024 (update was in sep 2024)
# this script is to test new features that may be more appropriate for SFREC data, and to concenrate analyses for tables and figs
# foundation for everything is in R script native_gllvm.R


# what this script does:
# 1) read in prepped data, review frequency and remove/lump rare
# 2) streamlines mixed model and gllvm analyses for ctw's chapter results

# legwork for all work in:
# native_mixedmodels.R
# native_gllvm.R
# native_ordinations.R



# -- SETUP -----
source("Native-analysis/native_prep_traits.R", print.eval = F, echo = F)
# ^note: warning message about NA coercion is okay. just converting notes about no sample avail

# load needed libraries
library(glmmTMB)
library(DHARMa)
library(performance)
library(broom.mixed)
library(bbmle) # for AICtab to compare models
library(gllvm)
library(cowplot)
library(modelsummary)
theme_ctw <- function(){
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 12))
}


# specify seed for reproducibility
myseed <- 1371
myseed2 <- 81615

# specify plotting colors for spp
plantcols <- c("orchid4", "seagreen4", "chocolate3", "orchid1", "seagreen1", "goldenrod")
fxnltargcols <- c("orchid4", "seagreen4", "chocolate3", "seagreen1", "orchid1", "chartreuse3", "purple2")

# 1. PREP Y, X, AND TR DATA  ----
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
# create hillpos/grazing legacy
xdat$hillpos <- with(xdat, ifelse(block %in% c(1,2), "downhill", "uphill")) # block 1, 2 = downhill (overgrazed)
xdat$hillpos <- factor(xdat$hillpos)


# already determined using negative binomial distribution for models. use 'count' form of veg data.
# species in count form
ydat_count <- merge(widesp_count, rowkey) # merge to be sure

# assign rownames and drop all else but species
rownames(ydat_count) <- ydat_count$rowid
ydat_count <- ydat_count[names(ydat_count)[grepl("^[A-Z]+", names(ydat_count))]]

ydat_summary <- group_by(tidysp0, plot, seedtrt, herbicide) %>%
  mutate(totcov = sum(count_cover)) %>%
  ungroup() %>%
  group_by(code4) %>%
  reframe(n_all = length(count_cover),
          n_seeded = length(count_cover[!grepl("Unse", seedtrt)]),
          seeded_pres = sum(count_cover > 0 & !grepl("Unse", seedtrt)),
          all_pres = sum(count_cover >0),
          prop_allpres = all_pres/n_all,
          prop_seedpres = seeded_pres/n_seeded,
          maxcov = max(count_cover),
          maxcov_seeded = max(count_cover[!grepl("Unse", seedtrt)]),
          meancov = mean(count_cover),
          sdcov = sd(count_cover),
          uniquecov_all = length(unique(count_cover[count_cover > 0])),
          uniquecov_seeded = length(unique(count_cover[count_cover > 0 & !grepl("Unse", seedtrt)]))) %>%
  mutate(has_traits = code4 %in% natamb_avgtraits_wide$code4) %>%
  left_join(distinct(spplist[c("code4", "unknown", "fxnl_grp", "nativity")])) %>%
  mutate(grp = ifelse(grepl("Unk", nativity), paste("Exotic", fxnl_grp),
                      paste(nativity, fxnl_grp)))

# review
plot_grid(
  ggplot(ydat_summary, aes(prop_allpres, uniquecov_all, col = grp)) +
    geom_vline(aes(xintercept = 0.1)) +
    geom_hline(aes(yintercept = 3)) +
    geom_jitter(aes(size = maxcov, shape = code4 %in% nats), width = 0.02, height = 0.1) +
    ggtitle("species distribution in all plots (n = 141)") +
    scale_color_manual(values = plantcols) +
    scale_shape_manual(name = "Seeded", values = c(1,19)) +
    scale_size_continuous(name = "Max cover", breaks = c(0, 10, 25,50, 75)) +
    guides(color = guide_legend(override.aes = list(size = 3)),
           shape = guide_legend(override.aes = list(size = 3))) +
    scale_y_continuous(limit = c(-.1,max(ydat_summary$uniquecov_all)+2)) +
    #theme_bw() +
    theme(legend.position = "none") +
    theme_ctw(),
  
  ggplot(ydat_summary, aes(prop_seedpres, uniquecov_seeded, col = grp)) +
    geom_vline(aes(xintercept = 0.1)) +
    geom_hline(aes(yintercept = 3)) +
    geom_jitter(aes(size = maxcov_seeded, shape = code4 %in% nats), width = 0.02, height = 0.1) +
    ggtitle("species distribution in seeded plots (n = 70)") +
    scale_color_manual(values = plantcols) +
    scale_shape_manual(name = "Seeded", values = c(1,19)) +
    scale_size_continuous(name = "Max cover", breaks = c(0, 10, 25,50, 75)) +
    guides(color = guide_legend(override.aes = list(size = 3)),
           shape = guide_legend(override.aes = list(size = 3))) +
    scale_y_continuous(limits = c(-.1,max(ydat_summary$uniquecov_all)+2)) +
    #theme_bw() +
    theme_ctw(),
  rel_widths = c(1,1.4)
)

# check presence by group
with(ydat_summary, lapply(split(prop_allpres, grp), function(x) table(x > 0.05)))
with(ydat_summary, lapply(split(prop_seedpres, grp), function(x) table(x > 0.05)))

# check colsums
checkcolSums <- sapply(ydat_count, function(x) sum(x>0))
sort(checkcolSums)/nrow(ydat_count)
# crosscheck presence in trait dataset
notraits <- names(ydat_count)[!names(ydat_count) %in% natamb_avgtraits_wide$code4] # 22 species not in traitmat
notraits
# what is their occurrence?
sort(checkcolSums[notraits]) # only 4 would be retained: UNBU (native), AST3 (non-native), NAPU (navarettia), GAPA (non-native forb, galium parisiense)
# how abundant are species (meancov)
sort(sapply(ydat_count, mean)) # most are under 1% cover on average across all plots
table(sort(sapply(ydat_count, max)))
# look at variation in abundance values
sort(sapply(ydat_count, sd))
table(sort(sapply(ydat_count, function(x) length(unique(x[x>0]))))) # most species (21) have only 2 values--one of which is 0..
# in gjam, lack of variation in data would make modeling more challenging. not sure how it is for gllvm, but generally it's not helpful in ordination either
plot(sapply(ydat_count, max) ~ sapply(ydat_count, mean), main = "compost spp mean cover vs max")
plot(sapply(ydat_count, sd) ~sapply(ydat_count, function(x) length(unique(x))), main = "compost spp unique count values vs sd counts")

# rules for lumping spp: 
# rare = in less than 5% of plots (7 or fewer)
# only 1 unique non-zero cover value or max cover of 1
# no trait data (unless frequently occurring and has variation in data -- at least 3 non zero values)

# drop anything not in at least 5% of plots
nrow(ydat_count)*.05 # at least 7 plots

# pull rare (<5% plots to group by fxnl group)
# > to not over-inflate count on rare, sum on pct cover then transform to count
ydat_raresp <- names(ydat_count)[checkcolSums < nrow(ydat_count)*.05]
# ID non-rare (common)
ydat_nonrare_sp <- names(ydat_count)[checkcolSums >= nrow(ydat_count)*.05]

# screen variation in data
ydat_uniquevals <- sapply(ydat_count, function(x) length(unique(x[x >0])))
sort(ydat_uniquevals[ydat_nonrare_sp])

# screen max abundance
ydat_maxcov <- sapply(ydat_count, max)
sort(ydat_maxcov[ydat_nonrare_sp])

# screen trait availability
ydat_nonrare_sp[!ydat_nonrare_sp %in% natamb_avgtraits_wide$code4] # add AST3 and UNBU to spp to lump (bc unknown and relatively low cov and low var)
# NAPU and GAPA don't have traits but were abundant, so useful for variation in data

# add to rare spp
ydat_raresp <- sort(c(ydat_raresp, "UNBU", "AST3"))
# drop from nonrare
ydat_nonrare_sp <- ydat_nonrare_sp[!ydat_nonrare_sp %in% c("UNBU", "AST3")]

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
# remove raresp from ydat_count 
ydat_count <- cbind(ydat_count[!names(ydat_count) %in% ydat_raresp], ydat_count_rare)
# check colsums
checkcolSums <- sapply(ydat_count, function(x) sum(x>0))
sort(checkcolSums)/nrow(ydat_count)
# drop anything not present in 5% of all plots
ydat_count <- ydat_count[checkcolSums > nrow(ydat_count)*.05]
# combine background native nfix and native forb, and exo forb and exo nfix (note each nfixers only present in 2 plots)
ydat_count$NativeForbNfix <- ydat_count$NativeForb + ydat_count_rare$NativeNfix
ydat_count$ExoticForbNfix <- ydat_count$ExoticForb + ydat_count_rare$ExoticNfix
# drop NativeForb and ExoticForb cats
ydat_count <- ydat_count[!names(ydat_count) %in% c("ExoticForb", "NativeForb")]
# review
summary(ydat_count) # ok
# make sure no rare sp in final
summary(ydat_raresp %in% names(ydat_count)) # none



# specify design structure
sDesign_full <- xdat[c("plot", "block", "xpos", "ypos", "wholeplotID", "subplotID")]


# check that ydat, sDesign, and xdat are in same order
summary(rownames(ydat_count) == rownames(xdat))
head(rownames(ydat_count)) 
head(xdat) # rows ordered by factor levels rather than alpha order
ydat_count <- ydat_count[rownames(xdat),]


# -- prep widefxnl with native seeded separate -----
ydat_fxnl_wide <- merge(widefxnl_count, rowkey) %>%
  left_join(widesp_count[c("plot", "seedtrt", "herbicide", nats)]) %>%
  dplyr::select(rowid, ExoticForb:NativeGrass, FEMI, BRCA, ESCA, NEMA)
# assign rownames
rownames(ydat_fxnl_wide) <- ydat_fxnl_wide$rowid
# arrange according to xdat
ydat_fxnl_wide <- ydat_fxnl_wide[rownames(xdat), !grepl("rowid", names(ydat_fxnl_wide))]


# -- prep fxnl neighbors with seeded spp long for mixed models -----
target_long <- subset(tidysp0, code4 %in% nats, select = -pct_cover) %>%
  rename(spp = code4, count = count_cover) %>%
  # drop TRCI to low recruitment
  subset(spp != "TRCI")
fxnl_long <- subset(tidyfxnl0, grepl("Back.*Exo.*[Gras|For|N-f]", coarse_fxnl)) %>%
  rename(spp = coarse_fxnl, count = totcov_count)
fxnltarget_long <- rbind(target_long, fxnl_long[names(target_long)]) %>%
  arrange(plot, wholeplotID, seedtrt, herbicide, spp) %>%
  # create hillpos
  mutate(hillpos = ifelse(block %in% c(1,2), "downhill", "uphill"),
         hillpos = factor(hillpos))

fxnltarget_long %>%
  mutate(fulltrt = recode(fulltrt, 'XCXC' = "Control", "XCD" = "D", "XCW" = "W",
                          'FXC' = 'F', 'CXC' = 'C')) %>%
  ggplot(aes(fulltrt, count, group = spp, fill = spp)) +
  stat_summary(position = position_dodge(width = 0.5), pch = 21, size = .7) +
  scale_y_continuous(breaks = seq(0,150,25)) +
  scale_fill_manual(name = NULL, 
                    values = fxnltargcols, labels = c("Non-native forb", "Non-native grass", "Non-native N-fixer", 
                                                      "B. carinatus", "E. californica", "F. microstachys", "N. maculata")) +
  facet_grid(herbicide ~ seedtrt, scales = "free_y") +
  theme_ctw() +
  theme(panel.grid.major.x = element_line(color = "grey", linetype =2)) +
  labs(x = "Soil x precipitation treatment", y = "Abundance")
# theme(legend.position = "bottom",
#       legend.justification.bottom = "right")


# boxplot to show spread
fxnltarget_long %>%
  mutate(fulltrt = recode(fulltrt, 'XCXC' = "Control", "XCD" = "D", "XCW" = "W",
                          'FXC' = 'F', 'CXC' = 'C')) %>%
  ggplot(aes(fulltrt, count, group = paste(spp, fulltrt), fill = spp)) +
  geom_line(aes(col = spp), position = position_dodge(width = 0.75)) +
  #stat_summary(position = position_dodge(width = 0.5), pch = 21, size = .7) +
  geom_point(position = position_dodge(width = 0.75), pch = 21, size = 2) +
  scale_color_manual(values = fxnltargcols, guide = "none") +
  scale_fill_manual(name = NULL, 
                    values = fxnltargcols, labels = c("Non-native forb", "Non-native grass", "Non-native N-fixer", 
                                                      "B. carinatus", "E. californica", "F. microstachys", "N. maculata")) +
  facet_grid(herbicide ~ seedtrt, scales = "free_y") +
  theme_ctw() +
  labs(x = "Soil x precipitation treatment", y = "Abundance")
theme(legend.position = "bottom",
      legend.justification.bottom = "right",
      legend.box.just = "right")


# -- prep seeded plot y data -----
# > can re-use xdat and traitdat, just need to subset to relevant spp
# need to split out seeded so no columns with all 0s -- this also means re-aggregating fxnl cover for rare spp
# what is 5% of plots?
sum(xdat$seedtrt != "Unseeded") * 0.05 # at least 4 plots of 70 plots total

# build count first; build beta only if need (will probably use count)
ydat_count_seeded <- merge(subset(widesp_count, seedtrt != "Unseeded"), rowkey) %>%
  # order in same way as xdat (not sure why got out of order)
  arrange(plot, seedtrt, herbicide)

# assign rownames and drop all else but species
rownames(ydat_count_seeded) <- ydat_count_seeded$rowid
ydat_count_seeded <- ydat_count_seeded[names(ydat_count_seeded)[grepl("^[A-Z]+", names(ydat_count_seeded))]]
# check colsums
checkcolSums_seeded <- sapply(ydat_count_seeded, function(x) sum(x>0))
sort(checkcolSums_seeded)/nrow(ydat_count_seeded)

# check max cover, unique cover vals, and overlap with trait dat
seeded_maxcov <- sapply(ydat_count_seeded, max)
seeded_uniquecov <- sapply(ydat_count_seeded, function(x) length(unique(x[x > 0])))
seeded_notraits <- names(ydat_count_seeded)[!names(ydat_count_seeded) %in% natamb_avgtraits_wide$code4] # 22 spp
# review
sort(checkcolSums_seeded); table(checkcolSums_seeded)
sort(seeded_maxcov); table(seeded_maxcov)
sort(seeded_uniquecov); table(seeded_uniquecov)

# pulseeded_maxcov# pull rare (<5% plots to group by fxnl group)
# > to not over-inflate count on rare, sum on pct cover then transform to count
ydat_raresp_seeded <- checkcolSums_seeded[checkcolSums_seeded <= nrow(ydat_count_seeded)*.05]
ydat_nonrare_sp_seeded <- checkcolSums_seeded[checkcolSums_seeded > nrow(ydat_count_seeded)*.05]

ydat_raresp_seeded # are any of these in trait dataset?
sort(ydat_raresp_seeded[names(ydat_raresp_seeded) %in% natamb_avgtraits_wide$code4]) # split 19 yes, 19 no; 13 rare present ARE in traits db
# what of nonrare is NOT in traits db
sort(ydat_nonrare_sp_seeded[!names(ydat_nonrare_sp_seeded) %in% natamb_avgtraits_wide$code4]) # AST3, NAPU, GAPA. similar to all plots (except UNBU was in allplots)
# follow same as for all plots: add AST3 to rare and keep napu and gapa as species
sort(seeded_maxcov[names(ydat_nonrare_sp_seeded)]) # AST3 was only ever max 1, NAPU max cover is 2 -- but NAPU is in 13 plots
sort(seeded_uniquecov[names(ydat_nonrare_sp_seeded)]) # APOC and NAPU only have two non-zero unique vals


# keep pct cover df separate to recycle for beta distribution
ydat_rare_seeded <- subset(tidysp0, seedtrt != "Unseeded" & code4 %in% c(names(ydat_raresp_seeded), "AST3")) %>%
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

# drop anything in raresp
ydat_count_seeded <- ydat_count_seeded[!names(ydat_count_seeded) %in% c(names(ydat_raresp_seeded), "AST3")] %>%
  # set rowid
  mutate(rowid = rownames(.)) %>%
  # join based on rowid then keep only plant cols
  merge(ydat_count_rare_seeded) %>%
  arrange(plot, seedtrt, herbicide)
# set rownames
rownames(ydat_count_seeded) <- ydat_count_seeded$rowid
# drop anything not plant
ydat_count_seeded <- ydat_count_seeded[names(ydat_count_seeded)[grepl("^[A-Z]+|Nat|Exo", names(ydat_count_seeded))]]

# check colsums
checkcolSums_seeded <- sapply(ydat_count_seeded, function(x) sum(x>0))
sort(checkcolSums_seeded)/nrow(ydat_count_seeded)
checkcolSums_seeded # sum nfix and forbs for ambient native plants only

# combine background native nfix and native forb, and exo forb and exo nfix (note each nfixers only present in 2 plots)
ydat_count_seeded <- ydat_count_seeded[checkcolSums_seeded > nrow(ydat_count_seeded)*.05]
# sum native forb and nfix lumped spp
ydat_count_seeded$NativeForbNfix <- ydat_count_seeded$NativeForb + ydat_count_rare_seeded$NativeNfix
# repeat for exotic -- exo nfix are in 4 plots, and would keep consistent with all plots methods
ydat_count_seeded$ExoticForbNfix <- ydat_count_seeded$ExoticForb + ydat_count_rare_seeded$ExoticNfix
# drop NativeForb and ExoticForb cats
ydat_count_seeded <- ydat_count_seeded[!names(ydat_count_seeded) %in% c("NativeForb", "ExoticForb", "ExoticNfix")]
# review
summary(ydat_count_seeded) # ok
# make sure no rare sp in final
summary(names(ydat_raresp_seeded) %in% names(ydat_count_seeded)) # none
# aggregated 38 rare spp (including AST3) into their functional groups

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
names(ydat_count)[!names(ydat_count) %in% traitmat_full$code4] # napu and aggregates
names(ydat_count_seeded)[!names(ydat_count_seeded) %in% rownames(traitmat_seeded)] # napu and aggregate groups

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


# scale traitmat

# -- GLMMs FOR TREATMENT EFFECTS -----
# try all spp model but may need to run separate models for each seeded spp

# mod 0 = fully additive
# mod 1 = spp x (herb + ppt + nut) (spp effect additive)
# mod 3 = herb + spp x env interactive
# mod 2 = spp x herb + spp x env interactive


# 1) models for seeded spp only -----
# mod 0 -- this is too simple and doesn't cpapture spp effects, but for comparison
nb_nats_seeded_add <- glmmTMB(count ~ spp + herbicide + nut_trt + ppt_trt + (1|block),
                                  data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                                  family = "nbinom1")

# mod 1
nb_nats_seeded_sppxadd <- glmmTMB(count ~ spp * (herbicide + nut_trt + ppt_trt) + (1|block),
                             data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                             family = "nbinom1")
nb_nats_seeded_sppxadd2 <- glmmTMB(count ~ spp * (herbicide + nut_trt + ppt_trt) + (1|block),
                                  data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                                  family = "nbinom2")
nb_nats_seeded_sppxadd2.sum <- glmmTMB(count ~ spp * (herbicide + nut_trt + ppt_trt) + (1|block),
                                   data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                                   family = "nbinom2", contrasts = list(spp = "contr.sum"))
# mod 2
nb_nats_seeded_sppxenv <- glmmTMB(count ~ herbicide + spp * (nut_trt * ppt_trt) + (1|block), 
                      data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                      family = "nbinom1")

nb_nats_seeded_sppxenv2 <- glmmTMB(count ~ herbicide + spp * (nut_trt * ppt_trt) + (1|block), 
                                  data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                                  family = "nbinom2")
# mod 3
nb_nats_seeded_sppx <- glmmTMB(count ~ (spp * herbicide) + (spp * nut_trt * ppt_trt) + (1|block),
                               data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                               family = "nbinom1")
nb_nats_seeded_sppx2 <- glmmTMB(count ~ (spp * herbicide) + (spp * nut_trt * ppt_trt) + (1|block),
                               data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                               family = "nbinom2")

# mod 4 -- additive env, but interaction between spp, herb and each env manip
nb_nats_seeded_splitx <- glmmTMB(count ~ (spp * herbicide) * (nut_trt + ppt_trt) + (1|block),
                                 data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                                 family = "nbinom2")
summary(nb_nats_seeded_splitx)

# mod 5
nb_nats_seeded_split3x <- glmmTMB(count ~ spp * herbicide + (herbicide* nut_trt *ppt_trt) + spp * (nut_trt + ppt_trt) + (1|block),
                                 data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                                 family = "nbinom2")
summary(nb_nats_seeded_split3x)
# to be sure, try full interactive. won't run, convergence problems
# nb_nats_seeded_x <- glmmTMB(count ~ spp * herbicide * nut_trt * ppt_trt + (1|block),
#                                data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
#                                family = "nbinom2") # doesn't work w nbinom1 or 2

data.frame(anova(nb_nats_seeded_add, 
      nb_nats_seeded_sppxadd,
      nb_nats_seeded_sppxenv,
      nb_nats_seeded_sppx))

performance::compare_performance(nb_nats_seeded_add, nb_nats_seeded_sppxadd, nb_nats_seeded_sppxenv,
                                 nb_nats_seeded_sppx)
performance::check_model(nb_nats_seeded_sppx2)
performance::test_lrt(nb_nats_seeded_sppxadd2, nb_nats_seeded_sppx2)

# test nbinom2
anova(nb_nats_seeded_sppx2, nb_nats_seeded_sppxadd2,  nb_nats_seeded_splitx, nb_nats_seeded_sppxenv2)
anova(nb_nats_seeded_sppx, nb_nats_seeded_sppx2)

performance::compare_performance(nb_nats_seeded_sppx2, nb_nats_seeded_sppxadd2,  nb_nats_seeded_splitx, nb_nats_seeded_sppxenv2, rank = T)

# mods 2 and 3 best
car::Anova(nb_nats_seeded_sppxadd) 
car::Anova(nb_nats_seeded_sppx) # herb x spp is an important term, spp x each term additive is what's most important
car::Anova(nb_nats_seeded_splitx) # don't need 3-way interaction between spp, herb, and env like with background
car::Anova(nb_nats_seeded_sppxadd2.sum, type = 3)
summary(nb_nats_seeded_sppxadd)
summary(nb_nats_seeded_splitx)
summary(nb_nats_seeded_sppxadd2)
summary(nb_nats_seeded_sppxadd2.sum)
DHARMa::plotQQunif(nb_nats_seeded_sppx)
DHARMa::plotQQunif(nb_nats_seeded_sppx2)
DHARMa::plotQQunif(nb_nats_seeded_splitx)

plotResiduals(nb_nats_seeded_sppx)
plotResiduals(nb_nats_seeded_sppx2)
plotResiduals(nb_nats_seeded_splitx)

DHARMa::testResiduals(nb_nats_seeded_sppx)
DHARMa::testResiduals(nb_nats_seeded_sppx2)
testResiduals(nb_nats_seeded_splitx)

testDispersion(nb_nats_seeded_sppxenv2) # nbinom2 is a better fit than nbinom1
anova(nb_nats_seeded_sppxadd2, nb_nats_seeded_splitx) # 3-way herb x each env is marginally better

summary(nb_nats_seeded_splitx) # none of the 3way effects are signif

# > proceed for native seeded spp with additive model
# > compare results in herb vs unherb plots to be sure model picking up effects
# herbicided plots only
# interactive effects
nb_nats_seeded_herbx <- glmmTMB(count ~ spp * (nut_trt * ppt_trt) + (1|block), 
                          data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt) & !grepl("Non", herbicide)),
                          family = "nbinom2")

summary(nb_nats_seeded_herbx)
summary(nb_nats_seeded_sppxadd2) # <- this mod best, just spp interaction w each manipulation 
summary(nb_nats_seeded_sppx2)
summary(nb_nats_seeded_sppxadd2.sum)


# additive effects only
nb_nats_seeded_herb <- glmmTMB(count ~ spp * (nut_trt + ppt_trt) + (1|block), 
                               data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt) & !grepl("Non", herbicide)),
                               family = "nbinom2")
summary(nb_nats_seeded_herb)
anova(nb_nats_seeded_herb, nb_nats_seeded_herbx) # full interactive not better, just like with all plots model

# test additive in unherbicided
nb_nats_seeded_unherb <- glmmTMB(count ~ spp * (nut_trt + ppt_trt) + (1|block), 
                                data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt) & grepl("Non", herbicide)),
                                family = "nbinom2")
summary(nb_nats_seeded_unherb)

# interactive -- can't run model with 3-way interaction between spp and nut x ppt, only spp + nut x ppt
nb_nats_seeded_unherbx <- glmmTMB(count ~ spp + (nut_trt * ppt_trt) + (1|block), 
                                 data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt) & grepl("Non", herbicide)),
                                 family = "nbinom2")

summary(nb_nats_seeded_unherbx)
anova(nb_nats_seeded_unherb, nb_nats_seeded_unherbx) # additive mode is better


# for curiosity, test nema only
nb_forbs_seeded <- glmmTMB(count ~ spp * (herbicide + nut_trt + ppt_trt) + (1|block), 
                          data = subset(fxnltarget_long, spp %in% c("NEMA", "ESCA") & !grepl("Unseed", seedtrt)),
                          family = "nbinom2")
summary(nb_forbs_seeded)

nb_grass_seeded <- glmmTMB(count ~  spp * (herbicide + nut_trt + ppt_trt) + (1|block), 
                          data = subset(fxnltarget_long, spp %in% c("FEMI", "BRCA") & !grepl("Unseed", seedtrt)),
                          family = "nbinom2")
summary(nb_grass_seeded)


# compare BRCA coeffs across all mods but forb-only
natsonly_coeffs <- rbind(
#cbind(tidy(nb_nats_seeded_sppx), mod = "all plots x"),
cbind(tidy(nb_nats_seeded_sppxadd2, conf.int = T, exponentiate = F), mod = "all plots"),
cbind(tidy(nb_nats_seeded_splitx, conf.int = T, exponentiate = F), mod = "all plots, split"),
cbind(tidy(nb_nats_seeded_herb,  conf.int = T, exponentiate = F), mod = "herb plots"),
cbind(tidy(nb_nats_seeded_unherb,  conf.int = T, exponentiate = F), mod = "unherb plots"),
cbind(tidy(nb_forbs_seeded,  conf.int = T, exponentiate = F), mod = "forbs only"),
cbind(tidy(nb_grass_seeded,  conf.int = T, exponentiate = F), mod = "grass only"))


ggplot(subset(natsonly_coeffs, is.na(group) & !grepl("spp", term) & !grepl("forb", mod)), aes(term, estimate, group = mod, col = mod)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, position = position_dodge(width = 0.25)) +
  geom_point(position = position_dodge(width = 0.25), size = 3) +
  coord_flip()
 
# stack resids for comparison
natsonly_preds <- rbind(
cbind(augment(nb_nats_seeded_sppxadd2), mod = "all plots"),
cbind(augment(nb_nats_seeded_splitx), mod = "all plots, split"),
cbind(augment(nb_forbs_seeded), mod = "forbs only"),
cbind(augment(nb_grass_seeded), mod = "grass only")
) %>%
  rbind(cbind(augment(nb_nats_seeded_herb), herbicide = "Herbicided", mod = "herb plots only")[names(.)],
    cbind(augment(nb_nats_seeded_unherb), herbicide = "Non-herbicided", mod = "unherb plots only")[names(.)])
# drop .rownames
natsonly_preds <- natsonly_preds[!grepl("rown", names(natsonly_preds))]

# plot to compare
ggplot(subset(natsonly_preds), aes(mod, `.fitted`)) +
  geom_violin() +
  geom_jitter() +
  facet_wrap(~spp)

ggplot(subset(natsonly_preds), aes(paste(ppt_trt, nut_trt), exp(`.fitted`), group = paste(ppt_trt, nut_trt, mod), color = mod)) +
  geom_boxplot(position = position_dodge(width = 0.5)) +
  geom_jitter(position = position_dodge(width = 0.5)) +
  geom_point(aes(paste(ppt_trt, nut_trt), count, group = paste(ppt_trt, nut_trt, mod)), col = "black", pch = 1, size = 3) +
  facet_grid(spp ~ herbicide, scales = "free_y")
# > it's better to keep herbicide in the model -- run w all seeded plots, all spp


ggplot(subset(natsonly_preds), aes(paste(ppt_trt, nut_trt), `.resid`, group = paste(ppt_trt, nut_trt, mod), color = mod)) +
  geom_hline(aes(yintercept = 0)) +
  geom_boxplot(position = position_dodge(width = 0.5)) +
  stat_summary(position = position_dodge(width = 0.6), pch = 4) +
  geom_jitter(position = position_dodge(width = 0.4)) +
  #geom_point(aes(paste(ppt_trt, nut_trt), count, group = paste(ppt_trt, nut_trt, mod)), col = "black", pch = 1, size = 3) +
  facet_grid(spp ~ herbicide, scales = "free_y")


# model for background response
# make grass reference group``
fxnlexo_long <- subset(fxnltarget_long,!spp %in% nats) %>%
  mutate(spp = factor(spp, levels = c("Background Exotic Grass", "Background Exotic Forb", "Background Exotic N-fixer")))


# species specific effects just for herbicide, + env manips
exo_allplots_seedadd <- glmmTMB(count ~ seedtrt + spp * (herbicide + nut_trt + ppt_trt) + (1|block/wholeplotID),
                            data = fxnlexo_long, family = "nbinom2")

# 3-way interactions between spp, ppt, and nut trt
exo_allplots_2envx <- glmmTMB(count ~  spp * (seedtrt + herbicide + nut_trt * ppt_trt) + (1|block/wholeplotID), 
                data = fxnlexo_long, family = "nbinom2")

# interactions with herbicide and env
exo_allplots_3envx <- glmmTMB(count ~ seedtrt + spp * (herbicide * nut_trt * ppt_trt)+ (1|block/wholeplotID), 
                             data = fxnlexo_long, family = "nbinom2")

# all interactions (not expecting this to run?)
exo_allplots_x <- glmmTMB(count ~  spp * seedtrt * herbicide * nut_trt * ppt_trt + (1|block/wholeplotID), 
                          data = fxnlexo_long, family = "nbinom2")

# build model that has 3-way interaction with nut and ppt separately 
exo_allplots_splitenvx <- glmmTMB(count ~  seedtrt + (spp * herbicide) * (nut_trt + ppt_trt) + (1|block/wholeplotID), 
                                  data = fxnlexo_long, family = "nbinom2")

# use sum to one for average effect
exo_allplots_splitenvx.sum <- glmmTMB(count ~  seedtrt + (spp * herbicide) * (nut_trt + ppt_trt) + (1|block/wholeplotID), 
                                  data = fxnlexo_long, family = "nbinom2", contrasts = list(spp = "contr.sum"))

summary(exo_allplots_splitenvx.sum)

anova(exo_allplots_seedadd, exo_allplots_2envx, exo_allplots_3envx, exo_allplots_splitenvx,exo_allplots_x)
# seed add with interactions between all or split interaction mods are best
car::Anova(exo_allplots_seedadd)
car::Anova(exo_allplots_2envx)
car::Anova(exo_allplots_3envx)
car::Anova(exo_allplots_splitenvx)
car::Anova(exo_allplots_x)

anova(exo_allplots_splitenvx, exo_allplots_3envx) # can still use the 3-way even if full interaction not signif

# review coeffs
summary(exo_allplots_seedadd)
summary(exo_allplots_splitenvx) # best mod
summary(exo_allplots_3envx)
summary(exo_allplots_2envx)


# create forb only model to compare
exoforb_allplots_splitenvx <- glmmTMB(count ~  seedtrt + (spp * herbicide) * (nut_trt + ppt_trt) + (1|block/wholeplotID), 
                               data = subset(fxnlexo_long, !grepl("Grass", spp)), family = "nbinom2")

summary(exoforb_allplots_splitenvx)
car::Anova(exoforb_allplots_splitenvx)

exoonly_coeffs <- rbind(
  cbind(tidy(exo_allplots_2envx, conf.int = T, exponentiate = T), mod = "exo, 2envx"),
  cbind(tidy(exo_allplots_3envx, conf.int = T, exponentiate = T), mod = "exo, 3envx"),
  cbind(tidy(exo_allplots_splitenvx,  conf.int = T, exponentiate = T), mod = "exo, split"),
  cbind(tidy(exo_allplots_x,  conf.int = T, exponentiate = T), mod = "exo x"),
  cbind(tidy(exo_allplots_seedadd,  conf.int = T, exponentiate = T), mod = "exo add"))


ggplot(subset(exoonly_coeffs, is.na(group) & !grepl("spp", term)), aes(term, estimate, group = mod, col = mod)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, position = position_dodge(width = 0.25)) +
  geom_point(position = position_dodge(width = 0.25), size = 3) +
  coord_flip()

# stack resids for comparison
exoonly_preds <- rbind(
  cbind(augment(exo_allplots_2envx), mod = "exo, 2envx"),
  cbind(augment(exo_allplots_3envx), mod = "exo, 3envx"),
  cbind(augment(exo_allplots_splitenvx), mod = "exo, split"),
  cbind(augment(exo_allplots_x), mod = "exo x"),
  cbind(augment(exo_allplots_seedadd), mod = "exo add"),
  cbind(augment(exoforb_allplots_splitenvx), mod = "forb exo"))

# drop .rownames
exoonly_preds <- exoonly_preds[!grepl("rown", names(exoonly_preds))]

ggplot(subset(exoonly_preds), aes(paste(ppt_trt, nut_trt), exp(`.fitted`), group = paste(seedtrt, herbicide, ppt_trt, nut_trt, mod), color = mod)) +
  geom_boxplot(position = position_dodge(width = 0.5)) +
  geom_jitter(position = position_dodge(width = 0.5)) +
  geom_point(aes(paste(ppt_trt, nut_trt), count, group = paste(ppt_trt, nut_trt, mod)), col = "black", pch = 1, size = 3) +
  facet_grid(spp ~ herbicide + seedtrt, scales = "free_y")

ggplot(subset(exoonly_preds), aes(paste(ppt_trt, nut_trt), `.resid`, group = paste(ppt_trt, nut_trt, mod), color = mod)) +
  geom_hline(aes(yintercept = 0)) +
  geom_boxplot(position = position_dodge(width = 0.5)) +
  stat_summary(position = position_dodge(width = 0.6), pch = 4) +
  geom_jitter(position = position_dodge(width = 0.4)) +
  #geom_point(aes(paste(ppt_trt, nut_trt), count, group = paste(ppt_trt, nut_trt, mod)), col = "black", pch = 1, size = 3) +
  facet_grid(spp ~ herbicide + seedtrt, scales = "free_y")



# > for native seeded plants, direct interactions with herbcide and env is best
# > for background exotic, spp x herbicide x each env treat is best, with seedtrt additive


# plot coefficients for best exo mod
exosplit_coeffs <- tidy(exo_allplots_splitenvx, effects = "fixed", conf.int = T, exponentiate = F) %>%
  mutate(spp = str_extract(term, "Forb|N-fixer"),
         spp = ifelse(is.na(spp), "Grass", spp),
         term2 = gsub("sppBackground.*Forb:", "", term),
         term2 = gsub("sppBack.*N-fixer:", "", term2),
         term2 = ifelse(grepl("sppBack|Interc", term2), "Intercept", term2))

ggplot(exosplit_coeffs, aes(term2, estimate, col = spp, group = spp)) +
  geom_hline(aes(yintercept = 0)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  coord_flip()


exosplit.sum_coeffs <- tidy(exo_allplots_splitenvx.sum, effects = "fixed", conf.int = T, exponentiate = F) %>%
  mutate(spp = str_extract(term, "spp[1-3]"),
         spp = ifelse(is.na(spp), "Grass", spp),
         term2 = gsub("spp[1-2]:", "", term),
         term2 = ifelse(grepl("spp|Interc", term2), "Intercept", term2))

ggplot(exosplit.sum_coeffs, aes(term2, estimate, col = spp, group = spp)) +
  geom_hline(aes(yintercept = 0)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  coord_flip()


exoadd_coeffs <- tidy(exo_allplots_seedadd, effects = "fixed", conf.int = T, exponentiate = F) %>%
  mutate(spp = str_extract(term, "Forb|N-fixer"),
         spp = ifelse(is.na(spp), "Grass", spp),
         term2 = gsub("sppBackground.*Forb:", "", term),
         term2 = gsub("sppBack.*N-fixer:", "", term2),
         term2 = ifelse(grepl("sppBack|Interc", term2), "Intercept", term2))

ggplot(exoadd_coeffs, aes(term2, estimate, col = spp, group = spp)) +
  geom_hline(aes(yintercept = 0)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  coord_flip()


# plot coefficients for best native model
natssppx2_coeffs <- tidy(nb_nats_seeded_sppxadd2, effects = "fixed", conf.int = T, exponentiate = F) %>%
  mutate(spp = str_extract(term, "spp[A-Z]+"),
         spp = ifelse(is.na(spp), "BRCA", gsub("spp", "", spp)),
         term2 = gsub("spp[A-Z]+:", "", term),
         term2 = ifelse(grepl("spp|Interc", term2), "Intercept", term2))

ggplot(subset(natssppx2_coeffs, !(grepl("C:.*D", term2) & spp == "NEMA")), aes(term2, estimate, col = spp, group = spp)) +
  geom_hline(aes(yintercept = 0)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  coord_flip()


natssppx2.sum_coeffs <- tidy(nb_nats_seeded_sppxadd2.sum, effects = "fixed", conf.int = T, exponentiate = F) %>%
  mutate(spp = str_extract(term, "spp[1-3]+"),
         spp = ifelse(is.na(spp), "BRCA", gsub("spp", "", spp)),
         term2 = gsub("spp[1-3]+:", "", term),
         term2 = ifelse(grepl("spp|Interc", term2), "Intercept", term2))

ggplot(subset(natssppx2.sum_coeffs, !(grepl("C:.*D", term2) & spp == "NEMA")), aes(term2, estimate, col = spp, group = spp)) +
  geom_hline(aes(yintercept = 0)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  coord_flip()


car::linearHypothesis(nb_nats_seeded_add, "0 = nut_trtF")
plot(effects::allEffects(nb_nats_seeded_sppxadd2))
emmeans::emmeans(nb_nats_seeded_add, ~ herbicide + nut_trt + ppt_trt | spp, type = "response")
emmeans::emmeans(exo_allplots_splitenvx, ~ seedtrt + herbicide + nut_trt + ppt_trt | spp, type = "response")
plot(effects::allEffects.default(exo_allplots_splitenvx))



# -- SCREEN BEST LV ------
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


# -- COMMUNITY ORDINATION -----
# > run with ydat_count (all plots)
# what are communiites like by block or hillslope?
# 1. no predictors to show ordination
gllvm_ordination_all <- gllvm(y = ydat_fxnl_wide,
                              X = xdat, 
                              family = "negative.binomial", 
                              #reltol.c = 1e-15,
                              #starting.val = "zero",
                              #num.lv = 2,
                              num.lv.c = 4, 
                              lv.formula = ~ seedtrt * herbicide + nut_trt * ppt_trt,
                              #lv.formula  = as.formula( ~ herbicide +  nut_trt * ppt_trt),
                              # formula = as.formula(~ herbicide + seedtrt), # if have seeedtrt and 3 latent vars, have zero in covar matrix
                              # nest subplot if can, but block/wholeplot captures most of the nested similarity
                              studyDesign = sDesign_full[c("block", "wholeplotID", "subplotID")], 
                              row.eff = ~ (1|block/wholeplotID), # with subplot throws singular fit error, also throws error with wholeplotID nested [if have 3 lvs]
                              seed = myseed)
ordiplot(gllvm_ordination_all, rotate = T, jitter = T, biplot = F)
summary(gllvm_ordination_all,rotate = T)
gllvm_ordination_all$convergence
plotVarPartitioning(varPartitioning(gllvm_ordination_all), mar = c(10,10,10,10), las = 2,
                    col = palette(hcl.colors(10, "Roma")))
plot(gllvm_ordination_all)
coef(gllvm_ordination_all)
corrplot::corrplot(getResidualCov(gllvm_ordination_all)$cov, is.corr = F,  method = "color")


# how did communities separate by plot condition? (herb x seeded x env)

gllvm::predict.gllvm(gllvm_ordination_all)
gllvm_ordination_all_envpred <- gllvm(y = ydat_fxnl_wide, #[,!names(ydat_fxnl_wide) %in% nats], 
                                      X = xdat, 
                                      family = "negative.binomial", 
                                      #num.lv = 1,
                                      num.lv = 2, 
                                      formula = ~ herbicide + seedtrt + nut_trt * ppt_trt,
                                      #lv.formula  = as.formula( ~  nut_trt * ppt_trt),
                                      # formula = as.formula(~ herbicide + seedtrt), # if have seeedtrt and 3 latent vars, have zero in covar matrix
                                      # nest subplot if can, but block/wholeplot captures most of the nested similarity
                                      studyDesign = sDesign_full[c("block", "wholeplotID", "subplotID")], 
                                      row.eff = ~ (1|block/wholeplotID), # with subplot throws singular fit error, also throws error with wholeplotID nested [if have 3 lvs]
                                      seed = myseed)
ordiplot(gllvm_ordination_all_envpred, biplot =T)
summary(gllvm_ordination_all_envpred)
plotVarPartitioning(varPartitioning(gllvm_ordination_all_envpred), mar = c(10,10,10,10), las = 2,
                    col = palette(hcl.colors(10, "Roma")))
plot(gllvm_ordination_all)




# -- SEEDED ONLY -----
sDesign_seeded$hillpos <- ifelse(sDesign_seeded$block %in% c(1,2), "downhill", "uphill")
# purpose: within the restored plots
# > what is the effect of treatment on native and neighbor abundances?
gllvm_seeded <- gllvm(y = ydat_count_seeded,
                      X = xdat_seeded, 
                      family = "negative.binomial", 
                      #reltol.c = 1e-15,
                      #starting.val = "zero",
                      num.lv = 2,
                      #num.RR = 1, 
                      #lv.formula = ~ block,
                      formula = ~ herbicide + (nut_trt + ppt_trt),
                      #lv.formula  = as.formula( ~ herbicide +  nut_trt * ppt_trt),
                      # formula = as.formula(~ herbicide + seedtrt), # if have seeedtrt and 3 latent vars, have zero in covar matrix
                      # nest subplot if can, but block/wholeplot captures most of the nested similarity
                      studyDesign = sDesign_seeded[c("hillpos", "block", "wholeplotID", "subplotID")], 
                      row.eff = ~ (1|block/wholeplotID), # with subplot throws singular fit error, also throws error with wholeplotID nested [if have 3 lvs]
                      seed = myseed)

gllvm::ordiplot(gllvm_seeded, jitter = T, biplot = T)
summary(gllvm_seeded)
coefplot(gllvm_seeded)
dev.off()
plotVarPartitioning(varPartitioning(gllvm_seeded), mar = c(10,10,10,10), las = 2,
                    col = palette(hcl.colors(10, "Roma")))
plot(gllvm_seeded)
corrplot::corrplot(getResidualCor(gllvm_seeded), is.corr = F, type = "lower",diag = F,
                   order = "hclust", method = "color")
corrplot::corrplot(getResidualCov(gllvm_seeded)$cov, is.corr = F, type = "lower",diag = F,
                   order = "hclust", method = "color")


gllvm_seeded_hill <- gllvm(y = ydat_count_seeded,
                      X = xdat_seeded, 
                      family = "negative.binomial", 
                      #reltol.c = 1e-15,
                      starting.val = "zero",
                      #num.lv = 2,
                      num.lv.c= 3, 
                      #lv.formula = ~ hillpos,
                      lv.formula = ~ herbicide + nut_trt * ppt_trt,
                      #lv.formula  = as.formula( ~ herbicide +  nut_trt * ppt_trt),
                      # formula = as.formula(~ herbicide + seedtrt), # if have seeedtrt and 3 latent vars, have zero in covar matrix
                      # nest subplot if can, but block/wholeplot captures most of the nested similarity
                      studyDesign = sDesign_seeded[c("block", "wholeplotID", "subplotID")], 
                      row.eff = ~ (1|block), # with wholeplot throws singular fit error, also throws error with wholeplotID nested [if have 3 lvs]
                      seed = myseed)

ordiplot(gllvm_seeded_hill, jitter = T, biplot = F)
summary(gllvm_seeded_hill)
coefplot(gllvm_seeded)
dev.off()
plotVarPartitioning(varPartitioning(gllvm_seeded_hill), mar = c(10,10,10,10), las = 2,
                    col = palette(hcl.colors(10, "Roma")))



# -- 4th corner model -------
traitmat_scaled <- mutate_at(traitmat_full_norm, .vars = names(traitmat_full_norm)[!grepl("code4|durat|fxnl_|nativ", names(traitmat_full_norm))], function(x) as.numeric(scale(x)))



# seeded plots, only spp with traits available
gllvm_seeded_4th <- gllvm(y = ydat_count_seeded[,names(ydat_count_seeded) %in% rownames(traitmat_scaled)], 
                                        X = xdat_seeded,
                                        TR = traitmat_scaled[rownames(traitmat_scaled) %in% names(ydat_count_seeded),c(reduced_traits)],
                                        family = "negative.binomial", 
                                        num.lv = 2, 
                                        formula = y ~ (herbicide + (nut_trt * ppt_trt)) * (Coarse.root.diameter.mm + Fine.root.specific.length.cm.g + Root.volume.cm3 + Height.cm + LDMC + RMF + SLA.cm2.g),
                                        #lv.formula = ~ block,
                                        studyDesign = sDesign_seeded[c("block", "wholeplotID", "subplotID")], 
                                        #corWithin = T,
                                        row.eff = ~ (1|block/wholeplotID),
                                        #  gradient.check = T,
                                        seed = myseed)

ordiplot(gllvm_seeded_4th)
summary(gllvm_seeded_4th)
coefplot(gllvm_seeded_4th, cex.ylab = 1, mar = c(5,20,2,2))
gllvm_seeded_4th$fourth.corner

# 4th but with herbicide interacting with env vars
gllvm_seeded_4th_herbx <- gllvm(y = ydat_count_seeded[,names(ydat_count_seeded) %in% rownames(traitmat_scaled)], 
                          X = xdat_seeded,
                          TR = traitmat_scaled[rownames(traitmat_scaled) %in% names(ydat_count_seeded),c(reduced_traits)],
                          family = "negative.binomial", 
                          num.lv = 2, 
                          formula = y ~ herbicide * (nut_trt + ppt_trt) * (Coarse.root.diameter.mm + Fine.root.specific.length.cm.g + Root.volume.cm3 + Height.cm + LDMC + RMF + SLA.cm2.g),
                          #lv.formula = ~ block,
                          studyDesign = sDesign_seeded[c("block", "wholeplotID", "subplotID")], 
                          corWithin = F,
                          row.eff = ~ (1|block/wholeplotID),
                          #  gradient.check = T,
                          seed = myseed)

ordiplot(gllvm_seeded_4th_herbx, biplot = T)
summary(gllvm_seeded_4th_herbx)
coefplot(gllvm_seeded_4th_herbx, cex.ylab = 1, mar = c(5,20,2,2))
gllvm_seeded_4th_herbx$X.design

gllvm_seeded_4th_herbx$fourth.corner %>%
  data.frame() %>%
  rownames_to_column(var = "treatment") %>%
  gather(trait, val, 2:ncol(.)) %>%
  # mutate(trait = factor(trait#, 
  #                       #levels = c("Coarse.root.diameter.mm", "Fine.root.specific.length.cm.g", "Root.volume.cm3", "RMF", "LDMC", "SLA.cm2.g", "Height.cm")
  #                       )) %>%
  ggplot(aes(treatment, trait)) +
  geom_tile(aes(fill = val), col = "black") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_distiller(name = "Mean effect", type = "div", direction = 1) +
  labs(x = NULL, y = "Functional trait") +
  theme_ctw()

# additive model onlyy
gllvm_seeded_4th_add <- gllvm(y = ydat_count_seeded[,names(ydat_count_seeded) %in% rownames(traitmat_scaled)], 
                                X = xdat_seeded,
                                TR = traitmat_scaled[rownames(traitmat_scaled) %in% names(ydat_count_seeded),c(reduced_traits)],
                                family = "negative.binomial", 
                                num.lv = 2, 
                                formula = y ~ (herbicide + nut_trt + ppt_trt) * (Coarse.root.diameter.mm + Fine.root.specific.length.cm.g + Root.volume.cm3 + Height.cm + LDMC + RMF + SLA.cm2.g),
                                #lv.formula = ~ block,
                                studyDesign = sDesign_seeded[c("block", "wholeplotID", "subplotID")], 
                                #corWithin = T,
                                row.eff = ~ (1|block/wholeplotID),
                                #  gradient.check = T,
                                seed = myseed)

anova(gllvm_seeded_4th, gllvm_seeded_4th_herbx, gllvm_seeded_4th_add, gllvm_seeded_4th_traitadd)
ordiplot(gllvm_seeded_4th_add)
coefplot(gllvm_seeded_4th_add, cex.ylab = 1, mar = c(5,18,2,2))
corrplot::corrplot(gllvm::getResidualCor(gllvm_seeded_4th_add), order = "hclust", method = "color")
corrplot::corrplot(gllvm_seeded_4th_add$fourth.corner, is.corr = F, method = "color")

gllvm_seeded_4th_add$fourth.corner %>%
  data.frame() %>%
  rownames_to_column(var = "treatment") %>%
  gather(trait, val, 2:ncol(.)) %>%
  mutate(trait = factor(trait, 
                        levels = c("Coarse.root.diameter.mm", "Fine.root.specific.length.cm.g", "Root.volume.cm3", "RMF", "LDMC", "SLA.cm2.g", "Height.cm"))) %>%
  ggplot(aes(treatment, trait)) +
  geom_tile(aes(fill = val), col = "black") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_distiller(name = "Mean effect", type = "div", direction = 1) +
  labs(x = NULL, y = "Functional trait") +
  theme_ctw()

  



gllvm_seeded_4th_traitadd <- gllvm(y = ydat_count_seeded[,names(ydat_count_seeded) %in% rownames(traitmat_scaled)], 
                              X = xdat_seeded,
                              TR = traitmat_scaled[rownames(traitmat_scaled) %in% names(ydat_count_seeded),c(reduced_traits)],
                              family = "negative.binomial", 
                              num.lv = 2, 
                              formula = ~ (herbicide + nut_trt + ppt_trt) + Coarse.root.diameter.mm + Fine.root.specific.length.cm.g + Root.volume.cm3 + Height.cm + LDMC + RMF + SLA.cm2.g,
                              studyDesign = sDesign_seeded[c("block", "wholeplotID", "subplotID")], 
                              #corWithin = T,
                              row.eff = ~ (1|block/wholeplotID),
                              #  gradient.check = T,
                              seed = myseed)
ordiplot(gllvm_seeded_4th_traitadd)
AIC(gllvm_seeded_4th_traitadd)
summary(gllvm_seeded_4th_traitadd)
coefplot(gllvm_seeded_4th_traitadd, cex.ylab = 1, mar = c(5,15,2,2))


AIC(gllvm_seeded_4th, gllvm_seeded_4th_add, gllvm_seeded_4th_herbx, gllvm_seeded_4th_traitadd)
anova(gllvm_seeded_4th, gllvm_seeded_4th_add, gllvm_seeded_4th_herbx, gllvm_seeded_4th_traitadd)
anova(gllvm_seeded_4th_add, gllvm_seeded_4th_traitadd)
BIC(gllvm_seeded_4th, gllvm_seeded_4th_add, gllvm_seeded_4th_herbx, gllvm_seeded_4th_traitadd)




# -- EXTRA CODE, not needed ----
# nb2_seedspp <- glmmTMB(count ~  spp + seedtrt + (1|block) + (1|nut_trt/ppt_trt/herbicide), 
#                        data = subset(fxnltarget_long, spp %in% nats),
#                        family = "nbinom1")
# zinb_seedspp <- glmmTMB(count ~ -1 + spp * herbicide + (1|block) + (1|nut_trt/ppt_trt), 
#                         data = subset(fxnltarget_long, spp %in% nats), 
#                         ziformula = ~ seedtrt,
#                         family = "nbinom1")
# hurdnb_seedspp <- glmmTMB(count ~ -1 + spp * herbicide + (1|block/nut_trt/ppt_trt), 
#                           data = subset(fxnltarget_long, spp %in% nats), 
#                           ziformula = ~ seedtrt,
#                           family = "truncated_nbinom1")
# summary(nb_seedspp)
# summary(nb2_seedspp)
# summary(zinb_seedspp)
# summary(hurdnb_seedspp)
# AICctab(nb2_seedspp, zinb_seedspp, hurdnb_seedspp)
# 
# 
# nb3_seedspp <- glmmTMB(count ~ herbicide + spp * (nut_trt + ppt_trt) + (1|block), 
#                        data = subset(fxnltarget_long, spp %in% nats & !grepl("Unse", seedtrt)),
#                        family = "nbinom1")
# summary(nb3_seedspp)
# 
# nb_all <- glmmTMB(count ~ -1 + spp * (nut_trt + ppt_trt + herbicide) + (1|block), 
#                   data = subset(fxnltarget_long, !grepl("Unse", seedtrt)),
#                   family = "nbinom1")
# summary(nb_all)
