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
library(multcomp)
library(emmeans)
library(broom.mixed)
library(gllvm)
library(cowplot)
library(modelsummary)
theme_ctw <- function(){
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
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
         hillpos = factor(hillpos),
         # make spp factor
         spp = factor(spp))

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


# points and lines to show spread
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
traitmat_scaled <- mutate_at(traitmat_full_norm, .vars = names(traitmat_full_norm)[!grepl("code4|durat|fxnl_|nativ", names(traitmat_full_norm))], function(x) as.numeric(scale(x)))


# -- CTW CUSTOM FUNCTIONS FOR MIXED MODS ------
# function to compare model stats and confidence intervals
getmod_effects <- function(mod, modlab = "anpp 2021"){
  coefdf <- tidy(mod, conf.int = T)
  coefdf$call <- str_flatten(as.character(mod$call))
  coefdf$mod <- modlab
  coefdf <- cbind(coefdf, glance(mod)) # redundant info but so it's there
  return(coefdf)
}

# function to compile emmeans, contrasts, with confidence intervals
getemmeans <- function(emobj, conobj, modlab = "anpp 2021"){
  marg <- tidy(emobj, conf.int = T)
  marg$mod <- modlab
  emdf <- tidy(conobj, conf.int = T)
  emdf$mod <- modlab
  tmplist <- list(marg, emdf)
  names(tmplist) <- c("marginal_means", "comparisons") 
  return(tmplist)  
}


# -- GLMMs FOR TREATMENT EFFECTS -----
# means of seeded spp presence and cover
subset(tidysp0_seeded, code4 %in% nats) %>%
  group_by(code4) %>%
  reframe(nobs = length(plot),
          pres = (sum(count_cover > 0)/length(plot))*100,
          meancov = mean(count_cover),
          secov = sd(count_cover)/sqrt(length(plot)))

# add plot size to widefxnl_count
# seeded comp plots are 0.5x0.5m, main comp is 1x1m, main comp is 4x greater
widefxnl_count <- mutate(widefxnl_count, plot_scale = ifelse(grepl("Unseed", seedtrt) & grepl("Non-", herbicide), 4, 1))

# ---- herbicide effect on total cover -----
# compare totcov in unseeded herbicide and no herbicide
hist(widefxnl_count$AllNeighbors) # pretty normal looking, all neighbors is seeded + neighbors
hist(widefxnl_count$AllNeighbors[widefxnl_count$seedtrt == "Unseeded"])
summary(glmmTMB(formula = as.formula(AllNeighbors ~ herbicide + nut_trt + ppt_trt + (1|block/wholeplotID)), 
                              data = subset(widefxnl_count, seedtrt == "Unseeded"))) # convergence issue with wholeplotID nested in blcok
summary(glmmTMB(formula = as.formula(AllNeighbors ~ herbicide + nut_trt + ppt_trt + (1|block/wholeplotID) + (1|seedtrt)), 
                data = widefxnl_count)) # still get a convergence error with seeded plots included

#report this model for average effects on total cover. can run it with wholeplot ID nested
totcov_trteffects <- glmmTMB(formula = as.formula(AllNeighbors ~ herbicide *seedtrt + nut_trt + ppt_trt + (1|block/wholeplotID) + (1|plot_scale)), 
                             data = widefxnl_count) 
plotResiduals(totcov_trteffects) # running with Gaussian is okay
check_overdispersion(totcov_trteffects) # no dispersion
summary(totcov_trteffects)
emtoteffects <- emmeans(totcov_trteffects, ~ herbicide + seedtrt| ppt_trt + nut_trt)
tidy(contrast(emtoteffects, "trt.vs.ctrl", simple = "herbicide", by = "seedtrt", combine = T), conf.int = T)
tidy(contrast(emtoteffects, "trt.vs.ctrl", simple = "seedtrt", by = "herbicide", combine = T), conf.int = T)
tidy(pairs(emtoteffects, simple = "seedtrt"), conf.int = T)
# compare totcover in seeded herbicide vs seeded ambient
tidy(contrast(emtoteffects, "pairwise", by = c("ppt_trt", "nut_trt")), conf.int = T)
pairs(emtoteffects, simple = "herbicide")
# try all spp model but may need to run separate models for each seeded spp

# mod 0 = fully additive
# mod 1 = spp x (herb + ppt + nut) (spp effect additive)
# mod 3 = herb + spp x env interactive
# mod 2 = spp x herb + spp x env interactive

mod0 <- as.formula(count ~ spp + herbicide + nut_trt + ppt_trt + (1|block)) 
mod1 <- as.formula(count ~ spp * (herbicide + nut_trt + ppt_trt) + (1|block))
mod2 <- as.formula(count ~ herbicide + (spp * nut_trt * ppt_trt) + (1|block))
mod3 <- as.formula(count ~ (spp * herbicide) + (spp * nut_trt * ppt_trt) + (1|block))
mod4 <- as.formula(count ~ (spp * herbicide) + (herbicide * nut_trt * ppt_trt) + (1|block))
mod5 <- as.formula(count ~ (spp * herbicide) * (nut_trt + ppt_trt) + (1|block))

# 1a) models for seeded spp only -----
# mod 0 -- this is too simple and doesn't cpapture spp effects, but for comparison
nb_nats_seeded_mod0nb1 <- glmmTMB(formula = mod0,
                                  data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                                  family = "nbinom1")

nb_nats_seeded_mod0nb2 <- glmmTMB(formula = mod0,
                                  data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                                  family = "nbinom2")

check_model(nb_nats_seeded_mod0nb1)
check_model(nb_nats_seeded_mod0nb2)
check_overdispersion(nb_nats_seeded_mod0nb1) # overdispersion detected for nb1
check_overdispersion(nb_nats_seeded_mod0nb2) # no overdispersion

compare_performance(nb_nats_seeded_mod0nb1, nb_nats_seeded_mod0nb2) # nb 2 better fit
par(mfrow = c(2,1))
DHARMa::testDispersion(nb_nats_seeded_mod0nb1)
DHARMa::testDispersion(nb_nats_seeded_mod0nb2) # nb2 better fit
DHARMa::testResiduals(nb_nats_seeded_mod0nb1) # dispersion test signif (same test is in performance package)
DHARMa::testResiduals(nb_nats_seeded_mod0nb2)

# > proceed with nbinom2 distribution


# mod 1
nb_nats_seeded_mod1nb2 <- glmmTMB(formula = mod1,
                                  data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                                  family = "nbinom2")

# mod 2
nb_nats_seeded_mod2nb2 <- glmmTMB(formula = mod2, 
                                  data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                                  family = "nbinom2")

# mod 3
nb_nats_seeded_mod3nb2 <- glmmTMB(formula = mod3,
                                  data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                                  family = "nbinom2")

# mod 4
nb_nats_seeded_mod4nb2 <- glmmTMB(formula = mod4,
                                  data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                                  family = "nbinom2")
# mod 5
nb_nats_seeded_mod5nb2 <- glmmTMB(formula = mod5,
                                  data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt)),
                                  family = "nbinom2")
# > fully interactive model does not converge

compare_performance(nb_nats_seeded_mod0nb2, nb_nats_seeded_mod1nb2, 
                    nb_nats_seeded_mod2nb2,  nb_nats_seeded_mod3nb2,
                    nb_nats_seeded_mod4nb2, nb_nats_seeded_mod5nb2) # mod 1 seems best, but consider 2, 3 and 5
anova(nb_nats_seeded_mod0nb2, nb_nats_seeded_mod1nb2, 
      nb_nats_seeded_mod2nb2,  nb_nats_seeded_mod3nb2,
      nb_nats_seeded_mod4nb2, nb_nats_seeded_mod5nb2) # 1, 2, 3, 5 = mods to consider. mod 1 and 5 esp, but esp 1
compare_performance(nb_nats_seeded_mod1nb2, 
                    nb_nats_seeded_mod2nb2,  nb_nats_seeded_mod3nb2,
                    nb_nats_seeded_mod5nb2, rank = T)

par(mfrow = c(2,2))
lapply(list(nb_nats_seeded_mod1nb2, nb_nats_seeded_mod2nb2, nb_nats_seeded_mod3nb2, nb_nats_seeded_mod5nb2), plotResiduals)
# > mod 1 is best ITO no trend in residuals
par(mfrow = c(1,1))
testZeroInflation(nb_nats_seeded_mod1nb2) # no zerof inflation

summary(nb_nats_seeded_mod1nb2)
summary(nb_nats_seeded_mod5nb2) # not many meaningful interactions
summary(nb_nats_seeded_mod2nb2) # terrible bounds for NEMA
summary(nb_nats_seeded_mod3nb2) # bad SE for NEMA
# > only mods 1 and 5 give realistic results
plot(emmeans::emmeans(nb_nats_seeded_mod1nb2, ~ nut_trt + ppt_trt  + herbicide| spp), type = "response")
plot(emmeans::emmeans(nb_nats_seeded_mod5nb2, ~ nut_trt + ppt_trt  + herbicide| spp), type = "response") # NE

# assess terms
car::Anova(nb_nats_seeded_mod1nb2) # all but nutrient (direct) important.
testResiduals(nb_nats_seeded_mod1nb2) # everything good


# > proceed for native seeded spp with additive model
# > compare results in herb vs unherb plots to be sure model picking up effects
# herbicided plots only
# interactive effects
nb_nats_seeded_herbx <- glmmTMB(count ~ spp * (nut_trt * ppt_trt) + (1|block), 
                                data = subset(fxnltarget_long, spp %in% nats & !grepl("Unseed", seedtrt) & !grepl("Non", herbicide)),
                                family = "nbinom2")

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
anova(nb_nats_seeded_unherb, nb_nats_seeded_unherbx) # additive model is better


# for curiosity, test nema only
nb_forbs_seeded_mod1nb2 <- glmmTMB(formula = mod1, 
                                   data = subset(fxnltarget_long, spp %in% c("NEMA", "ESCA") & !grepl("Unseed", seedtrt)),
                                   family = "nbinom2")
summary(nb_forbs_seeded_mod1nb2)

nb_grass_seeded_mod1nb2 <- glmmTMB(formula = mod1, 
                                   data = subset(fxnltarget_long, spp %in% c("FEMI", "BRCA") & !grepl("Unseed", seedtrt)),
                                   family = "nbinom2")
summary(nb_grass_seeded_mod1nb2)


# compare BRCA coeffs across all mods but forb-only
natsonly_coeffs <- rbind(
  #cbind(tidy(nb_nats_seeded_sppx), mod = "all plots x"),
  cbind(tidy(nb_nats_seeded_mod1nb2, conf.int = T, exponentiate = F), mod = "all plots"),
  cbind(tidy(nb_nats_seeded_herb,  conf.int = T, exponentiate = F), mod = "herb plots"),
  cbind(tidy(nb_nats_seeded_unherb,  conf.int = T, exponentiate = F), mod = "unherb plots"),
  cbind(tidy(nb_forbs_seeded_mod1nb2,  conf.int = T, exponentiate = F), mod = "forbs only"),
  cbind(tidy(nb_grass_seeded_mod1nb2,  conf.int = T, exponentiate = F), mod = "grass only"))


ggplot(subset(natsonly_coeffs, is.na(group) & !grepl("spp", term) & !grepl("forb", mod)), aes(term, estimate, group = mod, col = mod)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, position = position_dodge(width = 0.25)) +
  geom_point(position = position_dodge(width = 0.25), size = 3) +
  coord_flip()

# stack resids for comparison
natsonly_preds <- rbind(
  cbind(augment(nb_nats_seeded_mod1nb2), mod = "all plots"),
  cbind(augment(nb_forbs_seeded_mod1nb2), mod = "forbs only"),
  cbind(augment(nb_grass_seeded_mod1nb2), mod = "grass only")
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



# 1b) figure for nat effects ----- 

# gather marginal means to compare with means from raw data
summary(glht(nb_nats_seeded_mod1nb2,  emmeans::lsm(~ spp : herbicide : ppt_trt : nut_trt)))
mod_means <- emmeans(nb_nats_seeded_mod1nb2, ~ herbicide + ppt_trt + nut_trt | spp , type = "response")
emmeans(nb_nats_seeded_mod1nb2, "trt.vs.ctrl" ~ herbicide | spp , type = "response")
emmeans(nb_nats_seeded_mod1nb2, "trt.vs.ctrl" ~ nut_trt | spp , type = "response")
emmeans(nb_nats_seeded_mod1nb2, "trt.vs.ctrl" ~ ppt_trt | spp , type = "response")
contrast(mod_means, "trt.vs.ctrl")
contrast(mod_means, "trt.vs.ctrl", simple = "herbicide", by ="spp")
contrast(mod_means, "trt.vs.ctrl", simple = "nut_trt", by ="spp")
cld(mod_means)
nb_nats_margmeans <- data.frame(mod_means)
mod_means_log <- emmeans(nb_nats_seeded_mod1nb2, ~ herbicide + ppt_trt + nut_trt | spp)

nb_nats_logmeans <- data.frame(mod_means_log) %>%
  # denote which ones have clear neg or pos trend
  mutate(sig95 = (asymp.LCL * asymp.UCL) >0)
# kns wants to see within which terms and interactions are significant

# > may need to run separate model for each species for within differences until can figure out 
brca_mod <- glmmTMB(formula = count ~ herbicide + (nut_trt * ppt_trt) + (1 | block), 
                    data = subset(fxnltarget_long, spp %in% c("BRCA") & !grepl("Unseed", seedtrt)),
                    family = "nbinom2")
summary(brca_mod)
car::Anova(brca_mod) # BRCA responds to each and interaction between ppt and nut
femi_mod <- update(brca_mod, 
                   formula = count ~ herbicide * (nut_trt + ppt_trt) + (1 | block),
                   data = subset(fxnltarget_long, spp %in% c("FEMI") & !grepl("Unseed", seedtrt))) 
car::Anova(femi_mod) # FEMI only cares about herbicide and herbicide x ppt
esca_mod <- update(brca_mod, 
                   data = subset(fxnltarget_long, spp %in% c("ESCA") & !grepl("Unseed", seedtrt))) 
car::Anova(esca_mod) # ESCA doesn't care about anything
nema_mod <- update(brca_mod, 
                   data = subset(fxnltarget_long, spp %in% c("NEMA") & !grepl("Unseed", seedtrt))) 
car::Anova(nema_mod) # NEMA cares about ppt most, then nutrient, then marginally herbicide

# joint test from emmeans package
data.frame(joint_tests(nb_nats_seeded_mod1nb2, by = "spp")) # this test does what i'm looking for, no need for separate mods

# plot for seeded natives with marginal means
natfig <- fxnltarget_long %>%
  subset(spp %in% nats & grepl("Nat", seedtrt)) %>%
  left_join(nb_nats_margmeans) %>%
  mutate(fulltrt = recode(fulltrt, 'XCXC' = "Control", "XCD" = "D", "XCW" = "W",
                          'FXC' = 'F', 'CXC' = 'C'),
         ppt_trt = relevel(ppt_trt, ref = "D")) %>%
  ggplot(aes(fulltrt, count, group = paste(spp, fulltrt, herbicide))) +
  geom_vline(aes(xintercept = 3.5), lty = 2, col = "grey50") +
  geom_vline(aes(xintercept = 6.5), lty = 2, col = "grey50") +
  geom_line(position = position_nudge(x = -0.3)) +
  #stat_summary(aes(fill = ppt_trt), position = position_nudge(x = -0.3), pch = 21, size = .7) +
  geom_point(aes(fill = ppt_trt), position = position_nudge(x = -0.3), pch = 21, size = 2) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0, lwd = 1) +
  geom_point(aes(fulltrt, response, fill = ppt_trt), pch = 21, size = 3) +
  scale_color_brewer(name = "Precipitation\ntreatment", palette = "BrBG", 
                     labels = c("Drought (D)", "Control", "Wet (W)")) +
  scale_fill_brewer(name = "Precipitation\ntreatment", palette = "BrBG",
                    labels = c("Drought (D)", "Control", "Wet (W)")) +
  facet_grid(spp ~ herbicide, scales = "free_y") +
  theme_ctw() +
  theme(legend.justification.right = "top") +
  labs(x = "Soil x precipitation treatment", y = "Abundance") +
  theme(plot.margin = unit(c(1, 4, 1, 1), "lines"))

# get anova stats
nats_best_anova <- data.frame(car::Anova(nb_nats_seeded_mod1nb2)) %>%
  mutate(term = rownames(.)) %>%
  rename_all(function(x) gsub("[.]", "", x))


natfig
# manually add anova stats for now
grid::grid.text(label = "Type II Wald Chisq\n spp ***\n herbicide ***\n soil trt n.s.\n precip trt ***\n spp:herbicide **\n spp:soil trt ***\n spp:precip trt ***", 
                x = unit(0.8, "npc"), y = unit(0.65, "npc"), just = "left")
# brca
grid::grid.text(label = "herb. **\nsoil ***\nppt ***", 
                x = unit(0.07, "npc"), y = unit(0.9, "npc"), just = "left")
# esca
grid::grid.text(label = "herb. n.s.\nsoil n.s.\nppt n.s.", 
                x = unit(0.07, "npc"), y = unit(0.68, "npc"), just = "left")
# femi
grid::grid.text(label = "herb. ***\nsoil *\nppt n.s.", 
                x = unit(0.07, "npc"), y = unit(0.46, "npc"), just = "left")
# nema
grid::grid.text(label = "herb. n.s.\nsoil *\nppt ***", 
                x = unit(0.07, "npc"), y = unit(0.24, "npc"), just = "left")

# add significant terms to each spp row
cld(mod_means, adjust = "sidak", Letters = letters, alpha = 0.05, by = "spp")
summary(glht(nb_nats_seeded_mod1nb2))
summary(nb_nats_seeded_mod1nb2)


# 1c) nat seed tables, contrasts ------
nats_bestcoeffs <- getmod_effects(nb_nats_seeded_mod1nb2, "native seeded glmmTMB")
nats_emlist <- getemmeans(mod_means, contrast(mod_means, "trt.vs.ctrl"), "native seeded glmmTMB")
# best model joint test, arrange by species
write_csv(arrange(data.frame(joint_tests(nb_nats_seeded_mod1nb2, by = "spp")), spp), "Native-analysis/outputs/natglmmTMB_jointtest.csv")
# wald test for best model
write_csv(nats_best_anova, "Native-analysis/outputs/natglmmTMB_anova.csv")
# model coeffs
write_csv(nats_bestcoeffs, "Native-analysis/outputs/natglmmTMB_coeffs.csv")
# marg means
write_csv(nats_emlist$marginal_means, "Native-analysis/outputs/natglmmTMB_margmeans.csv")
# contrasts
write_csv(nats_emlist$comparisons, "Native-analysis/outputs/natglmmTMB_margmeans.csv")
write_csv(data.frame(emmeans(nb_nats_seeded_mod1nb2, "trt.vs.ctrl" ~ herbicide | spp , type = "response")), "Native-analysis/outputs/natglmmTMB_herbcontrast.csv")
write_csv(data.frame(emmeans(nb_nats_seeded_mod1nb2, "trt.vs.ctrl" ~ nut_trt | spp , type = "response")), "Native-analysis/outputs/natglmmTMB_nutcontrast.csv")
write_csv(data.frame(emmeans(nb_nats_seeded_mod1nb2, "trt.vs.ctrl" ~ ppt_trt | spp , type = "response")), "Native-analysis/outputs/natglmmTMB_pptcontrast.csv")





# 2b) model for background response ----
# make grass reference group
fxnlexo_long <- subset(fxnltarget_long,!spp %in% nats) %>%
  mutate(spp = factor(spp, levels = c("Background Exotic Grass", "Background Exotic Forb", "Background Exotic N-fixer")))

exomod0 <- as.formula(count ~ spp + seedtrt + herbicide + nut_trt + ppt_trt + (1|block/wholeplotID)) 
exomod1 <- as.formula(count ~ spp * (seedtrt + herbicide + nut_trt + ppt_trt) + (1|block/wholeplotID))
exomod2 <- as.formula(count ~ seedtrt + herbicide + (spp * nut_trt * ppt_trt) + (1|block/wholeplotID))
exomod3 <- as.formula(count ~ seedtrt + (spp * herbicide) + (spp * nut_trt * ppt_trt) + (1|block/wholeplotID))
exomod4 <- as.formula(count ~ seedtrt + (spp * herbicide) + (herbicide * nut_trt * ppt_trt) + (1|block/wholeplotID))
exomod5 <- as.formula(count ~ seedtrt + (spp * herbicide) * (nut_trt + ppt_trt) + (1|block/wholeplotID))
exomod6 <- as.formula(count ~ seedtrt + (spp * herbicide * nut_trt * ppt_trt) + (1|block/wholeplotID))
exomod7 <- as.formula(count ~ (seedtrt * spp * herbicide) + (spp *nut_trt * ppt_trt) + (1|block/wholeplotID))
exomod8 <- as.formula(count ~ (seedtrt * spp) + (spp * herbicide) * (nut_trt + ppt_trt) + (1|block/wholeplotID))
exomod9 <- as.formula(count ~ seedtrt + (spp * herbicide) + spp * (nut_trt + ppt_trt) + (1|block/wholeplotID))

# run mods
# > check which fit is best
exo_allplots_mod0nb1 <- glmmTMB(formula = exomod0, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod0nb2 <- glmmTMB(formula = exomod0, data = fxnlexo_long, family = "nbinom2")

compare_performance(exo_allplots_mod0nb1, exo_allplots_mod0nb2)
check_model(exo_allplots_mod0nb1)
check_model(exo_allplots_mod0nb2)
testResiduals(exo_allplots_mod0nb1)
check_zeroinflation(exo_allplots_mod0nb1) # overfits zeros with nbinom1
check_zeroinflation(exo_allplots_mod0nb2) # underfits zeros with nbinom2
check_overdispersion(exo_allplots_mod0nb1) # no overdispersion
check_overdispersion(exo_allplots_mod0nb2) # underdispersion
plotResiduals(exo_allplots_mod0nb1)
plotResiduals(exo_allplots_mod0nb2)

# try zinb or hurdle
exo_allplots_mod0zinb1 <-glmmTMB(formula = exomod0, ziformula = as.formula(~1), data = fxnlexo_long, family = "truncated_nbinom1")
plotResiduals(exo_allplots_mod0zinb1)
check_overdispersion(exo_allplots_mod0zinb1) # even more underdispersed

# nbinom2
exo_allplots_mod1nb2 <- glmmTMB(formula = exomod1, data = fxnlexo_long, family = "nbinom2")
exo_allplots_mod2nb2 <- glmmTMB(formula = exomod2, data = fxnlexo_long, family = "nbinom2")
exo_allplots_mod3nb2 <- glmmTMB(formula = exomod3, data = fxnlexo_long, family = "nbinom2")
exo_allplots_mod4nb2 <- glmmTMB(formula = exomod4, data = fxnlexo_long, family = "nbinom2")
exo_allplots_mod5nb2 <- glmmTMB(formula = exomod5, data = fxnlexo_long, family = "nbinom2")
exo_allplots_mod6nb2 <- glmmTMB(formula = exomod6, data = fxnlexo_long, family = "nbinom2")
exo_allplots_mod7nb2 <- glmmTMB(formula = exomod7, data = fxnlexo_long, family = "nbinom2")
exo_allplots_mod8nb2 <- glmmTMB(formula = exomod8, data = fxnlexo_long, family = "nbinom2")
exo_allplots_mod9nb2 <- glmmTMB(formula = exomod9, data = fxnlexo_long, family = "nbinom2")

# nbinom1
exo_allplots_mod1nb1 <- glmmTMB(formula = exomod1, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod2nb1 <- glmmTMB(formula = exomod2, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod3nb1 <- glmmTMB(formula = exomod3, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod4nb1 <- glmmTMB(formula = exomod4, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod5nb1 <- glmmTMB(formula = exomod5, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod6nb1 <- glmmTMB(formula = exomod6, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod7nb1 <- glmmTMB(formula = exomod7, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod8nb1 <- glmmTMB(formula = exomod8, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod9nb1 <- glmmTMB(formula = exomod9, data = fxnlexo_long, family = "nbinom1")


# test dispersion on all
par(mfrow = c(2,5))
lapply(list(exo_allplots_mod0nb2, exo_allplots_mod1nb2, exo_allplots_mod2nb2, exo_allplots_mod3nb2,
            exo_allplots_mod4nb2, exo_allplots_mod5nb2, exo_allplots_mod6nb2, 
            exo_allplots_mod7nb2, exo_allplots_mod8nb2,exo_allplots_mod9nb2),
       testDispersion) # underdispersed on all
lapply(list(exo_allplots_mod0nb2, exo_allplots_mod1nb2, exo_allplots_mod2nb2, exo_allplots_mod3nb2,
            exo_allplots_mod4nb2, exo_allplots_mod5nb2, exo_allplots_mod6nb2, 
            exo_allplots_mod7nb2, exo_allplots_mod8nb2,exo_allplots_mod9nb2),
       check_overdispersion) # underdispersed on all
lapply(list(exo_allplots_mod0nb2, exo_allplots_mod1nb2, exo_allplots_mod2nb2, exo_allplots_mod3nb2,
            exo_allplots_mod4nb2, exo_allplots_mod5nb2, exo_allplots_mod6nb2, 
            exo_allplots_mod7nb2, exo_allplots_mod8nb2,exo_allplots_mod9nb2),
       check_zeroinflation) # underpredicting 0s by 5-8 (mod 8 is least worst), probable zero inflation
lapply(list(exo_allplots_mod0nb2, exo_allplots_mod1nb2, exo_allplots_mod2nb2, exo_allplots_mod3nb2,
            exo_allplots_mod4nb2, exo_allplots_mod5nb2, exo_allplots_mod6nb2, 
            exo_allplots_mod7nb2, exo_allplots_mod8nb2,exo_allplots_mod9nb2),
       plotResiduals) 
# compare with nbinom1
lapply(list(exo_allplots_mod0nb1, exo_allplots_mod1nb1, exo_allplots_mod2nb1, exo_allplots_mod3nb1,
            exo_allplots_mod4nb1, exo_allplots_mod5nb1, exo_allplots_mod6nb1, 
            exo_allplots_mod7nb1, exo_allplots_mod8nb1,exo_allplots_mod9nb1),
       testDispersion)
lapply(list(exo_allplots_mod0nb1, exo_allplots_mod1nb1, exo_allplots_mod2nb1, exo_allplots_mod3nb1,
            exo_allplots_mod4nb1, exo_allplots_mod5nb1, exo_allplots_mod6nb1, 
            exo_allplots_mod7nb1, exo_allplots_mod8nb1,exo_allplots_mod9nb1),
       check_overdispersion) # no overdispersion
lapply(list(exo_allplots_mod0nb1, exo_allplots_mod1nb1, exo_allplots_mod2nb1, exo_allplots_mod3nb1,
            exo_allplots_mod4nb1, exo_allplots_mod5nb1, exo_allplots_mod6nb1, 
            exo_allplots_mod7nb1, exo_allplots_mod8nb1,exo_allplots_mod9nb1),
       check_zeroinflation) # overfitting zeros (e.g., predict 39-40 but only 27 observed)
lapply(list(exo_allplots_mod0nb1, exo_allplots_mod1nb1, exo_allplots_mod2nb1, exo_allplots_mod3nb1,
            exo_allplots_mod4nb1, exo_allplots_mod5nb1, exo_allplots_mod6nb1, 
            exo_allplots_mod7nb1, exo_allplots_mod8nb1,exo_allplots_mod9nb1),
       plotResiduals) 


# test zinb
exo_allplots_mod1zinb1 <- glmmTMB(formula = exomod1, ziformula = ~ herbicide, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod2zinb1 <- glmmTMB(formula = exomod2, ziformula = ~ herbicide, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod3zinb1 <- glmmTMB(formula = exomod3, ziformula = ~ herbicide, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod4zinb1 <- glmmTMB(formula = exomod4, ziformula = ~ herbicide, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod5zinb1 <- glmmTMB(formula = exomod5, ziformula = ~ herbicide, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod6zinb1 <- glmmTMB(formula = exomod6, ziformula = ~ herbicide, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod7zinb1 <- glmmTMB(formula = exomod7, ziformula = ~ herbicide, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod8zinb1 <- glmmTMB(formula = exomod8, ziformula = ~ herbicide, data = fxnlexo_long, family = "nbinom1")
exo_allplots_mod9zinb1 <- glmmTMB(formula = exomod9, ziformula = ~ herbicide, data = fxnlexo_long, family = "nbinom1")

# are zinb mods any better?
lapply(list(exo_allplots_mod0zinb1, exo_allplots_mod1zinb1, exo_allplots_mod2zinb1, exo_allplots_mod3zinb1,
            exo_allplots_mod4zinb1, exo_allplots_mod5zinb1, exo_allplots_mod6zinb1, 
            exo_allplots_mod7zinb1, exo_allplots_mod8zinb1,exo_allplots_mod9zinb1),
       check_overdispersion) # no overdispersion
lapply(list(exo_allplots_mod0zinb1, exo_allplots_mod1zinb1, exo_allplots_mod2zinb1, exo_allplots_mod3zinb1,
            exo_allplots_mod4zinb1, exo_allplots_mod5zinb1, exo_allplots_mod6zinb1, 
            exo_allplots_mod7zinb1, exo_allplots_mod8zinb1,exo_allplots_mod9zinb1),
       check_zeroinflation)
lapply(list(exo_allplots_mod0zinb1, exo_allplots_mod1zinb1, exo_allplots_mod2zinb1, exo_allplots_mod3zinb1,
            exo_allplots_mod4zinb1, exo_allplots_mod5zinb1, exo_allplots_mod6zinb1, 
            exo_allplots_mod7zinb1, exo_allplots_mod8zinb1,exo_allplots_mod9zinb1),
       plotResiduals) # zinb doesn't help any of these residuals

# test zip?
exo_allplots_mod5zip <- glmmTMB(formula = exomod5, ziformula = ~ herbicide, data = fxnlexo_long, family = "poisson")
# hurdle
exo_allplots_mod5hurd1 <- glmmTMB(formula = exomod5, ziformula = ~ herbicide, data = fxnlexo_long, family = "truncated_nbinom1")
# lognorm
exo_allplots_mod5lognorm <- glmmTMB(formula = exomod5, ziformula = ~ herbicide, data = fxnlexo_long, family = "lognormal")
# generalized poisson
exo_allplots_mod5p <- glmmTMB(formula = exomod5, data = fxnlexo_long, family = "genpois")

par(mfrow = c(2,3))
lapply(list(exo_allplots_mod5zip, exo_allplots_mod5nb1, exo_allplots_mod5nb2, exo_allplots_mod5zinb1, exo_allplots_mod5hurd1, exo_allplots_mod5p), plotResiduals)
summary(exo_allplots_mod5zinb1) # zinb doesn't help, hurdle mod is worse, lognorm no good. get error messages with zinb + hurd models

# compare nbinom vs zinb
compare_performance(exo_allplots_mod5nb1, exo_allplots_mod5zinb1, exo_allplots_mod5p, rank = T) # simple nbinom1 is least terrible
anova(exo_allplots_mod5p, exo_allplots_mod5nb1) # nb1 not an improvement

summary(exo_allplots_mod5nb1)
summary(exo_allplots_mod5p) # pretty similar

check_zeroinflation(exo_allplots_mod5nb1)
check_zeroinflation(exo_allplots_mod5p) # generalized predicts fewer zeros

# test genp distributions for all mods
exo_allplots_mod1genp <- glmmTMB(formula = exomod1, data = fxnlexo_long, family = "genpois")
exo_allplots_mod2genp <- glmmTMB(formula = exomod2, data = fxnlexo_long, family = "genpois")
exo_allplots_mod3genp <- glmmTMB(formula = exomod3, data = fxnlexo_long, family = "genpois")
exo_allplots_mod4genp <- glmmTMB(formula = exomod4, data = fxnlexo_long, family = "genpois")
exo_allplots_mod5genp <- glmmTMB(formula = exomod5, data = fxnlexo_long, family = "genpois")
exo_allplots_mod6genp <- glmmTMB(formula = exomod6, data = fxnlexo_long, family = "genpois")
exo_allplots_mod7genp <- glmmTMB(formula = exomod7, data = fxnlexo_long, family = "genpois")
exo_allplots_mod8genp <- glmmTMB(formula = exomod8, data = fxnlexo_long, family = "genpois")
exo_allplots_mod9genp <- glmmTMB(formula = exomod9, data = fxnlexo_long, family = "genpois")

# compare nbinom1 mods
compare_performance(exo_allplots_mod1genp, exo_allplots_mod2genp, exo_allplots_mod3genp,
                    exo_allplots_mod4genp, exo_allplots_mod5genp, exo_allplots_mod6genp, exo_allplots_mod7genp,
                    exo_allplots_mod8genp, exo_allplots_mod9genp) # 5 still best then 8
# heuristic rank performance
compare_performance(exo_allplots_mod0nb1, exo_allplots_mod1nb1, exo_allplots_mod2nb1, exo_allplots_mod3nb1,
                    exo_allplots_mod4nb1, exo_allplots_mod5nb1, exo_allplots_mod6nb1, exo_allplots_mod7nb1,
                    exo_allplots_mod8nb1, exo_allplots_mod9nb1, rank = T) # mod 5 by far, then 8, 9, 6

par(mfrow = c(3,3))
lapply(list(exo_allplots_mod1genp, exo_allplots_mod2genp, exo_allplots_mod3genp,
            exo_allplots_mod4genp, exo_allplots_mod5genp, exo_allplots_mod6genp, exo_allplots_mod7genp,
            exo_allplots_mod8genp, exo_allplots_mod9genp), plotResiduals)
compare_performance(exo_allplots_mod5genp, exo_allplots_mod6genp,
                    exo_allplots_mod8genp, exo_allplots_mod4genp, rank = T) # still 5
compare_performance(exo_allplots_mod5genp, exo_allplots_mod5nb1)
check_overdispersion(exo_allplots_mod5genp) # not overdispersed, slightly overfits zeros 
check_distribution(exo_allplots_mod5zinb1) # says most likely cauchy resids with zinb, but that model is not better

summary(exo_allplots_mod5genp)


# evaluate terms
car::Anova(exo_allplots_mod5nb1) 
car::Anova(exo_allplots_mod5genp) # all but herb x nut and nut 

plot(emmeans(exo_allplots_mod5nb1, ~  nut_trt + ppt_trt | spp + seedtrt + herbicide)) # means look slightly reduced in seed trt across board
plot(emmeans(exo_allplots_mod5genp, ~  nut_trt + ppt_trt | spp + seedtrt + herbicide)) # prettyy similar



# create forb only model to compare
exoforb_allplots_mod0nb2 <- glmmTMB(formula = exomod0, 
                                    data = subset(fxnlexo_long, !grepl("Grass", spp)), family = "nbinom2")
exoforb_allplots_mod1nb2 <- glmmTMB(formula = exomod1, 
                                    data = subset(fxnlexo_long, !grepl("Grass", spp)), family = "nbinom2")
exoforb_allplots_mod2nb2 <- glmmTMB(formula = exomod2, 
                                    data = subset(fxnlexo_long, !grepl("Grass", spp)), family = "nbinom2")
exoforb_allplots_mod3nb2 <- glmmTMB(formula = exomod3, 
                                    data = subset(fxnlexo_long, !grepl("Grass", spp)), family = "nbinom2")
exoforb_allplots_mod4nb2 <- glmmTMB(formula = exomod4, 
                                    data = subset(fxnlexo_long, !grepl("Grass", spp)), family = "nbinom2")
exoforb_allplots_mod5nb2 <- glmmTMB(formula = exomod5, 
                                    data = subset(fxnlexo_long, !grepl("Grass", spp)), family = "nbinom2")
car::Anova(exoforb_allplots_mod3nb2)
car::Anova(exoforb_allplots_mod4nb2)
car::Anova(exoforb_allplots_mod1nb2)
anova(exoforb_allplots_mod0nb2, exoforb_allplots_mod1nb2, exoforb_allplots_mod2nb2, 
      exoforb_allplots_mod3nb2, exoforb_allplots_mod4nb2, exoforb_allplots_mod5nb2) 
# mod 1 and 5 are better, but mods 2 and 3 too
compare_performance(exoforb_allplots_mod0nb2, exoforb_allplots_mod1nb2, exoforb_allplots_mod2nb2, 
                    exoforb_allplots_mod3nb2, exoforb_allplots_mod4nb2, exoforb_allplots_mod5nb2) # 5 seems best
compare_performance(exoforb_allplots_mod0nb2, exoforb_allplots_mod1nb2, exoforb_allplots_mod2nb2, 
                    exoforb_allplots_mod3nb2, exoforb_allplots_mod4nb2, exoforb_allplots_mod5nb2, rank = T) # 5 best in rank, followed by 4. mods 1-2 comparable
par(mfrow = c(2,3))
lapply(list(exoforb_allplots_mod0nb2, exoforb_allplots_mod1nb2, exoforb_allplots_mod2nb2, exoforb_allplots_mod3nb2, exoforb_allplots_mod4nb2, exoforb_allplots_mod5nb2), plotResiduals)
# mod 1,2, 3 looks best ITO resids

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
natssppx2_coeffs <- tidy(exo_allplots_sppxadd2, effects = "fixed", conf.int = T, exponentiate = F) %>%
  mutate(spp = str_extract(term, "spp[A-Z]+"),
         spp = ifelse(is.na(spp), "BRCA", gsub("spp", "", spp)),
         term2 = gsub("spp[A-Z]+:", "", term),
         term2 = ifelse(grepl("spp|Interc", term2), "Intercept", term2))

ggplot(subset(natssppx2_coeffs, !(grepl("C:.*D", term2) & spp == "NEMA")), aes(term2, estimate, col = spp, group = spp)) +
  geom_hline(aes(yintercept = 0)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  coord_flip()


natssppx2.sum_coeffs <- tidy(exo_allplots_sppxadd2.sum, effects = "fixed", conf.int = T, exponentiate = F) %>%
  mutate(spp = str_extract(term, "spp[1-3]+"),
         spp = ifelse(is.na(spp), "BRCA", gsub("spp", "", spp)),
         term2 = gsub("spp[1-3]+:", "", term),
         term2 = ifelse(grepl("spp|Interc", term2), "Intercept", term2))

ggplot(subset(natssppx2.sum_coeffs, !(grepl("C:.*D", term2) & spp == "NEMA")), aes(term2, estimate, col = spp, group = spp)) +
  geom_hline(aes(yintercept = 0)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  coord_flip()


car::linearHypothesis(exo_allplots_add, "0 = nut_trtF")
plot(effects::allEffects(exo_allplots_sppxadd2))
emmeans::emmeans(exo_allplots_add, ~ herbicide + nut_trt + ppt_trt | spp, type = "response")
emmeans::emmeans(exo_allplots_splitenvx, ~ seedtrt + herbicide + nut_trt + ppt_trt | spp, type = "response")
plot(effects::allEffects.default(exo_allplots_splitenvx))


# 2b) figure for background effects ------
# pull marginal means to plot
# gather marginal means to compare with means from raw data
summary(glht(exo_allplots_mod5genp,  emmeans::lsm(~ spp : seedtrt: herbicide : ppt_trt : nut_trt)))
exomod5  
exo_mod_means<- emmeans(exo_allplots_mod5genp, ~ seedtrt + herbicide + ppt_trt + nut_trt | spp , type = "response")
exo_margmeans <- data.frame(exo_mod_means)

emmeans(exo_allplots_mod5genp, trt.vs.ctrl ~ herbicide | spp + seedtrt, type = "response")

emmeans(exo_allplots_mod5genp, trt.vs.ctrl ~ ppt_trt | spp + seedtrt + nut_trt, type = "response")
emmeans(exo_allplots_mod5genp, trt.vs.ctrl ~ nut_trt | spp + seedtrt + ppt_trt + herbicide, type = "response")


exofig <- fxnltarget_long %>%
  subset(!spp %in% nats) %>%
  left_join(exo_margmeans) %>%
  mutate(fulltrt = recode(fulltrt, 'XCXC' = "Con.", "XCD" = "D", "XCW" = "W",
                          'FXC' = 'F', 'CXC' = 'C'),
         ppt_trt = relevel(ppt_trt, ref = "D"),
         spp = gsub("Background Exotic", "Non-native", spp),
         spp = factor(spp),
         spp = relevel(spp, ref = "Non-native Grass")) %>%
  ggplot(aes(fulltrt, count, group = paste(spp, fulltrt, herbicide))) +
  geom_vline(aes(xintercept = 3.5), lty = 2, col = "grey50") +
  geom_vline(aes(xintercept = 6.5), lty = 2, col = "grey50") +
  geom_line(position = position_nudge(x = -0.3)) +
  #stat_summary(aes(fill = ppt_trt), position = position_nudge(x = -0.3), pch = 21, size = .7) +
  geom_point(aes(fill = ppt_trt), position = position_nudge(x = -0.3), pch = 21, size = 2) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0, lwd = 1) +
  geom_point(aes(fulltrt, response, fill = ppt_trt), pch = 21, size = 3) +
  scale_color_brewer(name = "Precipitation\ntreatment", palette = "BrBG", labels = c("Drought (D)", "Control", "Wet (W)")) +
  scale_fill_brewer(name = "Precipitation\ntreatment", palette = "BrBG", labels = c("Drought (D)", "Control", "Wet (W)")) +
  #scale_color_manual(values = fxnltargcols, guide = "none") +
  # scale_fill_manual(name = NULL, 
  #                   values = fxnltargcols, labels = c("Non-native forb", "Non-native grass", "Non-native N-fixer", 
  #                                                     "B. carinatus", "E. californica", "F. microstachys", "N. maculata")) +
  facet_grid(spp ~ herbicide + seedtrt, scales = "free_y") +
  theme_ctw() +
  theme(legend.justification.right = "top",
        plot.margin = unit(c(1, 6, 1, 1), "lines")) +
  labs(x = "Soil x precipitation treatment", y = "Abundance")

# get anova stats
exo_best_anova <- data.frame(car::Anova(exo_allplots_mod5genp)) %>%
  mutate(term = rownames(.)) %>%
  rename_all(function(x) gsub("[.]", "", x))

joint_tests(exo_allplots_mod5genp, by = "spp")

exofig
# manually add anova stats for now
grid::grid.text(label = "Type II Wald Chisq\n seedtrt ***\n spp ***\n herbicide ***\n soil trt n.s.\n precip trt ***\n spp:herbicide ***\n spp:soil trt ***\n spp:precip trt **\n herbicide:soil trt n.s.\n herbicide:precip trt **\n spp:herbicide:soil trt **\n spp:herbicide:precip trt ***", 
                x = unit(0.80, "npc"), y = unit(0.56, "npc"), just = "left")
# add signif terms for spp -- in two columns per row
# non-nat grams
grid::grid.text(label = "seed ***\nherb. ***\nsoil n.s.", 
                x = unit(0.07, "npc"), y = unit(0.69, "npc"), just = "left")
grid::grid.text(label = "ppt ***\nherb:soil n.s.\nherb:ppt ***", 
                x = unit(0.13, "npc"), y = unit(0.69, "npc"), just = "left")
# non-nat forb
grid::grid.text(label = "seed ***\nherb. **\nsoil ***", 
                x = unit(0.07, "npc"), y = unit(0.58, "npc"), just = "left")
grid::grid.text(label = "ppt n.s.\nherb:soil **\nherb:ppt n.s.", 
                x = unit(0.13, "npc"), y = unit(0.58, "npc"), just = "left")
# non-nat n fixers
grid::grid.text(label = "seed ***\nherb. n.s.\nsoil **", 
                x = unit(0.07, "npc"), y = unit(0.3, "npc"), just = "left")
grid::grid.text(label = "ppt ***\nherb:soil n.s.\nherb:ppt n.s.", 
                x = unit(0.13, "npc"), y = unit(0.3, "npc"), just = "left")


# 2c) background tables, contrasts ------
exo_bestcoeffs <- getmod_effects(exo_allplots_mod5genp, "background fxnl glmmTMB")
exo_emlist <- getemmeans(mod_means, contrast(mod_means, "trt.vs.ctrl"), "native seeded glmmTMB")
# best model joint test, arrange by species
write_csv(arrange(data.frame(joint_tests(nb_nats_seeded_mod1nb2, by = "spp")), spp), "Native-analysis/outputs/natglmmTMB_jointtest.csv")
# wald test for best model
write_csv(nats_best_anova, "Native-analysis/outputs/natglmmTMB_anova.csv")
# model coeffs
write_csv(exo_bestcoeffs, "Native-analysis/outputs/exoglmmTMB_coeffs.csv")
# marg means
write_csv(nats_emlist$marginal_means, "Native-analysis/outputs/natglmmTMB_margmeans.csv")
# contrasts
write_csv(nats_emlist$comparisons, "Native-analysis/outputs/natglmmTMB_margmeans.csv")
write_csv(data.frame(emmeans(nb_nats_seeded_mod1nb2, "trt.vs.ctrl" ~ herbicide | spp , type = "response")), "Native-analysis/outputs/natglmmTMB_herbcontrast.csv")
write_csv(data.frame(emmeans(nb_nats_seeded_mod1nb2, "trt.vs.ctrl" ~ nut_trt | spp , type = "response")), "Native-analysis/outputs/natglmmTMB_nutcontrast.csv")
write_csv(data.frame(emmeans(nb_nats_seeded_mod1nb2, "trt.vs.ctrl" ~ ppt_trt | spp , type = "response")), "Native-analysis/outputs/natglmmTMB_pptcontrast.csv")


# -- CTW CUSTOM GLLVM FUNCTIONS ------
# above pcoa and rda suggest maybe 3, but 1 axis (diving hillslope) captures the most variation in data 
bestlv <- function(y = ydat_count, X = NULL, numseq = 1:4, yform = as.formula(~ herbicide + seedtrt), lvform = NA,
                   smat = sDesign_full, rowform = as.formula(~ (1|block/wholeplotID)), seed = myseed){
  criteria <- NULL
  if(is_formula(yform)){
    for(i in numseq){
      print(i) # print to know where you are in loop while waiting
      fiti <- gllvm(y, X, family = "negative.binomial", num.lv = i, sd.errors = FALSE,
                    formula = yform, studyDesign = smat, row.eff = rowform, seed = myseed)
      criteria[i] <- summary(fiti)$AICc
      names(criteria)[i] = i
    }
  }else{
    for(i in numseq){
      print(i) # print to know where you are in loop while waiting
      fiti <- gllvm(y, X, family = "negative.binomial", num.lv.c = i, sd.errors = FALSE,
                    lv.formula = lvform, studyDesign = smat, row.eff = rowform, seed = myseed)
      criteria[i] <- summary(fiti)$AICc
      names(criteria)[i] = i
    }}
  return(criteria)
}


# to pull loadings from gllvm
get_loadings <- function(mod){
  # spp loadings
  sppscores <- getLoadings(mod) %>%
    data.frame() %>%
    rownames_to_column("spp")
  
  # site loadings
  sitescores <- data.frame(getLV(mod)) %>%
    rownames_to_column("rowid")
  # bind xdat
  sitescores <- left_join(sitescores, rownames_to_column(data.frame(mod$lv.X),var = "rowid")) %>%
    left_join(rownames_to_column(sDesign_full, "rowid"))
  # get env scores
  varscores <- coef(mod)$Canonical.coefficients %>%
    data.frame() %>%
    rownames_to_column("vari") %>%
    # clean up labels
    mutate(plot_label = gsub("Herbicided|nut_trt|ppt_trt|seedtrt", "", vari))
  
  # return dat
  templist <- list(sppscores, sitescores, varscores)
  names(templist) <- c("spp", "sites", "varscores")
  return(templist)
}

# plot predicted vs observed to gauge how well fits data
# functions to grab predictions and residuals
getgll_predred <- function(mod){
  preds <- data.frame(predict.gllvm(mod))
  preds <- cbind(rowid = rownames(mod$data), preds) %>%
    data.frame() %>%
    gather(spp, pred, 2:ncol(.)) %>%
    mutate(pred_exp = round(exp(pred),4))
  modat <- mod$data %>%
    rownames_to_column("rowid") %>%
    gather(spp, count, grep("[A-Z]{4}|Forb|Grass", names(.)))
  modres <- data.frame(resid(mod)[["residuals"]]) %>%
    rownames_to_column("rowid") %>%
    gather(spp, res, grep("[A-Z]{4}|Grass|Forb", names(.)))
  
  # put it all together
  alltog <- left_join(modat, preds) %>%
    left_join(modres) %>%
    mutate(diff_ctpred = count - pred_exp)
  return(alltog)
}


# -- FXNL COMMUNITY ORDINATION -----
# > run with ydat_count (all plots)
# what are communiites like by block or hillslope?
# 1. no predictors to show ordination
gllvm_ordination_all.lv4 <- gllvm(y = ydat_fxnl_wide,
                                  #X = xdat, 
                                  family = "negative.binomial", 
                                  num.lv = 4,
                                  # nest subplot if can, but block/wholeplot captures most of the nested similarity
                                  studyDesign = sDesign_full[c("block", "wholeplotID", "subplotID")], 
                                  row.eff = ~ (1|block), # with subplot throws singular fit error, also throws error with wholeplotID nested [if have 3 lvs]
                                  seed = myseed)
# compare AIC for different number of LVs
gllvm_ordination_all.lv2 <- update(gllvm_ordination_all.lv4, num.lv = 2)
gllvm_ordination_all.lv3 <- update(gllvm_ordination_all.lv4, num.lv = 3)
gllvm_ordination_all.lv5 <- update(gllvm_ordination_all.lv4, num.lv = 5)
#gllvm_ordination_all.lv6 <- update(gllvm_ordination_all.lv4, num.lv = 6) # error for this one
AICc(gllvm_ordination_all.lv2, gllvm_ordination_all.lv3, gllvm_ordination_all.lv4, gllvm_ordination_all.lv5)
# [1] 5636.675 5640.367 5653.286 5664.468
# suggests 2 LVs is best.. but we 4 levels of manipulations. it could be treatments don't come in well at functional level
par(mfrow = c(2,2))
lapply(list(gllvm_ordination_all.lv2,gllvm_ordination_all.lv3,gllvm_ordination_all.lv4,gllvm_ordination_all.lv5), function(x) ordiplot(x, rotate = T, jitter = T, biplot = T))
lapply(list(gllvm_ordination_all.lv3,gllvm_ordination_all.lv4,gllvm_ordination_all.lv5), function(x) ordiplot(x, which = c(2,3), rotate = T, jitter = T, biplot = T))
plot.new()
# plot without spp
lapply(list(gllvm_ordination_all.lv2,gllvm_ordination_all.lv3,gllvm_ordination_all.lv4,gllvm_ordination_all.lv5), function(x) ordiplot(x, rotate = T, jitter = T, biplot = F))
lapply(list(gllvm_ordination_all.lv3,gllvm_ordination_all.lv4,gllvm_ordination_all.lv5), function(x) ordiplot(x, which = c(2,3), rotate = T, jitter = T, biplot = F))
# they each are picking up similar things
par(mfrow = c(1,1))

# pull scores for plotting
scores.lv2 <- get_loadings(gllvm_ordination_all.lv2)
scores.lv2$sites <- left_join(scores.lv2$sites, rownames_to_column(xdat, "rowid"))
scores.lv3 <- get_loadings(gllvm_ordination_all.lv3)
scores.lv3$sites <- left_join(scores.lv3$sites, rownames_to_column(xdat, "rowid"))
scores.lv4 <- get_loadings(gllvm_ordination_all.lv4)
scores.lv4$sites <- left_join(scores.lv4$sites, rownames_to_column(xdat, "rowid"))
# plot
ggplot(scores.lv2$sites, aes(LV1, LV2)) +
  geom_jitter( aes(LV1, LV2, col = paste(herbicide,seedtrt)),# shape = paste(herbicide, seedtrt)), 
               width = 0.0, size = 3) +
  scale_color_brewer(palette = "Paired") +
  facet_grid(ppt_trt ~ nut_trt) +
  ggtitle("Unconstrained ordination: 2 LVs")
ggplot(scores.lv3$sites, aes(LV1, LV2)) +
  geom_jitter( aes(LV1, LV2, col = paste(herbicide,seedtrt)),# shape = paste(herbicide, seedtrt)), 
               width = 0.0, size = 3) +
  scale_color_brewer(palette = "Paired") +
  facet_grid(ppt_trt ~ nut_trt) +
  ggtitle("Unconstrained ordination: 3 LVs")
ggplot(scores.lv4$sites, aes(LV1, LV2)) +
  geom_jitter( aes(LV1, LV2, col = paste(herbicide,seedtrt)),# shape = paste(herbicide, seedtrt)), 
               width = 0.0, size = 3) +
  scale_color_brewer(palette = "Paired") +
  facet_grid(ppt_trt ~ nut_trt) +
  ggtitle("Unconstrained ordination: 4 LVs")
# nothing moves on first two axes between lv 3 and lv 4 mods
# compare variance partition
par(mfrow = c(2,2))
lapply(list(gllvm_ordination_all.lv2, gllvm_ordination_all.lv3, gllvm_ordination_all.lv4, gllvm_ordination_all.lv5),
       function(x) plotVarPartitioning(varPartitioning(x), mar = c(7,5,7,5), adj = 0.25, xlab = "", main = "",
                                       las = 2, args.legend = list(horiz = F, cex = 0.75, xjust = 0.9, yjust = -0.4),
                                       col = palette(hcl.colors(10, "Roma"))))
# > 3 axes is the way for unconstrained ordination. Even when have more axes, it's mostly 3 that are picking up variation.
par(mfrow = c(1,1))
# covariance between spp and functional groups
corrplot::corrplot(getResidualCov(gllvm_ordination_all.lv3)$cov, is.corr = F, type = "upper", method = "color")
# > native grass (stipa pulchra mostly) does its own thing

site_ord <- ggplot(scores.lv3$sites) +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey") +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey") +
  geom_point( aes(LV1, LV2, col = paste(herbicide,seedtrt, sep = ", ")),# shape = paste(herbicide, seedtrt)), 
              size = 3) +
  #geom_text(data = subset(scores.lv3$spp, !grepl("NativeGrass", spp)), aes(LV1, LV2, label = spp)) +
  scale_color_brewer(name = ("Herbicide, seeding"), palette = "Paired") +
  labs(x = "LV1 (35.1% spp variance)", y = "LV2 (35.1% spp variance)") +
  facet_grid(ppt_trt ~ nut_trt) +
  theme_ctw() +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.justification.right = "top")

site_ord_legend <- cowplot::get_legend(site_ord)

spp_ord <- ggplot(scores.lv3$sites) +
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey") +
  geom_vline(aes(xintercept = 0), lty = 2, col = "grey") +
  geom_point( aes(LV1, LV2,), color = "transparent",
              size = 3) +
  geom_point(data = subset(scores.lv3$spp), aes(LV1, LV2)) +
  ggrepel::geom_text_repel(data = subset(scores.lv3$spp), aes(LV1, LV2, label = spp),  size = 5) +
  labs(x = "LV1 (35.1% spp variance)", y = "LV2 (35.1% spp variance)") +
  scale_color_brewer(name = ("Herbicide, seeding"), palette = "Paired") +
  theme_ctw()
spp_ord


ord_panel <- plot_grid(site_ord_legend, spp_ord, nrow = 2,
                       rel_heights = c(0.25, .75),labels = list("", "B"),
                       axis = "t",
                       align = "vh")
# put it together
plot_grid(site_ord + theme(legend.position = "none"), ord_panel,
          labels = c("A"),
          rel_widths = c(1.25,0.75),
          nrow = 1,
          align = "vh", axis = "t"
)



# how did communities separate by plot condition? (herb x seeded x env)
gllvm_ordination_all_envpred <- gllvm(y = ydat_fxnl_wide, #[,!names(ydat_fxnl_wide) %in% nats], 
                                      X = xdat, 
                                      family = "negative.binomial", 
                                      #num.lv = 1,
                                      num.lv.c = 4, 
                                      lv.formula = ~ seedtrt + (herbicide * nut_trt * ppt_trt),
                                      #lv.formula  = as.formula( ~  nut_trt * ppt_trt),
                                      # formula = as.formula(~ herbicide + seedtrt), # if have seeedtrt and 3 latent vars, have zero in covar matrix
                                      # nest subplot if can, but block/wholeplot captures most of the nested similarity
                                      studyDesign = sDesign_full[c("block", "wholeplotID", "subplotID")], 
                                      row.eff = ~ (1|block/wholeplotID), # with subplot throws singular fit error, also throws error with wholeplotID nested [if have 3 lvs]
                                      seed = myseed)
par(mfrow = c(1,1))
ordiplot(gllvm_ordination_all_envpred, biplot =T, symbols = T, jitter = T, jitter.amount = 0.5)
summary(gllvm_ordination_all_envpred)
plotVarPartitioning(varPartitioning(gllvm_ordination_all_envpred), mar = c(8,5,2,20), las = 2,
                    adj = 0.25,
                    xlab = "", main = "",
                    args.legend = list(horiz = F, cex = 1, xjust = -.05, yjust = 1),
                    col = palette(hcl.colors(10, "Roma")))
coefplot(gllvm_ordination_all_envpred, which = 1:4, cex.ylab = 2)

bestlv(y = ydat_fxnl_wide, X = xdat, numseq = 1:4, yform = NA, 
       lvform = as.formula(~ herbicide + seedtrt + ppt_trt * nut_trt))
#     1        2        3        4 
# 5365.453 5253.358 5218.548 5184.911 
# get message aboutl selected optimizer not available for model. may need to run manually
# try 5-6 constrained LVs since 4 best
gllvm_ordination_all_envpred_clv5 <- update(gllvm_ordination_all_envpred, num.lv.c = 5)
AICc(gllvm_ordination_all_envpred, gllvm_ordination_all_envpred_clv5) # 4 is better than 5
# confirm 3 higher # <-- get singular fit error?
gllvm_ordination_all_envpred_clv3 <- update(gllvm_ordination_all_envpred, num.lv.c = 3)
gllvm_ordination_all_envpred_clv2 <- update(gllvm_ordination_all_envpred, num.lv.c = 2) # doesn't work with 2 either
# 4 LVs is best
par(mfrow = c(2,2))
ordiplot(gllvm_ordination_all_envpred,which.lvs = c(1,2), biplot =F, symbols = T, jitter = T, jitter.amount = 0.5)
ordiplot(gllvm_ordination_all_envpred,which.lvs = c(1,3), biplot =F, symbols = T, jitter = T, jitter.amount = 0.5)
ordiplot(gllvm_ordination_all_envpred,which.lvs = c(3,2), biplot =F, symbols = T, jitter = T, jitter.amount = 0.5)
ordiplot(gllvm_ordination_all_envpred,which.lvs = c(3,4), biplot =F, symbols = T, jitter = T, jitter.amount = 0.5)
par(mfrow = c(1,1))

varPartitioning(gllvm_ordination_all_envpred)


# use informed lvs
gllvm_ordination_all_envpred_rr <- gllvm(y = ydat_fxnl_wide, #[!grepl("NativeGra", names(ydat_fxnl_wide))], #[,!names(ydat_fxnl_wide) %in% nats], 
                                         X = xdat, 
                                         family = "negative.binomial", 
                                         #num.lv = 1,
                                         num.RR = 4, 
                                         lv.formula = ~ seedtrt + (herbicide * nut_trt * ppt_trt),
                                         #lv.formula  = as.formula( ~  nut_trt * ppt_trt),
                                         # formula = as.formula(~ herbicide + seedtrt), # if have seeedtrt and 3 latent vars, have zero in covar matrix
                                         # nest subplot if can, but block/wholeplot captures most of the nested similarity
                                         studyDesign = sDesign_full[c("block", "wholeplotID", "subplotID")], 
                                         row.eff = ~ (1|block/wholeplotID), # with subplot throws singular fit error, also throws error with wholeplotID nested [if have 3 lvs]
                                         seed = myseed)

par(mfrow = c(2,2))
ordiplot(gllvm_ordination_all_envpred_rr,which.lvs = c(1,2), biplot =F, symbols = T, jitter = T, jitter.amount = 0.5)
ordiplot(gllvm_ordination_all_envpred_rr,which.lvs = c(1,3), biplot =F, symbols = T, jitter = T, jitter.amount = 0.5)
ordiplot(gllvm_ordination_all_envpred_rr,which.lvs = c(3,2), biplot =F, symbols = T, jitter = T, jitter.amount = 0.5)
ordiplot(gllvm_ordination_all_envpred_rr,which.lvs = c(3,4), biplot =F, symbols = T, jitter = T, jitter.amount = 0.5)
par(mfrow = c(1,1))
# > seems similar to when variables are used to inform axes
summary(gllvm_ordination_all_envpred_rr)
summary(gllvm_ordination_all_envpred)



gllvm_ordrr_scores <- get_loadings(gllvm_ordination_all_envpred_rr)
gllvm_ordc_scores <- get_loadings(gllvm_ordination_all_envpred)


plot_grid(
  ggplot(gllvm_ordrr_scores$sites) +
    geom_point( aes(CLV1, CLV2, col = paste(nut_trt,ppt_trt), shape = paste(herbicide, seedtrt)),size = 3) +
    #geom_text(data = subset(gllvm_ordrr_scores$spp, !grepl("NativeGrass", spp)), aes(CLV1, CLV2, label = spp)) +
    geom_segment(data = gllvm_ordrr_scores$varscores, aes(x = 0, xend = CLV1, y = 0, yend = CLV2, group = vari), 
                 arrow = arrow(length = unit(0.25, "cm"))) +
    ggrepel::geom_text_repel(data = gllvm_ordrr_scores$varscores, aes(x = CLV1, y = CLV2, label = plot_label)) +
    scale_shape_manual(values = c(1:8)) +
    ggtitle("Reduced rank ordination"),
  ggplot(gllvm_ordc_scores$sites) +
    geom_jitter( aes(CLV2, CLV3, col = paste(nut_trt, ppt_trt), shape = paste(herbicide, seedtrt)), width = 0, size = 3) +
    geom_segment(data = gllvm_ordc_scores$varscores, aes(x = 0, xend = CLV2, y = 0, yend = CLV3, group = vari), 
                 arrow = arrow(length = unit(0.25, "cm"))) +
    ggrepel::geom_text_repel(data = gllvm_ordc_scores$varscores, aes(x = CLV2, y = CLV3, label = plot_label)) +
    #geom_text(data = subset(gllvm_ordc_scores$spp, !grepl("NativeGrass", spp)), aes(CLV1, CLV2, label = spp)) +
    scale_shape_manual(values = c(1:8)) +
    ggtitle("Constrained LV ordination"),
  nrow = 2
  
  
)
dev.off()
gllvm::plotVarPartitioning(gllvm::varPartitioning(gllvm_ordination_all_envpred_rr))
coefplot(gllvm_ordination_all_envpred_rr,mfrow = c(4,5), xlim.list = rep(list(c(-100,100)),20))

# add fulltrt to ordination scores
gllvm_ordc_scores$sites <- left_join(gllvm_ordc_scores$sites, xdat)
gllvm_ordrr_scores$sites <- left_join(gllvm_ordrr_scores$sites, xdat)
# create cols for fulltrt

fulltrt_cols <- c("darkorange", "darkorange3", "darkorange4", "orchid", "slateblue1", "slateblue4", "palegreen", "turquoise1", "turquoise4")
names(fulltrt_cols ) <- levels(xdat$fulltrt)

# plot arrows, sites, and spp separately
ggplot(gllvm_ordrr_scores$sites) +
  geom_point( aes(CLV1, CLV2, col = fulltrt, shape = paste(herbicide, seedtrt)), size = 3) +
  geom_segment(data = gllvm_ordrr_scores$varscores, aes(x = 0, xend = CLV1, y = 0, yend = CLV2, group = vari), 
               arrow = arrow(length = unit(0.25, "cm"))) +
  scale_color_manual(values = fulltrt_cols) +
  ggrepel::geom_text_repel(data = gllvm_ordrr_scores$varscores, aes(x = CLV1, y = CLV2, label = plot_label)) +
  #geom_text(data = subset(gllvm_ordc_scores$spp, !grepl("NativeGrass", spp)), aes(CLV1, CLV2, label = spp)) +
  scale_shape_manual(values = c(16, 10, 15, 12))

ggplot(gllvm_ordc_scores$sites) +
  geom_hline(aes(yintercept =  0), lty = 2, col = "grey") +
  geom_vline(aes(xintercept =  0), lty = 2, col = "grey") +
  # geom_jitter( aes(CLV1, CLV2, col = fulltrt, shape = paste(herbicide, seedtrt)), width = 0, size = 3) +
  geom_segment(data = gllvm_ordc_scores$varscores, aes(x = 0, xend = CLV1, y = 0, yend = CLV2, group = vari),
               arrow = arrow(length = unit(0.25, "cm"))) +
  #scale_color_manual(values = fulltrt_cols) +
  #scale_shape_manual(values = c(16, 10, 15, 12)) +
  ggrepel::geom_text_repel(data = gllvm_ordc_scores$varscores, aes(x = CLV1, y = CLV2, label = plot_label))  +
  #geom_text(data = subset(gllvm_ordc_scores$spp, !grepl("NativeGrass", spp)), aes(CLV1, CLV2, label = spp))
  
  ggplot(gllvm_ordc_scores$sites) +
  geom_hline(aes(yintercept =  0), lty = 2, col = "grey") +
  geom_vline(aes(xintercept =  0), lty = 2, col = "grey") +
  geom_jitter( aes(CLV1, CLV2, col = fulltrt, shape = paste(herbicide, seedtrt)), width = 0, size = 3) +
  # geom_segment(data = gllvm_ordc_scores$varscores, aes(x = 0, xend = CLV3, y = 0, yend = CLV4, group = vari), 
  #              arrow = arrow(length = unit(0.25, "cm"))) +
  scale_color_manual(values = fulltrt_cols) +
  scale_shape_manual(values = c(16, 10, 15, 12)) +
  # ggrepel::geom_text_repel(data = gllvm_ordc_scores$varscores, aes(x = CLV3, y = CLV4, label = plot_label)) 
  #geom_text(data = subset(gllvm_ordc_scores$spp, !grepl("NativeGrass", spp)), aes(CLV1, CLV2, label = spp))
  
  ggplot(gllvm_ordrr_scores$sites) +
  geom_hline(aes(yintercept =  0), lty = 2, col = "grey") +
  geom_vline(aes(xintercept =  0), lty = 2, col = "grey") +
  geom_point(data = subset(gllvm_ordrr_scores$spp), aes(CLV1, CLV2)) + 
  ggrepel::geom_text_repel(data = subset(gllvm_ordrr_scores$spp), aes(CLV1, CLV2, label = spp), max.overlaps = 20)


ggplot(gllvm_ordrr_scores$sites) +
  geom_hline(aes(yintercept =  0), lty = 2, col = "grey") +
  geom_vline(aes(xintercept =  0), lty = 2, col = "grey") +
  geom_jitter( aes(CLV1, CLV2, col = fulltrt, shape = paste(herbicide, seedtrt)), width = 0.05, height = .05, size = 3) +
  # geom_segment(data = gllvm_ordc_scores$varscores, aes(x = 0, xend = CLV3, y = 0, yend = CLV4, group = vari), 
  #              arrow = arrow(length = unit(0.25, "cm"))) +
  scale_color_manual(values = fulltrt_cols) +
  scale_shape_manual(values = c(16, 10, 15, 12)) +
  facet_grid(ppt_trt ~ nut_trt) +
  theme_ctw()

# -- SEEDED ONLY -----
sDesign_seeded$hillpos <- ifelse(sDesign_seeded$block %in% c(1,2), "downhill", "uphill")
# check distribution of species again (out of 70 plots)
sort(sapply(ydat_count_seeded, function(x) sum(x>0))) # lump anything present in fewer than 7 plots. lump hypochaeris (into group that has greater count)
# > note: APOC is noted in trait dataset as exotic and it is native
# > drop exotic grass since not present enough
ydat_count_seeded_common <- ydat_count_seeded
# add HYGL counts to HYRA
# add CAPY, LASE to ExoticForbNFix
# add APOC to NativeForbNFix
#ydat_count_seeded_common$Hypochaeris <- ydat_count_seeded_common$HYGL + ydat_count_seeded_common$HYRA
ydat_count_seeded_common$NativeForbNfix <- with(ydat_count_seeded_common, NativeForbNfix + APOC)
ydat_count_seeded_common$ExoticForbNfix <- with(ydat_count_seeded_common, ExoticForbNfix + LASE + CAPY + HYRA + HYGL)
ydat_count_seeded_common <- dplyr::select(ydat_count_seeded_common, -c(CAPY, LASE, APOC, HYRA, HYGL, ExoticGrass))
sort(sapply(ydat_count_seeded_common, function(x) sum(x>0))) # looks okay

# -- unconstrained species ordination -----
gllvm_seeded_unconstrained <- gllvm(y = ydat_count_seeded_common,#[!grepl("Hypo", names(ydat_count_seeded_common))],
                                    #X = xdat_seeded, 
                                    family = "negative.binomial", 
                                    #num.lv.c = 2,
                                    num.lv = 4, 
                                    #lv.formula = ~ herbicide * (nut_trt + ppt_trt),
                                    # formula = as.formula(~ herbicide + seedtrt), # if have seeedtrt and 3 latent vars, have zero in covar matrix
                                    # nest subplot if can, but block/wholeplot captures most of the nested similarity
                                    # studyDesign = sDesign_seeded[c("block", "wholeplotID", "subplotID")], 
                                    # row.eff = ~ (1|block/wholeplotID), # with subplot throws singular fit error, also throws error with wholeplotID nested [if have 3 lvs]
                                    seed = myseed)

summary(gllvm_seeded_unconstrained)
ordiplot(gllvm_seeded_unconstrained, biplot = T, symbols = T)

# -- constrained ordination ------
# purpose: within the restored plots
# > what is the effect of treatment on native and neighbor abundances?
gllvm_seeded_r.lv2 <- gllvm(y = ydat_count_seeded_common,#[!grepl("Hypo", names(ydat_count_seeded_common))],
                            X = xdat_seeded, 
                            family = "negative.binomial", 
                            #num.lv.c = 2,
                            num.RR = 2, 
                            lv.formula = ~ herbicide * (nut_trt + ppt_trt),
                            # formula = as.formula(~ herbicide + seedtrt), # if have seeedtrt and 3 latent vars, have zero in covar matrix
                            # nest subplot if can, but block/wholeplot captures most of the nested similarity
                            studyDesign = sDesign_seeded[c("block", "wholeplotID", "subplotID")], 
                            row.eff = ~ (1|block/wholeplotID), # with subplot throws singular fit error, also throws error with wholeplotID nested [if have 3 lvs]
                            seed = myseed)

# try 3 axes
gllvm_seeded_r.lv3 <- update(gllvm_seeded_r.lv2, num.RR = 3)
gllvm_seeded_r.lv4 <- update(gllvm_seeded_r.lv2, num.RR = 4)
AICc(gllvm_seeded_r.lv2, gllvm_seeded_r.lv3, gllvm_seeded_r.lv4) # two seems best but 3 ok
#[1] 7463.544 7487.065 7530.632 # with hypochaeris
# 7385.271 7406.619 7441.368 without hypochaeris
# [1] 7417.701 7439.855 7472.930 # with hypochaeris added to exotic forbs
lmtest::lrtest(gllvm_seeded_r.lv2, gllvm_seeded_r.lv3, gllvm_seeded_r.lv4) # 3 axes different (but not sure what test is doing w/out predictors)
par(mfrow = c(1,3))
lapply(list(gllvm_seeded_r.lv2, gllvm_seeded_r.lv3, gllvm_seeded_r.lv4), function(x) ordiplot(x, jitter = T, symbols = T, biplot = T))
# not 4 axes
ordiplot(gllvm_seeded_r.lv2, which = c(1,2), jitter = T, symbols = T, biplot = T)
ordiplot(gllvm_seeded_r.lv3, which = c(1,2), jitter = T, symbols = T, biplot = T)
ordiplot(gllvm_seeded_r.lv3, which = c(2,3), jitter = T, symbols = T, biplot = T)
# 3 axes, without hypochaeris, seems best
par(mfrow = c(1,2))
lapply(list(gllvm_seeded_r.lv2, gllvm_seeded_r.lv3), function(x) plotVarPartitioning(varPartitioning(x)))
# 3 axes seems best
coefplot(gllvm_seeded_r.lv3, which = 1:9, mfrow = c(2,5))
par(mfrow = c(2,3))
plot(gllvm_seeded_r.lv3)
dev.off()
plotVarPartitioning(varPartitioning(gllvm_seeded_r.lv3), mar = c(10,10,10,10), las = 2,
                    adj = 0.25,
                    xlab = "", main = "",
                    args.legend = list(horiz = F, cex = 0.75, xjust = 0.9, yjust = -0.2),
                    col = palette(hcl.colors(10, "Roma")))



# informed ordination
gllvm_seeded_c.lv2 <- gllvm(y = ydat_count_seeded_common,#[!grepl("Hypo", names(ydat_count_seeded_common))],
                            X = xdat_seeded, 
                            family = "negative.binomial", 
                            num.lv.c = 2,
                            #num.RR = 2, 
                            lv.formula = ~ herbicide * (nut_trt + ppt_trt),
                            # formula = as.formula(~ herbicide + seedtrt), # if have seeedtrt and 3 latent vars, have zero in covar matrix
                            # nest subplot if can, but block/wholeplot captures most of the nested similarity
                            studyDesign = sDesign_seeded[c("block", "wholeplotID", "subplotID")], 
                            row.eff = ~ (1|block/wholeplotID), # with subplot throws singular fit error, also throws error with wholeplotID nested [if have 3 lvs]
                            seed = myseed)

ordiplot(gllvm_seeded_c.lv2, biplot = T)

gllvm_seeded_c.lv3 <- update(gllvm_seeded_c.lv2, num.lv.c = 3) # does not like 3 with lv.formula ("canoncial coeffs are not orthogonal")
gllvm_seeded_c.lv3_mod0 <- update(gllvm_seeded_c.lv2,  lv.formula = ~ herbicide + nut_trt + ppt_trt, num.lv.c = 3) # did not converge with 3
gllvm_seeded_c.lv4_mod0 <- update(gllvm_seeded_c.lv2,  lv.formula = ~ herbicide + nut_trt + ppt_trt, num.lv.c = 4)
gllvm_seeded_c.lv4_mod2 <- update(gllvm_seeded_c.lv2,  lv.formula = ~ herbicide + nut_trt * ppt_trt, num.lv.c = 4)
gllvm_seeded_c.lv2_mod2 <- update(gllvm_seeded_c.lv2,  lv.formula = ~ herbicide + nut_trt * ppt_trt)

par(mfrow = c(2,2))
ordiplot(gllvm_seeded_c.lv4_mod0, biplot = T, which = c(1,2))
ordiplot(gllvm_seeded_c.lv4_mod0, biplot = T, which = c(2,3))
ordiplot(gllvm_seeded_c.lv4_mod0, biplot = T, which = c(3,4)) # odd how species separate from site so much
ordiplot(gllvm_seeded_c.lv4_mod0, biplot = T, which = c(1,4))
# additive model but with 2 lvs
gllvm_seeded_c.lv2_mod0 <- update(gllvm_seeded_c.lv2,  lv.formula = ~ herbicide + nut_trt + ppt_trt)
ordiplot(gllvm_seeded_c.lv2_mod0)
AIC(gllvm_seeded_c.lv2_mod0) # additive better than mod with interaction (lower AIC: 7318.321)
AIC(gllvm_seeded_c.lv2) # AIC 7328.315
AIC(gllvm_seeded_c.lv4_mod0) # 4 lvs better than 2, AIC: 7305.585
summary(gllvm_seeded_c.lv4_mod0)
summary(gllvm_seeded_c.lv2_mod0)
summary(gllvm_seeded_c.lv4_mod2)
summary(gllvm_seeded_c.lv2_mod2)
AIC(gllvm_seeded_c.lv2_mod2) # 7314.922 worse
AIC(gllvm_seeded_c.lv4_mod2) # 7370.138 worst




# ---- variables as predictors ---------
glmod0 <- as.formula(~ herbicide + nut_trt + ppt_trt)
glmod1 <- as.formula(~ herbicide * (nut_trt + ppt_trt))
glmod2 <- as.formula(~ herbicide + (nut_trt * ppt_trt))
glmod3 <- as.formula(~ herbicide * ppt_trt + nut_trt)
glmod4 <- as.formula(~ herbicide * nut_trt + ppt_trt)
glmod5 <- as.formula(~ herbicide * nut_trt * ppt_trt)
glmod6 <- as.formula(~ herbicide * ppt_trt + (nut_trt * ppt_trt))

# test best lv for seeded dat mods
bestlv(y = ydat_count_seeded_common, lvform = NA, X = xdat_seeded, numseq = 1:4, smat = sDesign_seeded, yform = glmod2)
#    1        2        3        4 
# 7783.786 7876.785 7982.025 8164.347 # with less common spp
# 7499.236 7571.034 7653.345 7834.354 # with common only
# > 1 lv is best
bestlv(y = ydat_count_seeded_common,  X = xdat_seeded, numseq = 1:4, smat = sDesign_seeded, yform = glmod1)
#     1        2        3        4 
# 7526.828 7596.674 7675.429 7760.961 
# > 1 best
bestlv(y = ydat_count_seeded_common,  X = xdat_seeded, numseq = 1:4, smat = sDesign_seeded, yform = glmod0)
#    1        2        3        4 
# 7608.197 7680.600 7771.695 7845.798 # with less common spp
# 7409.269 7468.211 7536.093 7850.801 
# > 1 lv is best for simple additive model
bestlv(y = ydat_count_seeded_common,  X = xdat_seeded, numseq = 1:4, smat = sDesign_seeded, yform = glmod5)
#    1        2        3        4 
# 8451.329 8585.101 8725.222 8864.909 # with less common spp
# 8026.290 8129.948 8243.309 8364.310 # with common only
# > 1 lv by far


# additive model only
gllvm_seeded_mod0 <- gllvm(y = ydat_count_seeded_common,
                           X = xdat_seeded, 
                           family = "negative.binomial",
                           num.lv = 1,
                           formula = glmod0,
                           # nest subplot if can, but block/wholeplot captures most of the nested similarity
                           studyDesign = sDesign_seeded[c("block", "wholeplotID", "subplotID")], 
                           row.eff = ~ (1|block/wholeplotID),
                           seed = myseed)
# herb interacts with each env
gllvm_seeded_mod1 <- update(gllvm_seeded_mod0, formula = glmod1)
# hypothesis model
gllvm_seeded_mod2 <- update(gllvm_seeded_mod0, formula = glmod2)
# herb x ppt
gllvm_seeded_mod3 <- update(gllvm_seeded_mod0, formula = glmod3)
# herb x nut
gllvm_seeded_mod4 <- update(gllvm_seeded_mod0, formula = glmod4)
# interactions between all
gllvm_seeded_mod5 <- update(gllvm_seeded_mod0, formula = glmod5)
# env interaction and interaction btwn herbicide and ppt
gllvm_seeded_mod6 <- update(gllvm_seeded_mod0, formula = glmod6)

# compare models
AIC(gllvm_seeded_mod0, gllvm_seeded_mod1, gllvm_seeded_mod2, gllvm_seeded_mod3, gllvm_seeded_mod4, gllvm_seeded_mod5, gllvm_seeded_mod6)
# df      AIC
# gllvm_seeded_mod0 242 7302.798
# gllvm_seeded_mod1 362 7375.526
# gllvm_seeded_mod2 362 7347.933
# gllvm_seeded_mod3 302 7332.909
# gllvm_seeded_mod4 302 7343.230
# gllvm_seeded_mod5 602 7541.312
# gllvm_seeded_mod6 422 7380.426
# > suggests mod 0, 3, 4, 2, then 1. not 5.6 is not too much worse than mod 1
AICc(gllvm_seeded_mod0, gllvm_seeded_mod1, gllvm_seeded_mod2, gllvm_seeded_mod3, gllvm_seeded_mod4, gllvm_seeded_mod5, gllvm_seeded_mod6)
# [1] 7366.133 7526.828 7499.236 7434.752 7445.073 8026.290 7593.314
# > mods 0, 3, 4, 2 should consider. this suggests 0 best ITO AIC penality for small sample size
lmtest::lrtest(gllvm_seeded_mod0, gllvm_seeded_mod1, gllvm_seeded_mod2,
               gllvm_seeded_mod3, gllvm_seeded_mod4, gllvm_seeded_mod6, gllvm_seeded_mod5) # they are all different
# compare non-additive model
lmtest::lrtest(gllvm_seeded_mod1, gllvm_seeded_mod2,
               gllvm_seeded_mod3, gllvm_seeded_mod4, gllvm_seeded_mod6, gllvm_seeded_mod5) # they are all different from 
anova(gllvm_seeded_mod0, gllvm_seeded_mod1, gllvm_seeded_mod2,
      gllvm_seeded_mod3, gllvm_seeded_mod4, gllvm_seeded_mod5, gllvm_seeded_mod6)
# this suggests models 3, 1, and 6 are different. 5 marginal
anova(gllvm_seeded_mod1, gllvm_seeded_mod2, gllvm_seeded_mod3, gllvm_seeded_mod6) # 3 is not necessarily better than 1
# mod 6 is maybe
dev.off()
par(mfrow= c(2,4))
lapply(list(gllvm_seeded_mod0, gllvm_seeded_mod1, gllvm_seeded_mod2, 
            gllvm_seeded_mod3, gllvm_seeded_mod4, gllvm_seeded_mod5, gllvm_seeded_mod6), function(x) plot(x, which =1))
plot.new()
# > mods 1,2,3 are most on the line; mod 5 resids are clumped
lapply(list(gllvm_seeded_mod0, gllvm_seeded_mod1, gllvm_seeded_mod2, 
            gllvm_seeded_mod3, gllvm_seeded_mod4, gllvm_seeded_mod5,gllvm_seeded_mod6), function(x) plot(x, which =2))
plot.new()
# > mods 0, 3, 2, 1 seems okay
lapply(list(gllvm_seeded_mod0, gllvm_seeded_mod1, gllvm_seeded_mod2, 
            gllvm_seeded_mod3, gllvm_seeded_mod4, gllvm_seeded_mod5, gllvm_seeded_mod6), function(x) plot(x, which =3))
plot.new()
# > residuals all look pretty okay, mod 3 has more equal residual range
lapply(list(gllvm_seeded_mod0, gllvm_seeded_mod1, gllvm_seeded_mod2, 
            gllvm_seeded_mod3, gllvm_seeded_mod4, gllvm_seeded_mod5, gllvm_seeded_mod6), function(x) plot(x, which =4))
plot.new()
# 1-3 good. mod 1 might have smaller range overall
lapply(list(gllvm_seeded_mod0, gllvm_seeded_mod1, gllvm_seeded_mod2, 
            gllvm_seeded_mod3, gllvm_seeded_mod4, gllvm_seeded_mod5, gllvm_seeded_mod6), function(x) plot(x, which =5))
plot.new()
# mods 0, 1, 3
lapply(list(gllvm_seeded_mod0, gllvm_seeded_mod1, gllvm_seeded_mod2, 
            gllvm_seeded_mod3, gllvm_seeded_mod4, gllvm_seeded_mod5, gllvm_seeded_mod6), function(x) ordiplot(x, symbols = T, biplot = T))

# compare other distributions for mod 2 to be sure
gllvm_seeded_mod2zinb <- update(gllvm_seeded_mod2,family = "ZINB")
gllvm_seeded_mod2zip <- update(gllvm_seeded_mod2,family = "ZIP")
# plot to compare
par(mfrow = c(2,3))
lapply(list(gllvm_seeded_mod2, gllvm_seeded_mod2zinb, gllvm_seeded_mod2zip), function(x) plot(x, which = c(1,2)))
# not zip model; QQ plot is more in line with zinb but there is a trend in residuals
par(mfrow = c(2,2))
lapply(list(gllvm_seeded_mod2, gllvm_seeded_mod2zinb), function(x) plot(x, which = c(3,4)))
# zero infl is an improvement enough over nb, which is consistent with results elsewhere
par(mfrow = c(1,1))


#creating 2lv mods here so code runs (originally below)
gllvm_seeded_mod0.2lv <- update(gllvm_seeded_mod0, num.lv = 2)
gllvm_seeded_mod1.2lv <- update(gllvm_seeded_mod1, num.lv = 2)
gllvm_seeded_mod2.2lv <- update(gllvm_seeded_mod2, num.lv = 2)
gllvm_seeded_mod3.2lv <- update(gllvm_seeded_mod3, num.lv = 2)
gllvm_seeded_mod6.2lv <- update(gllvm_seeded_mod6, num.lv = 2)

# get predictions and residuals for all models in consideration
# > note_ mod 4 has high AIC
predred_mod0nb <- getgll_predred(gllvm_seeded_mod0)
predred_mod1nb <- getgll_predred(gllvm_seeded_mod1)
predred_mod1nb.2lv <- getgll_predred(gllvm_seeded_mod1.2lv)
predred_mod2nb <- getgll_predred(gllvm_seeded_mod2)
predred_mod2zinb <- getgll_predred(gllvm_seeded_mod2zinb)
predred_mod2nb.2lv <- getgll_predred(gllvm_seeded_mod2.2lv)
predred_mod3nb <- getgll_predred(gllvm_seeded_mod3)
predred_mod4nb <- getgll_predred(gllvm_seeded_mod4)
predred_mod5nb <- getgll_predred(gllvm_seeded_mod5)
predred_mod6nb <- getgll_predred(gllvm_seeded_mod6)
predred_mod6nb.2lv <- getgll_predred(gllvm_seeded_mod6.2lv)


# join then plot
alltog_predred <- rbind(cbind(predred_mod0nb, mod = "mod0"),
                        cbind(predred_mod1nb, mod = "mod1"),
                        cbind(predred_mod1nb.2lv, mod = "mod1.2lv"),
                        cbind(predred_mod2nb, mod = "mod2"),
                        cbind(predred_mod2zinb, mod = "mod2 zinb"),
                        cbind(predred_mod2nb.2lv, mod = "mod2.2lv"),
                        cbind(predred_mod3nb, mod = "mod3"),
                        cbind(predred_mod4nb, mod = "mod4"),
                        cbind(predred_mod5nb, mod = "mod5"),
                        cbind(predred_mod6nb, mod = "mod6"),
                        cbind(predred_mod6nb.2lv, mod = "mod6.2lv"))

ggplot(alltog_predred, aes(spp, res, group = spp)) +
  #geom_boxplot(width = 0.25) +
  #geom_point(pch = 1) +
  geom_hline(aes(yintercept = 0)) +
  stat_summary(aes(col = spp %in% nats), position = position_nudge(x = -0.25)) +
  coord_flip() +
  facet_wrap(~mod, ncol = 4) +
  theme_ctw()
# mod 2 zinb looks worse compared to nb (average resid is more off)
# mod 3 looks best ITO no directionality in residuals for most spp, including natives seeded
# mod 2 and 5 are also ok. it looks like there are tradeoffs for predicting spp no matter what model is chosen bc spp are heterogeneous in response

# plot predicted exponentiated values vs counts: means and SEs
ggplot(subset(alltog_predred, !grepl("zinb", mod)), aes(spp, pred_exp, group = spp)) +
  geom_hline(aes(yintercept = 0)) +
  stat_summary(aes(col = spp %in% nats), position = position_nudge(x = -0.25)) +
  stat_summary(aes(spp, count, group = spp), position = position_nudge(x = 0)) +
  coord_flip() +
  facet_wrap(~mod) +
  theme_ctw()
# > mod 3 might be best?, mod 1 maybe better overall

# plot observations and predictions for all plots: data points
ggplot(alltog_predred, aes(spp, pred_exp, group = spp)) +
  geom_boxplot(width = 0.25) +
  geom_point(aes(col = spp %in% nats)) +
  geom_hline(aes(yintercept = 0)) +
  geom_line(aes(spp, count, group = spp), position = position_nudge(x = -0.5)) +
  geom_point(aes(spp, count, group = spp), position = position_nudge(x = -0.5)) +
  coord_flip() +
  facet_wrap(~mod) +
  theme_ctw()
# 2 or 3

# plot difference
ggplot(subset(alltog_predred, !grepl("zinb", mod)), aes(spp, diff_ctpred, group = spp)) +
  geom_boxplot(width = 0.25) +
  geom_point(aes(col = spp %in% nats), position = position_nudge(x = - 0.5)) +
  geom_hline(aes(yintercept = 0)) +
  coord_flip() +
  facet_wrap(~mod, scales = "free_x") +
  theme_ctw()
# look at overall mean error
sort(sapply(split(alltog_predred$diff_ctpred, alltog_predred$mod), function(x) mean(abs(x)))) 
# mean pulls mod 2 error greater than mod 1 and 3
sort(sapply(split(alltog_predred$diff_ctpred, alltog_predred$mod), function(x) median(abs(x)))) 
# median error lower for mod 2
# mod 6 has relative lower error in mean and median

# group it by spp to compare
sppdiffs <- group_by(alltog_predred, spp, mod) %>%
  reframe(mean_absdiff = mean(abs(diff_ctpred)),
          med_absdiff = median(abs(diff_ctpred))) %>%
  group_by(spp) %>%
  mutate(rankmean = rank(mean_absdiff)) %>%
  ungroup() %>%
  mutate(seeded = spp %in% nats)
# plot ranks for smallest mean abs resid
ggplot(sppdiffs, aes(spp, rankmean, fill = mod)) +
  geom_col(color = "black") +
  coord_flip()
# > mod 5 is best for most BUT not for abundant spp. mod 1 and 2 aren't terrible
# FEMI: mod 1 best, then mod 5, 6, then mod 3
# BRCA: mod 5, then 6, 2, then 3 (ignore zinb)
# ESCA: mod 5, then 6, 4, then 1
# NEMA: mod 5, then mod2 zinb, then mod 2
# > mod 5 is what comes up most conssitently
ggplot(subset(sppdiffs, !grepl("zinb",mod)), aes(spp, mean_absdiff, fill = mod)) +
  geom_jitter(pch = 21, size = 5, width = 0.2, height = 0.1, alpha = .75) +
  coord_flip() +
  scale_fill_brewer(palette = "Spectral") +
  facet_wrap(~spp %in% nats) +
  theme_ctw()
# > mod 6 is relatively low
ggplot(subset(sppdiffs, !grepl("zinb",mod)), aes(spp, med_absdiff, fill = mod)) +
  geom_jitter(pch = 21, size = 5, width = 0.2, height = 0.1, alpha = .75) +
  coord_flip() +
  scale_fill_brewer(palette = "Spectral") +
  facet_wrap(~spp %in% nats) +
  theme_ctw()

# rmspe
sort(sapply(split(alltog_predred$diff_ctpred, alltog_predred$mod), function(x) mean(sqrt(x^2))))
# > mod 5 has lowest error, followed by 1, 3 (but big outliers in residuals pull those up)

AIC(gllvm_seeded_mod0, gllvm_seeded_mod1, gllvm_seeded_mod2, gllvm_seeded_mod2zinb,
    gllvm_seeded_mod3, gllvm_seeded_mod4, gllvm_seeded_mod5, gllvm_seeded_mod6)
# mod 0, then 2
AICc(gllvm_seeded_mod0, gllvm_seeded_mod1, gllvm_seeded_mod2, gllvm_seeded_mod2zinb,
     gllvm_seeded_mod3, gllvm_seeded_mod4, gllvm_seeded_mod5)
# mod 0, then 3, 4 .. 5 has highest score
BIC(gllvm_seeded_mod0, gllvm_seeded_mod1, gllvm_seeded_mod2, gllvm_seeded_mod2zinb,
    gllvm_seeded_mod3, gllvm_seeded_mod4, gllvm_seeded_mod5)
# mod 0 then 3, 4

# plot coefficients
listx <- list(c(-20,20))
coefplot(gllvm_seeded_mod0, mfrow = c(2,3), cex.ylab = 1, xlim.list = rep(listx, 6))
coefplot(gllvm_seeded_mod1, mfrow = c(2,5), cex.ylab = 1, xlim.list = rep(list(c(-30,30)), 10))
coefplot(gllvm_seeded_mod2, mfrow = c(2,5), cex.ylab = 1, xlim.list = rep(list(c(-30,30)), 10))
coefplot(gllvm_seeded_mod3, mfrow = c(2,4), cex.ylab = 1, xlim.list = rep(list(c(-20,20)), 8))
coefplot(gllvm_seeded_mod4, mfrow = c(2,4), cex.ylab = 1, xlim.list = rep(listx, 8))
coefplot(gllvm_seeded_mod5, mfrow = c(3,7), cex.ylab = 1, xlim.list = rep(list(c(-20,20)), 21))
coefplot(gllvm_seeded_mod6, mfrow = c(2,6), cex.ylab = 1, xlim.list = rep(list(c(-30,30)), 12))


# plot variance
par(mfrow = c(2,3))
plotVarPartitioning(varPartitioning(gllvm_seeded_mod0), mar = c(8,3,5,1), las = 2,
                    col = palette(hcl.colors(10, "Roma")))
plotVarPartitioning(varPartitioning(gllvm_seeded_mod1), mar = c(8,3,5,1), las = 2,
                    col = palette(hcl.colors(10, "Roma")))
plotVarPartitioning(varPartitioning(gllvm_seeded_mod2), mar = c(8,3,5,1), las = 2,
                    col = palette(hcl.colors(10, "Roma")))
plotVarPartitioning(varPartitioning(gllvm_seeded_mod3), mar = c(8,3,5,1), las = 2,
                    col = palette(hcl.colors(10, "Roma")))
plotVarPartitioning(varPartitioning(gllvm_seeded_mod4), mar = c(8,3,5,1), las = 2,
                    col = palette(hcl.colors(10, "Roma")))
plotVarPartitioning(varPartitioning(gllvm_seeded_mod5), mar = c(8,3,5,1), las = 2,
                    col = palette(hcl.colors(10, "Roma")))

par(mfrow = c(1,1))
plotVarPartitioning(varPartitioning(gllvm_seeded_mod6), mar = c(8,3,5,1), las = 2,
                    col = palette(hcl.colors(10, "Roma")))
plotVarPartitioning(varPartitioning(gllvm_seeded_mod1), mar = c(8,3,10,1), las = 2,
                    adj = 0.25, xlab = "", main = "",
                    las = 2, args.legend = list(horiz = F, cex = 0.75, xjust = 1, yjust = -0.4),
                    col = palette(hcl.colors(10, "Roma")))
plotVarPartitioning(varPartitioning(gllvm_seeded_mod2), mar = c(8,3,10,1), las = 2,
                    adj = 0.25, xlab = "", main = "",
                    las = 2, args.legend = list(horiz = F, cex = 0.75, xjust = 1, yjust = -0.4),
                    col = palette(hcl.colors(10, "Roma")))
dev.off()


# try with two 2lvs for biplot
gllvm_seeded_mod2.2lv <- update(gllvm_seeded_mod2, num.lv = 2)
gllvm_seeded_mod1.2lv <- update(gllvm_seeded_mod1, num.lv = 2)
gllvm_seeded_mod6.2lv <- update(gllvm_seeded_mod6, num.lv = 2)
par(mfrow = c(1,3))
ordiplot(gllvm_seeded_mod2.2lv, biplot = T, symbols = T)
ordiplot(gllvm_seeded_mod1.2lv, biplot = T, symbols = T) # pretty similar
ordiplot(gllvm_seeded_mod6.2lv, biplot = T, symbols = T) # pretty similar
plotVarPartitioning(varPartitioning(gllvm_seeded_mod2.2lv), mar = c(8,3,10,1), las = 2,
                    adj = 0.25, xlab = "", main = "",
                    las = 2, args.legend = list(horiz = F, cex = 0.75, xjust = 1, yjust = -0.4),
                    col = palette(hcl.colors(10, "Roma")))
plotVarPartitioning(varPartitioning(gllvm_seeded_mod1.2lv), mar = c(8,3,10,1), las = 2,
                    adj = 0.25, xlab = "", main = "",
                    las = 2, args.legend = list(horiz = F, cex = 0.75, xjust = 1, yjust = -0.4),
                    col = palette(hcl.colors(10, "Roma")))
plotVarPartitioning(varPartitioning(gllvm_seeded_mod6.2lv), mar = c(8,3,10,1), las = 2,
                    adj = 0.25, xlab = "", main = "",
                    las = 2, args.legend = list(horiz = F, cex = 0.75, xjust = 1, yjust = -0.4),
                    col = palette(hcl.colors(10, "Roma")))

names(summary(gllvm_seeded_mod6.2lv))$Call
summary(gllvm_seeded_mod6.2lv)$Call
logLik.gllvm(gllvm_seeded_mod0)

# pull model performance stats
get_gllvmstats <- function(mod){
  tempstats <- data.frame(summary(mod)[c("family", "AIC", "AICc", "BIC", "log-likelihood", "df", "num.lv")])
  tempstats$formula <- as.character(mod$formula)[2]
  tempstats$row.eff <- as.character(mod$row.eff)[2]
  tempstats$sigma.lv1 <- summary(mod)$sigma.lv[1]
  if(unique(tempstats$num.lv) == 1){
    tempstats$sigma.lv2 <- NA
  }else{
    tempstats$sigma.lv2 <- summary(mod)$sigma.lv[2]
  }
  # rearrange cols
  tempstats <- dplyr::select(tempstats, family, formula, row.eff, AIC:num.lv, sigma.lv1, sigma.lv2)
  return(tempstats)
}

# notes mods in the environment to summarize
seedmods <- ls()[grep("gllvm_seeded_mod[0-9]",ls())]
gllvm_seeded_modstats <- data.frame()
for(i in seedmods){
  mod <- get(i)
  gllvm_seeded_modstats <- rbind(gllvm_seeded_modstats, cbind(mod = i, get_gllvmstats(mod)))
}

# add prediction error to model stats
sppdiff_sum <- group_by(sppdiffs, mod) %>%
  reframe(rmspe = mean(mean_absdiff),
          median_abserror = median(med_absdiff))


# -- FINAL GLLVM FIGS -----
# order spp by their fxnl grps
sppfac <- data.frame(spp = names(ydat_count_seeded_common)) %>%
  left_join(spplist[c("code4", "species", "fxnl_grp", "nativity")], by = c("spp" = "code4")) %>%
  mutate(nativity = gsub("Exotic", "Non-native", nativity),
         grp = paste(nativity, fxnl_grp),
         grp = ifelse(spp == "NativeForbNfix", "Native Forb/Nfix",
                      ifelse(spp == "ExoticForbNfix", "Non-native Forb/Nfix", grp)),
         grp2 = ifelse(spp %in% nats, "Native seeded", 
                       gsub("Forb/Nfix", "Forb", grp)),
         species = ifelse(is.na(species), grp, species),
         species = gsub("fix", "-fixer", species),
         sppcol = "orchid4")
# manually assign colors
sppfac$sppcol[sppfac$grp2 == "Non-native Grass"] <- "seagreen4"
sppfac$sppcol[sppfac$grp2 == "Non-native N-Fixer"] <- "chocolate3"
sppfac$sppcol[sppfac$grp2 == "Native Forb"] <- "purple2"
sppfac$sppcol[sppfac$grp == "Native Grass"] <- "chartreuse3"
sppfac <- mutate(sppfac, grp3 = ifelse(spp %in% c("ESCA", "NEMA"), "Native Forb (seeded)",
                                       ifelse(spp %in% c("BRCA", "FEMI"), "Native Grass (seeded)", 
                                              gsub("Forb/Nfix", "Forb", grp))))
fxnltargcols <- c("orchid4", "seagreen4", "chocolate3", "seagreen1", "orchid1", "chartreuse3", "purple2")

# go with mod 2 with 2 lv

# -- variance partition plot ----

# make variance plot with spp x sorted by functional group

modvp <- varPartitioning(gllvm_seeded_mod2.2lv)
modvp_exp <- data.frame(modvp$PropExplainedVarSp) %>%
  rename_all(function(x) gsub("Random.effect..", "", x)) %>%
  rownames_to_column("spp") %>%
  gather(vari, varexp, herbicide:ncol(.)) %>%
  # set factors
  mutate(vari = factor(vari, 
                       levels = c("herbicide", "nut_trt", "ppt_trt", "nut_trt.ppt_trt", "LV1", "LV2", "block", "wholeplotID.block"),
                       labels = c("herbicide", "soil amendment", "precipitation", "soil:precipitation", "LV1", "LV2", "RE: block", "RE: plot/block"))) %>%
  # left join functional info
  left_join(dplyr::select(spplist, code4, fxnl_grp, nativity), by = c("spp" = "code4")) %>%
  mutate(nativity = gsub("Exotic", "Non-native", nativity),
         grp = paste(nativity, fxnl_grp),
         grp = ifelse(spp == "NativeForbNfix", "Native Forb/Nfix",
                      ifelse(spp == "ExoticForbNfix", "Non-native Forb/Nfix", grp)),
         grp = ifelse(spp %in% nats, "Seeded", grp),
         grp = factor(grp, 
                      level = c("Seeded", "Native Forb", "Native Forb/Nfix", "Non-native Forb", "Non-native Forb/Nfix", "Non-native N-fixer", "Non-native Grass"),
                      labels = c("Native Seeded", "Native Forb", "", "Non-native Forb", " ", "Non-native N-fixer", "Non-native Grass")),
  )

# reverse levels of vari so herbicide on bottom
modvp_exp$vari <- factor(modvp_exp$vari, levels = rev(levels(modvp_exp$vari)))

# get mean variance to include in legend like in gllvm
meanvar <- group_by(modvp_exp, vari) %>%
  reframe(meanvar_full = mean(varexp),
    meanvar = round(mean(varexp),3),
          medvar = round(median(varexp),3),
          sevar = round(sd(varexp)/sqrt(length(varexp)), 3))

varleg_labels <- paste0(meanvar$vari, "\n(mean prop:", meanvar$meanvar, " ± ", meanvar$sevar, " se)")
names(varleg_labels) <- meanvar$vari
varleg_labels
# variance plot
ggplot(modvp_exp, aes(spp, varexp, fill = vari)) +
  geom_col(col = "black") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_brewer(name = NULL, palette = "BrBG", labels = varleg_labels) + #guide = guide_legend(reverse = T)) +
  facet_grid(~grp, scales = "free", space = "free_x", labeller = label_wrap_gen(width = 10))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.spacing.y = unit(5, "pt"),
        axis.title.y = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold")) +
  labs(x = NULL, y = "Proportion variance explained")


# -- variance for best mod (additive) -----

modvp0 <- varPartitioning(gllvm_seeded_mod0.2lv)
modvp0_exp <- data.frame(modvp0$PropExplainedVarSp) %>%
  rename_all(function(x) gsub("Random.effect..", "", x)) %>%
  rownames_to_column("spp") %>%
  gather(vari, varexp, herbicide:ncol(.)) %>%
  # set factors
  mutate(vari = factor(vari, 
                       levels = c("herbicide", "nut_trt", "ppt_trt", "nut_trt.ppt_trt", "LV1", "LV2", "block", "wholeplotID.block"),
                       labels = c("herbicide", "soil amendment", "precipitation", "soil:precipitation", "LV1", "LV2", "RE: block", "RE: plot/block"))) %>%
  # left join functional info
  left_join(dplyr::select(spplist, code4, fxnl_grp, nativity), by = c("spp" = "code4")) %>%
  mutate(nativity = gsub("Exotic", "Non-native", nativity),
         grp = paste(nativity, fxnl_grp),
         grp = ifelse(spp == "NativeForbNfix", "Native Forb/Nfix",
                      ifelse(spp == "ExoticForbNfix", "Non-native Forb/Nfix", grp)),
         grp = ifelse(spp %in% nats, "Seeded", grp),
         grp = factor(grp, 
                      level = c("Seeded", "Native Forb", "Native Forb/Nfix", "Non-native Forb", "Non-native Forb/Nfix", "Non-native N-fixer", "Non-native Grass"),
                      labels = c("Native Seeded", "Native Forb", "", "Non-native Forb", " ", "Non-native N-fixer", "Non-native Grass")),
  )

# reverse levels of vari so herbicide on bottom
modvp0_exp$vari <- factor(modvp0_exp$vari, levels = rev(levels(modvp0_exp$vari)))

# get mean variance to include in legend like in gllvm
meanvar0 <- group_by(modvp0_exp, vari) %>%
  reframe(meanvar_full = mean(varexp),
          meanvar = round(mean(varexp),3),
          medvar = round(median(varexp),3),
          sevar = round(sd(varexp)/sqrt(length(varexp)), 3))

varleg0_labels <- paste0(meanvar0$vari, "\n(mean prop:", meanvar0$meanvar, " ± ", meanvar0$sevar, " se)")
names(varleg0_labels) <- meanvar0$vari
varleg0_labels
# variance plot
ggplot(modvp0_exp, aes(spp, varexp, fill = vari)) +
  geom_col(col = "black") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_brewer(name = NULL, palette = "BrBG", labels = varleg0_labels) + #guide = guide_legend(reverse = T)) +
  facet_grid(~grp, scales = "free", space = "free_x", labeller = label_wrap_gen(width = 10))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.spacing.y = unit(5, "pt"),
        axis.title.y = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold")) +
  labs(x = NULL, y = "Proportion variance explained")


# -- coeff plot -----
# plots gllvm_seeded_mod2.lv2 but also grab info from others models to compare

# get coeffs
get_gllvmcoeffs <- function(mod){ 
  coefftemp <- data.frame(summary(mod)$Coef.tableX) %>%
    rownames_to_column("term") %>%
    # clean up duplicate periods
    rename_all(function(x) casefold(gsub("[.]+", ".", x))) %>%
    # remove trailing period
    rename_all(function(x) gsub("z.", "z", x))
  intertemp <- data.frame(summary(mod)$Coefficients) %>%
    rownames_to_column("term") %>%
    mutate(term = paste0("Intercept.", term)) %>%
    rename(estimate = Intercept) %>%
    mutate(std.error = NA, zvalue = NA, pr.z = NA)
  # add columns to rbind with coefftemp
  
  # pull confidence interval
  conftemp <- data.frame(confint(mod)) %>% # use base R (confint.gllvm does the same)
    rownames_to_column("term") %>%
    mutate(term = gsub("Xcoef.", "", term))
  # clean up colnames
  names(conftemp) <- gsub("X", "ci", names(conftemp))
  names(conftemp) <- gsub("5..", "5", names(conftemp))
  
  # don't think i care about theta or inv phi so drop those
  # merge to drop
  tempcoef <- rbind(intertemp, coefftemp) %>%
    merge(conftemp) %>%
    mutate(spp = str_extract(term, "[A-Z]{4}|Native.*|Exotic.*"),
           term = gsub(":[A-Z]{4}|:Native.*|:Exotic.*", "", term),
           term = ifelse(grepl("Inter", term), "intercept", term),
           term = gsub("herbicideHerbicided", "herbicide", term))
  return(tempcoef)
}

coefs_mod2seed.2lv <- get_gllvmcoeffs(gllvm_seeded_mod2.2lv)
coefs_mod6seed.2lv <- get_gllvmcoeffs(gllvm_seeded_mod6.2lv)
coefs_mod1seed.2lv <- get_gllvmcoeffs(gllvm_seeded_mod1.2lv)
coefs_mod0seed.2lv <- get_gllvmcoeffs(gllvm_seeded_mod0.2lv)
coefs_mod1seed <- get_gllvmcoeffs(gllvm_seeded_mod1)
coefs_mod2seed <- get_gllvmcoeffs(gllvm_seeded_mod2)
coefs_mod3seed <- get_gllvmcoeffs(gllvm_seeded_mod3)

coefdf <- rbind(cbind(coefs_mod1seed.2lv, mod = "mod1 2lv"),
                cbind(coefs_mod0seed.2lv, mod = "mod0 2lv"),
                cbind(coefs_mod2seed.2lv, mod = "mod2 2lv"),
                cbind(coefs_mod6seed.2lv, mod = "mod6 2lv"),
                cbind(coefs_mod2seed, mod = "mod2 1lv"),
                cbind(coefs_mod1seed, mod = "mod1 1lv"),
                cbind(coefs_mod3seed, mod = "mod3 1lv")) %>%
  # terms got switched in label for mod 5 (ppt before nut), not sure why. mod formula is correct
  mutate(term = recode(term, "nut_trtF:ppt_trtW" = "ppt_trtW:nut_trtF", "nut_trtC:ppt_trtW" = "ppt_trtW:nut_trtC", 
                       "nut_trtF:ppt_trtD" = "ppt_trtD:nut_trtF",  "nut_trtC:ppt_trtD" = "ppt_trtD:nut_trtC"))

ggplot(subset(coefdf, grepl("lv", mod)), aes(spp, estimate)) +
  geom_hline(aes(yintercept = 0)) +
  geom_line()+
  geom_point(aes(col = mod, size = (ci2.5 * ci97.5) >0), position = position_dodge(width = 0.5)) +
  facet_wrap(~term, nrow = 3, scales = "free_x") +
  coord_flip() +
  scale_color_brewer(palette = "Paired") +
  scale_size_discrete(range = c(3,5)) +
  theme_bw() +
  theme_ctw()

ggplot(subset(coefdf, spp %in% nats), aes(spp, estimate)) +
  geom_hline(aes(yintercept = 0)) +
  geom_line()+
  geom_point(aes(col = mod, size = (ci2.5 * ci97.5) >0), position = position_dodge(width = 0.5)) +
  facet_wrap(~term, nrow = 3, scales = "free_x") +
  coord_flip() +
  scale_color_brewer(palette = "Paired") +
  scale_size_discrete(range = c(3,5)) +
  theme_bw() +
  theme_ctw()

# assign spp factor levels
coefs_mod2seed.2lv <- left_join(coefs_mod2seed.2lv, sppfac) %>%
  mutate(spp = as.character(spp))

natforb <- sppfac$spp[sppfac$grp3 == "Native Forb"]
exoforb <- sppfac$spp[sppfac$grp3 == "Non-native Forb"]
exonfix <- sppfac$spp[sppfac$grp3 == "Non-native N-fixer"]
exogram <- sppfac$spp[sppfac$grp3 == "Non-native Grass"]
coefs_mod2seed.2lv$term <- factor(coefs_mod2seed.2lv$term, levels = c("herbicide", "nut_trtF", "nut_trtC", "ppt_trtD", "ppt_trtW",
                                                   "nut_trtF:ppt_trtD", "nut_trtF:ppt_trtW", "nut_trtC:ppt_trtD", "nut_trtC:ppt_trtW"))

coefs_mod2seed.2lv$spp <- factor(coefs_mod2seed.2lv$spp, levels = c(c("BRCA", "ESCA", "FEMI", "NEMA"), natforb, exoforb, exonfix, exogram))
# check levels
levels(coefs_mod2seed.2lv$spp) # looks good
coefs_mod2seed.2lv$spp <-factor(coefs_mod2seed.2lv$spp, levels= rev(levels(coefs_mod2seed.2lv$spp)))
# factor fxnl grps for display
coefs_mod2seed.2lv$grp3 <- factor(coefs_mod2seed.2lv$grp3, levels = c("Native Forb", "Native Forb (seeded)", "Native Grass (seeded)",
                                                                      "Non-native Forb", "Non-native N-fixer", "Non-native Grass"))

ggplot(data = subset(coefs_mod2seed.2lv, term != "intercept"), aes(spp, estimate)) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = (0.5 + length(exogram))), lty = 2, col = "grey") +
  geom_vline(aes(xintercept = (0.5 + length(c(exogram, exonfix)))), lty = 2, col = "grey") +
  geom_vline(aes(xintercept = (0.5 + length(c(exogram, exonfix, exoforb)))), lty = 2, col = "grey") +
  geom_vline(aes(xintercept = (0.5 + length(c(exogram, exonfix, exoforb, natforb)))), lty = 2, col = "grey") +
  # have to plot points before error bar, otherwise 
  geom_point(aes(shape = pr.z < 0.05, col = grp3), size = 3) +
  geom_errorbar(data = subset(coefs_mod2seed.2lv, term != "intercept" & pr.z < 0.05), aes(ymin = ci2.5, ymax = ci97.5), width = 0)+
  #geom_errorbar(data = subset(coefdf, grepl("mod2 2", mod) & term != "intercept" & pr.z < 0.05), aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0)+
  geom_point(aes(shape = pr.z < 0.05, col = grp3), size = 3) +
  labs(y = "Mean effect (log scale), with 95% CI when p<0.05", x = NULL) +
  scale_color_manual(name = NULL, values = c("purple1", "mediumpurple1", "seagreen2", "orchid1", "orange","chartreuse3")) +
  #scale_size_discrete(range = c(3,5)) +
  scale_shape_manual(values = c(1,19)) +
  coord_flip() +
  facet_wrap(~term, nrow = 3) +
  theme_ctw() +
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 14))



# make best mod coefplot (best was simple mod 0 with 1 lv)
coefs_mod0seed.2lv <- left_join(coefs_mod0seed.2lv, sppfac) %>%
  mutate(spp = as.character(spp))

coefs_mod0seed.2lv$term <- factor(coefs_mod0seed.2lv$term, levels = c("herbicide", "nut_trtF", "nut_trtC", "ppt_trtD", "ppt_trtW"),
                                  labels = c("Herbicide", "Fertilizer", "Compost", "Drought", "Wet"))

coefs_mod0seed.2lv$spp <- factor(coefs_mod0seed.2lv$spp, levels = c(c("BRCA", "ESCA", "FEMI", "NEMA"), natforb, exoforb, exonfix, exogram))
# check levels
levels(coefs_mod0seed.2lv$spp) # looks good
coefs_mod0seed.2lv$spp <-factor(coefs_mod0seed.2lv$spp, levels= rev(levels(coefs_mod0seed.2lv$spp)))
# factor fxnl grps for display
coefs_mod0seed.2lv$grp3 <- factor(coefs_mod0seed.2lv$grp3, levels = c("Native Forb", "Native Forb (seeded)", "Native Grass (seeded)",
                                                                      "Non-native Forb", "Non-native N-fixer", "Non-native Grass"))
ggplot(data = subset(coefs_mod0seed.2lv, term != "intercept"), aes(spp, estimate)) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = (0.5 + length(exogram))), lty = 2, col = "grey") +
  geom_vline(aes(xintercept = (0.5 + length(c(exogram, exonfix)))), lty = 2, col = "grey") +
  geom_vline(aes(xintercept = (0.5 + length(c(exogram, exonfix, exoforb)))), lty = 2, col = "grey") +
  geom_vline(aes(xintercept = (0.5 + length(c(exogram, exonfix, exoforb, natforb)))), lty = 2, col = "grey") +
  # have to plot points before error bar, otherwise 
  geom_point(aes(shape = pr.z < 0.05, col = grp3), size = 3) +
  geom_errorbar(data = subset(coefs_mod0seed.2lv, term != "intercept" & pr.z < 0.05), aes(ymin = ci2.5, ymax = ci97.5), width = 0)+
  #geom_errorbar(data = subset(coefdf, grepl("mod2 2", mod) & term != "intercept" & pr.z < 0.05), aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0)+
  geom_point(aes(shape = pr.z < 0.05, col = grp3), size = 3) +
  labs(y = "Mean effect (log scale), with 95% CI when p<0.05", x = NULL) +
  scale_color_manual(name = "Functional group", values = c("purple1", "mediumpurple1", "seagreen2", "orchid1", "orange","chartreuse3")) +
  #scale_size_discrete(range = c(3,5)) +
  scale_shape_manual(values = c(1,19)) +
  coord_flip() +
  facet_wrap(~term, nrow = 2, scales = "free_x") +
  theme_ctw() +
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 14),
        legend.position = "inside",
        legend.justification.inside = c("right", "bottom"))


# for quick comparison
par(mfrow = c(3,1))
plotVarPartitioning(varPartitioning(gllvm_seeded_mod0.2lv), mar = c(8,3,2,19), las = 2,
                    adj = 0.25, xlab = "", main = "",
                    las = 2, args.legend = list(horiz = F, cex = 1.2, xjust = 0, yjust = 1),
                    col = palette(hcl.colors(10, "Roma")))
title(main = "Variance partitioning: Additive model", cex.main = 2)
plotVarPartitioning(varPartitioning(gllvm_seeded_mod3.2lv), mar = c(8,3,2,19), las = 2,
                    adj = 0.25, xlab = "", main = "",
                    las = 2, args.legend = list(horiz = F, cex = 1.2, xjust = 0, yjust = 1),
                    col = palette(hcl.colors(10, "Roma")))
title(main = "Variance partitioning: with herbicide x precipitation interaction", cex.main = 2)
plotVarPartitioning(varPartitioning(gllvm_seeded_mod2.2lv), mar = c(8,3,2,19), las = 2,
                    adj = 0.25, xlab = "", main = "",
                    las = 2, args.legend = list(horiz = F, cex = 1.2, xjust = 0, yjust = 1),
                    col = palette(hcl.colors(10, "Roma")))
title(main = "Variance partitioning: with soil x precipitation interaction", cex.main = 2)



# -- residuals plot ----
# residual ordination plot
mod2lv2_loading <- get_loadings(gllvm_seeded_mod2.2lv)
# join spp info
mod2lv2_loading$spp <- left_join(mod2lv2_loading$spp, sppfac)
mod2lv2_loading$sites <- left_join(mod2lv2_loading$sites, rownames_to_column(xdat_seeded, "rowid"))

# panel plot showing sites then spp
plot_grid(
  ggplot() +
    geom_hline(aes(yintercept = 0), lty= 2, col = "grey") +
    geom_vline(aes(xintercept = 0), lty= 2, col = "grey") +
    geom_point(data = mod2lv2_loading$sites, aes(LV1, LV2, fill = block), pch = 21, size = 4) +
    geom_point(data = mod2lv2_loading$spp, aes(LV1, LV2, col = grp), size = 0.5, color = "transparent") +
    #scale_fill_brewer(name = "Block", palette = "Paired") +
    scale_fill_manual(name = "Block", values = c("white", "grey80", "grey55", "black")) +
    theme_ctw() +
    theme(axis.text = element_text(size = 14), 
          axis.title =element_text(size = 16),
          legend.background = element_rect(fill = "transparent"),
          legend.position = "inside",
          legend.justification.inside = c(0.01, 0.01),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)),
  ggplot() +
    geom_hline(aes(yintercept = 0), lty= 2, col = "grey") +
    geom_vline(aes(xintercept = 0), lty= 2, col = "grey") +
    geom_point(data = mod2lv2_loading$sites, aes(LV1, LV2), size = 1, color = "transparent") +
    #geom_point(data = mod2lv2_loading$spp, aes(LV1, LV2), size = 1) +
    ggrepel::geom_label_repel(data = mod2lv2_loading$spp, aes(LV1, LV2, fill = grp3, label = spp, size = spp %in% nats),
                              #color = "floralwhite",
                              min.segment.length = unit(3, "pt"),
                              # point.padding = unit(3, "pt"), 
                              max.overlaps = 100,
                              #label.padding = unit(3, "pt"), 
                              fontface = "bold"
                              ) +
    # ggrepel::geom_text_repel(data = subset(mod2lv2_loading$spp), aes(LV1, LV2, label = spp, size = spp %in% nats), seed = 13, 
    #                          min.segment.length = unit(10, "pt"),
    #                          box.padding = unit(5, "pt")) +
    # scale_fill_brewer(name = "Block", palette = "Paired") +
    scale_fill_manual(name = NULL, values = c("purple1", "mediumpurple1", "seagreen2", "orchid1", "chartreuse3",  "orange")) +
    scale_size_manual(values = c(`TRUE` = 5, `FALSE`=4), guide = "none") +
    theme_ctw() +
    theme(axis.text = element_text(size = 14), 
          axis.title =element_text(size = 16),
          legend.key.size = unit(1, "pt"),
          legend.background = element_rect(fill = "transparent"),
          legend.position = "inside",
          legend.justification.inside = c(0.01,0.01),
          legend.text = element_text(size = 14)),
  align = "vh", axis = "top", nrow = 2
)




# alternate point plot
ggplot() +
  geom_hline(aes(yintercept = 0), lty= 2, col = "grey") +
  geom_vline(aes(xintercept = 0), lty= 2, col = "grey") +
  geom_point(data = mod2lv2_loading$sites, aes(LV1, LV2), size = 1, color = "transparent") +
  geom_point(data = mod2lv2_loading$spp, aes(LV1, LV2, col = grp3), size = 3) +
  ggrepel::geom_text_repel(data = subset(mod2lv2_loading$spp), aes(LV1, LV2, label = spp, size = spp %in% nats), seed = 13, 
                           min.segment.length = unit(10, "pt"),
                           box.padding = unit(5, "pt")) +
  scale_fill_brewer(name = "Block", palette = "Paired") +
  scale_color_manual(name = NULL, values = c("purple2", "orchid4", "seagreen2", "orchid1", "chartreuse3",  "orange")) +
  scale_size_manual(values = c(`TRUE` = 6, `FALSE`=5), guide = "none") +
  theme_ctw() +
  theme(axis.text = element_text(size = 14), 
        axis.title =element_text(size = 16),
        legend.key.size = unit(1, "pt"),
        legend.background = element_rect(fill = "transparent"),
        legend.position = "inside",
        legend.justification.inside = c(0.01,0.01),
        legend.text = element_text(size = 14))


# pull residual spp correlation
# > KNS says just show seeded spp cor
modrescor_tidy <- getResidualCov(gllvm_seeded_mod2.2lv)$cov %>%
  data.frame() %>%
  rownames_to_column("spp") %>%
  gather(spp_pair, corval, AICA:ncol(.))

ggplot(subset(modrescor_tidy, spp %in% nats), aes(spp_pair, spp, fill = corval)) +
  geom_tile(col = "black") +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  coord_flip()

modrescor <- getResidualCor(gllvm_seeded_mod2.2lv)
modrescor <- modrescor[rownames(modrescor) %in% nats,!colnames(modrescor) %in% nats]

modrescov <- getResidualCov(gllvm_seeded_mod2.2lv)$cov
modrescov <- modrescov[rownames(modrescov) %in% nats,!colnames(modrescov) %in% nats]

library(ComplexHeatmap)
hmap <- ComplexHeatmap::Heatmap(t(modrescor), row_dend_width = unit(2, "cm"))
orderspp <- colnames(modrescor)[row_order(draw(hmap))]
sppfac
labeltest <- with(subset(sppfac, !spp %in% nats), grp[order(factor(spp, levels = colnames(modrescor)))])
labeltest <- gsub("Forb/Nfix", "Forb", labeltest)
coltest <- with(subset(sppfac, !spp %in% nats), sppcol[order(factor(spp, levels = colnames(modrescor)))])
row_ha <- rowAnnotation(border = T, simple_anno_size = unit(1, "cm"),
                        annotation_legend_param = list(title_gp = gpar(fontsize = 14),
                                                       labels_gp = gpar(fontsize = 12)),
                        Group = labeltest, 
                        col = list(Group = c("Non-native Forb" = "orchid1", "Non-native N-fixer" = "orange", "Non-native Grass" = "chartreuse3", "Native Grass" = "seagreen1", "Native Forb"= "purple2"))
)
Heatmap(t(modrescor), row_dend_width = unit(2, "cm"), name = "Residual correlation", left_annotation = row_ha,
        border = T, row_title = "Neighbors", column_title = "Native seeded", row_km = 2, column_km = 2,
        heatmap_legend_param = list(direction = "horizontal", 
                                    title_gp = gpar(fontsize = 14),
                                    labels_gp = gpar(fontsize = 12)))

# repeat for covariance
Heatmap(t(modrescov), row_dend_width = unit(2, "cm"), name = "Residual covariance", left_annotation = row_ha, 
        border = T, row_title = "Neighbors", column_title = "Native seeded", row_km = 2, column_km = 2, 
        heatmap_legend_param = list(direction = "horizontal", 
                                    title_gp = gpar(fontsize = 14),
                                    labels_gp = gpar(fontsize = 12)))






# --residuals plot with additive model ----
# residual ordination plot
mod0lv2_loading <- get_loadings(gllvm_seeded_mod0.2lv)
# join spp info
mod0lv2_loading$spp <- left_join(mod0lv2_loading$spp, sppfac)
mod0lv2_loading$sites <- left_join(mod0lv2_loading$sites, rownames_to_column(xdat_seeded, "rowid"))

# panel plot showing sites then spp
plot_grid(
  ggplot() +
    geom_hline(aes(yintercept = 0), lty= 2, col = "grey") +
    geom_vline(aes(xintercept = 0), lty= 2, col = "grey") +
    geom_point(data = mod0lv2_loading$sites, aes(LV1, LV2, fill = block), pch = 21, size = 4) +
    geom_point(data = mod0lv2_loading$spp, aes(LV1, LV2, col = grp), size = 0.5, color = "transparent") +
    #scale_fill_brewer(name = "Block", palette = "Paired") +
    scale_fill_manual(name = "Block", values = c("white", "grey80", "grey55", "black")) +
    theme_ctw() +
    theme(axis.text = element_text(size = 14), 
          axis.title =element_text(size = 16),
          legend.background = element_rect(fill = "transparent"),
          legend.position = "inside",
          legend.justification.inside = c(0.01, 0.01),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)),
  ggplot() +
    geom_hline(aes(yintercept = 0), lty= 2, col = "grey") +
    geom_vline(aes(xintercept = 0), lty= 2, col = "grey") +
    geom_point(data = mod0lv2_loading$sites, aes(LV1, LV2), size = 1, color = "transparent") +
    #geom_point(data = mod0lv2_loading$spp, aes(LV1, LV2, color = grp3), size = 3) +
    ggrepel::geom_label_repel(data = mod0lv2_loading$spp, aes(LV1, LV2, label = spp, size = spp %in% nats),seed = 9, fill = "transparent", #alpha = 0.75,
                              color = "black",
                              min.segment.length = unit(3, "pt"),
                              # point.padding = unit(3, "pt"),
                              max.overlaps = 100,
                              #label.padding = unit(3, "pt"),
                              fontface = "bold"
    ) +
    ggrepel::geom_label_repel(data = mod0lv2_loading$spp, aes(LV1, LV2, fill = grp3, label = spp, size = spp %in% nats), seed = 9,
                              alpha = 0.75,
                              #color = "floralwhite",
                              min.segment.length = unit(3, "pt"),
                              # point.padding = unit(3, "pt"),
                              max.overlaps = 100,
                              #label.padding = unit(3, "pt"),
                              fontface = "bold"
    ) +
    # ggrepel::geom_text_repel(data = subset(mod0lv2_loading$spp), aes(LV1, LV2, label = spp, size = spp %in% nats), seed = 13,
    #                          min.segment.length = unit(10, "pt"),
    #                          box.padding = unit(5, "pt")) +
    # scale_fill_brewer(name = "Block", palette = "Paired") +
    scale_fill_manual(name = NULL, values = c("purple1", "mediumpurple1", "seagreen2", "orchid1", "chartreuse3",  "orange")) +
    scale_size_manual(values = c(`TRUE` = 5, `FALSE`=4), guide = "none") +
    theme_ctw() +
    theme(axis.text = element_text(size = 14), 
          axis.title =element_text(size = 16),
          legend.key.size = unit(1, "pt"),
          legend.background = element_rect(fill = "transparent"),
          legend.position = "inside",
          legend.justification.inside = c(0.01,0.01),
          legend.text = element_text(size = 14)),
  align = "vh", axis = "top", nrow = 2
)

modrescor0 <- getResidualCor(gllvm_seeded_mod0.2lv)
modrescor0 <- modrescor0[rownames(modrescor0) %in% nats,!colnames(modrescor0) %in% nats]

modrescov0 <- getResidualCov(gllvm_seeded_mod0.2lv)$cov
modrescov0 <- modrescov0[rownames(modrescov0) %in% nats,!colnames(modrescov0) %in% nats]

library(ComplexHeatmap)
hmap0 <- ComplexHeatmap::Heatmap(t(modrescor0), row_dend_width = unit(2, "cm"))
orderspp0 <- colnames(modrescor0)[row_order(draw(hmap0))]
sppfac
labeltest0 <- with(subset(sppfac, !spp %in% nats), grp[order(factor(spp, levels = colnames(modrescor0)))])
labeltest0 <- gsub("Forb/Nfix", "Forb", labeltest0)
coltest0 <- with(subset(sppfac, !spp %in% nats), sppcol[order(factor(spp, levels = colnames(modrescor0)))])
row_ha0 <- rowAnnotation(border = T, simple_anno_size = unit(1, "cm"),
                        annotation_legend_param = list(title_gp = gpar(fontsize = 14),
                                                       labels_gp = gpar(fontsize = 12),
                                                       labels = c("Native Forb", "Non-native Forb", "Non-native N-fixer","Non-native Grass"),
                                                       at = c("Native Forb", "Non-native Forb","Non-native N-fixer","Non-native Grass")
                                                       ),
                        Group = labeltest0, 
                        col = list(Group = c("Non-native Forb" = "orchid1", "Non-native N-fixer" = "orange", "Non-native Grass" = "chartreuse3", "Native Grass" = "seagreen1", "Native Forb"= "purple2"))
)
Heatmap(t(modrescor0), row_dend_width = unit(2, "cm"), name = "Residual correlation", left_annotation = row_ha0,
        border = T, row_title = "Neighbors", column_title = "Native seeded", row_km = 2, column_km = 1,
        heatmap_legend_param = list(direction = "horizontal", 
                                    title_gp = gpar(fontsize = 14),
                                    labels_gp = gpar(fontsize = 12)))

# repeat for covariance
Heatmap(t(modrescov0), row_dend_width = unit(2, "cm"), name = "Residual covariance", left_annotation = row_ha0, 
        border = T, row_title = "Neighbors", column_title = "Native seeded", row_km = 2, column_km = 1, 
        heatmap_legend_param = list(direction = "horizontal",
                                    title_gp = gpar(fontsize = 14),
                                    labels_gp = gpar(fontsize = 12)))



# -- unconstrained ordination pretty plot -----
unconstrained_seeded_loading <- get_loadings(gllvm_seeded_unconstrained) 
unconstrained_seeded_loading$spp <- left_join(unconstrained_seeded_loading$spp, distinct(coefs_mod0seed.2lv, spp, grp2, grp3)) %>%
  mutate(grp3 = factor(grp3, levels = c("Native Forb (seeded)", "Native Grass (seeded)", "Native Forb", "Non-native Forb", "Non-native N-fixer", "Non-native Grass")))
unconstrained_seeded_loading$sites <- left_join(unconstrained_seeded_loading$sites, cbind(rowid = rownames(xdat_seeded), xdat_seeded))

# spp plot
uncon_spp <- 
  ggplot(unconstrained_seeded_loading$spp, aes(LV1, LV2)) +
  geom_hline(aes(yintercept = 0), col = "grey50", lty = 2) +
  geom_vline(aes(xintercept = 0), col = "grey50", lty = 2) +
  geom_point(aes(col = grp3, size = spp %in% nats)) +
  # add nats text to plot
  geom_text(data = subset(unconstrained_seeded_loading$spp, spp %in% nats), aes(label = spp), 
            vjust = -1, hjust = -0.2, nudge_y = -0.01) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(name = NULL, values = c("mediumpurple1", "seagreen2", "purple1", "orchid1", "orange","chartreuse3")) +
  scale_size_manual(name = "Native seeded", values = c(`TRUE` = 5, `FALSE` = 4), guide = "none") +
  theme_ctw() +
  theme(legend.position = "top")

# site plot
uncon_site <- unconstrained_seeded_loading$sites %>%
  mutate(nut_trt = factor(nut_trt, labels = c("Control", "Fertilizer", "Compost")),
         ppt_trt = factor(ppt_trt, labels = c("Control", "Drought", "Wet"))) %>%
  ggplot(aes(LV1, LV2)) +
  geom_hline(aes(yintercept = 0), col = "grey50", lty = 2) +
  geom_vline(aes(xintercept = 0), col = "grey50", lty = 2) +
  geom_jitter(aes(shape = herbicide, fill = factor(block)), size = 4) +# width = 0.005, height = 0.005) +
  # geom_segment(data = informed_ord$varscores, aes(x = 0, xend = CLV1, y = 0, yend = CLV2, group = vari), 
  #              arrow = arrow(length = unit(0.25, "cm"))) +
  # ggrepel::geom_text_repel(data = informed_ord$varscores, aes(label = plot_label)) +
  scale_fill_manual(name = "Block", values = c("white", "grey80", "grey55", "black")) +
  scale_shape_manual(name = NULL, values = c(21, 22)) +
  guides(fill = guide_legend(override.aes = list(pch = 21))) +
  theme_ctw() +
  theme(legend.position = "top") +
  facet_grid(nut_trt ~ ppt_trt)

plot_grid(uncon_site, uncon_spp, align = "v", axis = "t", labels = "AUTO", label_size = 24)

# best informed mod is additive only (just like best predictive model, and 4 LVs)
informed_ord <- get_loadings(gllvm_seeded_c.lv4_mod0)
# pair spp info to spp table
informed_ord$spp <- left_join(informed_ord$spp, distinct(coefs_mod0seed.2lv, spp, grp2, grp3))  %>%
  mutate(grp3 = factor(grp3, levels = c("Native Forb (seeded)", "Native Grass (seeded)", "Native Forb", "Non-native Forb", "Non-native N-fixer", "Non-native Grass")))

informed_spp <- ggplot(informed_ord$spp, aes(CLV1, CLV2)) +
  #geom_point(data = informed_ord$sites, size = 0.25, col = "transparent", pch = 1) +
  geom_hline(aes(yintercept = 0), col = "grey50", lty = 2) +
  geom_vline(aes(xintercept = 0), col = "grey50", lty = 2) +
  geom_point(aes(col = grp3, size = spp %in% nats)) +
  geom_text(data = subset(informed_ord$spp, spp %in% c("BRCA", "ESCA", "FEMI")), aes(label = spp), 
            nudge_y = -0.07, vjust = 1, hjust = 0.35) +
  geom_text(data = subset(informed_ord$spp, spp == "NEMA"), aes(label = spp), 
            nudge_y = 0.07, vjust = 0, hjust = 0.35) +
  # ggrepel::geom_text_repel(data = subset(informed_ord$spp, spp %in% nats), aes(label = spp), seed  = 1,
  #                          point.padding = unit(5, "pt"),
  #                          box.padding = unit(5, "pt")) +
  # geom_segment(data = informed_ord$varscores, aes(x = 0, xend = CLV1, y = 0, yend = CLV2, group = vari), 
  #              arrow = arrow(length = unit(0.25, "cm"))) +
  # ggrepel::geom_text_repel(data = informed_ord$varscores, aes(label = plot_label)) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(name = NULL, values = c("mediumpurple1", "seagreen2", "purple1", "orchid1", "orange","chartreuse3")) +
  scale_size_manual(name = "Native seeded", values = c(`TRUE` = 5, `FALSE` = 4), guide = "none") +
  theme_ctw() +
  theme(legend.position = "top" #,
        #legend.justification.inside = c("right", "top")
        )

ggplot(informed_ord$spp, aes(CLV3, CLV4)) +
  geom_hline(aes(yintercept = 0), col = "grey50", lty = 2) +
  geom_vline(aes(xintercept = 0), col = "grey50", lty = 2) +
  geom_point(aes(col = grp3, size = spp %in% nats)) +
  geom_segment(data = informed_ord$varscores, aes(x = 0, xend = CLV3, y = 0, yend = CLV4, group = vari), 
               arrow = arrow(length = unit(0.25, "cm"))) +
  ggrepel::geom_text_repel(data = informed_ord$varscores, aes(label = plot_label)) +
  scale_color_manual(name = "Functional group", values = c("purple1", "mediumpurple1", "seagreen2", "orchid1", "orange","chartreuse3")) +
  scale_size_manual(name = "Native seeded", values = c(`TRUE` = 5, `FALSE` = 3)) +
  theme_ctw()

informed_sites <-  informed_ord$sites %>%
  mutate(nut_trt = factor(nut_trt, labels = c("Control", "Fertilizer", "Compost")),
         ppt_trt = factor(ppt_trt, labels = c("Control", "Drought", "Wet"))) %>%
  ggplot(aes(CLV1, CLV2)) +
  geom_hline(aes(yintercept = 0), col = "grey50", lty = 2) +
  geom_vline(aes(xintercept = 0), col = "grey50", lty = 2) +
  geom_jitter(aes(shape = herbicide, fill = factor(block)), size = 4) +# width = 0.005, height = 0.005) +
  # geom_segment(data = informed_ord$varscores, aes(x = 0, xend = CLV1, y = 0, yend = CLV2, group = vari), 
  #              arrow = arrow(length = unit(0.25, "cm"))) +
  # ggrepel::geom_text_repel(data = informed_ord$varscores, aes(label = plot_label)) +
  scale_fill_manual(name = "Block", values = c("white", "grey80", "grey55", "black")) +
  scale_shape_manual(name = NULL, values = c(21, 22)) +
  guides(fill = guide_legend(override.aes = list(pch = 21))) +
  theme_ctw() +
  facet_grid(nut_trt ~ ppt_trt) +
  theme(legend.position = "top")


ggplot(informed_ord$sites, aes(CLV3, CLV4)) +
  geom_hline(aes(yintercept = 0), col = "grey50", lty = 2) +
  geom_vline(aes(xintercept = 0), col = "grey50", lty = 2) +
  geom_jitter(aes(shape = herbicide, fill = factor(block)), size = 3, width = 0.005, height = 0.005) +
  scale_fill_manual(name = "Block", values = c("white", "grey80", "grey55", "black")) +
  guides(fill = guide_legend(override.aes = list(pch = 21))) +
  theme_ctw() +
  scale_shape_manual(values = c(21, 22)) +
  facet_grid(nut_trt ~ ppt_trt) +
  theme(legend.position = "top")
# 2 and 3, and 3 and 4 isn't revealing much else

# show gradients
informed_grad <- 
informed_ord$varscores %>%
  mutate(plot_label = recode(plot_label, "C" = "Compost", "D" = "Drought", "W" = "Wet", "F" = "Fertilizer", "herbicide" = "Herbicide")) %>%
  ggplot(aes(CLV1, CLV2)) +
  geom_hline(aes(yintercept = 0), col = "grey50", lty = 2) +
  geom_vline(aes(xintercept = 0), col = "grey50", lty = 2) +
  geom_segment(aes(x = 0, xend = CLV1, y = 0, yend = CLV2, group = vari), 
               arrow = arrow(length = unit(0.25, "cm"))) +
  ggrepel::geom_text_repel(aes(label = plot_label),vjust = 1, nudge_y = 0.005, point.padding = unit(6, "pt"), box.padding = unit(6, "pt")) +
  #geom_text(data = subset(informed_ord$varscores, plot_label ), aes(label = plot_label), nudge_y = 0.005, hjust = 1) +
  theme_ctw()


plot_grid(informed_sites,
          plot_grid(informed_spp, informed_grad, 
                    #align = "h", 
                    #axis = "l",
                    rel_heights = c(1,0.5),
                    nrow = 2,
                    labels = c("B", "C"), 
                    label_y = c(0.85, 1), label_x = -0.025,
                    label_size = 20),
          align = "v", axis = "t",
          labels = c("A"),
          label_size = 20, label_y = 0.9,
          rel_widths = c(1,0.85))

# -- EXTRA CODE, not needed ----
# # -- 4th corner model 
# 
# 
# 
# # seeded plots, only spp with traits available
# gllvm_seeded_4th <- gllvm(y = ydat_count_seeded[,names(ydat_count_seeded) %in% rownames(traitmat_scaled)], 
#                                         X = xdat_seeded,
#                                         TR = traitmat_scaled[rownames(traitmat_scaled) %in% names(ydat_count_seeded),c(reduced_traits)],
#                                         family = "negative.binomial", 
#                                         num.lv = 2, 
#                                         formula = y ~ (herbicide + (nut_trt * ppt_trt)) * (Coarse.root.diameter.mm + Fine.root.specific.length.cm.g + Root.volume.cm3 + Height.cm + LDMC + RMF + SLA.cm2.g),
#                                         #lv.formula = ~ block,
#                                         studyDesign = sDesign_seeded[c("block", "wholeplotID", "subplotID")], 
#                                         #corWithin = T,
#                                         row.eff = ~ (1|block/wholeplotID),
#                                         #  gradient.check = T,
#                                         seed = myseed)
# 
# ordiplot(gllvm_seeded_4th)
# summary(gllvm_seeded_4th)
# coefplot(gllvm_seeded_4th, cex.ylab = 1, mar = c(5,20,2,2))
# gllvm_seeded_4th$fourth.corner
# 
# # 4th but with herbicide interacting with env vars
# gllvm_seeded_4th_herbx <- gllvm(y = ydat_count_seeded[,names(ydat_count_seeded) %in% rownames(traitmat_scaled)], 
#                           X = xdat_seeded,
#                           TR = traitmat_scaled[rownames(traitmat_scaled) %in% names(ydat_count_seeded),c(reduced_traits)],
#                           family = "negative.binomial", 
#                           num.lv = 2, 
#                           formula = y ~ herbicide * (nut_trt + ppt_trt) * (Coarse.root.diameter.mm + Fine.root.specific.length.cm.g + Root.volume.cm3 + Height.cm + LDMC + RMF + SLA.cm2.g),
#                           #lv.formula = ~ block,
#                           studyDesign = sDesign_seeded[c("block", "wholeplotID", "subplotID")], 
#                           corWithin = F,
#                           row.eff = ~ (1|block/wholeplotID),
#                           #  gradient.check = T,
#                           seed = myseed)
# 
# ordiplot(gllvm_seeded_4th_herbx, biplot = T)
# summary(gllvm_seeded_4th_herbx)
# coefplot(gllvm_seeded_4th_herbx, cex.ylab = 1, mar = c(5,20,2,2))
# gllvm_seeded_4th_herbx$X.design
# 
# gllvm_seeded_4th_herbx$fourth.corner %>%
#   data.frame() %>%
#   rownames_to_column(var = "treatment") %>%
#   gather(trait, val, 2:ncol(.)) %>%
#   # mutate(trait = factor(trait#, 
#   #                       #levels = c("Coarse.root.diameter.mm", "Fine.root.specific.length.cm.g", "Root.volume.cm3", "RMF", "LDMC", "SLA.cm2.g", "Height.cm")
#   #                       )) %>%
#   ggplot(aes(treatment, trait)) +
#   geom_tile(aes(fill = val), col = "black") +
#   scale_x_discrete(expand = c(0,0)) +
#   scale_y_discrete(expand = c(0,0)) +
#   scale_fill_distiller(name = "Mean effect", type = "div", direction = 1) +
#   labs(x = NULL, y = "Functional trait") +
#   theme_ctw()
# 
# # additive model onlyy
# gllvm_seeded_4th_add <- gllvm(y = ydat_count_seeded[,names(ydat_count_seeded) %in% rownames(traitmat_scaled)], 
#                                 X = xdat_seeded,
#                                 TR = traitmat_scaled[rownames(traitmat_scaled) %in% names(ydat_count_seeded),c(reduced_traits)],
#                                 family = "negative.binomial", 
#                                 num.lv = 2, 
#                                 formula = y ~ (herbicide + nut_trt + ppt_trt) * (Coarse.root.diameter.mm + Fine.root.specific.length.cm.g + Root.volume.cm3 + Height.cm + LDMC + RMF + SLA.cm2.g),
#                                 #lv.formula = ~ block,
#                                 studyDesign = sDesign_seeded[c("block", "wholeplotID", "subplotID")], 
#                                 #corWithin = T,
#                                 row.eff = ~ (1|block/wholeplotID),
#                                 #  gradient.check = T,
#                                 seed = myseed)
# 
# anova(gllvm_seeded_4th, gllvm_seeded_4th_herbx, gllvm_seeded_4th_add, gllvm_seeded_4th_traitadd)
# ordiplot(gllvm_seeded_4th_add)
# coefplot(gllvm_seeded_4th_add, cex.ylab = 1, mar = c(5,18,2,2))
# corrplot::corrplot(gllvm::getResidualCor(gllvm_seeded_4th_add), order = "hclust", method = "color")
# corrplot::corrplot(gllvm_seeded_4th_add$fourth.corner, is.corr = F, method = "color")
# 
# gllvm_seeded_4th_add$fourth.corner %>%
#   data.frame() %>%
#   rownames_to_column(var = "treatment") %>%
#   gather(trait, val, 2:ncol(.)) %>%
#   mutate(trait = factor(trait, 
#                         levels = c("Coarse.root.diameter.mm", "Fine.root.specific.length.cm.g", "Root.volume.cm3", "RMF", "LDMC", "SLA.cm2.g", "Height.cm"))) %>%
#   ggplot(aes(treatment, trait)) +
#   geom_tile(aes(fill = val), col = "black") +
#   scale_x_discrete(expand = c(0,0)) +
#   scale_y_discrete(expand = c(0,0)) +
#   scale_fill_distiller(name = "Mean effect", type = "div", direction = 1) +
#   labs(x = NULL, y = "Functional trait") +
#   theme_ctw()
# 
#   
# 
# 
# 
# gllvm_seeded_4th_traitadd <- gllvm(y = ydat_count_seeded[,names(ydat_count_seeded) %in% rownames(traitmat_scaled)], 
#                               X = xdat_seeded,
#                               TR = traitmat_scaled[rownames(traitmat_scaled) %in% names(ydat_count_seeded),c(reduced_traits)],
#                               family = "negative.binomial", 
#                               num.lv = 2, 
#                               formula = ~ (herbicide + nut_trt + ppt_trt) + Coarse.root.diameter.mm + Fine.root.specific.length.cm.g + Root.volume.cm3 + Height.cm + LDMC + RMF + SLA.cm2.g,
#                               studyDesign = sDesign_seeded[c("block", "wholeplotID", "subplotID")], 
#                               #corWithin = T,
#                               row.eff = ~ (1|block/wholeplotID),
#                               #  gradient.check = T,
#                               seed = myseed)
# ordiplot(gllvm_seeded_4th_traitadd)
# AIC(gllvm_seeded_4th_traitadd)
# summary(gllvm_seeded_4th_traitadd)
# coefplot(gllvm_seeded_4th_traitadd, cex.ylab = 1, mar = c(5,15,2,2))
# 
# 
# AIC(gllvm_seeded_4th, gllvm_seeded_4th_add, gllvm_seeded_4th_herbx, gllvm_seeded_4th_traitadd)
# anova(gllvm_seeded_4th, gllvm_seeded_4th_add, gllvm_seeded_4th_herbx, gllvm_seeded_4th_traitadd)
# anova(gllvm_seeded_4th_add, gllvm_seeded_4th_traitadd)
# BIC(gllvm_seeded_4th, gllvm_seeded_4th_add, gllvm_seeded_4th_herbx, gllvm_seeded_4th_traitadd)





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


# # species specific effects just for herbicide, + env manips
# exo_allplots_seedadd <- glmmTMB(count ~ seedtrt + spp * (herbicide + nut_trt + ppt_trt) + (1|block/wholeplotID),
#                                 data = fxnlexo_long, family = "nbinom2")
# 
# # 3-way interactions between spp, ppt, and nut trt
# exo_allplots_2envx <- glmmTMB(count ~  spp * (seedtrt + herbicide + nut_trt * ppt_trt) + (1|block/wholeplotID), 
#                               data = fxnlexo_long, family = "nbinom2")
# 
# # interactions with herbicide and env
# exo_allplots_3envx <- glmmTMB(count ~ seedtrt + spp * (herbicide * nut_trt * ppt_trt)+ (1|block/wholeplotID), 
#                               data = fxnlexo_long, family = "nbinom2")
# 
# # all interactions (not expecting this to run?)
# exo_allplots_x <- glmmTMB(count ~  spp * seedtrt * herbicide * nut_trt * ppt_trt + (1|block/wholeplotID), 
#                           data = fxnlexo_long, family = "nbinom2")
# 
# # build model that has 3-way interaction with nut and ppt separately 
# exo_allplots_splitenvx <- glmmTMB(count ~  seedtrt + (spp * herbicide) * (nut_trt + ppt_trt) + (1|block/wholeplotID), 
#                                   data = fxnlexo_long, family = "nbinom2")
# 
# # use sum to one for average effect
# exo_allplots_splitenvx.sum <- glmmTMB(count ~  seedtrt + (spp * herbicide) * (nut_trt + ppt_trt) + (1|block/wholeplotID), 
#                                       data = fxnlexo_long, family = "nbinom2", contrasts = list(spp = "contr.sum"))
# 
# summary(exo_allplots_splitenvx.sum)
# 
# anova(exo_allplots_seedadd, exo_allplots_2envx, exo_allplots_3envx, exo_allplots_splitenvx,exo_allplots_x)
# # seed add with interactions between all or split interaction mods are best
# car::Anova(exo_allplots_seedadd)
# car::Anova(exo_allplots_2envx)
# car::Anova(exo_allplots_3envx)
# car::Anova(exo_allplots_splitenvx)
# car::Anova(exo_allplots_x)
# 
# anova(exo_allplots_splitenvx, exo_allplots_3envx) # can still use the 3-way even if full interaction not signif
# 
# # review coeffs
# summary(exo_allplots_seedadd)
# summary(exo_allplots_splitenvx) # best mod
# summary(exo_allplots_3envx)
# summary(exo_allplots_2envx)
