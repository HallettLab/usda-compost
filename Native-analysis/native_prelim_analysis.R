# prelim analysis on native recruitment data
# author(s): CTW
# date create: 2021-05-15

# purpose:
# quickplot for CTW's comps/to check treatment effect

# notes:
# check for treatment again once all may 2021 data in


# -- SETUP ------
rm(list = ls())
# load needed libraries
library(tidyverse)
library(cowplot)
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


# read in cleaned native data
natlong <- read.csv(paste0(datpath, "Native/Native_CleanedData/Compost_Native_LongClean.csv"), strip.white = T, na.strings = na_vals)
natwide <- read.csv(paste0(datpath, "Native/Native_CleanedData/Compost_Native_WideClean.csv"), strip.white = T, na.strings = na_vals)
# read in treatment key
trtkey <- read.csv(paste0(datpath, "Compost_TreatmentKey.csv"),na.strings = na_vals, strip.white = T)
# read in master spp list (lookup table)
spplist <- read.csv(paste0(datpath, "Compost_SppList.csv"), na.strings = na_vals, strip.white = T)

# check that all read in as expected
str(natlong)
str(natwide)
#fix date
natlong$date <- as.Date(natlong$date, format = "%Y-%m-%d")
natwide$date <- as.Date(natwide$date, format = "%Y-%m-%d")


# -- QUICKPLOTS -----
# check for treatment effect of seeding natives
# > differentiate seeded natives from unseeded natives
# > check if herbicide "headstart" advantaged natives
# > native grams v. forbs may respond differentially to ppt trts
# > expect exotics to perform better in amended plots.. unless if amendments helped seeded natives in herbicided plots

# try coarse cover first (exotic grass, exotic forb, native grass unseed, native forb unseed, native grass seeded, native forb seeded)
# > most likely unknown plants are exotic (they are all asters.. could be some agoseris or madia in there, but will lean towards exotics)
natcoarse <- mutate(natlong, nativity = ifelse(nativity == "Unknown", "Exotic", nativity),
                    native_seeded = ifelse(native_seeded == "No", "Background", "Seeded"),
                    coarse_fxnl = paste(native_seeded, nativity, fxnl_grp,sep = " ")) %>%
  group_by(plot, herbicide, fulltrt, block, ppt_trt, nut_trt, fxnl_grp, nativity, native_seeded, coarse_fxnl) %>%
  summarise(totcov = sum(pct_cover),
            spp = str_flatten(unique(code4), collapse = ", "),
            S = length(unique(code4))) %>%
  ungroup() %>%
  # change nut_trt control to XC
  mutate(nut_trt = recode(nut_trt, N = "XC"),
         fulltrt = gsub("N", "XC", fulltrt))

# specify plot cols
seedcols <- c("Background" = "bisque1", "Seeded" = "darkgoldenrod4")
fxncols <- c("Forb" ="mediumpurple2", "Grass" = "mediumseagreen", "N-fixer" = "red4")

# check cover
# > plot for legend
coarse_legend <- ggplot(subset(natcoarse, herbicide == "Non-herbicided"), aes(nut_trt, totcov, group = paste(native_seeded, nut_trt, nativity))) +
  geom_boxplot(aes(fill = native_seeded), outlier.shape = NA, alpha = 0.5) +
  #geom_boxplot(data = subset(natcoarse, herbicide == "Non-herbicided" & nativity == "Native"), aes(nut_trt, totcov, group = paste(native_seeded, nut_trt, nativity), fill = native_seeded), outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(col = fxnl_grp, group = paste(native_seeded, nut_trt, nativity)), position = position_jitterdodge(dodge.width = 0.75, jitter.height = 0, jitter.width = 0.2)) +
  scale_fill_manual(name = NULL, values = seedcols) +
  scale_color_manual(name = NULL, values = fxncols) +
  labs(x = NULL, y = "Total cover") +
  theme(strip.background = element_rect(fill = "transparent"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.spacing = unit(5, "pt"),
        #legend.box.background = element_rect(color = "black", fill = "transparent"),
        #legend.key.size = unit(3, "pt"),
        #strip.text.y = element_blank(),
        legend.box = "horizontal") +
  guides(fill = guide_legend(order =2, keyheight = unit(2, "pt")), color = guide_legend(order = 1, keyheight = unit(1.5, "pt"))) +
  facet_grid(ppt_trt ~ nativity, scales = "free", space = "free")

nh_nat_box <- ggplot(subset(natcoarse, herbicide == "Non-herbicided" & nativity == "Native"), aes(nut_trt, totcov, group = paste(native_seeded, nut_trt, nativity))) +
  geom_boxplot(aes(fill = native_seeded), outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(col = fxnl_grp), position = position_jitterdodge(dodge.width = 0.75, jitter.height = 0, jitter.width = 0.2)) +
  scale_fill_manual(name = NULL, values = seedcols) +
  scale_color_manual(name = NULL, values = fxncols, drop = F) +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(limits = c(0,105), expand = c(0.01,0))+
  theme(strip.background = element_rect(fill = "transparent"),
        legend.position = "none") +
  facet_grid(ppt_trt ~ nativity)

nh_exo_box <- ggplot(subset(natcoarse, herbicide == "Non-herbicided" & nativity == "Exotic"), aes(nut_trt, totcov, group = paste(native_seeded, nut_trt, nativity))) +
  geom_boxplot(aes(fill = native_seeded), outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(col = fxnl_grp), position = position_jitterdodge(dodge.width = 0.75, jitter.height = 0, jitter.width = 0.2)) +
  scale_fill_manual(name = NULL, values = seedcols) +
  scale_color_manual(name = NULL, values = fxncols, drop = F) +
  labs(x = NULL, y = "Total cover") +
  scale_y_continuous(limits = c(0,105), expand = c(0.01,0))+
  #scale_x_discrete(expand = c(0,0.5))+
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text.y = element_blank(),
        legend.position = "none") +
  facet_grid(ppt_trt ~ nativity)

nh_panel <- plot_grid(nh_exo_box, nh_nat_box, rel_widths = c(0.5,1))
leg_template <- get_legend(coarse_legend)

nh_cover <- ggdraw() + draw_plot(nh_panel) + draw_grob(leg_template, x = 0.28, y = 0, hjust = 0, vjust = 0.22)
ggsave(plot = nh_cover, filename = paste0(datpath, "Native/Native_Figures/nonherb_trends.pdf"), width = 6, height = 4)

# remake for herbcided
# > can use same legend
herb_nat_box <- ggplot(subset(natcoarse, herbicide != "Non-herbicided" & nativity == "Native"), aes(nut_trt, totcov, group = paste(native_seeded, nut_trt, nativity))) +
  geom_boxplot(aes(fill = native_seeded), outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(col = fxnl_grp), position = position_jitterdodge(dodge.width = 0.75, jitter.height = 0, jitter.width = 0.2)) +
  scale_fill_manual(name = NULL, values = seedcols) +
  scale_color_manual(name = NULL, values = fxncols, drop = F) +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(limits = c(0,105), expand = c(0.01,0))+
  theme(strip.background = element_rect(fill = "transparent"),
        legend.position = "none") +
  facet_grid(ppt_trt ~ nativity)

herb_exo_box <- ggplot(subset(natcoarse, herbicide != "Non-herbicided" & nativity == "Exotic"), aes(nut_trt, totcov, group = paste(native_seeded, nut_trt, nativity))) +
  geom_boxplot(aes(fill = native_seeded), outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(col = fxnl_grp), position = position_jitterdodge(dodge.width = 0.75, jitter.height = 0, jitter.width = 0.2)) +
  scale_fill_manual(name = NULL, values = seedcols) +
  scale_color_manual(name = NULL, values = fxncols, drop = F) +
  labs(x = NULL, y = "Total cover") +
  scale_y_continuous(limits = c(0,105), expand = c(0.01,0))+
  #scale_x_discrete(expand = c(0,0.5))+
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text.y = element_blank(),
        legend.position = "none") +
  facet_grid(ppt_trt ~ nativity)

herb_panel <- plot_grid(herb_exo_box, herb_nat_box, rel_widths = c(0.5,1))
herb_cover <- ggdraw() + draw_plot(herb_panel) + draw_grob(leg_template, x = 0.28, y = 0, hjust = 0, vjust = 0.22)
herb_cover
ggsave(plot = herb_cover, filename = paste0(datpath, "Native/Native_Figures/herb_trends.pdf"), width = 6, height = 4)


# look at richness
ggplot(subset(natcoarse, herbicide == "Non-herbicided"), aes(nut_trt, S, group = paste(native_seeded, nut_trt, nativity))) +
  geom_boxplot(aes(fill = native_seeded), outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(col = fxnl_grp), position = position_jitterdodge(dodge.width = 0.75, jitter.height = 0, jitter.width = 0.2)) +
  scale_fill_manual(name = NULL, values = seedcols) +
  scale_color_manual(name = NULL, values = fxncols, drop = F) +
  labs(x = NULL, y = NULL) +
  #scale_y_continuous(limits = c(0,105), expand = c(0.01,0))+
  theme(strip.background = element_rect(fill = "transparent"),
        legend.position = "none") +
  facet_grid(ppt_trt ~ nativity)

ggplot(subset(natcoarse, herbicide != "Non-herbicided"), aes(nut_trt, S, group = paste(native_seeded, nut_trt, nativity))) +
  geom_boxplot(aes(fill = native_seeded), outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(col = fxnl_grp), position = position_jitterdodge(dodge.width = 0.75, jitter.height = 0, jitter.width = 0.2)) +
  scale_fill_manual(name = NULL, values = seedcols) +
  scale_color_manual(name = NULL, values = fxncols, drop = F) +
  labs(x = NULL, y = NULL) +
  #scale_y_continuous(limits = c(0,105), expand = c(0.01,0))+
  theme(strip.background = element_rect(fill = "transparent"),
        legend.position = "none") +
  facet_grid(ppt_trt ~ nativity)

# > no strong differences in richness

# plot just exotic richness to be sure
ggplot(subset(natcoarse, nativity == "Exotic"), aes(nut_trt, S, group = paste(native_seeded, nut_trt, fxnl_grp))) +
  stat_summary(aes(col = fxnl_grp)) +
  #geom_boxplot(aes(fill = native_seeded), outlier.shape = NA, alpha = 0.5) +
  #geom_point(aes(col = fxnl_grp), position = position_jitterdodge(dodge.width = 0.75, jitter.height = 0, jitter.width = 0.2)) +
  scale_fill_manual(name = NULL, values = seedcols) +
  scale_color_manual(name = NULL, values = fxncols, drop = F) +
  labs(x = NULL, y = NULL) +
  #scale_y_continuous(limits = c(0,105), expand = c(0.01,0))+
  theme(strip.background = element_rect(fill = "transparent")) +
        #legend.position = "none") +
  facet_grid(ppt_trt ~ herbicide) # nope!

