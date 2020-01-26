library(tidyverse)

# specify dropbox pathway (varies by user -- ctw can tweak this later)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/Competition/Competition_EnteredData")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/Competition/Competition_EnteredData"
}else{
  ## LMH and AS
  datpath <- "~/Dropbox/USDA-compost/Data/Competition/Competition_EnteredData"
}

cover_master_wide <- read_csv(paste0(datpath,"Cover_CleanedData/Compost_Cover_WideClean.csv"))
cover_master_long <- read_csv(paste0(datpath,"Cover_CleanedData/Compost_Cover_LongClean.csv"))

dat <- cover_master_wide %>%
  tbl_df() %>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site)) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W"))) 
  

# variations on a theme with grass vs forb visualization
# compost seems to help grass at low site, marginally hurt at high site
ggplot(dat, aes(y = pct_grass, x=interaction(ppt_trt))) + geom_boxplot() + facet_wrap(~nut_trt)
ggplot(dat, aes(y = pct_grass, x=ppt_trt, color = site)) + geom_point() + facet_wrap(~nut_trt)
ggplot(dat, aes(y = pct_grass, x=ppt_trt, color = site)) + geom_boxplot() + facet_grid(site~nut_trt)
ggplot(dat, aes(y = pct_grass, x=nut_trt, color = site)) + geom_boxplot() + facet_grid(~site)
ggplot(dat, aes(y = pct_grass, x=ppt_trt, color = site)) + geom_boxplot() + facet_grid(~site)

# compost seems to help forb at high site, *marginally* hurt at low (makes sense, grass inverse)
# overall seems to like wet (maybe because of N-fixers?)
ggplot(dat, aes(y =  pct_forb, x=interaction(ppt_trt))) + geom_boxplot() + facet_wrap(~nut_trt)
ggplot(dat, aes(y = pct_forb, x=ppt_trt, color = site)) + geom_point() + facet_wrap(~nut_trt)
ggplot(dat, aes(y = pct_forb, x=ppt_trt, color = site)) + geom_boxplot() + facet_grid(site~nut_trt)
ggplot(dat, aes(y = pct_forb, x=nut_trt, color = site)) + geom_boxplot() + facet_grid(~site)
ggplot(dat, aes(y = pct_forb, x=ppt_trt, color = site)) + geom_boxplot() + facet_grid(~site)


# hordeum seems to <3 compost, but sample size is low because not in the lower site
ggplot(dat, aes(y = HOMU, x=ppt_trt)) + geom_boxplot() + facet_grid(~nut_trt)
ggplot(dat, aes(y = HOMU, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = HOMU, x=nut_trt)) + geom_boxplot() + facet_grid(~site)

# lomu seems to like compost, and is in both locations
ggplot(dat, aes(y = LOMU, x=interaction(ppt_trt))) + geom_boxplot() + facet_wrap(~nut_trt)
ggplot(dat, aes(y = LOMU, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = LOMU, x=nut_trt)) + geom_boxplot() + facet_grid(~site)


# avba doesn't seem to like the compost much, especially under wet conditions
# not much at the lower sites
ggplot(dat, aes(y = AVBA, x=interaction(ppt_trt))) + geom_boxplot() + facet_wrap(~nut_trt)
ggplot(dat, aes(y = AVBA, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = AVBA, x=nut_trt)) + geom_boxplot() + facet_grid(~site)
ggplot(dat, aes(y = AVBA, x=ppt_trt)) + geom_boxplot() + facet_grid(~site)



# erbo doesn't seem fussed really
ggplot(dat, aes(y = ERBO, x=interaction(ppt_trt))) + geom_boxplot() + facet_wrap(~nut_trt)
ggplot(dat, aes(y = ERBO, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = ERBO, x=nut_trt)) + geom_boxplot() + facet_grid(~site)
ggplot(dat, aes(y = ERBO, x=ppt_trt)) + geom_boxplot() + facet_grid(~site)


# trhi likes amendments - of both kinds, and especially when not in drought
ggplot(dat, aes(y = TRHI, x=interaction(ppt_trt))) + geom_boxplot() + facet_wrap(~nut_trt)
ggplot(dat, aes(y = TRHI, x=interaction(ppt_trt))) + geom_boxplot() + facet_grid(site~nut_trt)
ggplot(dat, aes(y = TRHI, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = TRHI, x=nut_trt)) + geom_boxplot() + facet_grid(~site)
ggplot(dat, aes(y = TRHI, x=ppt_trt)) + geom_boxplot() + facet_grid(~site)


# trsu likes amendments - of both kinds
ggplot(dat, aes(y = TRSU, x=interaction(ppt_trt))) + geom_boxplot() + facet_wrap(~nut_trt)
ggplot(dat, aes(y = TRSU, x=interaction(ppt_trt))) + geom_boxplot() + facet_grid(site~nut_trt)
ggplot(dat, aes(y = TRSU, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = TRSU, x=nut_trt)) + geom_boxplot() + facet_grid(~site)
ggplot(dat, aes(y = TRSU, x=ppt_trt)) + geom_boxplot() + facet_grid(~site)


# together trifoliums like amendments (especially straight fertilizer, especially at low site)
ggplot(dat, aes(y = (TRSU + TRHI), x=interaction(ppt_trt))) + geom_boxplot() + facet_wrap(~nut_trt)
ggplot(dat, aes(y = (TRSU + TRHI), x=interaction(ppt_trt))) + geom_boxplot() + facet_grid(site~nut_trt)
ggplot(dat, aes(y = (TRSU + TRHI), x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = (TRSU + TRHI), x=nut_trt)) + geom_boxplot() + facet_grid(~site)
ggplot(dat, aes(y = (TRSU + TRHI), x=ppt_trt)) + geom_boxplot() + facet_grid(~site)

# other grasses pretty low numbers
ggplot(dat, aes(y = BRHO, x=interaction(ppt_trt))) + geom_boxplot() + facet_wrap(~nut_trt)
ggplot(dat, aes(y = VUBR, x=interaction(ppt_trt))) + geom_boxplot() + facet_wrap(~nut_trt)


# look at diversity per trt (does compost * wet supress diversity?)
# richness
cover_master_long %>%
  subset(unknown == 0) %>%
  select(block:ppt_trt, code4) %>%
  group_by(block, ppt_trt, nut_trt) %>%
  summarise(S = length(unique(code4))) %>%
  ggplot(aes(ppt_trt, S)) +
  geom_boxplot() +
  geom_jitter(aes(col = as.factor(block)), width = 0.1) +
  labs(x = "Precipitation treatment", y = "Number of unique species", 
       title = "Species richness: more variability under drought, similar range within compost",
       subtitle = "Note: block 1-2 = forb dominant; blocks 3-4 = grass dominant") +
  facet_wrap(~nut_trt, labeller = labeller(nut_trt = c(C = "Compost", F = "Fertilizer", N = "Control"))) +
  scale_color_viridis_d(name = "Block") +
  theme_bw() +
  theme(plot.title = element_text(size = 12))

ggsave(paste0(datpath,"Cover_Figures/spprichness_bytrts.pdf"), 
       width = 6, height = 5, units = "in", scale = 1.1)

forage_outcome_grass <- subset(cover_master_long, fxnl_grp == "Grass" & lubridate::month(date) == 5) %>%
  mutate(desire = ifelse(grepl("madriten|diandr|Taen|Hord", species), "Weedy", "Desirable")) %>%
  group_by(block, nut_trt, ppt_trt, desire) %>%
  summarise(totcov = sum(pct_cover),
            S = length(unique(code4))) %>%
  ungroup()

ggplot(forage_outcome_grass, aes(ppt_trt, S)) +
  geom_boxplot(aes(fill = ppt_trt)) +
  #geom_point(aes(col = block)) +
  facet_grid(desire~nut_trt, scales = "free")

ggplot(forage_outcome_grass, aes(ppt_trt, totcov)) +
  geom_boxplot() +
  geom_jitter(aes(col = as.factor(block)), size = 2, width = 0.2) +
  labs(y = "Total cover (%)", x = "Precipitation treatment",
       title = "May 2019 grass cover, arrayed by nutrient treatment and forage desirability",
       subtitle = "Note: blocks 1-2 = forb dominant; blocks 3-4 = grass dominant") +
  scale_color_viridis_d(name = "Block") +
  facet_grid(desire~nut_trt, scales = "free", labeller = labeller(nut_trt = c(C = "Compost", F = "Fertilizer", N = "Control"))) +
  theme_bw() +
  theme(plot.title = element_text(size = 12))
  
ggsave(paste0(datpath,"Cover_Figures/peak_gramcover_desirability.pdf"), 
       width = 6, height = 5, units = "in", scale = 1.1)
