library(tidyverse)

# set path to compost data (dependent on user)
#datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/Cover/" #ctw's path
datpath <- "~/Dropbox/USDA-compost/Data/Cover/" # should work for LMH and AS

cover_master_wide <- read_csv(paste0(datpath,"Cover_CleanedData/Compost_Cover_WideClean.csv"))

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
