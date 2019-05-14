library(tidyverse)
library(readr)
library(stringr)

# read in data 
dat <- read_csv("~/Dropbox/USDA-compost/Data/ANPP/ANPP_EnteredData/Compost_ANPP_20190424.csv") %>%
  mutate(block = substr(plot, 1,1),
         treatment = substr(plot, 2,2)) %>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site)) %>%
  mutate(subplot=ordered(subplot, levels = c(D="D", XC="XC", W="W"))) 

# initial visualizations in relation to treatments and by functional group
ggplot(dat, aes(x=treatment, y = dry_wgt_g, fill = fxnl_grp)) + geom_boxplot() + facet_wrap(~subplot)

ggplot(dat, aes(x=treatment, y = dry_wgt_g, color = fxnl_grp)) + geom_point() + facet_grid(site~subplot)

ggplot(dat, aes(x=treatment, y = dry_wgt_g, fill = fxnl_grp)) + geom_boxplot() + facet_grid(~site)

ggplot(dat, aes(x=subplot, y = dry_wgt_g, fill = fxnl_grp)) + geom_boxplot() + facet_grid(~site)

# visualization of total biomass 
dat2 <- dat %>%
  group_by(plot, subplot, block, treatment, site) %>%
  summarize(dry_wgt_g = sum(dry_wgt_g))

ggplot(dat2, aes(x=treatment, y = dry_wgt_g)) + geom_boxplot() + facet_wrap(~subplot)

ggplot(dat2, aes(x=treatment, y = dry_wgt_g, fill = subplot)) + geom_boxplot() 

ggplot(dat2, aes(x=treatment, y = dry_wgt_g)) + geom_boxplot() + facet_grid(~site)

ggplot(dat2, aes(x=subplot, y = dry_wgt_g)) + geom_boxplot() + facet_grid(~site)

# visualization of forb:grass ratios
dat3 <- dat %>%
  select(plot, subplot, fxnl_grp, dry_wgt_g, block, treatment, site) %>%
  spread(fxnl_grp, dry_wgt_g) %>%
  mutate(FG = F/G)

ggplot(dat3, aes(x=treatment, y = FG)) + geom_boxplot() + facet_grid(~site)


ggplot(dat3, aes(x=subplot, y = FG)) + geom_boxplot() + facet_grid(~site)
