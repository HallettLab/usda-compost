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

# for bar graphs 
# visualization of total biomass
dat4 <- dat %>%
  group_by(subplot, treatment) %>%
  summarize(mean_dry_wgt_g = mean(dry_wgt_g), sedry=sd(dry_wgt_g)/sqrt(length(dry_wgt_g)))

ggplot(data=dat4, aes(x=treatment, y=mean_dry_wgt_g, fill=treatment))+
  theme_bw()+
  theme(strip.background = element_blank(), 
        text = element_text(size = 20), 
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  facet_wrap(~subplot)+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None")) +
  geom_errorbar(aes(ymin=mean_dry_wgt_g-sedry, ymax=mean_dry_wgt_g+sedry))+
  labs(x="Amendment", y="Total forage (grams of dry weight)") 

#bar graph visualization of forb and grass biomass
dat5 <- dat %>%
  group_by(treatment, fxnl_grp, site) %>%
  summarize(mean_dry_wgt_g = mean(dry_wgt_g), sedry=sd(dry_wgt_g)/sqrt(length(dry_wgt_g)))

ggplot(data=dat5, aes(x=treatment, y=mean_dry_wgt_g, fill=fxnl_grp))+
  theme_bw()+
  labs(x="Amendment", y="Total forage (grams of dry weight)") +
  theme(strip.background = element_blank(), 
        text = element_text(size = 20), 
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  facet_wrap(~site)+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values = c("darkgoldenrod1",  "darkolivegreen"), guide = guide_legend(title = "Functional \n Group"), labels=c("Forb", "Grass")) +
  geom_errorbar(aes(ymin=mean_dry_wgt_g-sedry, ymax=mean_dry_wgt_g+sedry), position=position_dodge(width=0.9))
