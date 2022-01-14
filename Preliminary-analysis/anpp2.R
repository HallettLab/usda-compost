theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 14),
              strip.text= element_text(size = 14), 
              axis.text = element_text(size = 12))

# read in data
dat <- read_csv("~/Dropbox/USDA-compost/Data/ANPP/ANPP_CleanedData/Compost_ANPP_clean.csv") %>%
  mutate(site = "Moderate",
         site = ifelse(block == 1 | block == 2, "Heavy", site))
# initial visualizations in relation to nut_trts and by functional group
ggplot(dat, aes(x=nut_trt, y = dry_wgt_g, fill = fxnl_grp)) + geom_boxplot() + facet_wrap(~ppt_trt)
ggplot(dat, aes(x=nut_trt, y = dry_wgt_g, color = fxnl_grp)) + geom_point() + facet_grid(site~ppt_trt)
ggplot(dat, aes(x=nut_trt, y = dry_wgt_g, fill = fxnl_grp)) + geom_boxplot() + facet_grid(~site)
ggplot(dat, aes(x=ppt_trt, y = dry_wgt_g, fill = fxnl_grp)) + geom_boxplot() + facet_grid(~site)
# visualization of total biomass
dat2 <- dat %>%
  group_by(plot, ppt_trt, block, nut_trt, site, yr, clip_event) %>%
  summarize(dry_wgt_g = sum(dry_wgt_g)) %>% mutate(ANPP=dry_wgt_g*16) %>% mutate(lbs=dry_wgt_g*16*8.92179)
f2a<-ggplot(subset(dat2, clip_event=2), aes(x=nut_trt, y = lbs, fill=nut_trt)) + geom_boxplot()+ # facet_wrap(~ppt_trt)+
  #theme_bw()+
  #facet_wrap(~clip_event)+
  labs(x="", y="Forage (lbs of dry weight per acre)") +
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))
f2a
f2b<-ggplot(subset(dat2, clip_event=2), aes(x=as.factor(yr), y = lbs, fill=nut_trt)) + geom_boxplot()+ # facet_wrap(~ppt_trt)+
  #theme_bw()+
  #facet_wrap(~clip_event)+
  labs(x="Year", y="") +
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))
f2b
f2c<-ggplot(subset(dat2, clip_event=2), aes(x=ppt_trt, y = lbs, fill = nut_trt)) + geom_boxplot() +
  #theme_bw()+
  labs(x="", y="") +
  #facet_wrap(~clip_event)+
  scale_x_discrete(labels=c("D" = "Dry", "XC" = "Ambient", "W" = "Wet"))+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))
f2c
library(cowplot)

p_col<-plot_grid(
  f2b + theme(legend.position="none"), f2c + theme(legend.position="none"),
  labels = c("B","C"), ncol = 1
)

legend <- get_legend(
  # create some space to the left of the legend
  f2a + theme(legend.box.margin = margin(0, 0, 0, 12))
)

f2<-plot_grid(
  f2a + theme(legend.position="none"), p_col, legend,
  labels = c("A", "", ""), ncol = 3)
f2

ggdraw(f2a) +
  draw_plot(f2b, .15, .6, .3, .3) +
  draw_plot(f2c, 0.5,0.6, 0.3,0.3)+
  draw_plot_label(
    c("A", "B", "C"),
    c(0, 0.15, 0.5),
    c(1, 0.92, 0.92),
    size = 14
  )

m1<-lme(log(lbs) ~nut_trt*ppt_trt, random=~1|site/yr, subset(dat2, clip_event=2), na.action=na.exclude)
summary(m1)
anova(m1)
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))

LS1<-lsmeans(m1, ~nut_trt, by="ppt_trt")
contrast(LS1, "pairwise")
#significant effect of nut_trt

ggplot(dat2, aes(x=nut_trt, y = dry_wgt_g, fill=nut_trt)) + geom_boxplot() + facet_grid(clip_event~site)+
  theme_bw()+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))

dat2$nut_trt_f = factor(dat2$nut_trt, levels=c('C','N','F'))
f3a<-ggplot(subset(dat2, clip_event=2), aes(x=nut_trt_f, y = lbs, fill=nut_trt)) + geom_boxplot() + facet_grid(~site)+
  #theme_bw()+
  labs(x="",y="Forage (lbs of dry weight per acre)")+
  scale_x_discrete(labels=c("C" = "Compost", "N" = "None", "F" = "Fertilizer"))+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))
f3a

ggplot(dat2, aes(x=ppt_trt, y = dry_wgt_g, fill=nut_trt)) + geom_boxplot() + facet_grid(clip_event~site)
# visualization of forb:grass ratios
dat3 <- dat %>%
  dplyr::select(plot, ppt_trt, fxnl_grp, dry_wgt_g, block, nut_trt, site, yr, clip_event) %>%
  spread(fxnl_grp, dry_wgt_g) %>%
  mutate(FG = Forb/Grass)
ggplot(dat3, aes(x=nut_trt, y = FG)) + geom_boxplot() + facet_grid(~site)
ggplot(dat3, aes(x=ppt_trt, y = FG)) + geom_boxplot() + facet_grid(~site)
# for bar graphs
# visualization of total biomass
dat4 <- dat2 %>%
  group_by(ppt_trt, nut_trt, site) %>%
  summarize(mean_dry_wgt_g = mean(dry_wgt_g), sedry=sd(dry_wgt_g)/sqrt(length(dry_wgt_g)))
ggplot(data=dat4, aes(x=nut_trt, y=mean_dry_wgt_g, fill=nut_trt))+
  theme_bw()+
  ylim(0,60)+
  theme(strip.background = element_blank(),
        text = element_text(size = 20),
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~ppt_trt*site)+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None")) +
  geom_errorbar(aes(ymin=mean_dry_wgt_g-sedry, ymax=mean_dry_wgt_g+sedry))+
  labs(x="Amendment", y="Biomass (g of dry weight per 0.25m^2)")
#bar graph visualization of forb and grass biomass
dat5 <- dat %>%
  group_by(nut_trt, fxnl_grp, site, ppt_trt) %>%
  summarize(mean_dry_wgt_g = mean(dry_wgt_g, na.rm = TRUE), sedry=(sd(dry_wgt_g, na.rm=T)/sqrt(length(dry_wgt_g))))
ggplot(data=dat5, aes(x=nut_trt, y=mean_dry_wgt_g, fill=fxnl_grp))+
  theme_bw()+
  ylim(0,50)+
  labs(x="Amendment", y="Total forage (grams of dry weight)") +
  theme(strip.background = element_blank(),
        text = element_text(size = 20),
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~site*ppt_trt)+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values = c("darkgoldenrod1",  "darkolivegreen"), guide = guide_legend(title = "Functional \n Group"), labels=c("Forb", "Grass")) +
  geom_errorbar(aes(ymin=mean_dry_wgt_g-sedry, ymax=mean_dry_wgt_g+sedry), position=position_dodge(width=0.9))
dat6 <- dat %>%
  group_by(nut_trt, fxnl_grp, site) %>%
  summarize(mean_dry_wgt_g = mean(dry_wgt_g, na.rm = TRUE), sedry=(sd(dry_wgt_g, na.rm=T)/sqrt(length(dry_wgt_g))))
ggplot(data=dat6, aes(x=nut_trt, y=mean_dry_wgt_g, fill=fxnl_grp))+
  theme_bw()+
  ylim(0,50)+
  labs(x="Amendment", y="Total forage (grams of dry weight)") +
  theme(strip.background = element_blank(),
        text = element_text(size = 20),
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~site)+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values = c("darkgoldenrod1",  "darkolivegreen"), guide = guide_legend(title = "Functional \n Group"), labels=c("Forb", "Grass")) +
  geom_errorbar(aes(ymin=mean_dry_wgt_g-sedry, ymax=mean_dry_wgt_g+sedry), position=position_dodge(width=0.9))
#test nut_trt effects on biomass
m1<-lme(dry_wgt_g ~nut_trt*site, random=~1|ppt_trt, dat2, na.action=na.exclude)
summary(m1)
anova(m1)
r.squaredGLMM(m1) #24% of variation explained by fixed effects, 24% by whole model (spatial variation?)
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))
#normally distributed, continue
LS1<-lsmeans(m1, ~nut_trt*site)
contrast(LS1, "pairwise")
#normally distributed, continue
LS1<-lsmeans(m1, ~nut_trt*site, by=site)
#normally distributed, continue
LS1<-lsmeans(m1, ~nut_trt*site, by="site")
contrast(LS1, "pairwise")
m2<-lme(dry_wgt_g ~nut_trt*fxnl_grp, random=~1|ppt_trt, dat, na.action=na.exclude)
summary(m2)
anova(m2)
r.squaredGLMM(m2) #24% of variation explained by fixed effects, 24% by whole model (spatial variation?)
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))
#normally distributed, continue
LS1<-lsmeans(m2, ~nut_trt*fxnl_grp)
#normally distributed, continue
LS1<-lsmeans(m2, ~nut_trt*fxnl_grp, by="fxnl_grp")
contrast(LS1, "pairwise")
m2<-lme(dry_wgt_g ~nut_trt*fxnl_grp, random=~1|yr, dat, na.action=na.exclude)
summary(m2)
anova(m2)
r.squaredGLMM(m2) #24% of variation explained by fixed effects, 24% by whole model (spatial variation?)
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))
#normally distributed, continue
LS1<-lsmeans(m2, ~nut_trt*fxnl_grp, by="fxnl_grp")
contrast(LS1, "pairwise")
#test nut_trt effects on biomass
m1<-lme(dry_wgt_g ~nut_trt*site, random=~1|yr, dat2, na.action=na.exclude)
summary(m1)
anova(m1)
r.squaredGLMM(m1) #24% of variation explained by fixed effects, 24% by whole model (spatial variation?)
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))
#normally distributed, continue
LS1<-lsmeans(m1, ~nut_trt*site, by="site")
contrast(LS1, "pairwise")
m2<-lme(dry_wgt_g ~nut_trt*fxnl_grp, random=~1|yr, dat, na.action=na.exclude)
summary(m2)
anova(m2)
r.squaredGLMM(m2) #24% of variation explained by fixed effects, 24% by whole model (spatial variation?)
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))
#normally distributed, continue
LS1<-lsmeans(m2, ~nut_trt*fxnl_grp, by="fxnl_grp")
contrast(LS1, "pairwise")
m2<-lme(dry_wgt_g ~nut_trt*ppt_trt, random=~1|site/yr, subset(dat2), na.action=na.exclude)
summary(m2)
anova(m2)
m2<-lme(dry_wgt_g ~nut_trt*ppt_trt, random=~1|site/yr, subset(dat2), na.action=na.exclude)
summary(m2)
anova(m2)
r.squaredGLMM(m2) #24% of variation explained by fixed effects, 24% by whole model (spatial variation?)
qqnorm(residuals(m2))
qqline(residuals(m2))
shapiro.test(residuals(m2))
#normally distributed, continue
LS1<-lsmeans(m2, ~nut_trt, by="ppt_trt")
contrast(LS1, "pairwise")
ggplot(dat2, aes(x=as.factor(yr), y = dry_wgt_g, fill=nut_trt)) + geom_boxplot()+ # facet_wrap(~ppt_trt)+
  theme_bw()+
  #facet_wrap(~clip_event)+
  labs(x="Year", y="Biomass (g of dry weight per 0.25m^2)") +
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))
ggplot(dat2, aes(x=ppt_trt, y = dry_wgt_g, fill = nut_trt)) + geom_boxplot() +
  theme_bw()+
  facet_wrap(~clip_event)+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))
ggplot(dat2, aes(x=ppt_trt, y = dry_wgt_g)) + geom_boxplot() +
  theme_bw()+
  facet_wrap(~clip_event)+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))
ggplot(dat2, aes(x=as.factor(yr), y = dry_wgt_g, fill=ppt_trt)) + geom_boxplot()+ # facet_wrap(~ppt_trt)+
  theme_bw()+
  #facet_wrap(~clip_event)+
  labs(x="Year", y="Biomass (g of dry weight per 0.25m^2)") #+
ggplot(dat2, aes(x=ppt_trt, y = dry_wgt_g, fill=ppt_trt)) + geom_boxplot()+ # facet_wrap(~ppt_trt)+
  theme_bw()+
  #facet_wrap(~clip_event)+
  labs(x="Year", y="Biomass (g of dry weight per 0.25m^2)") #+
m2<-lme(dry_wgt_g ~nut_trt*ppt_trt, random=~1|site/yr, subset(dat2), na.action=na.exclude)
summary(m2)
anova(m2)
r.squaredGLMM(m2) #24% of variation explained by fixed effects, 24% by whole model (spatial variation?)
qqnorm(residuals(m2))
qqline(residuals(m2))
shapiro.test(residuals(m2))
#normally distributed, continue
LS1<-lsmeans(m2, ~nut_trt, by="ppt_trt")
contrast(LS1, "pairwise")
#normally distributed, continue
LS1<-lsmeans(m2, ~ppt_trt, by="nut_trt")
contrast(LS1, "pairwise")
#normally distributed, continue
LS1<-lsmeans(m2, ~ppt_trt)
contrast(LS1, "pairwise")
m2<-lme(dry_wgt_g ~nut_trt*ppt_trt, random=~1|site/yr, subset(dat2), na.action=na.exclude)
summary(m2)
anova(m2)
r.squaredGLMM(m2) #24% of variation explained by fixed effects, 24% by whole model (spatial variation?)
qqnorm(residuals(m2))
qqline(residuals(m2))
shapiro.test(residuals(m2))
m2<-lme(dry_wgt_g ~nut_trt*ppt_trt, random=~1|yr, subset(dat2), na.action=na.exclude)
summary(m2)
anova(m2)
m2<-lme(dry_wgt_g ~nut_trt*ppt_trt, random=~1|yr, subset(dat2, sample_event==2), na.action=na.exclude)
m2<-lme(dry_wgt_g ~nut_trt*ppt_trt, random=~1|yr, subset(dat2, clip_event==2), na.action=na.exclude)
summary(m2)
anova(m2)
r.squaredGLMM(m2) #24% of variation explained by fixed effects, 24% by whole model (spatial variation?)
qqnorm(residuals(m2))
qqline(residuals(m2))
shapiro.test(residuals(m2))
#normally distributed, continue
LS1<-lsmeans(m2, ~ppt_trt, by="nut_trt")
contrast(LS1, "pairwise")
#normally distributed, continue
LS1<-lsmeans(m2, ~nut_trt, by="ppt_trt")
contrast(LS1, "pairwise")
m2<-lme(dry_wgt_g ~nut_trt*ppt_trt, random=~1|site/yr, subset(dat2, clip_event==2), na.action=na.exclude)
summary(m2)
anova(m2)
r.squaredGLMM(m2) #24% of variation explained by fixed effects, 24% by whole model (spatial variation?)
qqnorm(residuals(m2))
qqline(residuals(m2))
shapiro.test(residuals(m2))
#normally distributed, continue
LS1<-lsmeans(m2, ~nut_trt, by="ppt_trt")
contrast(LS1, "pairwise")
#normally distributed, continue
LS1<-lsmeans(m2, ~ppt_trt, by="nut_trt")
contrast(LS1, "pairwise")
contrast(LS1, "trt.vs.ctrl", ref = c("N","XC")
         contrast(LS1, "trt.vs.ctrl", ref = c("N","XC"))
         contrast(LS1, "trt.vs.ctrl", ref = c("N","XC"))
         contrast(LS1, "trt.vs.ctrl", ref = c("N"))
         contrast(LS1, "trt.vs.ctrl", ref = c("N","XC"))
         