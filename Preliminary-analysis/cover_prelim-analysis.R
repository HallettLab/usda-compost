library(tidyverse)
library(vegan)
library(lsmeans)
library(nlme)

# specify dropbox pathway (varies by user -- ctw can tweak this later)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/Competition/Competition_EnteredData")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/Competition/Competition_EnteredData"
}else{
  ## LMH and AS
  datpath <- "~/Dropbox/USDA-compost/Data/Competition/Competition_EnteredData"
}

cover_master_wide <- read_csv(paste0(datpath,"Cover/Cover_CleanedData/Compost_Cover_WideClean.csv"))

cover_master_long <- read_csv(paste0(datpath, "Cover/Cover_CleanedData/Compost_Cover_LongClean.csv"))


dat <- cover_master_wide %>%
  tbl_df() %>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site)) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W"))) 
  
cover_ord <- dat %>% dplyr::select(17:84)
dat$shannon <-diversity(cover_ord)

ggplot(dat, aes(y=shannon, x=nut_trt, fill=nut_trt))+geom_boxplot()+ #facet_wrap(~site)+
  theme_bw()+
  facet_grid(ppt_trt~yr)+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))

ggplot(dat, aes(y=shannon, x=nut_trt, fill=nut_trt))+geom_boxplot()+ #facet_wrap(~site)+
  theme_bw()+
  #facet_grid(ppt_trt~yr)+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), 
                    guide = guide_legend(title = "Amendment"), 
                    labels=c("Compost", "Fertilizer", "None"))+
  theme(strip.background = element_blank(), 
        text = element_text(size = 12), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  labs(y = "Shannon Diversity Index", x = "Amendment Treatment")+
  scale_x_discrete(labels=c("C" = "Compost", "F" = "Fertilizer", "N" = "None"))


m1<-lme(shannon ~nut_trt*ppt_trt, random=~1|site/yr, dat, na.action=na.exclude)
summary(m1)
anova(m1)
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))

LS1<-lsmeans(m1, ~nut_trt)
contrast(LS1, "pairwise")
#significant effect of nut_trt

# variations on a theme with grass vs forb visualization
# compost seems to help grass at low site, marginally hurt at high site
ggplot(dat, aes(y = pct_grass, x=interaction(ppt_trt))) + geom_boxplot() + facet_wrap(~nut_trt)
ggplot(dat, aes(y = pct_grass, x=ppt_trt, color = site, shape=as.factor(yr))) + geom_point() + facet_wrap(~nut_trt)
ggplot(dat, aes(y = pct_grass, x=ppt_trt, color = site)) + geom_boxplot() + facet_grid(site~nut_trt)
ggplot(dat, aes(y = pct_grass, x=nut_trt, color = site)) + geom_boxplot() + facet_grid(~site*yr)
ggplot(dat, aes(y = pct_grass, x=ppt_trt, color = site)) + geom_boxplot() + facet_grid(~site*yr)

# compost seems to help forb at high site, *marginally* hurt at low (makes sense, grass inverse)
# overall seems to like wet (maybe because of N-fixers?)
ggplot(dat, aes(y =  pct_forb, x=interaction(ppt_trt), fill=as.factor(yr))) + geom_boxplot() + facet_wrap(~nut_trt)
ggplot(dat, aes(y = pct_forb, x=ppt_trt, color = site, shape=as.factor(yr))) + geom_point() + facet_wrap(~nut_trt)
ggplot(dat, aes(y = pct_forb, x=ppt_trt, color = site)) + geom_boxplot() + facet_grid(site~nut_trt)
ggplot(dat, aes(y = pct_forb, x=nut_trt, color = site)) + geom_boxplot() + facet_grid(~site*yr)
ggplot(dat, aes(y = pct_forb, x=ppt_trt, color = site)) + geom_boxplot() + facet_grid(~site*yr)

#calculate sum of cover by fxnl-group to include N-fixers, may data only
cover.fxnl <- cover_master_long %>% filter (date=="2019-05-09"|date=="2019-05-08"|
                                              date=="2019-05-07"|date=="2020-04-28"|
                                              date=="2020-04-29"|date=="2020-04-30") %>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site)) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W")))%>%
  group_by(nut_trt, ppt_trt, block, date, site, fxnl_grp, yr) %>%
  summarise(cover=sum(pct_cover))%>%
  group_by(nut_trt, ppt_trt,fxnl_grp, yr, site)%>%
  summarise(meanc=mean(cover), se=sd(cover)/sqrt(length(cover)))

cover.fxnl.low<-cover.fxnl %>% filter(site=="low")
cover.fxnl.high<-cover.fxnl %>% filter(site=="high")

m1<-lme(meanc ~nut_trt*ppt_trt*fxnl_grp*yr, random=~1|site, cover.fxnl, na.action=na.exclude)
summary(m1)
anova(m1)
qqnorm(residuals(m1))
qqline(residuals(m1))
shapiro.test(residuals(m1))

LS1<-lsmeans(m1, ~fxnl_grp*nut_trt*ppt_trt, by=c("fxnl_grp","ppt_trt"))
contrast(LS1, "pairwise")
#significant effect of nut_trt

ggplot(cover.fxnl.high, aes(y=meanc, x=fxnl_grp, fill=ppt_trt))+
  geom_bar(stat="identity", position = position_dodge())+ 
  geom_errorbar(aes(ymax = meanc+se, ymin = meanc-se), position = position_dodge(.9), width=.25)+
  theme_bw()+
  facet_grid(nut_trt~yr)+
  scale_fill_manual(values = c("saddlebrown",  "gray", "lightblue"), 
                    guide = guide_legend(title = "Rainfall"), 
                    labels=c("Dry", "Ambient", "Wet"))+
  theme(strip.background = element_blank(), 
        text = element_text(size = 18), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  labs(y = "Percent Cover", x = "Functional Group")+
  scale_x_discrete(labels=c("C" = "Compost", "F" = "Fertilizer", "N" = "None"))


# hordeum seems to <3 compost, but sample size is low because not in the lower site
ggplot(dat, aes(y = HOMU, x=ppt_trt)) + geom_boxplot() + facet_grid(nut_trt~yr)
ggplot(dat, aes(y = HOMU, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = HOMU, x=ppt_trt)) + geom_boxplot() + facet_grid(~site*yr)
ggplot(dat, aes(y = HOMU, x=nut_trt)) + geom_boxplot() + facet_grid(~site*yr)

# bromus likes nutrients, takes off in 2020 in low site
ggplot(dat, aes(y = BRHO, x=ppt_trt)) + geom_boxplot() + facet_grid(nut_trt~yr)
ggplot(dat, aes(y = BRHO, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = BRHO, x=ppt_trt)) + geom_boxplot() + facet_grid(~site*yr)
ggplot(dat, aes(y = BRHO, x=nut_trt)) + geom_boxplot() + facet_grid(~site*yr)

# lomu seems to like compost, and is in both locations
ggplot(dat, aes(y = LOMU, x=interaction(ppt_trt))) + geom_boxplot() + facet_grid(nut_trt~yr)
ggplot(dat, aes(y = LOMU, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt*yr)
ggplot(dat, aes(y = LOMU, x=ppt_trt)) + geom_boxplot() + facet_grid(~site*yr)
ggplot(dat, aes(y = LOMU, x=nut_trt)) + geom_boxplot() + facet_grid(~yr)

# taca loves fertilizer
ggplot(dat, aes(y = TACA, x=interaction(ppt_trt))) + geom_boxplot() + facet_grid(nut_trt~yr)
ggplot(dat, aes(y = TACA, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt*yr)
ggplot(dat, aes(y = TACA, x=ppt_trt)) + geom_boxplot() + facet_grid(~site*yr)
ggplot(dat, aes(y = TACA, x=nut_trt)) + geom_boxplot() + facet_grid(~yr)


# avba doesn't seem to like the compost much, especially under wet conditions
# not much at the lower sites
ggplot(dat, aes(y = AVBA, x=interaction(ppt_trt))) + geom_boxplot() + facet_grid(nut_trt~yr)
ggplot(dat, aes(y = AVBA, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = AVBA, x=nut_trt)) + geom_boxplot() + facet_grid(~site*yr)
ggplot(dat, aes(y = AVBA, x=ppt_trt)) + geom_boxplot() + facet_grid(~site*yr)

# aira seems to dislike compost, but sample size is low because not in the lower site
ggplot(dat, aes(y = AICA, x=ppt_trt)) + geom_boxplot() + facet_grid(nut_trt~yr)
ggplot(dat, aes(y = AICA, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = AICA, x=ppt_trt)) + geom_boxplot() + facet_grid(~site*yr)
ggplot(dat, aes(y = AICA, x=nut_trt)) + geom_boxplot() + facet_grid(site~yr)


# erbo doesn't seem fussed really
ggplot(dat, aes(y = ERBO, x=interaction(ppt_trt))) + geom_boxplot() + facet_grid(nut_trt~yr)
ggplot(dat, aes(y = ERBO, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = ERBO, x=nut_trt)) + geom_boxplot() + facet_grid(~site*yr)
ggplot(dat, aes(y = ERBO, x=ppt_trt)) + geom_boxplot() + facet_grid(~site*yr)


# trhi likes amendments - of both kinds, and especially when not in drought
ggplot(dat, aes(y = TRHI, x=interaction(ppt_trt))) + geom_boxplot() + facet_wrap(yr~nut_trt)
ggplot(dat, aes(y = TRHI, x=interaction(ppt_trt))) + geom_boxplot() + facet_grid(site~nut_trt)
ggplot(dat, aes(y = TRHI, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = TRHI, x=nut_trt)) + geom_boxplot() + facet_grid(~yr)
ggplot(dat, aes(y = TRHI, x=ppt_trt)) + geom_boxplot() + facet_grid(~site*yr)

# trdu likes compost - and especially wet ppt
ggplot(dat, aes(y = TRDU, x=interaction(ppt_trt))) + geom_boxplot() + facet_wrap(yr~nut_trt)
ggplot(dat, aes(y = TRDU, x=interaction(ppt_trt))) + geom_boxplot() + facet_grid(site~nut_trt)
ggplot(dat, aes(y = TRDU, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = TRDU, x=nut_trt)) + geom_boxplot() + facet_grid(~yr)
ggplot(dat, aes(y = TRDU, x=ppt_trt)) + geom_boxplot() + facet_grid(~site*yr)

# trsu likes amendments - of both kinds
ggplot(dat, aes(y = TRSU, x=interaction(ppt_trt))) + geom_boxplot() + facet_grid(nut_trt~yr)
ggplot(dat, aes(y = TRSU, x=interaction(ppt_trt))) + geom_boxplot() + facet_grid(site~nut_trt)
ggplot(dat, aes(y = TRSU, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = TRSU, x=nut_trt)) + geom_boxplot() + facet_grid(yr~site)
ggplot(dat, aes(y = TRSU, x=ppt_trt)) + geom_boxplot() + facet_grid(yr~site)


# together trifoliums like amendments (especially straight fertilizer, especially at low site)
ggplot(dat, aes(y = (TRDU + TRHI + TRSU), x=interaction(ppt_trt))) + geom_boxplot() + facet_grid(nut_trt~yr)
ggplot(dat, aes(y = (TRDU + TRHI + TRSU), x=interaction(ppt_trt))) + geom_boxplot() + facet_grid(yr~nut_trt)
ggplot(dat, aes(y = (TRDU + TRHI + TRSU), x=ppt_trt, color = site, shape = as.factor(yr))) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = (TRDU + TRHI + TRSU), x=nut_trt)) + geom_boxplot() + facet_grid(~yr)
ggplot(dat, aes(y = (TRDU + TRHI + TRSU), x=ppt_trt)) + geom_boxplot() + facet_grid(site~yr)

# geranium seems to <3 compost
ggplot(dat, aes(y = GEMO, x=ppt_trt)) + geom_boxplot() + facet_grid(nut_trt~yr)
ggplot(dat, aes(y = GEMO, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = GEMO, x=ppt_trt)) + geom_boxplot() + facet_grid(~site*yr)
ggplot(dat, aes(y = GEMO, x=nut_trt)) + geom_boxplot() + facet_grid(site~yr)

# vetch seems to <3 compost
ggplot(dat, aes(y = VISA, x=ppt_trt)) + geom_boxplot() + facet_grid(nut_trt~yr)
ggplot(dat, aes(y = VISA, x=ppt_trt, color = site)) + geom_point() + facet_grid(~nut_trt)
ggplot(dat, aes(y = VISA, x=ppt_trt)) + geom_boxplot() + facet_grid(site~yr)
ggplot(dat, aes(y = VISA, x=nut_trt)) + geom_boxplot() + facet_grid(site~yr)

# other grasses pretty low numbers
ggplot(dat, aes(y = BRHO, x=interaction(ppt_trt))) + geom_boxplot() + facet_grid(nut_trt~yr)
ggplot(dat, aes(y = VUBR, x=interaction(ppt_trt))) + geom_boxplot() + facet_grid(nut_trt~yr)



dat_may <- cover_master_wide %>%
  tbl_df() %>%
  filter(date=="2019-05-09"|date=="2019-05-08"|date=="2019-05-07"|date=="2020-04-28"|date=="2020-04-29"|date=="2020-04-30"|date=="2021-05-13")%>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site)) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W"))) 


dat_may<-dat_may%>%mutate(site=ordered(site, levels = c(low="low", high="high")))%>%
  mutate(nut_trt=ordered(nut_trt, levels = c(C="C", F="F", N="N")))%>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W")))
  

cover_ord_may <- dat_may %>% dplyr::select(18:92)

cover_ord_all <- dat %>% dplyr::select(18:92)

dat<-dat%>%mutate(site=ordered(site, levels = c(low="low", high="high")))%>%
  mutate(nut_trt=ordered(nut_trt, levels = c(C="C", F="F", N="N")))



#make bray-curtis dissimilarity matrix
spp.bcd <- vegdist(cover_ord_may)
spp.mds<-metaMDS(cover_ord_may, trace = FALSE, autotransform=T, trymax=100, k=4) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
spp.mds #solution converged after 20 tries, stress = 9.04
summary(spp.mds)

stressplot(spp.mds, spp.bcd) #stressplot to show fit, fit is decent
ordiplot(spp.mds)
spscores1<-scores(spp.mds,display="sites",choices=1)
spscores2<-scores(spp.mds,display="sites",choices=2)
tplots<-dat_may[,2]
tplot_levels<-levels(tplots)
spscoresall<-data.frame(tplots,spscores1,spscores2)

legend2<-data.frame("date"=c("April19", "May19", "April20"))
legend2<-legend2 %>%mutate(date=ordered(date, levels = c(April19="April19", May19="May19", April20="April20" )))

bio.plot <- ordiplot(spp.mds,choices=c(1,2), type = "none")   #Set up the plot
colsT1 <- rep(c("indianred4",  "dodgerblue1", "darkgoldenrod"), each = 9) #color based on amendment
Lcols <- rep(c("indianred4",  "dodgerblue1", "darkgoldenrod"))
shapes <- rep(c(15, 3, 17), each=1) #shapes on rainfall
Lshapes <- rep(c(15,3,17))
points(spscoresall$NMDS1,spscoresall$NMDS2,col=colsT1,pch=shapes) 
text(spp.mds, display = "species", cex=0.5, col="grey30") #label species
legend("bottomright",legend=levels(as.factor(cover_master_wide$nut_trt)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
# add legend for treatment
legend("topright",legend=levels(as.factor(cover_master_wide$ppt_trt)), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
#help(ordiplot)

bio.plot2 <- ordiplot(spp.mds,choices=c(1,2), type = "none")   #Set up the plot
colsT2 <- rep(c("grey70",  "black"), each = 36) #color based on site
Lcols2 <- rep(c("grey70",  "black", "darkgoldenrod"))
shapes <- rep(c(15, 17), each=1) #shapes on year
Lshapes <- rep(c(15,17))
points(spscoresall$NMDS1,spscoresall$NMDS2,col=colsT2,pch=shapes) 
text(spp.mds, display = "species", cex=0.5, col="grey30") #label species
legend("bottomright",legend=levels(as.factor(dat_may$site)), col=Lcols2, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
# add legend for treatment
legend("topright",legend=levels(as.factor(dat_may$yr)), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
#help(ordiplot)

native_gf <- table(cover_master_long$nativity, cover_master_long$nut_trt)
  native_gf<-as.data.frame(native_gf)
  
ggplot(native_gf, aes(x=Var2, y=Freq, fill=Var1))+
  geom_bar(stat="identity", position=position_dodge())
  

dat_long <- cover_master_long %>%
  tbl_df() %>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site)) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W")))

dm <-dat_long%>%
  filter(date=="2019-05-09"|date=="2019-05-08"|date=="2019-05-07")%>%
  arrange_(~ desc(pct_cover)) %>%
  group_by_(~ plot) %>%
  slice(1:5)

da <-dat_long%>%
  filter(date=="2019-04-19"|date=="2019-04-20")%>%
  arrange_(~ desc(pct_cover)) %>%
  group_by_(~ plot) %>%
  slice(1:5)

dl <- cover_master_long %>%
  tbl_df() %>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site)) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W"))) 

summary <- dl %>% group_by(nut_trt, ppt_trt, site, yr, fxnl_grp, genus) %>%
  filter(pct_cover > 5)%>%
  summarize(meancov = mean(pct_cover))

library(RColorBrewer)
#stacked bar plot
ggplot(summary, aes(fill=genus, y=meancov, x=nut_trt)) +
  geom_bar( stat="identity")+
  facet_grid(site~ppt_trt*yr)+
  scale_fill_manual(values=c("#49ab98",
                             "#ab5bdb",
                             "#5bb839",
                             "#d659c8",
                             "#4ebd62",
                             "#674db9",
                             "#acb72e",
                             "#6d78f1",
                             "#e4a31f",
                             "#4691eb",
                             "#d1952d",
                             "#6373d2",
                             "#7faa3d",
                             "#b93f9a",
                             "#528632",
                             "#a770cf",
                             "#73b46c",
                             "#d14285",
                             "#4cae7c",
                             "#db3963",
                             "#4eaad4",
                             "#d15621",
                             "#617ec5",
                             "#c08938",
                             "#cc78c4",
                             "#969942",
                             "#907ebc",
                             "#d13f3c",
                             "#71874c",
                             "#955282",
                             "#cf995c",
                             "#d489b7",
                             "#825c24",
                             "#bd5c7d",
                             "#c19266",
                             "#c14b5c",
                             "#cc6347",
                             "#c26e78",
                             "#c67060"))

#rank shift
library(codyn)
dat_may_l <- cover_master_long %>%
  tbl_df() %>%
  filter(date=="2019-05-09"|date=="2019-05-08"|date=="2019-05-07"|date=="2020-04-28"|date=="2020-04-29"|date=="2020-04-30")%>%
  dplyr::select (plotid, yr, pct_cover, species) 

rs<-rank_shift (dat_may_l, time.var="yr", species.var="species", abundance.var="pct_cover", replicate.var = as.character("plotid"))
rs<-rs%>% separate(plotid, c("block","nut_trt","ppt_trt"), sep=cumsum(c(1,1,2)))%>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site)) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W")))
ggplot(rs, aes(x=ppt_trt, y=MRS, fill=nut_trt))+ geom_boxplot()

ggplot(data=rs, aes(x=ppt_trt, y=MRS, fill=nut_trt))+
  theme_bw()+
  facet_wrap(~site)+
  theme(strip.background = element_blank(), 
        text = element_text(size = 20), 
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  #facet_wrap(~ppt_trt)+
  geom_boxplot()+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None")) +
  #geom_errorbar(aes(ymin=mean_dry_wgt_g-sedry, ymax=mean_dry_wgt_g+sedry))+
  labs(x="Rainfall Treatment", y="Mean Rank Shift") +
  scale_x_discrete(labels=c("D" = "Dry", "XC" = "Ambient", "W" = "Wet"))

to<-turnover(
  dat_may_l,
  time.var="yr",
  species.var="species",
  abundance.var="pct_cover",
  replicate.var = "plotid",
  metric = "total")
to<-to%>% separate(plotid, c("block","nut_trt","ppt_trt"), sep=cumsum(c(1,1,2)))%>%
  mutate(site = "moderate_graze", 
         site = ifelse(block == 1 | block == 2, "heavy_graze", site)) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W")))
ggplot(to, aes(x=nut_trt, y=total))+ geom_boxplot() +facet_grid(~ppt_trt)
ggplot(to, aes(x=nut_trt, y=total))+ geom_boxplot()

ggplot(data=to, aes(x=ppt_trt, y=total, fill=nut_trt))+
  theme_bw()+
  #facet_wrap(~site)+
  theme(strip.background = element_blank(), 
        text = element_text(size = 20), 
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  #facet_wrap(~ppt_trt)+
  geom_boxplot()+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None")) +
  #geom_errorbar(aes(ymin=mean_dry_wgt_g-sedry, ymax=mean_dry_wgt_g+sedry))+
  labs(x="Rainfall Treatment", y="Species Turnover") +
  scale_x_discrete(labels=c("D" = "Dry", "XC" = "Ambient", "W" = "Wet"))


tol<-turnover(
  dat_may_l,
  time.var="yr",
  species.var="species",
  abundance.var="pct_cover",
  replicate.var = "plotid",
  metric = "disappearance")
tol<-tol%>% separate(plotid, c("block","nut_trt","ppt_trt"), sep=cumsum(c(1,1,2)))%>%
  mutate(site = "moderate_graze", 
         site = ifelse(block == 1 | block == 2, "heavy_graze", site)) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W")))
ggplot(tol, aes(x=ppt_trt, y=disappearance))+ geom_boxplot()+facet_grid(~nut_trt)
ggplot(tol, aes(x=ppt_trt, y=disappearance))+ geom_boxplot()

ggplot(data=tol, aes(x=ppt_trt, y=disappearance, fill=nut_trt))+
  theme_bw()+
  facet_wrap(~site)+
  theme(strip.background = element_blank(), 
        text = element_text(size = 20), 
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  #facet_wrap(~ppt_trt)+
  geom_boxplot()+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None")) +
  #geom_errorbar(aes(ymin=mean_dry_wgt_g-sedry, ymax=mean_dry_wgt_g+sedry))+
  labs(x="Rainfall Treatment", y="Species Losses") +
  scale_x_discrete(labels=c("D" = "Dry", "XC" = "Ambient", "W" = "Wet"))



tog<-turnover(
  dat_may_l,
  time.var="yr",
  species.var="species",
  abundance.var="pct_cover",
  replicate.var = "plotid",
  metric = "appearance")
tog<-tog%>% separate(plotid, c("block","nut_trt","ppt_trt"), sep=cumsum(c(1,1,2)))%>%
  mutate(site = "moderate_graze", 
         site = ifelse(block == 1 | block == 2, "heavy_graze", site)) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W")))
ggplot(tog, aes(x=ppt_trt, y=appearance))+ geom_boxplot()+facet_grid(~nut_trt)
ggplot(tog, aes(x=ppt_trt, y=appearance))+ geom_boxplot()

ggplot(data=tog, aes(x=ppt_trt, y=appearance, fill=nut_trt))+
  theme_bw()+
  facet_wrap(~site)+
  theme(strip.background = element_blank(), 
        text = element_text(size = 20), 
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  #facet_wrap(~ppt_trt)+
  geom_boxplot()+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None")) +
  #geom_errorbar(aes(ymin=mean_dry_wgt_g-sedry, ymax=mean_dry_wgt_g+sedry))+
  labs(x="Rainfall Treatment", y="Species Gains") +
  scale_x_discrete(labels=c("D" = "Dry", "XC" = "Ambient", "W" = "Wet"))

#richness and evenness
struc<-community_structure(
  dat_may_l,
  time.var = "yr",
  abundance.var="pct_cover",
  replicate.var = "plotid",
  metric = c("SimpsonEvenness"))

struc<-struc%>% separate(plotid, c("block","nut_trt","ppt_trt"), sep=cumsum(c(1,1,2))) %>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site))%>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W")))
ggplot(struc, aes(x=ppt_trt, y=SimpsonEvenness, fill=nut_trt))+ geom_boxplot()+facet_grid(yr~site)
ggplot(tog, aes(x=ppt_trt, y=appearance))+ geom_boxplot()

ggplot(data=struc, aes(x=ppt_trt, y=SimpsonEvenness, fill=nut_trt))+
  theme_bw()+
  facet_wrap(site~yr)+
  theme(strip.background = element_blank(), 
        text = element_text(size = 20), 
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  #facet_wrap(~ppt_trt)+
  geom_boxplot()+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None")) +
  #geom_errorbar(aes(ymin=mean_dry_wgt_g-sedry, ymax=mean_dry_wgt_g+sedry))+
  labs(x="Rainfall Treatment", y="Evenness (Simpson's)") +
  scale_x_discrete(labels=c("D" = "Dry", "XC" = "Ambient", "W" = "Wet"))

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

