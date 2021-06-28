setwd("~/Dropbox/USDA-compost/Data")
pheno<-read.csv("Phenology/Phenology_CleanedData/Compost_Phenology_Clean.csv")

#show phenology overtime by nut_trt, ppt_trt and elevation
pheno_gf<-pheno %>% filter(plot!="NA")%>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site)) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W"))) %>%
  group_by(date, site,nut_trt, ppt_trt, yr) %>%
  summarize(meanPG=mean(pct_green), sePG=sd(pct_green)/sqrt(length(pct_green)))

ggplot(data=pheno_gf, aes(x=date, y=meanPG, color=nut_trt))+
  facet_grid(site~ppt_trt*yr)+
  geom_point(aes(cex=1.5))+
  geom_errorbar(aes(ymax = meanPG+sePG, ymin = meanPG-sePG), width=.25)+
  geom_line()+
  labs(x="Date", y="Percent Green") +
  theme(text = element_text(size=15))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_color_manual(values = c( "indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Treatment"))

#pheno plot for low dry plots only
pheno_gf2<-pheno %>% filter(plot!="NA")%>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site)) %>%
  #filter(site=="low")%>%
  filter(ppt_trt=="D")%>%
  group_by(date, nut_trt, yr) %>%
  summarize(meanPG=mean(pct_green), sePG=sd(pct_green)/sqrt(length(pct_green)))

ggplot(data=pheno_gf2, aes(x=date, y=meanPG, color=nut_trt))+
  geom_point(aes())+
  geom_errorbar(aes(ymax = meanPG+sePG, ymin = meanPG-sePG), width=.25)+
  geom_line()+
  #facet_wrap(~yr)+
  labs(x="Date", y="Percent Green") +
  theme(text = element_text(size=20))+
  theme_bw()+
  scale_color_manual(values = c( "indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Treatment"), labels=c("Compost", "Fertilizer", "None"))


