#setwd("~/Dropbox/USDA-compost/Data")
# specify dropbox pathway (varies by user -- ctw can tweak this later)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/"
}else{
  ## LMH and AS
  datpath <- "~/Dropbox/USDA-compost/Data/"
}

pheno<-read.csv(paste0(datpath, "Phenology/Phenology_CleanedData/Compost_Phenology_Clean.csv"))

#show phenology overtime by nut_trt, ppt_trt and elevation
pheno_gf<-pheno %>% filter(plot!="NA")%>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site),
         date = as.Date(date, format = "%Y-%m-%d")) %>%
  mutate(ppt_trt=ordered(ppt_trt, levels = c(D="D", XC="XC", W="W"))) %>%
  group_by(date, site, nut_trt, ppt_trt) %>%
  summarize(meanPG=mean(pct_green), sePG=sd(pct_green)/sqrt(length(pct_green))) %>%
  ungroup()

ggplot(data=pheno_gf, aes(x=date, y=meanPG, color=nut_trt))+
  facet_wrap(~site*ppt_trt)+
  geom_point(cex=1.5)+
  geom_errorbar(aes(ymax = meanPG+sePG, ymin = meanPG-sePG), width=.25)+
  geom_line()+
  labs(x="Date", y="Percent Green") +
  theme(text = element_text(size=15))+
  theme_bw()+
  scale_color_manual(values = c( "indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Treatment"))

#pheno plot for low dry plots only
pheno_gf2<-pheno %>% filter(plot!="NA")%>%
  mutate(site = "high", 
         site = ifelse(block == 1 | block == 2, "low", site),
         date = as.Date(date, format = "%Y-%m-%d")) %>%
  filter(site=="low")%>%
  filter(ppt_trt=="D")%>%
  group_by(date, nut_trt) %>%
  summarize(meanPG=mean(pct_green), sePG=sd(pct_green)/sqrt(length(pct_green))) %>%
  ungroup()

ggplot(data=pheno_gf2, aes(x=date, y=meanPG, color=nut_trt))+
  geom_point(cex=1.5)+
  geom_errorbar(aes(ymax = meanPG+sePG, ymin = meanPG-sePG), width=.25)+
  geom_line()+
  labs(x="Date", y="Percent Green") +
  theme(text = element_text(size=20))+
  theme_bw()+
  scale_color_manual(values = c( "indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Treatment"), labels=c("Compost", "Fertilizer", "None"))


