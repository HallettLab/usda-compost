theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 14),
              strip.text= element_text(size = 14), 
              axis.text = element_text(size = 12))

dat <- read_csv("~/Dropbox/USDA-compost/Data/belowground_expt/ANPP_20210515.csv")
key <- read_csv("~/Dropbox/USDA-compost/Data/belowground_expt/compost_soil_treatmentkey.csv")
key2<- read_csv("~/Dropbox/USDA-compost/Data/Compost_Treatmentkey.csv")

dat<-merge(dat,key)
dat<-merge(dat,key2)

ggplot(subset(dat, ppt_trt=="D"&below_trt!="no_both"&below_trt!="no_nema"), aes(x=below_trt, y=dry_wgt_g))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))
  
ggplot(subset(dat, ppt_trt=="D"), aes(x=below_trt, y=(dry_wgt_g*100*8.92179), fill=below_trt))+
  geom_boxplot()+
  #theme_bw()+
  labs(x="Belowground Treatment", y="Forage (lbs per acre)")+
  scale_fill_manual(values= c("black","grey80","white","grey50"), guide = guide_legend(title = "Treatment"), labels=c("Control", "No AMF", "No Both", "No Nematodes"))

ggplot(subset(dat, ppt_trt=="D"&below_trt!="no_both"&below_trt!="no_amf"), aes(x=below_trt, y=(dry_wgt_g*100)))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))

ggplot(subset(dat, ppt_trt=="D"), aes(x=below_trt, y=(dry_wgt_g*100*8.9),  fill=nut_trt))+
  geom_boxplot()+
  #theme_bw()+
  labs(x="Belowground Treatment", y="Forage (lbs per acre)")+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))

