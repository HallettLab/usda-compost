nema <- read_csv("~/Dropbox/USDA-compost/Data/belowground_expt/compost_nematodes_subplot1_2021.csv")
micro<- read_csv("~/Dropbox/USDA-compost/Data/soil_community/microarthropods_2019.csv")
nema<-nema[-c(37, 38, 39), ]

nema$ppt_trt_f = factor(nema$Ppt_trt, levels=c('D','XC','W'))

ggplot(data=subset(nema, Block==3|Block==4), aes(x=ppt_trt_f, y = nema_kg_soil, fill=ppt_trt_f)) + 
  geom_boxplot()+ # facet_wrap(~ppt_trt)+
  #theme_bw()+
  #facet_wrap(~Ppt_trt)+
  labs(x="", y="Nematodes per kg dry soil") +
  scale_fill_manual(values = c("indianred1",  "lightgoldenrod2", "skyblue2"), guide = guide_legend(title = "Precipitation"), labels=c("Compost", "Fertilizer", "None"))

ggplot(data=subset(nema, Block==3|Block==4), aes(x=Nut_trt, y = nema_kg_soil, fill=ppt_trt_f)) + 
  geom_boxplot()+ # facet_wrap(~ppt_trt)+
  #theme_bw()+
  #facet_wrap(~Ppt_trt)+
  labs(x="", y="Nematodes per kg dry soil") +
  scale_fill_manual(values = c("indianred1",  "lightgoldenrod2", "skyblue2"), guide = guide_legend(title = "Precipitation"), labels=c("Dry", "Ambient", "Wet"))

ggplot(data=subset(nema, ppt_trt_f=="XC"&Block==1|ppt_trt_f=="XC"&Block==2), aes(x=Nut_trt, y = nema_kg_soil, fill=Nut_trt)) + 
  geom_boxplot()+ # facet_wrap(~ppt_trt)+
  #theme_bw()+
  #facet_wrap(~Ppt_trt)+
  labs(x="", y="Nematodes per kg dry soil") +
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))

ggplot(data=subset(micro, treatment_2=="w"&plot_num==3|treatment_2=="w"&plot_num==4), aes(x=treatment, y = col_kg, fill=treatment)) + 
  geom_boxplot()+ # facet_wrap(~ppt_trt)+
  #theme_bw()+
  #facet_wrap(~Ppt_trt)+
  labs(x="", y="Collembola per kg dry soil") +
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))

df <- data.frame (nut_trt  = c("C", "C", "C", "F","F","F","N","N","N", "C", "C", "C", "F","F","F","N","N","N"),
                  Group= c("Microbivores","Microbivores","Microbivores","Microbivores","Microbivores","Microbivores","Microbivores","Microbivores","Microbivores","Parasites","Parasites","Parasites","Parasites","Parasites","Parasites","Parasites","Parasites","Parasites"),
                  Percent = c("55", "65", "59", "35","50","25","29","60","45", "23","25","40","55","48","58", "45","25","30"))
      

ggplot(data=subset(df), aes(x=Group, y = as.numeric(Percent), fill=nut_trt)) + 
  geom_boxplot()+ # facet_wrap(~ppt_trt)+
  #theme_bw()+
  #facet_wrap(~Ppt_trt)+
  labs(x="", y="Nematode group (%)") +
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))

