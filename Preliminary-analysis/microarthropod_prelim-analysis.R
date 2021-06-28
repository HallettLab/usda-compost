library(nlme)
library(ggplot2)
library(multcomp)
library(lsmeans)
library(MuMIn)

# read in data 
ma <- read.csv("~/Dropbox/USDA-compost/Data/soil_community/microarthropods_2019.csv", header=T) %>%
  mutate(ppt_trt1=treatment_2) %>% dplyr::select(-treatment_2)%>%
  mutate(block=plot_num) %>% dplyr::select(-plot_num)%>%
  mutate(nut_trt=treatment) %>% dplyr::select(-treatment)%>%
  mutate(ppt_trt= ifelse(ppt_trt1=="d", "D", ifelse(ppt_trt1=="w", "W", ifelse(ppt_trt1=="xc", "XC", ppt_trt1)))) %>%
  dplyr::select(-ppt_trt1)

ggplot(ma, aes(x=nut_trt, y = col_kg, fill = ppt_trt)) + geom_boxplot() 
ggplot(ma, aes(x=nut_trt, y = mite_kg, fill = ppt_trt)) + geom_boxplot() 

#summarize collembola
ma2 <- ma %>%
  group_by(ppt_trt, nut_trt) %>%
  summarize(mean_col = mean(col_kg), secol=sd(col_kg)/sqrt(length(col_kg)))

ggplot(data=ma2, aes(x=nut_trt, y=mean_col, fill=nut_trt))+
  theme_bw()+
  theme(strip.background = element_blank(), 
        text = element_text(size = 20), 
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  facet_wrap(~ppt_trt)+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None")) +
  geom_errorbar(aes(ymin=mean_col-secol, ymax=mean_col+secol))+
  labs(x="Amendment", y="Collembola/kg dry soil") 

ma2b <- ma %>%
  group_by(nut_trt) %>%
  summarize(mean_col = mean(col_kg), secol=sd(col_kg)/sqrt(length(col_kg)))

ggplot(data=ma2b, aes(x=nut_trt, y=mean_col, fill=nut_trt))+
  theme_bw()+
  theme(strip.background = element_blank(), 
        text = element_text(size = 20), 
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  #facet_wrap(~ppt_trt)+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None")) +
  geom_errorbar(aes(ymin=mean_col-secol, ymax=mean_col+secol))+
  labs(x="Amendment", y="Collembolans per kg dry soil") 

#summarize mites
ma3 <- ma %>%
  group_by(nut_trt) %>%
  summarize(mean_mite = mean(mite_kg), secol=sd(mite_kg)/sqrt(length(mite_kg)))

ggplot(data=ma3, aes(x=nut_trt, y=mean_mite, fill=nut_trt))+
  theme_bw()+
  theme(strip.background = element_blank(), 
        text = element_text(size = 20), 
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  #facet_wrap(~ppt_trt)+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None")) +
  geom_errorbar(aes(ymin=mean_mite-secol, ymax=mean_mite+secol))+
  labs(x="Amendment", y="Mites/kg dry soil") 
