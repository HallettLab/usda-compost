
# -- SETUP ----
rm(list = ls()) # clean environment
library(tidyverse)
options(stringsAsFactors = F)
theme_set(theme_bw())

# specify dropbox pathway (varies by user)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/"
}else{
  ## LMH and AS
  datpath <- "~/Dropbox/USDA-compost/Data/"
}

# read in main soil nutrient datatset and treatment keys
soil_nut<-read.csv(paste0(datpath, "belowground_expt/Soil/Soil_RawData/Compost_soilnutrients_raw.csv"))
key <- read_csv("~/Dropbox/USDA-compost/Data/belowground_expt/compost_soil_treatmentkey.csv")
key2<- read_csv("~/Dropbox/USDA-compost/Data/Compost_Treatmentkey.csv")

soil_nut$plotid<-soil_nut$Plot.ID

dat<-merge(soil_nut,key)
dat<-merge(soil_nut,key2)

dat$ppt_trt_f = factor(dat$ppt_trt, levels=c('D','XC','W'))

ggplot(subset(dat, PO4.P!="BQL"), aes(x=nut_trt, y=as.numeric(PO4.P), fill=nut_trt))+
  facet_wrap(~ppt_trt_f)+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))

ggplot(subset(dat, NO3.N!="BQL"&Subplot==1), aes(x=nut_trt, y=as.numeric(NO3.N), fill=nut_trt))+
  facet_wrap(~ppt_trt_f)+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))

ggplot(subset(dat, NH4.N!="BQL"&Subplot==1), aes(x=nut_trt, y=as.numeric(NH4.N), fill=nut_trt))+
  facet_wrap(~ppt_trt_f)+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))

ggplot(subset(dat, C!="BQL"), aes(x=nut_trt, y=as.numeric(C), fill=nut_trt))+
  facet_wrap(~ppt_trt_f)+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))

ggplot(subset(dat, N!="BQL"), aes(x=nut_trt, y=as.numeric(N), fill=nut_trt))+
  facet_wrap(~ppt_trt_f)+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))

ggplot(subset(dat, OM!="BQL"), aes(x=nut_trt, y=as.numeric(OM), fill=nut_trt))+
  facet_wrap(~ppt_trt_f)+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = c("indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Amendment"), labels=c("Compost", "Fertilizer", "None"))



                    