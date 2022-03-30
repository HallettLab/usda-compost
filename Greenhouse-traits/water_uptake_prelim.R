# -- SETUP ----
rm(list = ls()) # clean enviroment
# libraries needed
library(tidyverse)
library(lubridate)

# change default settings
na_vals <- c("", " ", NA, "NA") #tells it what to call NA
options(stringsAsFactors = F) #setting for strings of characters 'abc'
theme_set(theme_bw()) #for plotting graphs later



# specify dropbox pathway (varies by user -- EAS can tweak this for others use)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/PlantTraits/greenhouse_traits/")){
  ## CTW pathway
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/PlantTraits/greenhouse_traits/"
}else{
  ## LMH and EAS
  datpath <- "~/Dropbox/USDA-compost/Data/PlantTraits/greenhouse_traits/"
}

# list files in entered data folder
datfiles <- dats <- list.files(paste0(datpath, "GH_EnteredData"), full.names = T)

# read in raw water uptake data
uptake<- read.csv(paste0(datpath, "GH_EnteredData/greenhouse_water-uptake.csv"), na.strings = na_vals, strip.white = T)
uptake_control <- read.csv(paste0(datpath, "GH_EnteredData/greenhouse_water-uptake-controls.csv"), na.strings = na_vals, strip.white = T)
root <- read.csv(paste0(datpath, "GH_EnteredData/Root_data_2.csv"), na.strings = na_vals, strip.white = T)

#calculate total water loss in the controls due to evaporation
dat_control<-uptake_control %>% mutate(starttime = ymd_hm(paste(uptake_control[,3], uptake_control[,4])), 
                                   endtime = ymd_hm(paste(uptake_control[,6], uptake_control[,7])),
                                   interval = difftime(endtime,starttime,units = "mins"),
                                   loss= (wet_wt_g-dry_wt_g), #calculate total loss
                                   loss_hr=(loss/as.numeric(interval))*60) #normalize to loss per hour

ggplot(dat_control, aes(x=date_dry, y=loss_hr, color=as.factor(Rep)))+
  #geom_errorbar(aes(ymax = per_germ+se, ymin = per_germ - se), position = position_dodge(width = 0.5), width = 0.1) +
  #geom_point( position = position_dodge(width = 0.5)) +
  geom_point()+
  labs(y = "Water Loss Per Hour", x="") +
  #scale_color_manual(values = c( "indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Treatment"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# average evaporation by day
sum_control <- dat_control %>% group_by(date_dry)%>% summarize(mean_evap_hr = mean(loss_hr)) 

#calculate total water loss in the samples due to evaporation + uptake
dat <- uptake %>% mutate(starttime = ymd_hm(paste(uptake[,5], uptake[,6])), 
                         endtime = ymd_hm(paste(uptake[,8], uptake[,9])),
                         interval = difftime(endtime,starttime,units = "mins"),
                         loss= (wet_wt_g-dry_wt_g),
                         loss_hr=(loss/as.numeric(interval))*60)

#join with evaporation data
dat <- left_join(dat, sum_control, by="date_dry")

#calculate total uptake by subtracting water loss due to evaporation
dat <- dat %>% mutate(uptake_hr=loss_hr-mean_evap_hr)

#join with prelim root scan data
dat <- left_join(dat, root, by=c("Code", "Rep"))




ggplot(dat, aes(x=Code, y=uptake_hr))+
  #geom_errorbar(aes(ymax = per_germ+se, ymin = per_germ - se), position = position_dodge(width = 0.5), width = 0.1) +
  #geom_point( position = position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, color = "red")+
  geom_boxplot()+
  labs(y = "Water Uptake (g/hr)", x="") +
  #scale_color_manual(values = c( "indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Treatment"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dat_sum<- dat %>% filter(Length.cm.>0&dry_wt_g>0)

ggplot(dat_sum, aes(x=SurfArea.cm2., y=uptake_hr, color=Code))+
  #geom_errorbar(aes(ymax = per_germ+se, ymin = per_germ - se), position = position_dodge(width = 0.5), width = 0.1) +
  #geom_point( position = position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, color = "red")+
  geom_point()+
  labs(y = "Water Uptake (g/hr)") +
  #scale_color_manual(values = c( "indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Treatment"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(dat_sum, aes(x=RootVolume.cm3., y=uptake_hr, color=Code))+
  #geom_errorbar(aes(ymax = per_germ+se, ymin = per_germ - se), position = position_dodge(width = 0.5), width = 0.1) +
  #geom_point( position = position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, color = "red")+
  geom_point()+
  labs(y = "Water Uptake (g/hr)") +
  #scale_color_manual(values = c( "indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Treatment"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(dat_sum, aes(x=AvgDiam.mm., y=uptake_hr, color=Code))+
  #geom_errorbar(aes(ymax = per_germ+se, ymin = per_germ - se), position = position_dodge(width = 0.5), width = 0.1) +
  #geom_point( position = position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, color = "red")+
  geom_point()+
  labs(y = "Water Uptake (g/hr)") +
  #scale_color_manual(values = c( "indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Treatment"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(dat_sum, aes(x=Tips, y=uptake_hr, color=Code))+
  #geom_errorbar(aes(ymax = per_germ+se, ymin = per_germ - se), position = position_dodge(width = 0.5), width = 0.1) +
  #geom_point( position = position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, color = "red")+
  geom_point()+
  labs(y = "Water Uptake (g/hr)") +
  #scale_color_manual(values = c( "indianred4",  "dodgerblue1", "darkgoldenrod"), guide = guide_legend(title = "Treatment"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#clean species list
spp<- read.csv(paste0(datpath, "GH_EnteredData/species_list_master.csv"), na.strings = na_vals, strip.white = T)

spp <- unique(spp[ , 1:6 ])
spp <- arrange(spp, code) 

#write_csv(spp, paste0(datpath, "GH_CleanedData/species_list_master.csv"))
