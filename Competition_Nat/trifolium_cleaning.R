# starter script to clean and format competition trifolium data 

# notes: 
# 1/27/22: quick code so NK and CA can begin preliminary data analysis 
# EAS will clean up and expand code later to integrate with the remaining species' competition data 
## MAJOR TO-DO: check and enter TRHI data from the datasheets for the backgrounds and other species phytos  
#From background datasheet: #note that TRHI background/competitor data STILL missing for plot 12
#From main competition datasheets: TRHI phytometer/invader data missing for plots 21-32

# Important note about seed weights and seeding weights from repo wiki ("2019 10_07_seeding expt setup" page):
# "All seeds were pre-weighed to 8g/m2, which is 2g per half meter squared of raw seed weight (not including husk or awns)"
#example change

# -- SETUP ----
rm(list = ls()) # clean enviroment
# libraries needed
library(tidyverse)
# change default settings
na_vals <- c("", " ", NA, "NA") #tells it what to call NA
options(stringsAsFactors = F) #setting for strings of characters 'abc'
theme_set(theme_bw()) #for plotting graphs later

# specify dropbox pathway (varies by user -- EAS can tweak this for CA or NK use)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/Competition/")){
  ## CTW pathway
 datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/Competition/"
}else{
  ## LMH and EAS
  datpath <- "~/Dropbox/USDA-compost/Data/Competition/"
}

# -- DATA IMPORT (NAT) ----
#1. set data pathway 
datpath <- "~/Desktop/USDA-compost/Data/Competition/"  ##NAT TO EDIT 

# list files in entered data folder
datfiles <- dats <- list.files(paste0(datpath, "Competition_EnteredData"), full.names = T)


# Data import (Carmen)
datpath <- "~/Research/Nat_Thesis/USDA-compost/Data/Competition/" 
datfiles <- dats <- list.files(paste0(datpath, "Competition_EnteredData"), full.names = T)



# read in raw trifolium data
trif<- read.csv(paste0(datpath, "Competition_EnteredData/trifolium_seeds_competition_summer2021.csv"), na.strings = na_vals, strip.white = T)
comp<-read.csv(paste0(datpath, "Competition_EnteredData/competition_seeds_summer2021_comp.csv"), na.strings = na_vals, strip.white = T)
phyto<-read.csv(paste0(datpath, "Competition_EnteredData/competition_seeds_summer2021_phyto.csv"), na.strings = na_vals, strip.white = T)
# read in treatment key
trtkey <- read.csv(paste0(datpath, "comp_trt_key.csv"), na.strings = na_vals, strip.white = T)
#read in field clip data to check total stems
clip_phyto<- read.csv(paste0(datpath, "Competition_EnteredData/competition_spring2020_phytometers_clipped.csv"), na.strings = na_vals, strip.white = T)
clip_back<- read.csv(paste0(datpath, "Competition_EnteredData/competition_spring2020_background_clipped.csv"), na.strings = na_vals, strip.white = T)
#read in background density
dens<- read.csv(paste0(datpath, "Competition_CleanedData/competitor_meandensity_jan2020.csv"), na.strings = na_vals, strip.white = T)


#quick check
glimpse(trif) #check if anything looks weird, check which data are numeric vs characters

#### PREP Trifolium data---
trif$phytonum <- as.numeric(trif$phytonum) #convert phytonum to numeric from character

#check range of values for plot (1-36), subplot (1-7), and phytonum (1-6)
range(trif$plot) #checks out, goes to 36
range(trif$subplot) #why does this go to 10? should be 7
trif$subplot[trif$subplot == 10] <- 6  #checking with the treatment key, this is entered wrong and should be subplot 6
trif$phytonum #NAs present for competitors (when TRHI is background, will fix/replace below with trtkey file) 

#remove notes and phytonum columns (with NA errors)
trif <- trif %>% select(-survey_date, -survey_init, -notes, -X, -X.1, -X.2, -phytonum) 
comp <- comp %>% select(-survey_date, -survey_init, -notes, -mature.flowers, -immature.flowers, -seed.mass..g., -biomass.w.o.stems.g.) 
#phyto <- phyto %>% select(-survey_date, -survey_init, -notes, -seed.mass..g., -biomass.w.o.stems.g., -Mature.Flowers..if.applicable.,-Immature.flowers..if.applicable., -X, -X.1, -X.2, -X.3) 


#bring in data recorded for how many stems were clipped in the field
clip_phyto<- clip_phyto %>% rename(background = comp., subplot=sub.plot) %>% 
  filter(phyto=="TRHI") %>% select(-page,-jan2020_stems,-jan2020_notes,-spr2020_notes,-clip_init,-date_clip)

clip_back<- clip_back %>% rename(background = competitor, spr2020_stems = count_clipped) %>% 
  filter(background=="TRHI") %>% select(-page,-clip_init,-clip_date, -order, -notes)

clip_check<-bind_rows(clip_phyto,clip_back)
clip_check <- clip_check %>% select(-phyto_num,-phyto) 
clip_check<-left_join(clip_check,trtkey)

trif <- left_join(trif, clip_check) #compare tot_stems data with spr2020_stems data
#use spr2020_stems data for future analyses?  this is what was recorded as clipped in the field

#get any trifolium data from full competition data sheets
comp <- comp %>% mutate(tot_seeds=as.numeric(seeds), tot_stems=as.numeric(tot_stems)) %>% 
  filter(background=="TRHI"&tot_seeds!="NA") %>% select(-seeds)


#summarize by summing the total seeds, mass, stems for each phytometer
trif_clean <- trif %>% group_by(nut_trt, ppt_trt, block, plot, subplot, phytonum, phyto, background ) %>% 
  summarize(tot_seeds=sum(as.numeric(seeds), na.rm = TRUE), tot_seed_mass=sum(as.numeric(seed_mass),na.rm = TRUE), 
            tot_stems=max(as.numeric(spr2020_stems), na.rm = TRUE), #change this to tot_stems from spr2020_stems?
            tot_stem_mass=sum(as.numeric(stem_mass),na.rm = TRUE), tot_mass=tot_seed_mass + tot_stem_mass)

#add missing trifolium background data
trif_clean2<-left_join(trif_clean, comp[, c("plot", "subplot", "tot_seeds","tot_stems")], by=c("plot", "subplot"))
trif_clean<-trif_clean2 %>% mutate(tot_seeds = coalesce(tot_seeds.x,tot_seeds.y), tot_stems=coalesce(tot_stems.x,tot_stems.y)) %>%
  select(-tot_seeds.x,-tot_seeds.y, -tot_stems.x,-tot_stems.y)

#### calculate quantity seeded for densities
# load seed mass data to crunch qty seeded
# > reading in others dats for now (do we have our own measurements?)
# read in J. Larson dry seed mass to screen for overcounts in density (more likely spp present in background seed bank so density enhanced)
seed_mass <- read.csv("Data-cleaning/Larson_CA_dryseedmass.csv")

# -- PREP LOOKUP TABLES ----
# make lookup table for subsample frame to scale to meter-square
scale_lt <- data.frame(subsample_cm = c("5x5", "10x10", "25x25", "50x50"),
                       area_cm2 = c(5*5, 10*10, 25*25, 50*50)) %>%
  # make half plot scale for reality check with density (i.e. does stems projected for 50x50cm exceed amount seeded?)
  mutate(scale_half = (50*50)/area_cm2,
         # full meter scale factor
         scale_m2 = (100*100)/area_cm2)

# join seed mass data to LUT to project max density possible per species
# > note: here is where to adjust based on whether seeds weighed out with or without awns/husks/attachments (looking at you ERBO)
seed_lt <- subset(seed_mass, grepl("avef|erob|lolm|taec|trih", species)) %>% #start with Julie's seed weights
  group_by(species) %>%
  summarise(Seed = mean(perseedwt)) %>%
  ungroup() %>%
  rename("background" = "species") %>%
  #mutate(source = "JL") %>%  # add jl for source
  mutate(background = casefold(background, upper = T),
         background = paste0(substr(background, 1,2), substr(background, 4,5)),
         # scale to half meter density -- seeded at 8g per m2 (2g per half m2)
         max_density_halfm2 = 2/Seed) 

#check for missing data
check<-trtkey%>% filter(phyto=="TRHI")
trif_clean<-left_join(check,trif_clean) 

trif_clean<- trif_clean %>%
  #sort by plot - check again if all are present (1-36, no 33)
  arrange(plot) %>%
  #specify if competitor or invader
  mutate(role = "invader", role = if_else(background =="TRHI", "competitor", 
                                          if_else(background=="Control","control",role))) 

#note that TRHI background/competitor data STILL missing for plot 12
#TRHI phytometer/invader data missing for plots 21-32

#add competitor/background density
trif_clean<-merge(trif_clean,dens) 
trif_clean<-left_join(trif_clean, seed_lt[, c("background", "max_density_halfm2")], by="background")
trif_clean <- trif_clean %>% select(-nobs,-mean_density_1m2)

#standardize data to seeds/stem
#NOTE: decide what to do with NAs 
trif_sum<-trif_clean %>% group_by(nut_trt, ppt_trt, block, plot, subplot, phytonum, phyto, background, role)%>%
  summarize(output=tot_seeds/tot_stems) #need to check NaNs for when there are 0 stems

#for Nat to do preliminary analyses (temporary file until we update data and later combine with full competition dataset)
write.csv(trif_clean, paste0(datpath, "Competition_CleanedData/competition_trifolium_seeds_2021.csv"))  


#VISUALIZE ----
# specify plotting cols
ppt_cols <- c(D = "brown", W = "dodgerblue", XC = "seagreen")
plant_cols <- c(AVBA = "darkgreen", HOMU = "lightgreen", TACA = "limegreen", LOMU = "blue", 
                ERBO  = "red", TRHI = "orchid", Control = "grey40")

#Visualize data, summarizing mean seed output by treatment and background competitor
trif_plot1<-trif_sum %>% group_by(nut_trt, ppt_trt, phyto, background, role)%>%
  summarize(se_output = sd(output, na.rm=T)/sqrt(length(output)),output=mean(output,na.rm=T))

ggplot(data=trif_plot1, aes(nut_trt, output, col = background)) +
  geom_errorbar(aes(ymax = output + se_output, ymin = output - se_output), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  scale_color_manual(values = plant_cols) +
  ggtitle("Mean Trifolium Seeds +- 1SE") +
  facet_grid(.~ppt_trt) # can also facet by background, trying alltog with colors to compare more directly

#Visualize data, summarizing mean seed output by treatment only (ignore competitor)
trif_plot2<-trif_sum %>% group_by(nut_trt, ppt_trt)%>%
  summarize(se_output = sd(output, na.rm=T)/sqrt(length(output)),output=mean(output,na.rm=T))

ggplot(data=trif_plot2, aes(nut_trt, output)) +
  geom_errorbar(aes(ymax = output + se_output, ymin = output - se_output), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  ggtitle("Mean Trifolium Seeds +- 1SE") +
  facet_grid(.~ppt_trt)


#### FOR LATER - calculate quantity seeded for densities
# load seed mass data to crunch qty seeded
# > reading in others dats for now (do we have our own measurements?)
# read in J. Larson dry seed mass to screen for overcounts in density (more likely spp present in background seed bank so density enhanced)
seed_mass <- read.csv("Data-cleaning/Larson_CA_dryseedmass.csv")


# -- PREP LOOKUP TABLES ----
# make lookup table for subsample frame to scale to meter-square
scale_lt <- data.frame(subsample_cm = c("5x5", "10x10", "25x25", "50x50"),
                       area_cm2 = c(5*5, 10*10, 25*25, 50*50)) %>%
  # make half plot scale for reality check with density (i.e. does stems projected for 50x50cm exceed amount seeded?)
  mutate(scale_half = (50*50)/area_cm2,
         # full meter scale factor
         scale_m2 = (100*100)/area_cm2)

# join seed mass data to LUT to project max density possible per species
# > note: here is where to adjust based on whether seeds weighed out with or without awns/husks/attachments (looking at you ERBO)
seed_lt <- subset(seed_mass, grepl("avef|erob|lolm|taec|trih", species)) %>% #start with Julie's seed weights
  group_by(species) %>%
  summarise(Seed = mean(perseedwt)) %>%
  ungroup() %>%
  rename("ID" = "species") %>%
  # add jl for source
  mutate(source = "JL") %>%
  mutate(ID = casefold(ID, upper = T),
         ID = paste0(substr(ID, 1,2), substr(ID, 4,5)),
         # scale to half meter density -- seeded at 8g per m2 (2g per half m2)
         max_density_halfm2 = 2/Seed) 


