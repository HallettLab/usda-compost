# adjust and infill usda-compost ongoing soil mositure
# author(s): CTW
# questions?: caitlin.t.white@colorado.edu


# script purpose:
# read in timestamp-corrected and QA'd ongoing soil mositure dataset
# ID timeseries that have low warnings (vwc slightly below 0)
# shift sub-0 periods so lowest point at 0 (these happen in the summer time when soil dry)

# (not added but one day or anyone else can): infill missing vwc


# notes:
# ctw only working on soil moisture through jan 2022 then needs to shift focus
# if no one on project is using soil moisture at high resolution and specific lines don't matter, infilling missing data is not a critical need
# e.g., there are several sensors whose data are missing june 2020 to end of data collection (sep 2021) and no soil moisture data was collected for block 4
# > i.e., infilling and interpolating would be a big effort
# shifting up low values to make summer data usable is relatively easy and perhaps preferable to NAing data
# > other option would be to set sub-0 values to 0, but then miss diurnal fluctation
# > values sub-0 are within accuracy of decagon logger, so shifting up seems preferable (i.e., pattern real, just a little off from sensors getting wacky in super dry environment)



# -- SETUP -----
# clear environment
rm(list=ls())
# load libraries needed
library(tidyverse) # for dplyr, tidyr, and ggplot
library(cowplot)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals <- c(" ", "", NA, "NA")

# set path to compost data (main level)
# > varies by user
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/"
}else{
  ## probs LMH and AS? (fix if not -- rproj is two levels down in github repo)
  datpath <- "~/Dropbox/USDA-compost/Data/"
}


# time-corrected and QA'd soil moisture data
smdat <- read.csv(paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/SoilMoisture_compiled2_time-corrected_flagged.csv"), na.strings = na_vals)
# cimis ppt aggregated to 2-hr soil moisture collection intervals
cimis <- read.csv(paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/CIMIS_084BrownsValley_pptforSoilMoisture.csv"), na.strings = na_vals)
# > note: timestamps and dates are read in as character by R but will try to shift using cleanorder (should work fine)
str(smdat)
str(cimis)


# -- PREP -----
# id portids that need adjustment
adjustid <- subset(smdat, grepl("adjust", flag_note), select = c(portid, fulltrt,ppt_trt, mon, waterYear)) %>%
  distinct()
length(unique(adjustid$portid)) # half of the sensors..

ggplot(subset(smdat, portid %in% unique(adjustid$portid)), aes(cleanorder, vwc, col = grepl("adjust", flag_note), group = portid)) +
  geom_line() +
  facet_grid(nut_trt~ ppt_trt)
# some just need first dats of series (at installation set to 0)

View(subset(smdat, portid %in% adjustid$portid & waterYear == 2019 & vwc < 0 & cleanorder < 1000))
# how many portids have sub-0 after install period
adjustseries <- with(subset(smdat, portid %in% adjustid$portid), unique(portid[(vwc<0 & !is.na(vwc) & cleanorder> 100)])) # 22 of 32
# pull out install adjustments corrections only
adjustearly <- unique(adjustid$portid[!adjustid$portid %in% adjustseries])


# -- TEST CASE -----
testcase <- adjustseries[21]
testorder <- 0000:15000
ggplot(subset(smdat, portid == testcase & cleanorder %in% testorder), aes(cleanorder, vwc, col = grepl("adjust", flag_note), col = portid)) +
  geom_line() +
  # are there any good references to help triage?
  geom_line(data = subset(smdat, !portid %in% adjustseries & fulltrt == unique(adjustid$fulltrt[adjustid$portid == testcase])),aes(col = portid), alpha = 0.5) #+
  facet_wrap(~block) #+
  theme(legend.position = "none")
