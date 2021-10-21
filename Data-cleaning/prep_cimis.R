# prep cimis data so can pair with compost data
# author(s): ctw (caitlin.t.white@colorado.edu)
# created: july 2020

# script purpose:
# read in latest downloaded CIMIS browns valley dataset (can't read in dynamically from web, must login to request data and download)
# treat dates and times so comparable with usda compost data dates and times
# write out to cleaned data


# notes:
# may also want to add in code later to summarize cimis data (e.g. collapse to every 2 hrs or daily.. could also do in an analysis script depending on needs)



# -- SETUP -----
# clear environment
rm(list=ls())
# load libraries needed
library(tidyverse) # for dplyr, tidyr, and ggplot
library(lubridate)
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

# read in CIMIS met data for QA checks after compilation
cimis <- read.csv(list.files(paste0(datpath, "CIMIS/Download_RawData"), full.names = T), na.strings = na_vals) 


# 1. Prep CIMIS data for comparison -----
# prep cimis
# remove row with all NAs
cimis_hrly <- cimis[!apply(cimis,1,function(x) all(is.na(x))),] %>%
  # clean up colnames
  rename_all(function(x) gsub("[.]{2}", ".",x)) %>%
  rename_all(function(x) gsub("[.]$", "", x)) %>%
  # clean up cols to pair with soilmoisture dat
  mutate(date = as.Date(Date, format = "%m/%d/%Y"),
         # adjust 24:00 date to next day (to match SFREC format)
         date = ifelse(Hour.PST == 2400, as.character(date +1), as.character(date)),
         # adjust Julian day
         # doy = ifelse(Hour.PST == 2400 & grepl("^12/31", Date), 1, 
         #              ifelse(Hour.PST == 2400, Jul+1, Jul)),
         time = ifelse(nchar(Hour.PST)==3, paste0(0, substr(Hour.PST,1,1), ":", substr(Hour.PST, 2,3), ":00"),
                       paste(substr(Hour.PST,1,2), substr(Hour.PST, 3, 4), "00", sep =":")),
         # change 24:00 to 00:00 to match SFRECdat
         time = gsub("24", "00", time),
         # clean_datetime = as.POSIXct(paste(date, time), format = "%Y-%m-%d %H:%M:%S"),
         # using as.character date-time because R doesn't know how to convert 2am timestamp in March when switch from PST to PDT
         # CIMIS is report in PST-converted times only, so will deal with time zone DST issue in QA or analysis later if needed
         clean_datetime = paste(date, time),
         date = as.Date(date))


# add info for water year and day of water year
# create lookup table for all dates, full water years
firstWYcimis <- as.Date(paste0(min(year(cimis_hrly$date))-1, "-10-1"))
lastWYcimis <- as.Date(paste0(max(year(cimis_hrly$date))+1, "-9-30"))
date_lookup_cimis <- data.frame(date = seq.Date(firstWYcimis, lastWYcimis, 1)) %>%
  # add in year, doy, join water year info, etc back in..
  mutate(mon = month(date), 
         doy = yday(date),
         waterYear = ifelse(mon %in% 10:12, year(date)+1, year(date))) %>%
  # create day of water year
  arrange(waterYear, date) %>%
  group_by(waterYear) %>%
  mutate(dowy = seq(1, length(date), 1)) %>%
  ungroup()

# add in water year then clean up cols
cimis_hrly <- left_join(cimis_hrly, date_lookup_cimis)
summary(cimis_hrly)
# check: any NAs present in cols where shouldn't be?
summary(is.na(cimis_hrly))



# -- WRITE OUT -----
# hourly dataset
write.csv(cimis_hrly, paste0(datpath, "CIMIS/CIMIS_084BrownsValley_hourly_prepped.csv"), row.names = F)

