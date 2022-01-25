# QAQC usda-compost ongoing soil moisture dataset
# author(s): CTW
# questions?: caitlin.t.white@colorado

# SCRIPT PURPOSE:
# read in timestamp-corrected soil moisture data (from soilmoisture1_correct_timestamps.R)
# process browns valley hourly cimis data to 2-hr timestamp intervals matching the soil moisture dataset 
# > include start of water year 2019 when project started (cleanorder for those are negative, count backwards from start of soil moisture dataset)
# > derive rain and dry events over lifetime of project and use to determine start and end of growing season and summer dry seasons 
# adjust sensor data where bulk of timeseries dropped below 0 so more flagging functions will apply correctly
# > majority of shifts upwards for summer drops will be done in next workflow script, just coarsely treating substantial drops here
# create flagging functions based on references (see notes)
# test flagging functions on a few sensors that have a complete timeseries for the whole project and visually appear reasonable
# run flagging functions on full dataset
# screen flagged values to determine rules for removal
# NA values that meet criteria for removal and annotate flag_note based on reason
# follow with more manual checks to flag values that were missed by flagging functions
# > e.g., flagging functions that assess timestep changes will not work for observations preceeded by missing data
# write out intermediate datasets (e.g., precip, raw flags) and QA'd dataset to DataQA folder
# > also write out QA'd dataset as cleaned dataset until full workflow developed (e.g., sub-0 adjustments and infilling)


# NOTES:
# there's a cimis R package! 'cimir'. CTW isn't using it, but good to know. can read in data dynamically instead of download static files.
# https://cran.r-project.org/web/packages/cimir/cimir.pdf

# in the process of reviewing flags, visually determined several sensors were labeled with incorrect treatments (it's in the raw field notes)
# this script corrects those and also writes out the corrected logger key dataset
# because previous two scripts built on treatments as they were assigned in the raw notes, to be on the safe side (not have to re-write code), will keep raw treatment assignments as is and make the fix here

# REFERENCES (for flagging):
# Dorigo et al. 2013:
# Dorigo, W., Xaver, A., Vreugdenhil, M., Gruber, A., Hegyiová, A., Sanchis-Dufau, A., Zamojski, D., Cordes, C., Wagner, W. and Drusch, M. (2013), 
# Global Automated Quality Control of In Situ Soil Moisture Data from the International Soil Moisture Network. Vadose Zone Journal, 12: 1-21 vzj2012.0097. 
# https://doi.org/10.2136/vzj2012.0097

# Quiring et al. 2016:
# Quiring, S. M., Ford, T. W., Wang, J. K., Khong, A., Harris, E., Lindgren, T., Goldberg, D. W., & Li, Z. (2016). 
# The North American Soil Moisture Database: Development and Applications, Bulletin of the American Meteorological Society, 97(8), 1441-1459. 
# Retrieved Jan 25, 2022, from https://journals.ametsoc.org/view/journals/bams/97/8/bams-d-13-00263.1.xml



# -- SETUP -----
# clear environment
rm(list=ls())
# load libraries needed
library(tidyverse) # for dplyr, tidyr, and ggplot
library(lubridate)
library(cowplot)
library(pracma)
library(zoo)
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


# read in time-corrected soil moisture data
smdat <- read.csv(paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/SoilMoisture_compiled_time-corrected.csv"), na.strings = na_vals)

## prepped hourly CIMIS met data for QA checks after compilation
cimis_hrly <- list.files(paste0(datpath, "CIMIS"), full.names = T) %>%
  subset(grepl(".csv",.)) %>%
  read.csv(na.strings = na_vals) 

# how did both read in?
str(smdat)
str(cimis_hrly)
# convert date cols in all
smdat$clean_datetime <- as.POSIXct(smdat$clean_datetime, tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
smdat$date <- as.Date(smdat$date, format = "%Y-%m-%d")

cimis_hrly$clean_datetime <- as.POSIXct(cimis_hrly$clean_datetime,  tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
cimis_hrly$date <- as.Date(cimis_hrly$date, format = "%Y-%m-%d")
cimis_hrly$Date <- as.Date(cimis_hrly$Date, format = "%m/%d/%Y")



# -- PREVIEW RAIN V VWC, PREP DATA -----
plot_grid(
  ggplot(smdat, aes(clean_datetime, vwc, col = vwc < 0)) +
    geom_point(alpha = 0.5) +
    theme(legend.position = "none") +
    facet_grid(ppt_trt~.),
  ggplot(subset(cimis_hrly, date %in% smdat$date), aes(clean_datetime, Precip.mm)) +
    geom_point(),
  align = "v", nrow = 2)

# what are the ppt flags?
summary(as.factor(cimis_hrly$qc.1))
# from manual: R = "data is far our of historical limits", Y = "data moderately out of historical limits"
with(cimis_hrly, sapply(split(Precip.mm, qc.1), summary))
# plot with flags colored
ggplot(cimis_hrly, aes(clean_datetime, Precip.mm, col = qc.1)) +
  geom_point(alpha = 0.5) +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b-%y")
# flagged vals look like tiny precip in the summer months.. but doesn't look like it will make a big difference anywhere

# where are rain data missing?
ggplot(cimis_hrly, aes(clean_datetime, is.na(Precip.mm))) +
  geom_jitter(alpha = 0.5) +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b-%y") # not too bad. comes mostly towards or after end of project

# collapse precip to same temporal resolution as soil moisture
# > project not installed until november 2019, but rain season started before that
gplot(subset(cimis_hrly, year(date) == 2018)) + geom_point(aes(clean_datetime, Precip.mm)) # okay to start oct 1 2018 (beginning of water year)
# create smdat clean_datetime back to start of oct 1
smdat_datetime <- distinct(smdat[c("clean_datetime", "date", "cleanorder")])
start2018 <- data.frame(clean_datetime = seq.POSIXt(min(smdat_datetime$clean_datetime), as.POSIXct("2018-10-01 00:00:00", tz = "UTC"), -(2*60*60))) %>%
  mutate(cleanorder = -1:-nrow(.)) %>%
  arrange(clean_datetime) %>%
  mutate(date = date(clean_datetime)) %>%
  dplyr::select(names(smdat_datetime))
# stack
smdat_datetime <- rbind(start2018, smdat_datetime) %>%
  # add waterYear info
  left_join(distinct(cimis_hrly[c("date", "mon","doy", "waterYear", "dowy")])) %>%
  mutate(yr = year(date)) %>%
  dplyr::select(clean_datetime:doy, yr, waterYear, dowy)
cimis_ppt2hr <- subset(cimis_hrly, date%in% smdat_datetime$date, select = c(clean_datetime, date, time, mon:dowy, Precip.mm)) %>%
  left_join(smdat_datetime) %>%
  #fill cleanorder upevent to pair hourly data
  fill(cleanorder, .direction = "up") %>%
  #subset to 1 hour before soil moisture dat starts
  #subset(clean_datetime >= min(smdat$clean_datetime) - (60*60)) %>%
  # sum precip by cleanorder
  group_by(cleanorder) %>%
  summarise(ppt_mm = sum(Precip.mm, na.rm = T),
            nobs = length(Precip.mm)) %>%
  # drop anything that's not 2 obs per step
  filter(nobs == 2) %>%
  # add clean datetime and water year info back in
  left_join(distinct(smdat_datetime)) %>%
  # track wetup events in rainfall
  mutate(rain = ppt_mm > 0 | (ppt_mm == 0 & dplyr::lag(ppt_mm) > 0 & dplyr::lead(ppt_mm > 0)), # if 0 sandwiched between non-zero, count it in same event
         # but if negligible rain sandwiched between zeroes don't count
         rain2 = ppt_mm == 0.1 & dplyr::lag(ppt_mm) == 0  & dplyr::lag(ppt_mm,2) == 0 & dplyr::lead(ppt_mm) == 0 & dplyr::lead(ppt_mm,2) == 0,
         # adjust whether rain event or not based on rain2
         rain = ifelse(rain2, FALSE, rain),
         # create empty event cols for tracking rain and dry events
         rainevent = as.numeric(NA),
         dryspell = as.numeric(NA)) %>%
  group_by(waterYear) %>%
  mutate(wY_accumulated_ppt = cumsum(ppt_mm)) %>%
  ungroup() %>%
  # drop nobs
  dplyr::select(-nobs)
# track rain events based on rain (T/F) col
trackrain <- list2DF(rle(cimis_ppt2hr$rain))

r <- 1
rainevent <- 1
dryspell <- 1
for(i in 1:nrow(trackrain)){
  # note where to stop
  tempend <- r+ (trackrain$lengths[i] -1)
  # note event
  if(!trackrain$values[i]){
    cimis_ppt2hr$dryspell[r:tempend] <- dryspell
    # increment r and dryspell
    r <- tempend+1
    dryspell <- dryspell + 1
  }else{
    cimis_ppt2hr$rainevent[r:tempend] <- rainevent
    # increment r and rainevent
    r <- tempend+1
    rainevent <- rainevent + 1
  }
}

# track rainy season start and end each water year (for screening out dry period wetup spikes)
# 0.5" = amount for first germination (12.7mm)
# per water year (look at sep too), find first event where total precip is >= 12.7mm
# find dry period with longest stretch after spring (should correspond to dry period?)
# > note: sometimes it rains a teensy bit (<5mm) for a day or two in the summer. ignore these events.
scan_rainseason <- function(dat, threshold = 12.7, groupvars = c("rainevent", "dryspell")){
  # season start occurs on the day accumulation reaches 12.7 within a given event
  tempdat <- grouped_df(dat, groupvars) %>%
    mutate(event_ppt = sum(ppt_mm),
           event_accum = cumsum(ppt_mm),
           # within germevent, there has to be another rain event within a certain time, otherwise assume germinates would die
           # > 2 weeks?
           n_days = length(unique(date)), # this is just days involved, not actual event length       
           # grab as character or else R converts to number
           date_threshold = ifelse(event_ppt >= threshold, as.character(min(clean_datetime[(event_accum >= threshold)])), NA),
           # convert back to date
           date_threshold = as.POSIXct(date_threshold, tz = "UTC"),
           # grab days to next event,
           start_event = min(clean_datetime),
           end_event = max(clean_datetime),
           event_time= end_event - start_event, #in seconds
           event_hrs = as.numeric(event_time/(60*60)),
           event_days = round((event_hrs/24),2)) %>% 
    ungroup()
  
  # NOTE TO SELF: can fuss with details later -- events that are short duration (dry or wet) should be folder into larger events
  # for now just focusing on season demarcation
  eventdat <- distinct(subset(tempdat, select = c(rainevent, dryspell, event_ppt, n_days:ncol(tempdat)))) %>%
    gather(event_type, event_num, rainevent, dryspell) %>%
    mutate(start_mon = month(start_event),
           start_yr = year(start_event),
           eco_yr = ifelse(start_mon %in% 9:12, start_yr+ 1, start_yr),
           eco_mon = ifelse(start_mon %in% 9:12, start_mon - 8, start_mon + 4)) %>%
    # re-organize df
    subset(!is.na(event_num), select= c(event_type, event_num, eco_yr, start_yr, eco_mon, start_mon, event_ppt:event_hrs)) %>%
    group_by(event_type) %>%
    mutate(next_event = dplyr::lead(start_event) - end_event,
           next_event_days = round(as.numeric(next_event/24),2),
           next_ppt = dplyr::lead(event_ppt)) %>% # in hours
    ungroup()
  
  # create growing season via for loop
  seasondat <- distinct(subset(eventdat, !is.na(date_threshold), select = eco_yr))
  seasondat$season_start <- NA
  seasondat$season_end <- NA
  
  for(y in unique(seasondat$eco_yr)){
    scandat <- subset(eventdat, eco_yr == y & event_type == "rainevent")
    startchoices <- as.character(with(scandat, date_threshold[!is.na(date_threshold) & next_event_days < 5]))
    for(u in startchoices){
      u <- as.POSIXct(u, tz = "UTC")
      # check if rain fell the next month also
      lookahead <- subset(scandat, start_event > u & start_event <= date(u)+30)
      # if nothing present or it didn't rain in the next month, next
      # > note: germ rain in sep 2019 and rained several times that month, but then didn't rain until nov 27 (assume germinates died)
      if(nrow(lookahead) == 0 | !(month(u) +1) %in% lookahead$start_mon){
        next
      }
      # if there are more germinating events start season
      if(any(!is.na(lookahead$date_threshold)) | sum(lookahead$event_ppt) > 10){
        seasondat$season_start[seasondat$eco_yr == y] <- as.character(u)
        break
      }else{next}
    }
    # search for the season end after march
    endchoices <- as.character(with(scandat, end_event[eco_mon > 7 & (next_event_days > 5 |  next_ppt < 1)]))
    for(e in endchoices){
      e <- as.POSIXct(e, tz = "UTC")
      # check if rain fell the next month also
      lookahead <- subset(scandat, start_event > e)
      # these are sort of arbitrary criteria -- just trying to match what I know should be the end date, can improve rules later
      if(any(!is.na(lookahead$date_threshold)) | any(lookahead$event_ppt > 10) | sum(lookahead$event_ppt) > 10){
        next
      }
      seasondat$season_end[seasondat$eco_yr == y] <- as.character(e)
      break
    }
  }
  # clean up
  datout <- subset(eventdat,!is.na(eco_yr), select = c(event_type:start_yr ,start_mon:event_ppt, date_threshold:end_event, event_hrs)) %>%
    left_join(subset(seasondat))
  
  return(datout)
}

gsdat <- scan_rainseason(cimis_ppt2hr)
# add end of watering treatments to growing season rain event dat
gsdat$stopwater <- as.Date(with(gsdat, ifelse(eco_yr == 2021, "2021-04-15",
                                              ifelse(eco_yr == 2020, "2020-04-15", "2019-04-30"))))
# extract season dat
seasondat <- distinct(subset(gsdat, eco_yr < 2022, select = c(eco_yr, season_start:stopwater)))
# fix 00:00:00 time
seasondat$season_start <- with(seasondat, ifelse(!grepl(":", season_start), paste(season_start, "00:00:00"), season_start))

raw_v_ppt <- plot_grid(
  ggplot(smdat, aes(clean_datetime, vwc)) +
    geom_line(aes(group = portid), alpha = 0.5) +
    geom_hline(aes(yintercept = 0), col = "red") +
    scale_x_datetime(date_breaks = "2 months", date_labels = "%b-%y") +
    ggtitle("USDA Compost QA: time-corrected soil moisture (by water treatment) vs. observed CIMIS Browns Valley precip") +
    labs(x = NULL, subtitle = Sys.Date()) +
    theme(axis.text.x = element_blank()) +
    facet_grid(ppt_trt~.),
  # start cimis data at compost project start
  ggplot(distinct(subset(cimis_ppt2hr, cleanorder >= 1, select = c("clean_datetime", "ppt_mm"))), aes(clean_datetime, ppt_mm)) +
    geom_point(alpha = 0.5) +
    scale_x_datetime(date_breaks = "2 months", date_labels = "%b-%y") +
    # create strip so lines up
    facet_grid("precip"~.),
  nrow = 2, align = "vh",
  rel_heights = c(1.25, 0.75)
)
raw_v_ppt

ggsave(filename = paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/PrelimQA_Figures/Compost_timecorrectedVWC_vPrecip.pdf"),
       plot = raw_v_ppt,
       width = 10, height = 6, units = "in")

# clean up environment
rm(trackrain, i, dryspell, rainevent, star2018, tempend, r)



# -- SHIFT SERIES DROPS BEFORE FLAGGING ----
# there are a few cases where the timeseries dropped and remained at that level
# to not lose data, try to shift dropped series and note adjustment
# ID cases first
ggplot(smdat, aes(clean_datetime, vwc, group = portid)) +
  geom_line() +
  facet_wrap(~logger)
ggplot(smdat, aes(clean_datetime, vwc, group = portid, col = factor(block))) +
  geom_line(alpha = 0.5) +
  facet_grid(ppt_trt~nut_trt)
# three cases in ND, NXC and FW (1 each) -- look at lowest vals in Sep 2021 for who
drops <- subset(smdat, fulltrt %in% c("ND", "NXC", "FW") & yr == 2021 & mon == 9, select = c(portid, fulltrt, vwc, yr, mon)) %>%
  distinct() %>%
  grouped_df("fulltrt") %>%
  filter(vwc == min(vwc, na.rm = T))
# confirm correct cases
ggplot(smdat, aes(clean_datetime, vwc, group = portid, col = factor(block))) +
  geom_line(alpha = 0.5) +
  geom_line(data = subset(smdat, portid %in% drops$portid), aes(clean_datetime, vwc, group = portid), col = "black", lty = 2) +
  facet_grid(ppt_trt~nut_trt) # yep

# go case by case. diff series to look for the breaks
drops_df <- subset(smdat, portid %in% drops$portid) %>%
  mutate(adjust_vwc = vwc) %>%
  # try to ID bigger breaks
  group_by(portid) %>%
  mutate(diffvwc = vwc - lag(vwc),
         breakvwc = vwc < 0 & diffvwc <0 & lag(vwc > 0),
         absdiff = abs(diffvwc)) %>%
  arrange(portid, diffvwc) %>%
  mutate(rankdiff = rank(diffvwc)) %>%
  ungroup()

# XC, B1L2_5
drop1 <- drops$portid[1]
ggplot(subset(smdat, fulltrt %in% drops$fulltrt[1]), aes(cleanorder, vwc, group = portid)) +
  ggtitle(drops$portid[1])+
  geom_line(alpha = 0.5) +
  geom_line(data = subset(smdat, portid == drop1), aes(cleanorder, vwc), col = "red", alpha =.5) +
  geom_text(data = subset(drops_df, portid == drop1 & rankdiff < 6), aes(label = rankdiff), col = "blue", alpha = 0.7) #diffvwc == min(diffvwc, na.rm = T)
# 1+2 diffs = adjustment for series drop starting at 2
adjust5 <- abs(with(drops_df, diffvwc[portid == drop1 & rankdiff == 5]))
adjust1 <- abs(with(drops_df, diffvwc[portid == drop1 & rankdiff == 1]))
adjust2 <- abs(with(drops_df, diffvwc[portid == drop1 & rankdiff == 2]))
break5 <- with(drops_df, cleanorder[portid == drop1 & rankdiff == 5])
break1 <-  with(drops_df, cleanorder[portid == drop1 & rankdiff == 1])
break2 <-  with(drops_df, cleanorder[portid == drop1 & rankdiff == 2])

# review 
suborder <- c(8000:12000)
ggplot(subset(smdat, fulltrt %in% drops$fulltrt[1] & cleanorder %in% suborder), aes(cleanorder, vwc, col = factor(block), group = portid)) +
  ggtitle(drops$portid[1])+
  geom_line(alpha = 0.5) +
  geom_line(data = subset(smdat, portid == drop1 & cleanorder %in% suborder), aes(cleanorder, vwc), col = "red", alpha =.5) +
  geom_line(data = subset(drops_df, portid == drop1 & cleanorder %in% suborder), aes(cleanorder, adjust_vwc), col = "orange") +
  #geom_text(data = subset(drops_df, portid == drop1 & cleanorder %in% suborder), aes(cleanorder, adjust_vwc, label = cleanorder), col = "purple") +
  geom_text(data = subset(drops_df, portid == drop1 & cleanorder %in% suborder & rankdiff < 6), aes(label = rankdiff), col = "blue", alpha = 0.7)
# 6358 and 6359 are first drops
# 6398 and 6399 are next drops
# 8928 and 8929 are next (ranks 1 and 2)
# 11000 and 11001 jump up
# 11028 and 11029 jump up
# 11243 is last jump up
adjust6358 <- abs(with(drops_df, diffvwc[portid == drop1 & cleanorder == 6358]))
adjust6359 <- abs(with(drops_df, diffvwc[portid == drop1 & cleanorder == 6359]))
adjust6398 <- abs(with(drops_df, diffvwc[portid == drop1 & cleanorder == 6398]))
adjust6399 <- abs(with(drops_df, diffvwc[portid == drop1 & cleanorder == 6399]))
adjust8928 <- abs(with(drops_df, diffvwc[portid == drop1 & cleanorder == 8928]))
adjust8929 <- abs(with(drops_df, diffvwc[portid == drop1 & cleanorder == 8929]))
adjust11000 <- abs(with(drops_df, diffvwc[portid == drop1 & cleanorder == 11000]))*1.5
adjust11001 <- abs(with(drops_df, diffvwc[portid == drop1 & cleanorder == 11001]))*1.5
# adjust11028 <- abs(with(drops_df, diffvwc[portid == drop1 & cleanorder == 11028]))
# adjust11029 <- abs(with(drops_df, diffvwc[portid == drop1 & cleanorder == 11029])) # to 11080
adjust11243 <- abs(with(drops_df, diffvwc[portid == drop1 & cleanorder == 11243]))
# need to adjust in order
drops_df$adjust_vwc[drops_df$portid == drop1 & drops_df$cleanorder >= 6358] <- (with(drops_df, vwc[portid == drop1 & cleanorder >= 6358]) + adjust6358)
drops_df$adjust_vwc[drops_df$portid == drop1 & drops_df$cleanorder >= 6359] <- (with(drops_df, adjust_vwc[portid == drop1 & cleanorder >= 6359]) + adjust6359)
drops_df$adjust_vwc[drops_df$portid == drop1 & drops_df$cleanorder >= 6398] <- (with(drops_df, adjust_vwc[portid == drop1 & cleanorder >= 6398]) + adjust6398)
drops_df$adjust_vwc[drops_df$portid == drop1 & drops_df$cleanorder >= 6399] <- (with(drops_df, adjust_vwc[portid == drop1 & cleanorder >= 6399]) + adjust6399)
drops_df$adjust_vwc[drops_df$portid == drop1 & drops_df$cleanorder >= 8928] <- (with(drops_df, adjust_vwc[portid == drop1 & cleanorder >= 8928]) + adjust8928)
drops_df$adjust_vwc[drops_df$portid == drop1 & drops_df$cleanorder >= 8929] <- (with(drops_df, adjust_vwc[portid == drop1 & cleanorder >= 8929]) + adjust8929)
drops_df$adjust_vwc[drops_df$portid == drop1 & drops_df$cleanorder >= 11000] <- (with(drops_df, adjust_vwc[portid == drop1 & cleanorder >= 11000]) - adjust11000)
drops_df$adjust_vwc[drops_df$portid == drop1 & drops_df$cleanorder >= 11001] <- (with(drops_df, adjust_vwc[portid == drop1 & cleanorder >= 11001]) - adjust11001)
#drops_df$adjust_vwc[drops_df$portid == drop1 & drops_df$cleanorder %in% c(11029:11080)] <- (with(drops_df, adjust_vwc[portid == drop1 & cleanorder %in% c(11029:11080)]) - adjust11029)
drops_df$adjust_vwc[drops_df$portid == drop1 & drops_df$cleanorder >= 11243] <- (with(drops_df, adjust_vwc[portid == drop1 & cleanorder >= 11243]) - adjust11243)
# as much as I can do with justification

# 2nd -- just need to bring end up
drop2 <- drops$portid[2]
suborder <- c(9000:16000)
ggplot(subset(smdat, fulltrt %in% drops$fulltrt[2] & cleanorder %in% suborder), aes(cleanorder, vwc, col = factor(block), group = portid)) +
  ggtitle(drops$portid[1])+
  geom_line(alpha = 0.5) +
  geom_line(data = subset(smdat, portid == drop2 & cleanorder %in% suborder), aes(cleanorder, vwc), col = "red", alpha =.5) +
  geom_line(data = subset(drops_df, portid == drop2 & cleanorder %in% suborder), aes(cleanorder, adjust_vwc), col = "orange") +
  #geom_text(data = subset(drops_df, portid == drop2 & cleanorder %in% suborder), aes(cleanorder, adjust_vwc, label = cleanorder), col = "purple") +
  geom_text(data = subset(drops_df, portid == drop2 & cleanorder %in% suborder & rankdiff < 6), aes(label = rankdiff), col = "blue", alpha = 0.7)

# how far is lowest point from 0?
drops_df$adjust_vwc[drops_df$portid == drop2]<- drops_df$vwc[drops_df$portid == drop2] 
adjustdrop2 <- 0- with(drops_df, min(vwc[cleanorder > 12000 & portid == drop2], na.rm = T))
drops_df$adjust_vwc[drops_df$portid == drop2 & drops_df$cleanorder >= 10658] <- (with(drops_df, adjust_vwc[portid == drop2 & cleanorder >= 10658]) + adjustdrop2)
# looks reasonable

# 3rd 
drop3 <- drops$portid[3]
suborder <- c(10000:16000)
ggplot(subset(smdat, fulltrt %in% drops$fulltrt[3] & cleanorder %in% suborder), aes(cleanorder, vwc,lty = comp_trt, col = factor(block), group = portid)) +
  ggtitle(drops$portid[1])+
  geom_line(alpha = 0.5) +
  geom_line(data = subset(smdat, portid == drop3 & cleanorder %in% suborder), aes(cleanorder, vwc, lty = comp_trt), col = "red", alpha =.5) +
  geom_line(data = subset(drops_df, portid == drop3 & cleanorder %in% suborder), aes(cleanorder, adjust_vwc, lty = comp_trt), col = "orange") +
  #geom_text(data = subset(drops_df, portid == drop3 & cleanorder %in% suborder), aes(cleanorder, adjust_vwc, label = cleanorder), col = "purple") +
  geom_text(data = subset(drops_df, portid == drop3 & cleanorder %in% suborder & rankdiff < 6), aes(label = rankdiff), col = "blue", alpha = 0.7)

# fix summer 2020 drop first
# 6765 and then 6766
# 6837 then 6838
drops_df$adjust_vwc[drops_df$portid == drop3]<- drops_df$vwc[drops_df$portid == drop3] 
adjust6765 <- abs(with(drops_df, diffvwc[portid == drop3 & cleanorder == 6765]))
adjust6766 <- abs(with(drops_df, diffvwc[portid == drop3 & cleanorder == 6766]))
# have to diff the last non-NA val for this next adjustment
adjust6837 <- abs(with(drops_df, vwc[portid == drop3 & cleanorder == 6835] - vwc[portid == drop3 & cleanorder == 6837]))
adjust6838 <- abs(with(drops_df, diffvwc[portid == drop3 & cleanorder == 6838]))
drops_df$adjust_vwc[drops_df$portid == drop3 & drops_df$cleanorder >= 6765] <- (with(drops_df, vwc[portid == drop3 & cleanorder >= 6765]) + adjust6765)
drops_df$adjust_vwc[drops_df$portid == drop3 & drops_df$cleanorder >= 6766] <- (with(drops_df, adjust_vwc[portid == drop3 & cleanorder >= 6766]) + adjust6766)
drops_df$adjust_vwc[drops_df$portid == drop3 & drops_df$cleanorder >= 6837] <- (with(drops_df, adjust_vwc[portid == drop3 & cleanorder >= 6837]) + adjust6837)
drops_df$adjust_vwc[drops_df$portid == drop3 & drops_df$cleanorder >= 6838] <- (with(drops_df, adjust_vwc[portid == drop3 & cleanorder >= 6838]) + adjust6838)

# adjust at the break (timestep at rank 4 and the step before or after (two-step drop))
# but like above so lowest val after break is at 0
break4 <- with(drops_df, cleanorder[portid == drop3 & rankdiff == 4])
adjust4 <- abs(with(drops_df, diffvwc[portid == drop3 & rankdiff == 4]))
#adjust4next <- abs(with(drops_df, diffvwc[portid == drop3 & cleanorder == (break4+1)]))
adjustdrop3 <- 0- with(drops_df, min(adjust_vwc[cleanorder > break4 & portid == drop3], na.rm = T))
drops_df$adjust_vwc[drops_df$portid == drop3 & drops_df$cleanorder == break4] <- (with(drops_df, adjust_vwc[portid == drop3 & cleanorder == break4]) + adjust4)
drops_df$adjust_vwc[drops_df$portid == drop3 & drops_df$cleanorder > break4] <- (with(drops_df, adjust_vwc[portid == drop3 & cleanorder > break4]) + adjustdrop3)
# ok. the end is weird, but not much to do. adjusted the diff from more neg point to 0 so position final vwc vals not negative
# if adjust by the actual drop, then ends almost as wet as wettest block 1 line


# subset out of master and add adjusted vals in with note on adjustment
smdat_raw <- smdat
smdat <- subset(smdat, !portid %in% drops$portid)
smdat$adjust_note <- NA
smdat$adjust_note <- as.character(smdat$adjust_note)

drops_df$adjustment <- drops_df$adjust_vwc - drops_df$vwc
drops_toadd <- arrange(drops_df, portid, cleanorder) %>%
  mutate(adjust_note = ifelse(!is.na(adjustment) & adjustment > 0, paste("vwc adjusted", round(adjustment,4), "for new break in data"), NA)) %>%
  dplyr::select(-vwc) %>%
  rename(vwc = adjust_vwc) %>%
  dplyr::select(names(smdat))
smdat <- rbind(smdat, drops_toadd)
smdat <- arrange(smdat, portid, cleanorder)

# all good?
ggplot(smdat, aes(clean_datetime, vwc)) + 
  geom_line(aes(group = portid, col = !is.na(adjust_note))) +
  facet_grid(ppt_trt~nut_trt)
# ok

# clean up -- keep drop data frames just in case need
rm(list = ls()[grepl("^adjust|^drop[0-9]|^break", ls())])



  # -- FLAGGING TEST CASE -----
# apply QA/QC methods from lit
# subset dat to a relatively clean/trustworthy line
# > who has the fewest qa notes?
group_by(smdat, portid) %>%
  mutate(totobs = length(vwc)) %>%
  subset(is.na(qa_note)) %>%
  group_by(portid, fulltrt, totobs) %>%
  summarise(prop.clean = (length(vwc) / totobs)) %>%
  distinct() %>%
  arrange(totobs, prop.clean) %>%
  data.frame()
# b2l1 = non-irrigated test case, b3l1 = irrigated fert and comp, and one fert drought
testcase <- subset(smdat, logger == "B2L1")
rownames(testcase) <- NULL  # clear rownames
testcase$rowid <- as.numeric(rownames(testcase))

testcase2 <- subset(smdat, logger %in% c("B2L1", "B3L1"))

# range check
check_range <- function(dat, id = c("portid", "clean_datetime", "cleanorder"), hi = 0.6, low = -0.04){
  tempdat <- subset(dat, select = c(id, "vwc"))
  tempdat$flag_range[tempdat$vwc > hi] <- "high warning"
  tempdat$flag_range[tempdat$vwc <= low] <- "outside low"
  tempdat$flag_range[tempdat$vwc < 0 & tempdat$vwc > low] <- "low warning"
  return(tempdat)
}

# apply
rangetest <- check_range(testcase)
# plot
ggplot(rangetest, aes(clean_datetime, vwc, group = portid)) +
  geom_line() +
  geom_point(data=subset(rangetest, !is.na(flag_range)), aes(clean_datetime, vwc, col = flag_range)) +
  facet_wrap(~portid)


# flatline of 10+ timesteps (don't think there are any but to be sure)
check_streak <- function(dat, id = c("clean_datetime", "cleanorder"), groupvars = c("portid"), steps = 10){
  tempdat <- subset(dat, select = c(id, groupvars, "vwc")) %>%
    grouped_df(groupvars) %>%
    mutate(diff = vwc == dplyr::lag(vwc))
  streakruns <- list2DF(rle(tempdat$diff))
  # calculate startpos for streakruns
  streakruns$endpos <- cumsum(streakruns$lengths)
  streakruns$startpos <- with(streakruns, endpos - lengths)
  # subset to streakruns that meet or exceed flatline threshold
  streakruns <- subset(streakruns, lengths >= steps & values)
  # create flag col
  tempdat$flag_streak <- NA
  tempdat$streak_num <- NA
  for(i in 1:nrow(streakruns)){
    tempdat$flag_streak[streakruns$startpos[i]:streakruns$endpos[i]] <- paste0("flatline (", streakruns$lengths[i]+1, " timesteps)")
    tempdat$streak_num[streakruns$startpos[i]:streakruns$endpos[i]] <- i
  }
  # clean up
  tempdat <- ungroup(subset(tempdat, select = -diff))
  return(tempdat)
}

# apply
streaktest <- check_streak(testcase)
# plot
ggplot(subset(streaktest, portid %in% portid[!is.na(flag_streak)] & cleanorder %in% cleanorder[!is.na(flag_streak)]),aes(clean_datetime, vwc, group = portid)) +
  geom_line() +
  geom_point(data=subset(streaktest, !is.na(flag_streak)), aes(clean_datetime, vwc, col = flag_streak)) +
  facet_wrap(~streak_num, scales = "free")


# target day exceeds +3 z-score in 30d window
check_deviance <- function(dat, id = c("clean_datetime", "cleanorder"), groupvars = c("portid", "fulltrt"), rollwin = (30*12)){
  tempdat <- subset(dat, select = c(id, groupvars, "vwc"))
  tempdat <- grouped_df(tempdat, groupvars)
  # calculate moving window mean, sd, and # non-NA observations
  tempdat$rollmean <- rollmean(tempdat$vwc, rollwin, fill = NA, align = "left", na.rm = T) # look from target forward in time rollwin steps
  tempdat$rollsd <- rollapply(tempdat$vwc, rollwin, align = "left", fill = NA, function(x) sd(x, na.rm = T))
  tempdat$rollnobs <- rollapply(tempdat$vwc, rollwin, align = "left", fill = NA, function(x) length(x[!is.na(x)]))
  tempdat <- ungroup(tempdat)
  # calculate z score (mean-centered vwc)/(sd)
  tempdat$rollz <- with(tempdat, (vwc-rollmean)/rollsd)
  return(tempdat)
}
test_deviance <- check_deviance(testcase)

ggplot(test_deviance, aes(clean_datetime, vwc)) +
  geom_line() +
  geom_point(data = subset(test_deviance, abs(rollz)>3), aes(size = abs(rollz)), col = "pink", alpha = 0.5) +
  facet_wrap(~portid)

# target day exceeds +3 z-score in same-day over all years step change
check_step <- function(dat, id = c("clean_datetime", "date", "cleanorder"), groupvars = c("portid", "fulltrt"), grouptrt = c("doy", "ppt_trt"), rollwin = (30)){
  tempdat <- subset(dat, select = c(id, groupvars, grouptrt, "vwc")) %>% #portid == "B2L1_1", 
    grouped_df(groupvars) %>%
    # calculate diff
    mutate(diff = dplyr::lead(vwc) - vwc) %>% # subtract target from next (forward) so is left-aligned
    #mutate(diff = vwc - dplyr::lag(vwc)) %>% # try other way
    ungroup %>%
    # calculate scaled val by day of year and precip treatment
    grouped_df(grouptrt) %>%
    mutate(meandiff = mean(diff, na.rm = T),
           sddiff = sd(diff, na.rm = T),
           nobsdiff = length(diff[!is.na(diff)])) %>%
    ungroup() %>%
    # calculate z score (mean centered step change)/(sd step change in window)
    mutate(zdiff = (diff - meandiff)/sddiff)
  return(tempdat)
}

step_test <- check_step(testcase)
ggplot(step_test, aes(clean_datetime, vwc)) +
  geom_line(aes(group = portid)) + 
  geom_point(data = subset(step_test, abs(zdiff) > 3), aes(size = abs(zdiff)), alpha = 0.5, col = "yellow") +
  facet_wrap(~portid)

# see what was flagged in both
combine_stepdev <- left_join(step_test, test_deviance) %>%
  mutate(flag_extreme = (abs(zdiff) > 5 & abs(rollz) > 5)) # any z greater than 3 in paper (but that flags too many points here)
# plot
ggplot(combine_stepdev, aes(clean_datetime, vwc)) +
  geom_line(aes(group = portid)) +
  geom_point(data = subset(combine_stepdev, flag_extreme),aes(clean_datetime, vwc), col = "orchid", size = 2) +
  facet_wrap(~portid)

# try it with two loggers
step_test2 <- check_step(testcase2)
dev_test2 <- check_deviance(testcase2)
stepdev2 <-  left_join(step_test2, dev_test2) %>%
  mutate(flag_extreme = (abs(zdiff) > 5 & abs(rollz) > 5))
# plot
ggplot(stepdev2, aes(clean_datetime, vwc)) +
  geom_line(aes(group = portid, col = substr(portid, 1,2)), alpha = 0.5) +
  geom_point(data = subset(stepdev2, flag_extreme),aes(clean_datetime, vwc), col = "orchid", size = 2) +
  facet_wrap(~portid)
#facet_wrap(~gsub("_[1-5]", "", portid))


# all 3 criteria must fail to flag
# update: this won't work when NAs present in time series .. need to infill through NA spline in zoo
check_spike <- function(dat, id = c("clean_datetime", "date", "cleanorder", "portid", "fulltrt", "ppt_trt"), groupvars = c("portid")){
  tempdat <- subset(dat, select = unique(c(id, groupvars, "vwc")))
  # > before grouping data to iterate by port, run zoo tasks (for step 2 check)
  
  # calculate Savitzky Golay 2nd order derivative (they used 3hrs for filter window.. but that seems small?)
  # need to make data zoo object and iterate through each portid for this to work..
  # > first initiate data frame for storing results
  zoodat <- data.frame()
  for(i in unique(tempdat$portid)){
    portdat <- ungroup(subset(tempdat, portid == i, select = unique(c(id, groupvars, "vwc"))))
    x <- zoo::zoo(portdat$vwc, portdat$clean_datetime)
    portdat$deriv2 <- pracma::savgol(zoo::na.spline(x), fl = 3, forder = 2, dorder = 2)
    zoodat <- rbind(zoodat, portdat)
  }
  # join 2nd deriv info back to tempdat
  tempdat <- left_join(tempdat, zoodat[c("clean_datetime", "cleanorder", "portid", "deriv2")])
  # proceed with checks (via grouping data)
  tempdat <- grouped_df(tempdat, groupvars)
  
  # 1. subsequent time step change of 15%
  # > ctw also adding abs change of at least 0.01m3^m-3 to avoid low moisture wobbles in summer/dry periods
  tempdat$stepchange <- abs(tempdat$vwc / dplyr::lag(tempdat$vwc))
  tempdat$abschange <- abs(tempdat$vwc - dplyr::lag(tempdat$vwc))
  # 2. second deriv outside 0.8 to 1.2 (based on calibrated data in Dorigo et al. 2013)
  # calculate rel change between lag and lead
  tempdat$step2derv <- abs(dplyr::lag(tempdat$deriv2) / dplyr::lead(tempdat$deriv2))
  # 3. 24hr CV exceeds 1
  # 12 hrs prior to target t, and 12 hours ahead, but *don't* include t
  tempdat$mu24 <- rollapply(tempdat$vwc, list(c(-6:-1, 1:6)), fill = NA, function(x) mean(x, na.rm = T))
  tempdat$sd24 <- rollapply(tempdat$vwc, list(c(-6:-1, 1:6)), fill = NA, function(x) sd(x, na.rm = T))
  tempdat$nobs24 <- rollapply(tempdat$vwc, list(c(-6:-1, 1:6)), fill = NA, function(x) length(x[!is.na(x)]))
  # calculate cv
  tempdat <- ungroup(tempdat)
  tempdat$cv24 <- with(tempdat, sd24/mu24)
  # flag if all three criteria violated
  tempdat$spike <- with(tempdat, (stepchange > 1.15 | stepchange < 0.85) & (step2derv < 1.2 & step2derv > 0.8) & cv24 < 1)
  # return simple tempdat
  return(tempdat[unique(c(id, groupvars, "vwc", "nobs24", "abschange", "spike"))])
}

test_spike <- check_spike(testcase2)
ggplot(test_spike, aes(clean_datetime, vwc))+
  geom_line(aes(group = portid)) +
  geom_point(data = subset(test_spike, spike & abschange > 0.015), col = "orchid", size = 2, alpha = 0.5)+
  facet_wrap(~portid)
# more helpful when add in abschange check
ggplot(subset(test_spike, portid == "B2L1_1" & cleanorder %in% 3000:4000), aes(cleanorder, vwc))+
  geom_line(aes(group = portid)) +
  geom_point(data = subset(test_spike, spike & portid == "B2L1_1" & cleanorder %in% 3000:4000), col = "orchid", alpha = 0.5) +
  geom_point(data = subset(test_spike, spike & portid == "B2L1_1" & cleanorder %in% 3000:4000 & abschange < 0.01), col = "yellow", alpha = 0.75)

# build the break check since easy enough..
check_break <- function(dat, id = c("clean_datetime", "date", "cleanorder", "portid", "fulltrt", "ppt_trt"), groupvars = c("portid")){
  tempdat <- subset(dat, select = unique(c(id, groupvars, "vwc")))
  # > before grouping data to iterate by port, run zoo tasks (for step 2 check)
  
  # calculate Savitzky Golay 1st and 2nd order derivatives (they used 3hrs for filter window.. but that seems small?)
  # need to make data zoo object and iterate through each portid for this to work..
  # > first initiate data frame for storing results
  zoodat <- data.frame()
  for(i in unique(tempdat$portid)){
    portdat <- ungroup(subset(tempdat, portid == i, select = unique(c(id, groupvars, "vwc"))))
    x <- zoo::zoo(portdat$vwc, portdat$clean_datetime)
    portdat$deriv1 <- pracma::savgol(zoo::na.spline(x), fl = 3, forder = 2, dorder = 1)
    portdat$deriv2 <- pracma::savgol(zoo::na.spline(x), fl = 3, forder = 2, dorder = 2)
    portdat$vwc_interp <- as.numeric(zoo::na.spline(x))
    zoodat <- rbind(zoodat, portdat)
  }
  # join 2nd deriv info back to tempdat
  tempdat <- left_join(tempdat, zoodat[c("clean_datetime", "cleanorder", "portid", "deriv1", "deriv2", "vwc_interp")])
  # proceed with checks (via grouping data)
  tempdat <- grouped_df(tempdat, groupvars)
  
  # 1. rel step change at least 10% and abs change at least 0.01m3^m-3
  tempdat$relchange <- abs((tempdat$vwc - dplyr::lag(tempdat$vwc)) / tempdat$vwc)
  tempdat$abschange <- abs(tempdat$vwc - dplyr::lag(tempdat$vwc))
  
  # 2. first deriv at target t not more than 10x the 24hr average 2nd deriv, centered at t
  # calculate 24 hr average of 2nd deriv
  tempdat$mu24_2derv <- rollmean(tempdat$deriv2, 12, align = "center", fill = NA)
  
  # 3. ratios of 2nd derivs (two to crunch)
  tempdat$rat1 <- with(tempdat, round(abs(deriv2 / dplyr::lead(deriv2)), 0))
  tempdat$rat2 <- with(tempdat, round(abs(dplyr::lead(deriv2) / dplyr::lead(deriv2,2)), 0))
  
  # ungroup
  tempdat <- ungroup(tempdat)
  # flag if all three criteria violated
  tempdat$cr1 <- with(tempdat, relchange > 0.1 & abschange > 0.01)
  tempdat$cr2 <- tempdat$deriv1 > (10*tempdat$deriv2)
  tempdat$cr3 <- with(tempdat, abs(rat1) == 1 & abs(rat2) > 10) # abs to be sure
  tempdat$flagbreak <- with(tempdat, cr1 & cr2 & cr3)
  # return simple tempdat
  return(tempdat[unique(c(id, groupvars, "vwc", "vwc_interp", "flagbreak"))])
}

test_break <- check_break(testcase2)
ggplot(test_break, aes(clean_datetime, vwc))+
  geom_line(aes(group = portid)) +
  geom_point(data = subset(test_break, flagbreak), col = "orchid", alpha = 0.5, size = 2)+
  facet_wrap(~portid)
# okay 
ggplot(subset(test_break, portid == "B2L1_1" & cleanorder %in% 0:1000), aes(cleanorder, vwc))+
  geom_line(aes(group = portid)) +
  geom_point(data = subset(test_break, flagbreak & portid == "B2L1_1" & cleanorder %in% 0:1000), col = "orchid", alpha = 0.75, size = 2)+
  facet_wrap(~portid)

# wetup event and no observed rainfall in the previous 24 hrs
# > either get irrigation schedule for W plots or only run this on XC and D ppt trts (W in summer okay)
# for this, need a prepped rainfall dataset, the soil moisture data, and will combine both then flag any wetup not preceded by a rainfall event in previous 24hrs

outside_precip <- function(soildat, raindat, keepsoil = c("clean_datetime", "cleanorder", "portid", "fulltrt", "ppt_trt", "nut_trt", "vwc"), groupvars = c("portid")){
  # subset soildat
  tempdat <- subset(soildat, select = keepsoil)
  
  # 1. track wetup in soilmoisture
  rundf <- data.frame() # initiate data frame
  # iterate through each portid in the soil moisture dataset
  for(p in unique(tempdat$portid)){
    portdat <- subset(tempdat, portid == p) %>%
      mutate(diffvwc = round(vwc - dplyr::lag(vwc),2),
             # increase of at least 15%, or anticipates a wetup of at least 15%, or is still wetting up in timestep after 0.15 increase
             wetup = diffvwc >= 0.15 | (diffvwc > 0.01 & dplyr::lead(diffvwc) >= 0.15) | (diffvwc >= 0.01 & dplyr::lag(diffvwc) > 0.15),
             # keep track of events per portid if ever want to plot
             wetup_event = as.numeric(NA))
    temprun <- list2DF(rle(portdat$wetup))
    
    # annotate events
    r <- 1
    wetupevent <- 0 # ignore the "first" wetup (installation)
    for(i in 1:nrow(temprun)){
      # note where to stop
      tempend <- r+ (temprun$lengths[i] -1)
      # note event
      # > if NA or FALSE ignore
      if((is.na(temprun$values[i]) | !temprun$values[i])){
        # increment r and dryspell
        r <- tempend+1
      }else{
        portdat$wetup_event[r:tempend] <- wetupevent
        # increment r and rainevent
        r <- tempend+1
        wetupevent <- wetupevent + 1
      }
    }
    rundf <- rbind(rundf, portdat)
  }
  
  # 2. join raindat
  rundf <- left_join(rundf, raindat)
  
  # 3. check for wetup events outside of precip events
  rundf <- grouped_df(rundf, groupvars)
  # screen previous 24 hours for *any* precip event
  rundf$rain24 <- rollapply(rundf$rainevent, list(c(-12:0)), fill = NA, function(x) any(!is.na(x)))
  rundf$flag_wetup <- rundf$wetup & !rundf$rain24
  
  # return
  return(ungroup(rundf))
}

test_precip <- outside_precip(smdat, cimis_ppt2hr)
# add seasoninfo
test_precip2 <- left_join(test_precip, seasondat, by = c("waterYear" = "eco_yr")) %>%
  mutate(season_start = as.POSIXct(season_start, tz = "UTC", season_end = as.POSIXct(season_end, tz = "UTC"))) %>%
  mutate(irrigated = ifelse(ppt_trt == "W" & flag_wetup & clean_datetime >= season_start & date <= stopwater, 1, NA))
ggplot(test_precip, aes(clean_datetime, vwc)) +
  geom_line(aes(group = portid)) +
  geom_point(data = subset(test_precip, flag_wetup), col = "red", size = 2, alpha = 0.75) +
  geom_point(data = subset(test_precip2, irrigated == 1), col = "yellow", size = 1, alpha = 0.75) +
  facet_wrap(nut_trt~ppt_trt)


# clean up environment before proceeding
rm(testcase, testcase2, streaktest, combine_stepdev, dev_test2, rangetest, step_test, 
   step_test2, stepdev2, test_break, test_deviance, test_precip, test_precip2, test_spike)



# -- QA/QC ALL DATA -----
# logic checks:
# plausible range check
smdat_qa <- check_range(smdat)
smdat_qa <- cbind(smdat, flag_range = smdat_qa$flag_range)
# retain original vwc as raw_vwc
smdat_qa$raw_vwc <- smdat_qa$vwc ## NOTE: will need to add raw unadjusted vwc for three portids with adjusted data
# NA vwc outside range check
smdat_qa$vwc[grepl("outside|high", smdat_qa$flag_range)] <- NA

# outside precipitation and irrigation events
precip_qa <- outside_precip(smdat_qa, cimis_ppt2hr)
# join irrigation info
precip_qa <- left_join(precip_qa, seasondat, by = c("waterYear" = "eco_yr")) %>%
  mutate(season_start = as.POSIXct(season_start, tz = "UTC"), season_end = as.POSIXct(season_end, tz = "UTC")) %>%
  mutate(irrigated = (ppt_trt == "W" & clean_datetime >= season_start & date <= stopwater))


# behavior checks:
# streak (flatlines) check
streak_qa <- check_streak(smdat_qa)

# deviance (extreme values) check
deviance_qa <- check_deviance(smdat_qa)
# step check
step_qa <- check_step(smdat_qa)
# pair deviance and step
devstep_qa <-  left_join(deviance_qa, step_qa) %>%
  mutate(flag_extreme = (abs(zdiff) > 5 & abs(rollz) > 5)) # go with +5 sd outside

# spike check
spike_qa <- check_spike(smdat_qa)

# data break check
break_qa <- check_break(smdat_qa)


# combine all flags
smdat_qa_all <- smdat_qa %>%
  left_join(subset(precip_qa, select = c(cleanorder:portid, flag_wetup, irrigated))) %>%
  left_join(subset(streak_qa, select = c(cleanorder:portid, flag_streak, streak_num))) %>%
  left_join(subset(devstep_qa, select = c(cleanorder:portid, flag_extreme))) %>%
  left_join(subset(spike_qa, select = c(cleanorder:portid, spike))) %>%
  left_join(subset(break_qa, select = c(cleanorder:portid, flagbreak)))

flagged_data <- subset(smdat_qa_all, ((flag_wetup & !irrigated) | !is.na(flag_streak) | flag_extreme) | spike | flagbreak) %>%
  mutate(flag_sensor = (flag_extreme & spike) | (flag_extreme & flagbreak) | (spike & flagbreak))



# -- REVIEW FLAGS -----
# start with 1 logger
unique(smdat_qa_all$logger)
testlog <- unique(smdat_qa_all$logger)
ggplot(subset(smdat_qa_all, logger %in% testlog), aes(clean_datetime, vwc, group = portid)) +
  geom_line(alpha =0.5) +
  # what happens if remove anything flagged?
  geom_point(data = subset(flagged_data, logger %in% testlog & flag_wetup), col = "blue", alpha = 0.5) +
  geom_point(data = subset(flagged_data, logger %in% testlog & !is.na(flag_streak)), col = "orchid", alpha = 0.5) +
  geom_point(data = subset(flagged_data, logger %in% testlog & flag_extreme), col = "orange", alpha = 0.5) +
  geom_point(data = subset(flagged_data, logger %in% testlog & spike), col = "green", alpha = 0.5) +
  geom_point(data = subset(flagged_data, logger %in% testlog & flagbreak), col = "red", alpha = 0.5) +
  facet_wrap(~fulltrt) #+portid

# look at outside precip flags and values that have multiple behavior flags only
ggplot(subset(smdat_qa_all, portid %in% with(flagged_data, portid[flag_wetup | flag_sensor])), aes(clean_datetime, vwc, group = portid)) +
  geom_line(alpha =0.5) +
  # what happens if remove anything outside precip event and if had more than 1 unsual behavior flag
  geom_point(data = subset(flagged_data, flag_wetup), col = "blue", size = 2, alpha = 0.5) +
  geom_point(data = subset(flagged_data, flag_sensor), col = "orange", size = 2, alpha = 0.5) +
  facet_wrap(~fulltrt)

# what would I be missing if only removed those?
ggplot(subset(smdat_qa_all, !portid %in% with(flagged_data, portid[flag_wetup | flag_sensor])), aes(clean_datetime, vwc, group = portid)) +
  geom_line(alpha =0.5) +
  # what happens if remove anything outside precip event and if had more than 1 unsual behavior flag
  #geom_point(data = subset(flagged_data, flag_wetup), col = "blue", size = 2, alpha = 0.5) +
  #geom_point(data = subset(flagged_data, flag_sensor), col = "orange", size = 2, alpha = 0.5) +
  facet_wrap(nut_trt~ppt_trt)
# not perfect but at least looks a little better and treated systematically

# plot streaks
ggplot(subset(streak_qa, !is.na(flag_streak)), aes(clean_datetime, vwc, col = substr(portid, 6,6))) +
  geom_point(alpha = 0.6) +
  facet_wrap(~gsub("_[0-9]", "", portid))



# -- TREATMENT CONGRUENCY (optional, not coded) ----
# some particular treatments stood out above for portids that don't follow most profiles
# first review identity to be sure treatments are labeled correctly (i've reviewed before, but to be quadrupley sure).. and to see where they are in space
sustrts <- c("ND", "FD", "CW", "FW")
colorder <- data.frame(plot = 1:36, colorder = rep(1:3, 12))
trtrev <- distinct(subset(smdat_qa_all, select = c(logger:comp_trt)))

# first treatment group (ND)
ggplot(subset(smdat_qa_all, fulltrt %in% sustrts[1]), aes(clean_datetime, vwc, group = portid)) +
  geom_line(aes(col = portid, lty = comp_trt)) +
  facet_wrap(~fulltrt+comp_trt, nrow= 2)

# could B1L3 have gotten downhill runoff?? (altho control doesn't respond as seeded does in same subplot)
ggplot(subset(smdat_qa_all, block == 1 & nut_trt == "N"), aes(clean_datetime, vwc, group = portid)) +
  geom_line(aes(col = portid, lty = comp_trt), alpha = 0.5) +
  facet_grid(paste(plot, ppt_trt)~nut_trt) # kind of looks like it.. altho both are in seeded and responding differently. i don't know. will document it for AS.

# second treatment group (FD)
ggplot(subset(smdat_qa_all, fulltrt %in% sustrts[2]), aes(clean_datetime, vwc, group = portid)) +
  geom_line(aes(col = portid, lty = comp_trt), alpha = 0.5) +
  facet_grid(comp_trt~block)

ggplot(subset(smdat_qa_all, fulltrt %in% sustrts[2] & block == 3), aes(clean_datetime, vwc, group = portid)) +
  geom_line(aes(col = portid, lty = comp_trt), alpha = 0.5) +
  facet_wrap(~plot) # there are 5 different control lines in block 3 FD. B3L1_3 is different than the rest all three years.. the other four are similar

# plot block
subset(smdat_qa_all, block == 3) %>%
  left_join(colorder) %>%
ggplot(aes(clean_datetime, vwc, group = portid)) +
  geom_line(aes(col = ppt_trt, lty = comp_trt), alpha = 0.5) +
  facet_grid(colorder~nut_trt)
# i feel like a line each in F W + D were switched when writing down what was what. they each consistently match the other ppt_trt visually
# including wetup in W plots when XC does not wetup.. i'm going to switch them and then confirm w AS.

# see who it is:
subset(smdat_qa_all, block == 3 & nut_trt == "F") %>%
  left_join(colorder) %>%
  ggplot(aes(clean_datetime, vwc, group = portid)) +
  geom_line(aes(col = paste(plot, port), lty = ppt_trt), alpha = 0.8) +
  facet_grid(colorder~nut_trt) # b3l1_3 (plot 22 port 3, currently FD) and b3l1_4 (plot 23 port 4, currently FW)

# third treatment group (CW)
subset(smdat_qa_all, fulltrt %in% sustrts[3]) %>%
  left_join(colorder) %>%
  ggplot(aes(clean_datetime, vwc, group = portid)) +
  geom_line(aes(col = portid, lty = comp_trt)) +
  facet_grid(paste(colorder, nut_trt)~block) # issue only in block 1, plot 3 (compost row)

subset(smdat_qa_all, block == 1) %>%
  left_join(colorder) %>%
  ggplot(aes(clean_datetime, vwc, group = portid)) +
  geom_line(aes(col = paste(plot, port), lty = ppt_trt), alpha = 0.6) +
  # array in spatial order
  facet_grid(colorder~factor(nut_trt, levels = c("C", "N", "F")))
# i think a line is missing for the XC ppt plot in compost.. there should be one and the one in wet that's different looks more like drought
# also what's coded as drought in compost looks more like XC..
# looking at Word doc key, it says block 2 logger 5 ports 4 and 5 go to compost control in block 1 (but mb wet composition-control bc other b2l5 are block 2 compost wet, which is nextdoor)

subset(smdat_qa_all, block == 1 & nut_trt == "C") %>%
  left_join(colorder) %>%
  left_join(seasondat, by = c("waterYear" = "eco_yr")) %>%
  mutate(season_start = as.POSIXct(season_start, tz = "UTC"), season_end = as.POSIXct(season_end, tx = "UTC")) %>%
  ggplot(aes(clean_datetime, vwc, group = portid)) +
  # plot ribbons for rain seasondat
  geom_ribbon(aes(xmin = season_start, xmax = season_end, group = waterYear), fill = "grey50", alpha = 0.25) +
  geom_line(aes(col = portid, lty = ppt_trt), alpha = 0.6) +
  # array in spatial order
  facet_grid(paste(plot, ppt_trt)~comp_trt) # 7 portids in one subplot seems like a lot..
# what if b1l1_1 and _2 go to XC, B1L1_4 and _5 go to D and B2L5_4 and _5 go to W?

# look at all again to confirm
subset(smdat_qa_all, block == 1) %>%
  left_join(colorder) %>%
  left_join(seasondat, by = c("waterYear" = "eco_yr")) %>%
  mutate(season_start = as.POSIXct(season_start, tz = "UTC"), season_end = as.POSIXct(season_end, tx = "UTC")) %>%
  ggplot(aes(clean_datetime, vwc, group = portid)) +
  # plot ribbons for rain seasondat
  geom_ribbon(aes(xmin = season_start, xmax = season_end, group = waterYear), fill = "grey50", alpha = 0.25) +
  geom_line(aes(col = portid, lty = ppt_trt), alpha = 0.6) +
  # array in spatial order
  facet_wrap(~nut_trt+ppt_trt +plot) 
# seems appropriate correction

# last treatment to review (also might as well review block 2 since i've looked at 1 and 3 now)
ggplot(subset(smdat_qa_all, fulltrt %in% sustrts[4]), aes(clean_datetime, vwc, group = portid)) +
  geom_line(aes(col = portid, lty = comp_trt)) +
  facet_wrap(~fulltrt+comp_trt, nrow= 2)

subset(smdat_qa_all, fulltrt %in% sustrts[4]) %>%
  left_join(colorder) %>%
  ggplot(aes(clean_datetime, vwc, group = portid)) +
  geom_line(aes(col = portid, lty = comp_trt)) +
  facet_grid(paste(colorder, nut_trt)~block) # might be one i've already figured out since in block 3.. b3l1_4. yes. should be switched with with b3l1_3

# review block 2
subset(smdat_qa_all, block == 2) %>%
  left_join(colorder) %>%
  left_join(seasondat, by = c("waterYear" = "eco_yr")) %>%
  mutate(season_start = as.POSIXct(season_start, tz = "UTC"), season_end = as.POSIXct(season_end, tx = "UTC")) %>%
  ggplot(aes(clean_datetime, vwc, group = portid)) +
  # plot ribbons for rain seasondat
  geom_ribbon(aes(xmin = season_start, xmax = season_end, group = waterYear), fill = "grey50", alpha = 0.25) +
  geom_line(aes(col = ppt_trt, lty = comp_trt), alpha = 0.6) +
  # array in spatial order
  facet_grid(colorder~factor(nut_trt, levels = c("F", "N", "C"))) 
# looks okay.. maybe some drought plots just got wet in spots in 2020
# drought shelters here don't seem quite as effective in year 1

# look at drought plots in block 2
subset(smdat_qa_all, block == 2 & ppt_trt == "D") %>%
  left_join(colorder) %>%
  left_join(seasondat, by = c("waterYear" = "eco_yr")) %>%
  mutate(season_start = as.POSIXct(season_start, tz = "UTC"), season_end = as.POSIXct(season_end, tx = "UTC")) %>%
  ggplot(aes(clean_datetime, vwc, group = portid)) +
  # plot ribbons for rain seasondat
  geom_ribbon(aes(xmin = season_start, xmax = season_end, group = waterYear), fill = "grey50", alpha = 0.25) +
  geom_line(aes(col = portid, lty = comp_trt), alpha = 0.6) +
  # array in spatial order
  facet_grid(paste(plot,nut_trt)~ppt_trt) # shelter doesn't seem effective here? or in year 1

# look at single nut trt
subset(smdat_qa_all, block == 2 & nut_trt == "F") %>%
  left_join(colorder) %>%
  left_join(seasondat, by = c("waterYear" = "eco_yr")) %>%
  mutate(season_start = as.POSIXct(season_start, tz = "UTC"), season_end = as.POSIXct(season_end, tx = "UTC")) %>%
  ggplot(aes(clean_datetime, vwc, group = portid)) +
  # plot ribbons for rain seasondat
  geom_ribbon(aes(xmin = season_start, xmax = season_end, group = waterYear), fill = "grey50", alpha = 0.25) +
  geom_line(aes(col = portid, lty = comp_trt), alpha = 0.6) +
  # array in spatial order
  facet_grid(paste(plot,ppt_trt)~nut_trt) 
# i guess fine, just not as effective in this block



# -- 1. correct mislabeled treatments -----
trts2correct_raw <- subset(trtrev, portid %in% c("B3L1_3", "B3L1_4", "B1L1_1", "B1L1_2", "B1L1_4", "B1L1_5"), select = c(logger:plot, fulltrt:comp_trt))
trts2correct_new <- trts2correct_raw

# 1. address ND and FD: switch plot 22 port 3 (current block 3 FD) and plot 23 port 4 (current block 3 FW) treatments (nut_trt, ppt_trt, and comp_trt)
trts2correct_new[which(trts2correct_raw$plot == 22 & trts2correct_raw$port == 3),c("logger", "port", "portid")] <- trts2correct_raw[which(trts2correct_raw$plot == 23 & trts2correct_raw$port == 4),c("logger", "port", "portid")] 
trts2correct_new[which(trts2correct_raw$plot == 23 & trts2correct_raw$port == 4),c("logger", "port", "portid")] <- trts2correct_raw[which(trts2correct_raw$plot == 22 & trts2correct_raw$port == 3),c("logger", "port", "portid")] 

# 2. address block 1 Compost:
# B1l1_1 and _2 (current drought) becomes Compost Control; B1L1_4 and _5 should be Compost Drought
# 2 lines
trts2correct_new[with(trts2correct_raw, which(nut_trt == "C" & ppt_trt == "D")), "plot"] <-7
trts2correct_new[with(trts2correct_raw, which(nut_trt == "C" & ppt_trt == "D")), "fulltrt"] <- "CXC"
trts2correct_new[with(trts2correct_raw, which(nut_trt == "C" & ppt_trt == "D")), "ppt_trt"] <- "XC"
trts2correct_new[with(trts2correct_raw, which(nut_trt == "C" & ppt_trt == "W")), "plot"] <- 8
trts2correct_new[with(trts2correct_raw, which(nut_trt == "C" & ppt_trt == "W")), "fulltrt"] <- "CD"
trts2correct_new[with(trts2correct_raw, which(nut_trt == "C" & ppt_trt == "W")), "ppt_trt"] <- "D"
trts2correct_new$plotid <- with(trts2correct_new, paste0(plot, fulltrt))
trts2correct_new$corrected <- TRUE

# look okay?
trts2correct_raw
trts2correct_new # yes

# review correction applied and for seeding treatment (do those seem fine)
updat <- subset(smdat_qa_all, portid %in% trts2correct_new$portid, select = -c(plot:comp_trt)) %>%
  left_join(trts2correct_new) %>%
  dplyr::select(c(names(smdat_qa_all), "corrected")) %>%
  rbind(cbind(subset(smdat_qa_all, !portid%in% trts2correct_new$portid), corrected = FALSE)) %>%
  left_join(colorder)
# go by block or treatment
ggplot(subset(updat, block == 3), aes(clean_datetime, vwc, group = portid))+
  geom_line(aes(lty = comp_trt, col = corrected), alpha = 0.5) +
  facet_grid(paste("block", block, "row", colorder) ~ nut_trt)
ggplot(subset(updat, block == 3), aes(clean_datetime, vwc, group = portid))+
  geom_line(aes(lty = ppt_trt, col = comp_trt), alpha = 0.75) +
  facet_grid(paste("block", block, "row", colorder) ~ nut_trt)
# seems like generally seeded is drier than control (which is what i'd expect)
# check with ashley if most plots should have at least one seeded and one control comp? corrected portids could easily be seeded where they got moved (both relatively drier)


# block 1
ggplot(subset(updat, block == 1), aes(clean_datetime, vwc, group = portid))+
  geom_line(aes(lty = comp_trt, col = corrected), alpha = 0.5) +
  facet_grid(colorder ~ nut_trt) # looks better
ggplot(subset(updat, block == 1), aes(clean_datetime, vwc, group = portid))+
  geom_line(aes(lty = ppt_trt, col = comp_trt), alpha = 0.75) +
  facet_grid(paste("block", block, "row", colorder) ~ nut_trt)

# how often does block 2 have only one type of seed treatment?
ggplot(subset(updat, block == 2), aes(clean_datetime, vwc, group = portid))+
  geom_line(aes(lty = ppt_trt, col = comp_trt), alpha = 0.75) +
  facet_grid(paste("block", block, "row", colorder) ~ nut_trt)
# ? seems like a mistake too if there are only sensors assigned to one comp treatment?
# re-checked word doc treatment key and leaving as is unless says otherwise  


# 2. re-run precip qa -----
# proceed with updated data, drop precip flag (will need to re-run with updated treatment info)
smdat_qa_all_goodtrts <- arrange(updat, logger, portid, clean_datetime) %>%
  subset(select = -c(flag_wetup, irrigated))

# outside precipitation and irrigation events
precip_qa <- outside_precip(smdat_qa_all_goodtrts, cimis_ppt2hr)
# join irrigation info
precip_qa <- left_join(precip_qa, seasondat, by = c("waterYear" = "eco_yr")) %>%
  mutate(season_start = as.POSIXct(season_start, tz = "UTC"), season_end = as.POSIXct(season_end, tz = "UTC")) %>%
  mutate(irrigated = (ppt_trt == "W" & clean_datetime >= season_start & date <= stopwater))

smdat_qa_all_goodtrts <- left_join(smdat_qa_all_goodtrts, subset(precip_qa, select = c(cleanorder:portid, flag_wetup, irrigated)))
# cleanup
rm(updat, sustrts)



# 3. deviation from same treatment sensors -----
# for each treatment, compare to homogenized (mean) others
# check daily deviance and daily delta deviance
# anything wildly outside of average relationship gets flagged (4 sd seems to catch values that probably are questionable)
congruency_check <- data.frame()
congruency_sd_threshold <- 5 # can mod this to play around with flagging outcomes
congruency_raw_threshold <- 0.2
# unique(smdat_qa_all_goodtrts$portid)
for(u in unique(smdat_qa_all_goodtrts$portid)){
  # look at it by full trt and just ppt trt, just in case there aren't enough fulltrt lines for last year because of missing data
  # pull ppt trt of interest
  temp_ppttrt <- unique(smdat_qa_all_goodtrts$ppt_trt[smdat_qa_all_goodtrts$portid == u])
  tempfulltrt <- unique(smdat_qa_all_goodtrts$fulltrt[smdat_qa_all_goodtrts$portid == u])
  
  # split the data
  refdata <- subset(smdat_qa_all_goodtrts, portid != u & ppt_trt == temp_ppttrt)
  # id rows that are flagged
  refflag <- with(refdata, which(!is.na(flag_streak) | flag_extreme | spike | flagbreak | (flag_wetup & !irrigated)))
  # NA anything that's been flagged to be conservative
  refdata$ref_vwc <- refdata$vwc
  refdata$ref_vwc[refflag] <- NA
  # since each line still has error, perhaps use average relationship of the target to each line, if deviates from both then flag
  refdata <- arrange(refdata, portid, cleanorder) %>% group_by(portid) %>%
    mutate(ref_lag = dplyr::lag(ref_vwc)) %>%
    ungroup()
  refdata$ref_abschange <- with(refdata, round(ref_vwc - ref_lag,6))
  refdata$ref_relchange <- with(refdata, round((ref_vwc - ref_lag)/ref_lag, 4))
  
  targetdat <- subset(smdat_qa_all_goodtrts, portid == u)
  # NA anything that violates precip flag
  targetdat$target_vwc <- targetdat$vwc
  targetdat$target_vwc[which(targetdat$flag_wetup & ! targetdat$irrigated)] <- NA
  targetdat$target_lag <- dplyr::lag(targetdat$target_vwc)
  targetdat$target_abschange <- with(targetdat, round(target_vwc - target_lag,6))
  targetdat$target_relchange <- with(targetdat, round((target_vwc - target_lag)/target_lag,4)) 
  
  compdata <- subset(refdata, select = c(portid, fulltrt:clean_datetime, ref_vwc:ref_relchange)) %>%
    left_join(subset(targetdat, select = c(clean_datetime, cleanorder, target_vwc:target_relchange))) %>%
    # check deviance from vwc val, abschange, and relchange by portid
    mutate(diff_vwc = target_vwc - ref_vwc,
           diff_abschange = target_abschange - ref_abschange,
           diff_relchange = target_relchange - ref_relchange,
           comp_vwc = as.numeric(scale(diff_vwc)),
           # to screen for raw absolute difference
           comp_vwc_rawdiff = diff_vwc,
           comp_abschange = as.numeric(scale(diff_abschange)),
           comp_relchange = as.numeric(scale(diff_relchange))) %>% # as.numeric to remove scale attributes
    # count number of full trt portids and ppt trt portids available each day for flagging threshold %>%
    group_by(cleanorder) %>%
    mutate(ports_ppt = length(unique(portid[!is.na(ref_vwc)])),
           ports_fulltrt = length(unique(portid[fulltrt == tempfulltrt & !is.na(ref_vwc)]))) %>%
    ungroup()
  
  
  # extract flags and return
  compflags <- subset(compdata, select = c(portid:clean_datetime, comp_vwc:ports_fulltrt)) %>%
    gather(check, val, comp_vwc:comp_relchange) %>%
    # keep only values that exceed threshold
    subset(abs(val) > congruency_sd_threshold| (check == "comp_vwc_rawdiff" & abs(val) >= congruency_raw_threshold)) %>%
    mutate(sametrt = fulltrt == tempfulltrt) %>%
    # tally how deviations per step
    group_by(check, cleanorder) %>%
    mutate(nports = length(unique(portid)),
           nports_fulltrt = length(unique(portid[sametrt]))) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(flag_fulltrt = (nports_fulltrt >= (ports_fulltrt -1)) & (ports_fulltrt >= 3), # there have to be at least 2 full trtment ports in the comparison
           flag_ppttrt = nports >= (ports_ppt *0.5),
           flag_both = flag_fulltrt & flag_ppttrt) %>%
    ungroup() %>%
    subset(flag_fulltrt | flag_ppttrt | flag_both, select = c(cleanorder, clean_datetime, check, flag_fulltrt:flag_both)) %>%
    mutate(portid = u) %>%
    left_join(targetdat[c("cleanorder", "target_vwc")]) %>%
    distinct()
  
  # figure out if the absolute diff is negative or positive for flagging purposes
  compdirection <- subset(compdata, cleanorder %in% compflags$cleanorder, select = c(portid, cleanorder, comp_vwc_rawdiff)) %>%
    group_by(cleanorder) %>%
    summarise(mean_vwc_rawdiff = mean(comp_vwc_rawdiff, na.rm = T)) %>%
    ungroup() #%>%
    #mutate(check = "comp_vwc_abs")
  
  # join to compflags
  compflags <- left_join(compflags, compdirection) %>%
    #NA any mean_vwcdiff that's not for the raw diff check
    mutate(mean_vwc_rawdiff = ifelse(check == "comp_vwc_rawdiff", mean_vwc_rawdiff, NA))
  
  # append to master
  congruency_check <- rbind(congruency_check, compflags)
  
  
}

# testcase plot (run on one target logger only)
compflags2 <- left_join(compflags, distinct(smdat_qa_all_goodtrts[c("portid", "fulltrt", "ppt_trt", "nut_trt")])) 
# will only plot the last ppt trt run in the loop above
ggplot(compdata) +
  geom_line(aes(cleanorder, ref_vwc, group = portid), alpha = 0.5, col = "grey80") +
  geom_line(aes(cleanorder, target_vwc), alpha = 0.5, col = "blue") +
  geom_point(data = compflags, aes(cleanorder, target_vwc, col = check), alpha = 0.5) +
  #geom_point(data = subset(compflags, flag_ppttrt), aes(cleanorder, target_vwc), col = "pink", alpha = 0.5) +
  #geom_point(data = subset(compflags, flag_fulltrt), aes(cleanorder, target_vwc), col = "green", alpha = 0.5) +
  geom_point(data = subset(compflags2, flag_both), aes(cleanorder, target_vwc, fill = check), pch= 21, alpha = 0.5) +
  facet_wrap(~fulltrt, nrow= 3)
  
# review each flag type case (3: comparative vwc deviance, comparative absolute change deviance, comparative relative change deviance)
unique(congruency_check$check)
comp2check <- "comp_relchange" # "comp_vwc_rawdiff" "comp_abschange"   "comp_relchange"   "comp_vwc"
test <- left_join(smdat_qa_all_goodtrts, subset(congruency_check, check == comp2check))
ggplot(test, aes(clean_datetime, vwc)) +
  geom_line(aes(group = portid, col = factor(block)), alpha = 0.5) +
  #geom_point(data = subset(test, check == comp2check & flag_fulltrt), col = "pink", alpha = 0.7) +
  #geom_point(data = subset(test, check == comp2check & flag_ppttrt), col = "orange", alpha = 0.7) +
  geom_point(data = subset(test, check == comp2check & flag_both), pch = 21, alpha = 0.7) + #aes(fill = mean_vwc_rawdiff <0), 
  scale_color_viridis_d() +
  facet_grid(nut_trt ~ ppt_trt)
 
# both for vwc deviance seem like it correctly IDs odd values (maybe also comparing against general ppt trt)
# abschange seems useful, but not conservative enough? (too many values flagged)
# whatever in relchange that should be flagged for being an odd value is probably getting flagged for abschange or deviance too
# count number of flags per obs for *both* flag only

congruency_check2 <- group_by(congruency_check, portid, cleanorder) %>%
  mutate(nflags_both = length(unique(check[flag_both]))) %>%
  ungroup() %>%
  # add treatment info for facetting
  left_join(distinct(smdat_qa_all_goodtrts[c("portid", "fulltrt", "block", "ppt_trt", "nut_trt")]))

# who are the flagged values?
ggplot(test, aes(clean_datetime, vwc)) +
  geom_line(aes(group = portid, col = factor(block)), alpha = 0.5) +
  #geom_point(data = subset(test, check == comp2check & flag_fulltrt), col = "pink", alpha = 0.7) +
  #geom_point(data = subset(test, check == comp2check & flag_ppttrt), col = "orange", alpha = 0.7) +
  #geom_point(data = subset(congruency_check2, flag_both & nflags_both == 1), aes(clean_datetime, target_vwc), col = "red", alpha = 0.7) +
  geom_point(data = subset(congruency_check2, flag_both & nflags_both > 1), aes(clean_datetime, target_vwc, fill = factor(nflags_both)), pch = 21, size = 2, alpha = 0.7) +
  scale_color_viridis_d() +
  facet_grid(ppt_trt ~ nut_trt)

#closer look at one treatment
ggplot(subset(test, fulltrt == "ND"), aes(clean_datetime, vwc)) +
  geom_line(aes(group = portid), alpha = 0.5) +
  #geom_point(data = subset(congruency_check2, fulltrt == "CXC" & flag_both), aes(clean_datetime, target_vwc, fill = check), pch = 21, size = 2, alpha = 0.7) +
  geom_point(data = subset(congruency_check2, fulltrt == "ND" & check == "comp_vwc_rawdiff" & flag_both), aes(clean_datetime, target_vwc, fill = mean_vwc_rawdiff), pch = 21, size = 2, alpha = 0.7) +
  #scale_color_viridis_d() +
  facet_wrap(~ portid)

# zoom in
ggplot(subset(test, portid == "B2L3_1" & cleanorder > 10000), aes(cleanorder, vwc)) +
  geom_line(aes(group = portid), alpha = 0.5) +
  #geom_point(data = subset(congruency_check2, fulltrt == "CXC" & flag_both), aes(clean_datetime, target_vwc, fill = check), pch = 21, size = 2, alpha = 0.7) +
  geom_point(data = subset(congruency_check2, portid == "B2L3_1" & check == "comp_vwc_rawdiff" & flag_ppttrt & cleanorder > 10000), aes(cleanorder, target_vwc, fill = mean_vwc_rawdiff), pch = 21, alpha = 0.7) +
  #scale_color_viridis_d() +
  facet_wrap(~ portid)

# how many flags per logger? [how many data points potentially removed]
sort(sapply(split(congruency_check[c("flag_fulltrt", "flag_ppttrt", "flag_both")], congruency_check$portid), function(x) sum(x == TRUE)))

# think it seems reasonable to move forward with using flag_both for comp_vwc
# maybe combine other congruency flags with other flags (e.g., if three flagged, review)
# keep congruency check df for clean up at the end (ctw added comp_vwc_rawdiff after going through flagging to try to catch spikes during dry periods)
# make wideformat to join with smdat_qa_all_goodtrts
comparative_flags <- subset(congruency_check, flag_both, select = -c(flag_ppttrt, flag_fulltrt, mean_vwc_rawdiff)) %>% #& !grepl("relchange|rawdiff", check)
  distinct() %>%
  spread(check, flag_both) %>%
  rename(flag_abschange_congruency = comp_abschange, 
         flag_vwc_congruency = comp_vwc,
         flag_relchange_congruency = comp_relchange,
         flag_rawdiff_congruency = comp_vwc_rawdiff)

# add to flagged and treatment-corrected dataset
smdat_qa_all_goodtrts <- left_join(smdat_qa_all_goodtrts, comparative_flags)

# clean up environment -- keep congruency check for later
rm(test, comp2check, congruency_check2, temp_ppttrt, tempfulltrt,compdata, compdirection, 
   refflag, targetdat, refdata, u, testlog, compflags, compflags2)

# remake flagged dataset with congruency flags and new outside precip flags
flagged_data <- subset(smdat_qa_all_goodtrts, !is.na(flag_range) | (flag_wetup & !irrigated) | !is.na(flag_streak) | flag_extreme | spike | flagbreak | flag_abschange_congruency | flag_vwc_congruency | flag_relchange_congruency | flag_abschange_congruency) %>%
  mutate(flag_sensor = (flag_extreme & spike) | (flag_extreme & flagbreak) | (spike & flagbreak))
# be sure at least one flag col as a true or not NA value (sometimes NAs get pulled in T/F subsetting)
flagged_data$nflags <- apply(flagged_data[names(flagged_data)[grepl("flag|spike", names(flagged_data))]], 1, function(x) sum(x != FALSE & !is.na(x)))
summary(factor(flagged_data$nflags)) # good, all have at least one. most have just one flag

# review
ggplot(smdat_qa_all_goodtrts, aes(clean_datetime, vwc)) +
  geom_line(aes(group = portid), alpha = 0.5) +
  geom_point(data = subset(flagged_data, nflags > 2), aes(clean_datetime, vwc, col = factor(nflags))) +
  facet_grid(ppt_trt ~ nut_trt)

# which points have behavior flag + congruency flag?
ggplot(smdat_qa_all_goodtrts, aes(clean_datetime, vwc)) +
  geom_line(aes(group = portid), alpha = 0.5) +
  geom_point(data = subset(flagged_data, flag_extreme & flag_abschange_congruency), aes(clean_datetime, vwc, col = factor(nflags))) +
  # like what this check pulls
  #geom_point(data = subset(flagged_data, (flagbreak|flag_extreme) & flag_rawdiff_congruency), aes(clean_datetime, vwc, col = factor(nflags))) + 
  #geom_point(data = subset(flagged_data, flag_vwc_congruency), aes(clean_datetime, vwc), col = "orange", alpha = 0.5) + 
  #geom_point(data = subset(flagged_data, flag_sensor), aes(clean_datetime, vwc), col = "orchid", alpha = 0.5) + 
  facet_grid(ppt_trt ~ nut_trt)
# don't like using relchange, flags points where sensor drops below 0



# -- REMOVE VALUES AND ANNOTATE -----
# to NA:
# 1. outside plausible range (already done)
# 2. wetup outside precip event
# 3. flatlines
# 4. has flag_sensor (at least two of: flag extreme [flagged stepchange + deviance from quiring et al. 2016], flag break, and flag spike from doriga et al. 2013 )
# 5. congruency flags (outside 5 sd of at least half of the same ppt trt sensors and all but 1 of the same full treatment sensors)
# > any flag_vwc_congruency
# combo of (flagbreak|flagextreme) & flag_rawdiff_congruency
# write some sort of manual code to catch spikes in dry period of summer that isn't caught by outside precip
# > or can re-run outside precip after NA flagged vals

# id which rows have one of the above conditions
# > flagging done on vwc with plausible range flagged vals already removed, so annotate those before the rest
needs_NA <- with(flagged_data, which( (flag_wetup & !irrigated) | flag_sensor | ((flagbreak|flag_extreme) & flag_rawdiff_congruency) | flag_vwc_congruency))
View(flagged_data[needs_NA,])

# copy data before proceeding
copy <- smdat_qa_all_goodtrts

# start with flag_range col as flag note col(qa_note col has timestamp adjustments annotated currently)
smdat_qa_all_goodtrts$flag_note <- smdat_qa_all_goodtrts$flag_range
smdat_qa_all_goodtrts$flag_note <- with(smdat_qa_all_goodtrts, gsub("high warning", "vwc removed: above 0.6 (unlikely high value)", flag_note))
smdat_qa_all_goodtrts$flag_note <- with(smdat_qa_all_goodtrts, gsub("outside low", "vwc removed: below -0.05 (impossible value)", flag_note))
smdat_qa_all_goodtrts$flag_note <- with(smdat_qa_all_goodtrts, gsub("low warning", "vwc needs adjustment: below 0", flag_note))

# add caution for flatline but don't remove data points (too many)
# do any streaks have notes already? (don't think so but to be sure)
with(smdat_qa_all_goodtrts, summary(is.na(flag_note[!is.na(streak_num)]))) # all NA
smdat_qa_all_goodtrts$flag_note[!is.na(smdat_qa_all_goodtrts$streak_num)] <- paste("caution: vwc", smdat_qa_all_goodtrts$flag_streak[!is.na(smdat_qa_all_goodtrts$streak_num)])

# make dataframe with notes to use
#notes_df <- data.frame(flagname = names(flagged_data)[grep("^flag|spike|congruency", names(flagged_data))])
# drop flag_range
notes_df <- flagged_data[needs_NA,]
notes_df$flag_precip <- with(notes_df, flag_wetup & !irrigated)
# keep only identifying and flag cols for annotations
notes_df <- subset(notes_df, select = names(notes_df)[grep("clean|portid|ppt_trt|^flag|spike|congruency", names(notes_df))])
notes_df <- subset(notes_df, select = -c(flag_range, flag_streak, flag_wetup, flag_sensor, flag_relchange_congruency))
summary(notes_df)
# convert flags to notes
notes_df$flag_extreme <- ifelse(notes_df$flag_extreme, "vwc stepchange and spike (Quiring et al. 2016)", NA)
notes_df$spike <- ifelse(notes_df$spike, "vwc spike (Dorigo et al. 2013)", NA)
notes_df$flagbreak <- ifelse(notes_df$flagbreak, "break in vwc (Dorigo et al. 2013)", NA)
notes_df$flag_abschange_congruency <- ifelse(notes_df$flag_abschange_congruency, "absolute stepchange exceeds 5 sd average difference between sensor and comparable treatment sensors and all but one or all same full treatment sensors", NA)
notes_df$flag_vwc_congruency <- ifelse(notes_df$flag_vwc_congruency, "vwc exceeds 5 sd average difference between sensor and >50% same ppt treatment sensors and all but one or all same full treatment sensors", NA)
notes_df$flag_rawdiff_congruency <- ifelse(notes_df$flag_vwc_congruency, "vwc exceeds 0.15 difference between sensor and >50% same ppt treatment sensors and all but one or all same full treatment sensors", NA)
notes_df$flag_precip <- with(notes_df, ifelse(flag_precip & ppt_trt == "W", "wetup without rainfall or possible irrigation in previous 24 hours", 
                                              ifelse(flag_precip & ppt_trt != "W", "wetup without rainfall in previous 24 hours", NA)))

# arrange flag cols in order of clear error (logic error to questionable)
notes_df <- subset(notes_df, select = c(portid:clean_datetime, flag_precip, flag_extreme:flag_rawdiff_congruency))
copy <- smdat_qa_all_goodtrts

# go row by row for other flags (this can take a minute to run)
for(i in 1:nrow(notes_df)){
  # grab notes to use
  tempnotes <- notes_df[i,]
  # remove any col that's NA
  tempnotes <- tempnotes[,sapply(tempnotes, function(x) !is.na(x))]
  # make note
  tempnotes$note <- paste("vwc removed, flag(s):", str_flatten(tempnotes[,5:ncol(tempnotes)], collapse = "; "))
  # id which row in soilmoisture dataset needs annotation
  temprow <- with(smdat_qa_all_goodtrts, which(portid == tempnotes$portid & cleanorder == tempnotes$cleanorder))
  # apply
  smdat_qa_all_goodtrts$vwc[temprow] <- NA
  smdat_qa_all_goodtrts$flag_note[temprow] <- tempnotes$note
}

# look at data
ggplot(smdat_qa_all_goodtrts) +
  geom_point(data = subset(smdat_qa_all_goodtrts, !is.na(flag_note) & !grepl("adjustment|caution", flag_note)), aes(clean_datetime, raw_vwc, group = portid), col = "orchid", alpha = 0.5) +
  geom_line(aes(clean_datetime, vwc, group = portid, col = factor(block)), alpha = 0.5) +
  scale_color_viridis_d() +
  facet_wrap(ppt_trt ~ nut_trt)
# point removed look okay

# look at CXC, FW, and DN to see how to write check to grab orphan outlier points
ggplot(subset(smdat_qa_all_goodtrts, fulltrt == "CXC")) +
  geom_point(aes(clean_datetime, vwc, group = portid), alpha = 0.5) +
  # overlay raw in pink to see what changed
  geom_line(aes(clean_datetime, raw_vwc, group = portid), col = "orchid", alpha = 0.5) +
  facet_wrap(~ portid) # think I'm okay with what's here now that I've tweaked flagging
# b2l4_3 needs clean up in dry period

ggplot(subset(smdat_qa_all_goodtrts, fulltrt == "FW")) +
  geom_point(aes(clean_datetime, vwc, group = portid), alpha = 0.5) +
  geom_line(aes(clean_datetime, raw_vwc, group = portid), col = "orchid", alpha = 0.5) +
  facet_wrap(~ portid)
# b1l4_4 needs spot clean up at end of series

ggplot(subset(smdat_qa_all_goodtrts, fulltrt == "ND")) +
  geom_point(aes(clean_datetime, vwc, group = portid), alpha = 0.5) +
  geom_line(aes(clean_datetime, raw_vwc, group = portid), col = "orchid", alpha = 0.5) +
  facet_wrap(~ portid)
# could maybe catch some points in dry periods (e.g., b2l3_1, b3l4_3, b3l4_5)



# dry/rain-season flagging clean up ------
# write check for long stretches without rain
dryseason <- seasondat %>% #[c("eco_yr", "season_end")]
  # add 24 hrs to the end of rain season
  mutate(season_end = as.POSIXct(season_end),
         drystart = season_end + (24*60*60))
# add dry end
dryseason$dryend <- c(seasondat$season_start[seasondat$eco_yr>2019], as.character(max(smdat_qa_all_goodtrts$clean_datetime)))
# change dry season end to first germ rain of 2019 in sep (even tho growing season starts in nov)
dryseason$dryend[dryseason$eco_yr == 2019] <- as.character(min(gsdat$date_threshold[gsdat$eco_yr == 2020], na.rm = T))
dryseason$dryend <- as.POSIXct(dryseason$dryend)
# move dry end earlier for 2019 and 2020 by 24hrs to be conservative (since season started at when germain threshold was met)
# > 2021 ends before start of water year
dryseason$dryend[dryseason$eco_yr < 2021] <- dryseason$dryend[dryseason$eco_yr < 2021] - (24*60*60)
dryseason$season_start <- as.POSIXct(dryseason$season_start, tz = "UTC")

# pull flag_vwc_rawdiff when flag_both true
flag_dryseason <- flagged_data[-needs_NA,] %>% # remove points that have already been flagged
  subset(flag_rawdiff_congruency | flag_relchange_congruency | flag_vwc_congruency | flag_abschange_congruency) %>%
  # left_join(distinct(smdat_qa_all_goodtrts[c("portid", "fulltrt", "ppt_trt", "nut_trt")])) %>%
  # mutate(yr = year(clean_datetime)) %>%
  left_join(dryseason[c("eco_yr", "drystart", "dryend", "stopwater")], by = c("yr"= "eco_yr")) %>%
  mutate(dryseason = clean_datetime >= drystart & clean_datetime < dryend) %>%
  left_join(cimis_ppt2hr[c("cleanorder", "dryspell")]) #%>%
  #subset(dryseason & !is.na(dryspell))

# see if these will catch bad values in dry periods
ggplot(subset(smdat_qa_all_goodtrts, portid %in% unique(flag_dryseason$portid))) +
  geom_line(aes(clean_datetime, raw_vwc, group = portid), col = "orchid", alpha = 0.5) +
  geom_line(aes(clean_datetime, vwc, group = portid), alpha = 0.5) +
  geom_point(data= subset(flag_dryseason, flag_rawdiff_congruency & dryseason & !is.na(dryspell)), aes(clean_datetime, target_vwc), alpha = 0.5, col = "green") +
  geom_point(data= subset(flag_dryseason, nflags == 3), aes(clean_datetime, target_vwc), alpha = 0.5, col = "orange") +
  #geom_point(data= subset(flag_rawdiff, flag_both & dryseason & !is.na(dryspell)), aes(clean_datetime, target_vwc), alpha = 0.5, col = "orange") +
  facet_wrap(~ portid) 
# > just go with flag_rawdiff and dry season 

# does it work for extreme low values in rainy season? 
ggplot(smdat_qa_all_goodtrts) +
  geom_line(aes(clean_datetime, raw_vwc, group = portid), col = "orchid", alpha = 0.5) +
  geom_line(aes(clean_datetime, vwc, group = portid), alpha = 0.5) +
  geom_point(data= subset(flag_dryseason, flag_rawdiff_congruency & !dryseason), aes(clean_datetime, target_vwc), alpha = 0.5, col = "orange") +
  facet_wrap(~ portid)
# yes for b1l4_4 and b2l4_3 but not others. will need manual clean up on those. not too many

# apply dry season non-rain-event raw diff flags, just flag_both only. manual spot clean up to follow for finishing.
# > also include points for b1l4_4 and b2l4_3 during non-dry season since catches those points
flag_dryseason <- subset(flag_dryseason, flag_rawdiff_congruency & dryseason & !is.na(dryspell)| 
                           flag_rawdiff_congruency & !dryseason & !is.na(dryspell) & portid %in% c("B1L4_4", "B2L4_3"))

# review one more time
ggplot(subset(smdat_qa_all_goodtrts, portid %in% unique(flag_dryseason$portid))) +
  geom_line(aes(clean_datetime, raw_vwc, group = portid), col = "orchid", alpha = 0.5) +
  geom_line(aes(clean_datetime, vwc, group = portid), alpha = 0.5) +
  geom_point(data= flag_dryseason, aes(clean_datetime, target_vwc, col = dryseason, shape = clean_datetime < stopwater), alpha = 0.5) +
  facet_wrap(~ paste(portid, fulltrt)) # ok to proceed

# just in case mess up
copy <- smdat_qa_all_goodtrts
drynote <- "vwc differs from comparable treatment sensors and all but one or all same full treatment sensors by >0.15 during dry season"
wetnote <- "vwc differs from comparable treatment sensors and all but one or all same full treatment sensors by >0.15 during dry spell"
# annotate
for(i in 1:nrow(flag_dryseason)){
  # id which row in soilmoisture dataset needs annotation
  temprow <- with(smdat_qa_all_goodtrts, which(portid == flag_dryseason$portid[i] & cleanorder == flag_dryseason$cleanorder[i]))
  vwc_present <- !is.na(smdat_qa_all_goodtrts$vwc[temprow])
  flagnote_empty <- is.na(smdat_qa_all_goodtrts$flag_note[temprow])
  # assign note whether dry season or outside
  if(flag_dryseason$dryseason[i]){
    tempnote <- drynote
  }else{tempnote <- wetnote}
  #if vwc still there, NA
  if(vwc_present){
    smdat_qa_all_goodtrts$vwc[temprow] <- NA
  }
  # if flag_note empty add new note, otherwise append this flag to list of reasons
  if(flagnote_empty){
    smdat_qa_all_goodtrts$flag_note[temprow] <- paste("vwc removed, flag(s):", tempnote)
  }else{
    # append
    smdat_qa_all_goodtrts$flag_note[temprow] <- paste(smdat_qa_all_goodtrts$flag_note[temprow], tempnote, sep = "; ")
  }
}


# clean up environment before proceeding
rm(flagnote_empty, i, temprow, tempnote, tempnotes, vwc_present, suborder, drynote, wetnote)




# manual spot clean up -----
# which portids still have values that need NA or caution?
# > go by ppt treatment
## control ppt sensors
left_join(subset(smdat_qa_all_goodtrts, ppt_trt == "XC"), dryseason, by = c("waterYear" = "eco_yr")) %>%
  ggplot(aes(clean_datetime, vwc)) +
  geom_ribbon(aes(xmin = drystart, xmax = dryend, group = waterYear), fill = "chocolate1", alpha = 0.25) +
  geom_line(aes(clean_datetime, raw_vwc, group = portid), color = "orchid", alpha = 0.5) +
  geom_point(alpha = 0.2) +
  facet_wrap(~paste(nut_trt, portid))
# > b2l4_1 in summer 2020, 
# > b2l4_3 in summer 2019 (and maybe before start of growing season 2020)
manual_xc <- casefold(c("b2l4_1", "b2l4_3", "b1l2_5", "b1l3_4"), upper = T)

## drought ppt sensors
left_join(subset(smdat_qa_all_goodtrts, ppt_trt == "D"), dryseason, by = c("waterYear" = "eco_yr")) %>%
  ggplot(aes(clean_datetime, vwc)) +
  geom_ribbon(aes(xmin = drystart, xmax = dryend, group = waterYear), fill = "chocolate1", alpha = 0.25) +
  geom_line(aes(clean_datetime, raw_vwc, group = portid), color = "orchid", alpha = 0.5) +
  geom_point(alpha = 0.2) +
  facet_wrap(~paste(nut_trt, portid))
# > b1l3_3 should be removed -- all unreliable data
# > b2l3_1 in summer 2021
# > b3l4_3 summer 2020 and summer 2021
manual_d <- casefold(c("b1L3_3", "b2l3_1","b3l4_3", "b3l4_5"),upper = T)

## wet ppt sensors
left_join(subset(smdat_qa_all_goodtrts, ppt_trt == "W"), dryseason, by = c("waterYear" = "eco_yr")) %>%
  ggplot(aes(clean_datetime, vwc)) +
  geom_ribbon(aes(xmin = drystart, xmax = dryend, group = waterYear), fill = "chocolate1", alpha = 0.25) +
  geom_line(aes(clean_datetime, raw_vwc, group = portid), color = "orchid", alpha = 0.5) +
  geom_point(alpha = 0.2) +
  facet_wrap(~paste(nut_trt, portid))
# b1l4_4 summer 2021
# b3l3_1 some high points in rainy season winter 2019/2020
# b2l2_1 summer 2019
manual_w <- casefold(c("b1l4_4", "b3l3_1", "b2l2_1", "b2l2_4", "b1l2_3"), upper = T)

# create dat of just portids to correct with dry and rain seasons
dryseason2 <- data.frame(dryseason) %>% 
  mutate_at(.vars = names(.)[grep("start|end", names(.))], as.character) %>%
  gather(met, val, season_start:season_end, drystart:dryend) %>%
  mutate(val = as.POSIXct(val, tz = "UTC")) %>%
  left_join(smdat_datetime[c("clean_datetime", "cleanorder")], by = c("val"="clean_datetime"))%>%
  mutate_at(c("val", "cleanorder"), as.character) %>%
  gather(thing, val, val, cleanorder) %>%
  mutate(thing = ifelse(thing == "cleanorder", "order", "time")) %>%
  unite(met, met, thing) %>%
  spread(met, val) %>%
  mutate_at(names(.)[grep("_order", names(.))], as.numeric) %>%
  mutate_at(names(.)[grep("_time", names(.))], function(x) as.POSIXct(x, tz = "UTC"))

# make manual dat with season dat attached
manualdat <- subset(smdat_qa_all_goodtrts, portid %in% c(manual_xc, manual_d, manual_w)) %>%
  left_join(dryseason2, by = c("waterYear" = "eco_yr")) %>%
  mutate(dryseason = clean_datetime > drystart_time & clean_datetime <= dryend_time,
         dryseason_plus = clean_datetime > (drystart_time + (60*60*24*21)) & clean_datetime <= dryend_time,
         growseason = clean_datetime >= season_start_time & clean_datetime <= season_end_time)

ggplot(manualdat, aes(cleanorder, vwc)) +
  geom_ribbon(aes(xmin = season_start_order, xmax = season_end_order, group = waterYear), fill = "lightblue", alpha = 0.25) +
  geom_ribbon(aes(xmin = drystart_order, xmax = dryend_order, group = waterYear), fill = "chocolate1", alpha = 0.25) +
  geom_line(aes(cleanorder, raw_vwc, group = portid), color = "orchid", alpha = 0.5) +
  geom_point(alpha = 0.5) +
  geom_point(data = subset(manualdat, dryseason_plus & vwc > 0.1), col = "red") +
  facet_wrap(~paste(fulltrt, portid))

# go one by one... no quick logic code to catch all

# 1. CXC manual ---- 
ggplot(subset(smdat_qa_all_goodtrts, fulltrt == "CXC"), aes(cleanorder, vwc)) +
  geom_ribbon(data = subset(manualdat, fulltrt == "CXC"), aes(xmin = season_start_order, xmax = season_end_order, group = waterYear), fill = "lightblue", alpha = 0.25) +
  geom_ribbon(data = subset(manualdat, fulltrt == "CXC"), aes(xmin = drystart_order, xmax = dryend_order, group = waterYear), fill = "chocolate1", alpha = 0.25) +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_line() +
  geom_point(data = subset(smdat_qa_all_goodtrts, fulltrt == "CXC" & flag_extreme), aes(col = flag_extreme), alpha = 0.5) + 
  scale_x_continuous(breaks = seq(0,12000,1000)) +
  facet_wrap(~portid)
# b2l4_1 summer 2020 and high point in rain season 2021
# b2l4_3 NA summer 2019

#b2l4_1
b2l4order <- 4000:11000
b2l4_1_removal2020 <- with(smdat_qa_all_goodtrts, which(portid == "B2L4_1" & cleanorder %in% c(7200:8700) & is.na(flag_note) & vwc > 0.07))
b2l4_1_removal2021_1 <- with(smdat_qa_all_goodtrts, which(portid == "B2L4_1" & cleanorder %in% c(9000:9250) & is.na(flag_note) & vwc > 0.35))
b2l4_1_removal2021_2 <- with(smdat_qa_all_goodtrts, which(portid == "B2L4_1" & cleanorder %in% c(9400:9500) & is.na(flag_note) & vwc > 0.29))
b2l4_1_removal2021_3 <- with(smdat_qa_all_goodtrts, which(portid == "B2L4_1" & (cleanorder > 10500 & is.na(flag_note) & vwc > 0.13)))
b2l4_1_removal2021_4 <- with(smdat_qa_all_goodtrts, which(portid == "B2L4_1" & (cleanorder > 10800 & is.na(flag_note) & !is.na(vwc))))
b2l4_1_removal <- c(b2l4_1_removal2020, b2l4_1_removal2021_1, b2l4_1_removal2021_2, b2l4_1_removal2021_3, b2l4_1_removal2021_4)
ggplot(subset(manualdat, portid == "B2L4_1" & cleanorder %in% b2l4order), aes(cleanorder, vwc)) +
  #geom_hline(aes(yintercept = .065), alpha = 0.5, lty = 2, col= "blue") +
  #geom_vline(aes(xintercept = 10830), lty = 2, col= "purple") +
  scale_x_continuous(breaks = seq(7000, 12000, 500))+
  geom_ribbon(aes(xmin = drystart_order, xmax = dryend_order, group = waterYear), fill = "chocolate1", alpha = 0.25) +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_point(alpha =0.5) +
  # add good line
  geom_line(data = subset(smdat_qa_all_goodtrts, portid == "B2L4_2" & cleanorder %in% b2l4order), col = "blue", alpha =0.5) +
  geom_point(data = smdat_qa_all_goodtrts[b2l4_1_removal,], col= "red") # ok


#b2l4_3
b2l4_3order <- 4000:5000
b2l4_3_removal <- with(smdat_qa_all_goodtrts, which(portid == "B2L4_3" & cleanorder %in% c(2500:4500) & is.na(flag_note) & vwc > 0.1))
ggplot(subset(manualdat, portid == "B2L4_3" & cleanorder %in% b2l4_3order), aes(cleanorder, vwc)) +
  #geom_hline(aes(yintercept = .065), alpha = 0.5, lty = 2, col= "blue") +
  #geom_vline(aes(xintercept = 10830), lty = 2, col= "purple") +
  scale_x_continuous(breaks = seq(1000, 12000, 500))+
  #geom_ribbon(aes(xmin = season_start_order, xmax = season_end_order, group = waterYear), fill = "lightblue", alpha = 0.25) +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_point(alpha =0.5) +
  # add good line
  geom_line(data = subset(smdat_qa_all_goodtrts, portid == "B2L4_2" & cleanorder %in% b2l4_3order), col = "blue", alpha =0.5) +
  geom_point(data = smdat_qa_all_goodtrts[b2l4_3_removal,], col= "red") # ok



# 2. FW manual ----  
# B1L4_4
ggplot(subset(smdat_qa_all_goodtrts, fulltrt == "FW"), aes(cleanorder, vwc)) +
  geom_ribbon(data = subset(manualdat, fulltrt == "FW"), aes(xmin = season_start_order, xmax = season_end_order, group = waterYear), fill = "lightblue", alpha = 0.25) +
  geom_ribbon(data = subset(manualdat, fulltrt == "FW"), aes(xmin = drystart_order, xmax = dryend_order, group = waterYear), fill = "chocolate1", alpha = 0.25) +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_line() +
  geom_point(data = subset(smdat_qa_all_goodtrts, fulltrt == "FW" & flag_extreme), aes(col = flag_extreme), alpha = 0.5) +
  scale_y_continuous(breaks = seq(0,0.6,.05)) +
  scale_x_continuous(breaks = seq(0,12000,1000)) +
  facet_wrap(~portid)

# zoom in to b1l4_4
b1l4_4_removal <- with(smdat_qa_all_goodtrts, which(portid == "B1L4_4" & cleanorder >= 10250 & is.na(flag_note) & !is.na(vwc)))
b1l4order <- 10000:11000
ggplot(subset(manualdat, portid == "B1L4_4" & cleanorder %in%  b1l4order), aes(cleanorder, vwc)) +
  #geom_ribbon(aes(xmin = season_start_order, xmax = season_end_order, group = waterYear), fill = "lightblue", alpha = 0.25) +
  #geom_ribbon(aes(xmin = drystart_order, xmax = dryend_order, group = waterYear), fill = "chocolate1", alpha = 0.25) +
  geom_vline(aes(xintercept = 10250), alpha = 0.5, lty = 2, col= "blue") +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_line()  +
  # add good line
  geom_line(data = subset(smdat_qa_all_goodtrts, portid == "B1L4_5" & cleanorder %in% b1l4order), col = "blue", alpha =0.5) +
  geom_point(data = smdat_qa_all_goodtrts[b1l4_4_removal,], col= "red") #>=10315 looks good for threshold to exclude data

# b2l2_1 (summer 2019)
b2l2order <- 0:4000
b2l2_1_removal_1 <- with(smdat_qa_all_goodtrts, which(portid == "B2L2_1" & cleanorder %in% 2000:2100 & is.na(flag_note) & vwc > 0.2))
b2l2_1_removal_2 <- with(smdat_qa_all_goodtrts, which(portid == "B2L2_1" & cleanorder %in% 3000:3500 & is.na(flag_note) & vwc > 0.03))
b2l2_1_removal <- c(b2l2_1_removal_1, b2l2_1_removal_2)
ggplot(subset(manualdat, portid == "B2L2_1" & cleanorder %in%  b2l2order), aes(cleanorder, vwc)) +
  geom_ribbon(aes(xmin = season_start_order, xmax = season_end_order, group = waterYear), fill = "lightblue", alpha = 0.25) +
  #geom_ribbon(aes(xmin = drystart_order, xmax = dryend_order, group = waterYear), fill = "chocolate1", alpha = 0.25) +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_line()  +
  # add good line
  geom_line(data = subset(smdat_qa_all_goodtrts, portid == "B2L2_2" & cleanorder %in% b2l2order), col = "blue", alpha =0.5) +
  geom_point(data = smdat_qa_all_goodtrts[b2l2_1_removal,], col= "red")


# 3. ND manual -----
# B1L3_3
b1l33_removal <- with(smdat_qa_all_goodtrts, which(portid == "B1L3_3" & !is.na(vwc)  & (is.na(flag_note) | grepl("warning", flag_range))))
View(smdat_qa_all_goodtrts[b1l33,])
# > note: most of these don't have flags probably because they require a value for the timestep before

# ND: B2L3_1, B3L4_3, B3L4_5
ggplot(subset(smdat_qa_all_goodtrts, fulltrt == "ND"), aes(cleanorder, vwc)) + 
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_line() +
  geom_point(data = subset(smdat_qa_all_goodtrts, fulltrt == "ND" & flag_extreme), aes(col = flag_extreme), alpha = 0.5) +
  geom_ribbon(data = subset(manualdat, fulltrt == "ND"), aes(xmin = season_start_order, xmax = season_end_order, group = waterYear), fill = "lightblue", alpha = 0.25) +
  geom_ribbon(data = subset(manualdat, fulltrt == "ND"), aes(xmin = drystart_order, xmax = dryend_order, group = waterYear), fill = "chocolate1", alpha = 0.25) +
  scale_y_continuous(breaks = seq(0,0.6,.05)) +
  scale_x_continuous(breaks = seq(0,12000,1000)) +
  facet_wrap(~portid)

# zoom in to b2l3_1
b2l3_1_removal <- with(smdat_qa_all_goodtrts, which(portid == "B2L3_1" & cleanorder >= 11100 & !flag_extreme & vwc > .10))
ggplot(subset(manualdat, portid == "B2L3_1" & cleanorder > 11000), aes(cleanorder, vwc)) +
  geom_hline(aes(yintercept = .10), alpha = 0.5, lty = 2, col= "blue") +
  geom_vline(aes(xintercept = 11100), alpha = 0.5, lty = 2, col= "blue") +
  #geom_ribbon(aes(xmin = season_start_order, xmax = season_end_order, group = waterYear), fill = "lightblue", alpha = 0.25) +
  geom_ribbon(aes(xmin = drystart_order, xmax = dryend_order, group = waterYear), fill = "chocolate1", alpha = 0.25) +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_point(data = smdat_qa_all_goodtrts[b2l3_1_removal,], col= "red") +
  geom_line() # .10 is fine

# zoom in to b3l4_3
b3l4_3_removal2020_1 <- with(smdat_qa_all_goodtrts, which(portid == "B3L4_3" & cleanorder %in% c(6258:6497) & is.na(flag_note) & !is.na(vwc)))
b3l4_3_removal2020_2 <- with(smdat_qa_all_goodtrts, which(portid == "B3L4_3" & cleanorder %in% c(6900:8000) & is.na(flag_note) & !is.na(vwc)))
b3l4_3_removal2021_1 <- with(smdat_qa_all_goodtrts, which(portid == "B3L4_3" & cleanorder %in% c(dryseason2$drystart_order[dryseason2$eco_yr == 2021]:10800) & vwc > 0.072 & is.na(flag_note) & !is.na(vwc)))
b3l4_3_removal2021_2 <- with(smdat_qa_all_goodtrts, which(portid == "B3L4_3" & cleanorder >=10830 & vwc > 0.06 & is.na(flag_note) & !is.na(vwc)))
b3l4_3_removal <- c(b3l4_3_removal2020_1, b3l4_3_removal2020_2, b3l4_3_removal2021_1, b3l4_3_removal2021_2)
b3l4order <- 0000:15000
ggplot(subset(manualdat, portid == "B3L4_3" & cleanorder %in% b3l4order), aes(cleanorder, vwc)) +
  #geom_hline(aes(yintercept = .065), alpha = 0.5, lty = 2, col= "blue") +
  #geom_vline(aes(xintercept = dryseason2$drystart_order[dryseason2$eco_yr == 2021]), lty = 2, col= "chocolate") +
  #geom_vline(aes(xintercept = 10830), lty = 2, col= "purple") +
  geom_ribbon(aes(xmin = season_start_order, xmax = season_end_order, group = waterYear), fill = "lightblue", alpha = 0.25) +
  geom_ribbon(aes(xmin = drystart_order, xmax = dryend_order, group = waterYear), fill = "chocolate1", alpha = 0.25) +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  # add good line
  geom_line(data = subset(smdat_qa_all_goodtrts, portid == "B2L3_2" & cleanorder %in% b3l4order), col = "blue", alpha =0.5) +
  geom_line(alpha =0.5) +
  geom_point(data = smdat_qa_all_goodtrts[b3l4_3_removal,], col= "red") # looks fine


#b3l4_5
b3l4_5_removal <- with(smdat_qa_all_goodtrts, which(portid == "B3L4_5" & cleanorder %in% 8000:8500 & vwc > 0.2 & is.na(flag_note)))
b3l4order <- 7000:9000
ggplot(subset(manualdat, portid == "B3L4_5" & cleanorder %in% b3l4order), aes(cleanorder, vwc)) +
  #geom_hline(aes(yintercept = .065), alpha = 0.5, lty = 2, col= "blue") +
  #geom_vline(aes(xintercept = dryseason2$drystart_order[dryseason2$eco_yr == 2021]), lty = 2, col= "chocolate") +
  #geom_vline(aes(xintercept = 10830), lty = 2, col= "purple") +
  #geom_ribbon(aes(xmin = season_start_order, xmax = season_end_order, group = waterYear), fill = "lightblue", alpha = 0.25) +
  #geom_ribbon(aes(xmin = drystart_order, xmax = dryend_order, group = waterYear), fill = "chocolate1", alpha = 0.25) +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  # add good line
  geom_line(data = subset(smdat_qa_all_goodtrts, portid == "B3L4_4" & cleanorder %in% b3l4order), col = "blue", alpha =0.5) +
  geom_line(alpha =0.5) +
  geom_point(data = smdat_qa_all_goodtrts[b3l4_5_removal,], col= "red") # looks fine



# 4. NW manual -----
# B3L3_1, B2L2_4, B1L2_3
ggplot(subset(smdat_qa_all_goodtrts, fulltrt == "NW"), aes(cleanorder, vwc)) +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_line() +
  geom_point(data = subset(smdat_qa_all_goodtrts, fulltrt == "NW" & flag_extreme), aes(col = flag_extreme), alpha = 0.5) + 
  scale_y_continuous(breaks = seq(0,0.6,.05)) +
  scale_x_continuous(breaks = seq(0,12000,1000)) +
  facet_wrap(~portid)
# zoom in
ggplot(subset(smdat_qa_all_goodtrts,portid %in% c("B2L2_4", "B3L3_1")), aes(cleanorder, vwc)) + 
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_line() +
  geom_point(data = subset(smdat_qa_all_goodtrts, portid %in% c("B2L2_4", "B3L3_1") & flag_extreme), aes(col = flag_extreme), alpha = 0.5) + 
  scale_y_continuous(breaks = seq(0,0.6,.05)) +
  scale_x_continuous(breaks = seq(0,12000,1000)) +
  facet_wrap(~portid, scales = "free")

## b2l2_4
ggplot(subset(smdat_qa_all_goodtrts, portid == "B2L2_4" & cleanorder %in% c(2500:3000, 4000:4250)), aes(cleanorder, vwc)) +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_line() +
  geom_hline(aes(yintercept = .13), lty = 2)
ggplot(subset(smdat_qa_all_goodtrts, portid == "B2L2_4" & cleanorder %in% c(3800:3850)), aes(cleanorder, vwc)) +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_line() +
  geom_hline(aes(yintercept = .13))
#  manual pull
b2l2_4_removal <- with(smdat_qa_all_goodtrts, which(portid == "B2L2_4" & is.na(flag_note) & (
                                                        (cleanorder %in% c(2750:3500, 4100:4200) & vwc >= 0.1) | (cleanorder %in% c(3820:3830) & vwc < 0.13)
                                                        )))
# grab dip in both block 2 lines also
b2l2_5drop <- with(smdat_qa_all_goodtrts, which(portid == c("B2L2_5") & is.na(flag_note) & (cleanorder %in% c(3820:3828) & vwc < 0.13)))
b2l2_4order <- 2750:4250
ggplot(subset(smdat_qa_all_goodtrts, portid == "B2L2_4" & cleanorder %in% b2l2_4order), aes(cleanorder, vwc)) +
  geom_vline(aes(xintercept = 4100), lty = 2) +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_line() +
  # add good line
  geom_line(data = subset(smdat_qa_all_goodtrts, portid == "B2L2_5" & cleanorder %in% b2l2_4order), col = "blue", alpha = 0.6) +
  geom_point(data = smdat_qa_all_goodtrts[b2l2_4_removal,], aes(col = portid)) +
  geom_point(data = smdat_qa_all_goodtrts[b2l2_5drop,], aes(col = portid)) +
  geom_hline(aes(yintercept = .13), lty = 2)  #correct points

# b3l3_1
b3l3_1_removal <- with(smdat_qa_all_goodtrts, which(portid == "B3L3_1" & is.na(flag_note) & cleanorder %in% 592:676 & !is.na(vwc)))
ggplot(subset(manualdat, portid == "B3L3_1" & cleanorder %in% 500:1000), aes(cleanorder, vwc, group = portid)) +
  scale_x_continuous(breaks = seq(100,12000,50)) +
  geom_vline(aes(xintercept = 592), lty = 2) +
  geom_vline(aes(xintercept = 676), lty = 2) +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_line() +
  geom_line(dat = subset(smdat_qa_all_goodtrts, logger == "B3L3" & ppt_trt == "W" & cleanorder %in% 500:1000), col = "blue", alpha = 0.5) +
  geom_point(data = smdat_qa_all_goodtrts[b3l3_1_removal,], aes(col = portid)) # good

# b1l2_3
b1l2_3_removal <- with(smdat_qa_all_goodtrts, which(portid == "B1L2_3" & is.na(flag_note) & !is.na(vwc) & 
                                                      ((waterYear == 2020  & flag_extreme) | (cleanorder > 12300 & vwc > 0.043))
                                                    ))
b1l2_3order <- 12000:15000
ggplot(subset(manualdat, portid == "B1L2_3" & cleanorder %in% b1l2_3order), aes(cleanorder, vwc, group = portid)) +
  #scale_x_continuous(breaks = seq(1000,12000,1000)) +
  #geom_vline(aes(xintercept = 592), lty = 2) +
  geom_hline(aes(yintercept = 0.043), lty = 2) +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_point() +
  geom_line(dat = subset(smdat_qa_all_goodtrts, logger == "B3L3" & ppt_trt == "W" & cleanorder %in% b1l2_3order), col = "blue", alpha = 0.5) +
  geom_point(data = smdat_qa_all_goodtrts[b1l2_3_removal[-c(1:2)],], aes(col = portid)) # good

# 5. NXC manual -----
# need to remove spike in drydown for block 1 adjuted data (B1L2_5)
# b1l2_5
b1l2_5_removal <- with(smdat_qa_all_goodtrts, which(portid == "B1L2_5" & is.na(flag_note) & !is.na(vwc) & cleanorder %in% 10705:11500 & vwc > 0.085))
b1l2_5order <- 10000:15000
ggplot(subset(manualdat, portid == "B1L2_5" & cleanorder %in% b1l2_5order), aes(cleanorder, vwc, group = portid)) +
  #scale_x_continuous(breaks = seq(1000,12000,1000)) +
  geom_vline(aes(xintercept = 10705), lty = 2) +
  geom_vline(aes(xintercept = 11500), lty = 2) +
  geom_hline(aes(yintercept = 0.085), lty = 2) +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_line() +
  geom_line(dat = subset(smdat_qa_all_goodtrts, portid == "B1L2_4" & cleanorder %in% b1l2_5order), col = "blue", alpha = 0.5) +
  geom_point(data = smdat_qa_all_goodtrts[b1l2_5_removal[-c(1:2)],], aes(col = portid)) # good


# 6. FXC manual -----
# block 1 spike (b1l3_4)
b1l3_4_removal <- with(smdat_qa_all_goodtrts, which(portid == "B1L3_4" & is.na(flag_note) & !is.na(vwc) & cleanorder >12200 & vwc > 0.07))
b1l3_4order <- 10000:15000
ggplot(subset(manualdat, portid == "B1L3_4" & cleanorder %in% b1l3_4order), aes(cleanorder, vwc, group = portid)) +
  #scale_x_continuous(breaks = seq(1000,12000,1000)) +
  geom_line(aes(cleanorder, raw_vwc), col = "orchid") +
  geom_line() +
  geom_line(dat = subset(smdat_qa_all_goodtrts, portid == "B1L2_4" & cleanorder %in% b1l3_4order), col = "blue", alpha = 0.5) +
  geom_point(data = smdat_qa_all_goodtrts[b1l3_4_removal[-c(1:2)],], aes(col = portid)) # good


# 7. combine all manual -----
# points removed purely due to manual review (not flagged, but in cluster of flagged points)
unreliable_data <- c(b2l4_1_removal, b2l4_3_removal, #CXC
                     b1l4_4_removal, b2l2_1_removal, #FW
                     b1l33_removal, b2l3_1_removal, b3l4_3_removal, b3l4_5_removal, #ND
                     b2l2_4_removal, b2l2_5drop, b3l3_1_removal, b1l2_3_removal, #NW
                     b1l2_5_removal, #NXC,
                     b1l3_4_removal #XC
)

# does this seem right?
ggplot(manualdat, aes(cleanorder, vwc)) +
  geom_ribbon(aes(xmin = season_start_order, xmax = season_end_order, group = waterYear), fill = "lightblue", alpha = 0.25) +
  geom_ribbon(aes(xmin = drystart_order, xmax = dryend_order, group = waterYear), fill = "chocolate1", alpha = 0.25) +
  geom_line(alpha = 0.5) +
  geom_point(data = smdat_qa_all_goodtrts[unreliable_data,], col = "red", alpha = 0.5) +
  facet_wrap(~paste(fulltrt, portid)) # looks fine -- only removed 1 point for b2l2_5

# see if they have any flags
View(smdat_qa_all_goodtrts[unreliable_data,]) # yes
# save copy in case need to redo
copy <- smdat_qa_all_goodtrts

# pull just the flag annotations
flags <- distinct(subset(notes_df, select = flag_precip:flag_rawdiff_congruency)) %>%
  gather(flag_name, flag_descrip, flag_precip:ncol(.)) %>%
  distinct()
# add flag descrip for rawdiff_congruency
flags$flag_descrip[grep("rawdiff", flags$flag_name)] <- "vwc differs from comparable treatment sensors and all but one or all same full treatment sensors by >0.15"
flags <- rbind(flags, data.frame(flag_name = names(flagged_data)[grep("relchange", names(flagged_data))],
                                 flag_descrip = gsub("absolute", "relative", with(flags, flag_descrip[grepl("abschange", flag_name) & !is.na(flag_descrip)]))))
# remove NAs and wetup flag descrip for irrigated sensors
flags <- subset(flags, !is.na(flag_descrip) & !grepl("irrigat", flag_descrip))
# change flag_precip to flag wetup to match
flags$flag_name[grepl("precip", flags$flag_name)] <- "flag_wetup"

# annotate manual -- can a little time to run
for(i in unreliable_data){
  # NA vwc
  smdat_qa_all_goodtrts$vwc[i] <- NA
  # check if other flags present
  tempflags <- smdat_qa_all_goodtrts[i, flags$flag_name] %>%
    gather(flag_name, flag_val, 1:ncol(.)) %>%
    left_join(flags, by = "flag_name") %>%
    subset(flag_val)
  # if other flags present, append to manual review reason
  if(nrow(tempflags)>0){
    manual_note <- paste0("vwc removed, flag(s): unreliable data (manual review), additional flags: ",
                          str_flatten(tempflags$flag_descrip, collapse = "; "))
  }else{
    manual_note <- "vwc removed, flag(s): unreliable data (manual review)"
  }
  # check if flag annotation already present (e.g., flatline or low value warning)
  flagnote_empty <- is.na(smdat_qa_all_goodtrts$flag_note[i])
  # if flag_note empty add new note, prefix this flag to below 0 warning
  if(flagnote_empty){
    smdat_qa_all_goodtrts$flag_note[i] <- manual_note
  }else{
    # append
    smdat_qa_all_goodtrts$flag_note[i] <- paste(manual_note, smdat_qa_all_goodtrts$flag_note[i], sep = "; ")
  }
}

# review
ggplot(smdat_qa_all_goodtrts, aes(clean_datetime, vwc)) +
  geom_hline(aes(yintercept = 0), col = "red") +
  geom_point(data = subset(smdat_qa_all_goodtrts, grepl("removed", flag_note)), aes(clean_datetime, raw_vwc, group = portid, col = paste(block, comp_trt)), alpha = 0.5) +
  geom_line(aes(group = portid), alpha = 0.5) +
  labs(title = "USDA Compost QA: VWC data removed in flagging workflow (points, colored by block-seeding trt to distinguish sensors)",
       subtitle = Sys.Date()) +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.direction = "horizontal",
        legend.position = c(0.95,1.07),
        legend.justification = "right") +
  scale_color_viridis_d() +
  facet_grid(nut_trt~ppt_trt)

# write out
ggsave(filename = paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/PrelimQA_Figures/Compost_timecorrectedVWC_flagged.pdf"),
       width = 7, height = 4.5, units = "in", scale = 1.5)

# plot data kept, colored by portid to see better
ggplot(smdat_qa_all_goodtrts, aes(clean_datetime, vwc)) +
  # highlight low warning points
  geom_line(aes(group = portid, col = paste(block, comp_trt)), alpha = 0.5) +
  geom_point(data = subset(smdat_qa_all_goodtrts, !grepl("removed|flatline", flag_note) & !is.na(flag_note)), aes(clean_datetime, raw_vwc, group = portid, col = paste(block, comp_trt)), pch = 1, alpha = 0.5) +
  geom_hline(aes(yintercept = 0), col = "red", lty = 3) +
  labs(title = "USDA Compost QA: VWC data kept post-flagging workflow (points == retained but <0 warning)",
       subtitle = Sys.Date()) +
  scale_color_viridis_d(name = "block-seed trt") +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.direction = "horizontal",
        legend.position = c(0.95,1.07),
        legend.justification = "right") +
  facet_grid(nut_trt~ppt_trt)

# write out
ggsave(filename = paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/PrelimQA_Figures/Compost_timecorrectedVWC_QAd.pdf"),
       width = 7, height = 4.5, units = "in", scale = 1.5)


# -- review values removed due to sensor flag (potentially add back in) ----
# go by ppt trt
trt <- "CD"
ggplot(subset(smdat_qa_all_goodtrts, fulltrt == trt), aes(clean_datetime, vwc)) +
  geom_line(data = subset(smdat_qa_all_goodtrts, fulltrt == trt), aes(clean_datetime, raw_vwc), col = "orange", alpha = 0.5) +
  geom_line(aes(group = portid), alpha = 0.5) +
  geom_point(data = subset(smdat_qa_all_goodtrts, fulltrt == trt & !is.na(streak_num) & !is.na(flag_note)), aes(clean_datetime, raw_vwc), col = "orange") +
  geom_point(data = subset(flagged_data,fulltrt == trt & flag_abschange_congruency), aes(clean_datetime, target_vwc), col = "orchid", alpha = 0.5) +
  geom_point(data = subset(flagged_data,fulltrt == trt & flag_vwc_congruency), aes(clean_datetime, target_vwc), col = "blue") +
  facet_wrap(~portid)

# just a handful of portids ID'd to add points back in:
# FW: B2L2_2, B3L1_5 (rainy season 2020)
# NW: B1L2_1 (onset of rainy season 2020) 
plot_grid(
  ggplot(subset(cimis_ppt2hr, cleanorder >= 0), aes(cleanorder, wY_accumulated_ppt, col = waterYear, group = waterYear)) +
    geom_line(lwd = 2) +
    #scale_color_viridis_c() +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = seq(0,12000,1000)),
  ggplot(subset(smdat_qa_all_goodtrts, portid %in% c("B2L2_2", "B3L1_5", "B1L2_1")), aes(cleanorder, raw_vwc, group = portid)) +
  geom_line(col = "orchid") +
  geom_line(aes(cleanorder, vwc)) +
  scale_x_continuous(breaks = seq(0,12000,2000))+
  facet_wrap(~paste(portid, fulltrt), nrow = 3),
  ggplot(subset(cimis_ppt2hr, cleanorder >= 0), aes(cleanorder, ppt_mm, col = waterYear, group = waterYear)) +
    geom_line() +
    #scale_color_viridis_c() +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = seq(0,12000,1000)),
  rel_heights = c(0.5, 1, 0.5),
  nrow = 3, align= "vh"
)

ggplot(subset(smdat_qa_all_goodtrts, portid == "B2L2_2" & cleanorder %in% 20:30), aes(cleanorder, raw_vwc, group = portid)) +
  geom_line(col = "orchid") +
  geom_text(data = subset(smdat_qa_all_goodtrts, portid == "B2L2_2" & cleanorder %in% 20:30 & !is.na(flag_note)), aes(label = cleanorder), col = "orchid") +
  #geom_text(aes(label = cleanorder), col = "orchid") +
  geom_line(aes(cleanorder, vwc)) # put back 26

ggplot(subset(smdat_qa_all_goodtrts, portid == "B3L1_5" & cleanorder %in% 6525:6550), aes(cleanorder, raw_vwc, group = portid)) +
  geom_line(col = "orchid") +
  geom_text(data = subset(smdat_qa_all_goodtrts, portid == "B3L1_5" & cleanorder %in% 6525:6550 & !is.na(flag_note)), aes(label = cleanorder), col = "orchid") +
  #geom_text(aes(label = cleanorder), col = "orchid") +
  geom_line(aes(cleanorder, vwc)) #6538

ggplot(subset(smdat_qa_all_goodtrts, portid == "B1L2_1" & cleanorder %in% 3700:4470), aes(cleanorder, raw_vwc, group = portid)) +
  geom_line(col = "orchid") +
  geom_text(data = subset(smdat_qa_all_goodtrts, portid == "B1L2_1" & cleanorder %in% 3700:4470 & !is.na(flag_note)), aes(label = cleanorder), col = "orchid") +
  geom_line(aes(cleanorder, vwc)) #3760, 4367, 4368
addback <- with(smdat_qa_all_goodtrts, which((portid == "B2L2_2" & cleanorder %in% 20:30 & !is.na(flag_note)) |
                                             (portid == "B3L1_5" & cleanorder %in% 6525:6550 & !is.na(flag_note)) |
                                             (portid == "B1L2_1" & cleanorder %in% 3700:4470 & !is.na(flag_note))
)) 
# review
ggplot(subset(smdat_qa_all_goodtrts, portid %in% c("B2L2_2", "B3L1_5", "B1L2_1")), aes(cleanorder, raw_vwc, group = portid)) +
  geom_line(col = "orchid") +
  geom_line(aes(cleanorder, vwc)) +
  geom_point(data = smdat_qa_all_goodtrts[addback,], aes(col = substr(flag_note, 0, 50))) +
  scale_x_continuous(breaks = seq(0,12000,2000))+
  theme(legend.position = "bottom")+
  facet_wrap(~paste(portid, fulltrt), nrow = 3) #
smdat_qa_all_goodtrts$flag_note[addback] <- gsub("removed,", "allowed (manual review),", smdat_qa_all_goodtrts$flag_note[addback])
smdat_qa_all_goodtrts$vwc[addback] <- smdat_qa_all_goodtrts$raw_vwc[addback]


# annotate sensor treatment corrections -----
# > note what it was and what it is now
for(i in 1:nrow(trts2correct_new)){
  oldrow <- which(trts2correct_raw$portid == trts2correct_new$portid[i])
  trts2correct_new$logger_note[i] <- paste("treatment for portID", trts2correct_new$portid[i],"corrected from block",trts2correct_raw$block[oldrow], "plot", 
        trts2correct_raw$plot[oldrow], trts2correct_raw$fulltrt[oldrow], trts2correct_raw$comp_trt[oldrow], "composition to block",
        trts2correct_new$block[i], "plot",trts2correct_new$plot[i], trts2correct_new$fulltrt[i], trts2correct_new$comp_trt[i], 
        "composition after soil moisture review")
}
# > keep as its own column instead of combining with qa_note
smdat_qa_out <- left_join(smdat_qa_all_goodtrts, trts2correct_new[c("portid", "corrected", "logger_note")]) %>%
  # clean up for writing out
  subset(select= c(logger:raw_datetime, raw_vwc, sourcefile, qa_note, logger_note, flag_note)) %>%
  arrange(portid, cleanorder)

# last check all looks good
ggplot(smdat_qa_out, aes(clean_datetime, vwc, col = factor(block), group = portid, lty = comp_trt)) +
  geom_point(data = subset(smdat_qa_out, grepl("removed", flag_note)), aes(clean_datetime, raw_vwc), pch = 1, alpha = 0.5) +
  geom_line(alpha = 0.5) +
  facet_grid(nut_trt ~ ppt_trt) # okay



# -- FINISHING -----
# timestamp-corrected, flagged + NA'd data out to DataQA folder as intermediate dataset
write.csv(smdat_qa_out, paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/SoilMoisture_compiled2_time-corrected_flagged.csv"), row.names = F)
# also as all_clean (until infill/adjustment script created, if that ever happens)
write.csv(smdat_qa_out, paste0(datpath, "SoilMoisture/SoilMoisture_CleanedData/SoilMoisture_all_clean.csv"), row.names = F)

# >> write out intermediate datasets if need/want to revisit
# growing season (from first germ threshold to last spring rain, dry season [period between], irrigation stop dates)
# rearrange cols
dryseason2 <- dplyr::select(dryseason2, eco_yr, season_start_time, season_start_order, season_end_time, season_end_order,
                            stopwater, drystart_time, drystart_order, dryend_time, dryend_order)
write.csv(dryseason2, paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/growing_season.csv"), row.names = F)

# corrected sensor treatments
# include spatial order in this as well
loggerkey_corrected <- distinct(subset(smdat_qa_all_goodtrts, select = c(logger:comp_trt, corrected, colorder))) %>%
  left_join(trts2correct_new[c("portid", "logger_note")])
# manually add roworder (1-3 within block)
loggerkey_corrected$roworder <- 1 # to start numeric col (block 1 C = 1; block 2 F = 1; block 3 C = 1)
loggerkey_corrected$roworder[loggerkey_corrected$block == 1 & loggerkey_corrected$nut_trt == "N"] <- 2
loggerkey_corrected$roworder[loggerkey_corrected$block == 1 & loggerkey_corrected$nut_trt == "F"] <- 3
loggerkey_corrected$roworder[loggerkey_corrected$block == 2 & loggerkey_corrected$nut_trt == "N"] <- 2
loggerkey_corrected$roworder[loggerkey_corrected$block == 2 & loggerkey_corrected$nut_trt == "C"] <- 3
loggerkey_corrected$roworder[loggerkey_corrected$block == 3 & loggerkey_corrected$nut_trt == "F"] <- 2
loggerkey_corrected$roworder[loggerkey_corrected$block == 3 & loggerkey_corrected$nut_trt == "N"] <- 3
# rearrange cols
loggerkey_corrected <- subset(loggerkey_corrected, select = c(logger:comp_trt, roworder, colorder, logger_note))
# plot to be sure it looks correct
ggplot(loggerkey_corrected, aes(roworder, colorder)) +
  geom_text(aes(label = plot)) +
  facet_wrap(~block) # when facetting, this will be arrayed in the correct order

# annotate position so it's clear
loggerkey_corrected$roworder[loggerkey_corrected$roworder == 1] <- "1 east"
loggerkey_corrected$roworder[loggerkey_corrected$roworder == 2] <- "2 center"
loggerkey_corrected$roworder[loggerkey_corrected$roworder == 3] <- "3 west"
loggerkey_corrected$colorder[loggerkey_corrected$colorder == 1] <- "1 uphill"
loggerkey_corrected$colorder[loggerkey_corrected$colorder == 2] <- "2 middle"
loggerkey_corrected$colorder[loggerkey_corrected$colorder == 3] <- "3 downhill"
ggplot(loggerkey_corrected, aes(1, 1)) +
  geom_text(aes(label = plot, col = fulltrt)) +
  facet_grid(colorder~paste(block, roworder)) # ok

write.csv(loggerkey_corrected, paste0(datpath, "SoilMoisture/SoilMoisture_CleanedData/decagon_loggerkey_corrected.csv"), row.names = F)

# flags (all, not just removed)
# include manual flags in this dataset
flagged_data_out <- dplyr::select(flagged_data, portid, cleanorder, clean_datetime, raw_vwc, flag_range, flag_wetup, irrigated, flag_streak:flagbreak, flag_abschange_congruency:flag_rawdiff_congruency) %>%
  full_join(subset(smdat_qa_all_goodtrts, !is.na(flag_note), select = c(names(.),"flag_note"))) %>%
  left_join(flagged_data[c("portid", "cleanorder", "flag_sensor")]) %>%
  left_join(smdat_qa_out[c("portid", "cleanorder", "vwc")]) %>%
  left_join(flag_dryseason[c("portid", "cleanorder", "dryseason", "dryspell")]) %>%
  left_join(notes_df[c("portid", "cleanorder", "flag_precip")]) %>%
  distinct() %>%
  mutate(vwc_allowed = !is.na(vwc)) %>%
  dplyr::select(portid:irrigated, flag_precip, dryseason:dryspell, flag_streak:flagbreak, flag_sensor, flag_abschange_congruency:flag_note, vwc_allowed) %>%
  arrange(portid, cleanorder)
summary(flagged_data_out)
summary(is.na(flagged_data_out$flag_note))
write.csv(flagged_data_out, paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/SoilMoisture_all_flags.csv"), row.names = F)
write.csv(flags, paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/SoilMoisture_flagkey.csv"), row.names = F)

## derived raindats to CIMIS subfolder 
# clean up: ignore rain and rain2 logic cols, add location/station info so provenance clear
cimis_out <- cbind(distinct(cimis_hrly[c("Stn.Id", "Stn.Name", "CIMIS.Region")]), 
                   subset(cimis_ppt2hr, select = -c(rain, rain2))) %>%
  rename(stnId = Stn.Id, stnName = Stn.Name, cimis_region = CIMIS.Region) %>%
  arrange(clean_datetime)
write.csv(cimis_out, paste0(datpath,"SoilMoisture/SoilMoisture_DataQA/CIMIS_084BrownsValley_pptforSoilMoisture.csv"), row.names = F)
