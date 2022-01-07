# flag sus vwc in usda-compost soil moisture
# author(s): CTW
# questions?: caitlin.t.white@colorado


# notes
# there's a cimis R package! 'cimir'. CTW isn't using it, but good to know. can read in data dynamically instead of download static files.
# https://cran.r-project.org/web/packages/cimir/cimir.pdf

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
cimis_ppt2hr <- subset(cimis_hrly, date%in% smdat$date, select = c(clean_datetime, date, time, mon:dowy, Precip.mm)) %>%
  left_join(distinct(smdat[c("clean_datetime", "date", "cleanorder")])) %>%
  #fill cleanorder upevent to pair hourly data
  fill(cleanorder, .direction = "up") %>%
  #subset to 1 hour before soil moisture dat starts
  subset(clean_datetime >= min(smdat$clean_datetime) - (60*60)) %>%
  # sum precip by cleanorder
  group_by(cleanorder) %>%
  summarise(ppt_mm = sum(Precip.mm, na.rm = T)) %>%
  # add clean datetime and water year info back in
  left_join(distinct(subset(smdat, select = c(cleanorder:dowy)))) %>%
  # track wetup events in rainfall
  mutate(rain = ppt_mm > 0 | (ppt_mm == 0 & lag(ppt_mm) > 0 & lead(ppt_mm > 0)), # if 0 sandwiched between non-zero, count it in same event
         # but if negligible rain sandwiched between zeroes don't count
         rain2 = ppt_mm == 0.1 & lag(ppt_mm) == 0  & lag(ppt_mm,2) == 0 & lead(ppt_mm) == 0 & lead(ppt_mm,2) == 0,
         rain = ifelse(rain2, FALSE, rain),
         #rain2 = ifelse()
         rainevent = as.numeric(NA),
         dryspell = as.numeric(NA)) %>%
  group_by(waterYear) %>%
  mutate(wY_accumulated_ppt = cumsum(ppt_mm)) %>%
  ungroup()
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
    mutate(next_event = lead(start_event) - end_event,
           next_event_days = round(as.numeric(next_event/24),2),
           next_ppt = lead(event_ppt)) %>% # in hours
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
    labs(x = NULL) +
    theme(axis.text.x = element_blank()) +
    facet_grid(ppt_trt~.),
  ggplot(distinct(cimis_ppt2hr[c("clean_datetime", "ppt_mm")]), aes(clean_datetime, ppt_mm)) +
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
rm(trackrain, i, dryspell, rainevent)


# -- TEST CASE -----
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
    mutate(diff = vwc == lag(vwc))
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

allstreak <- check_streak(smdat)

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
    mutate(diff = lead(vwc) - vwc) %>% # subtract target from next (forward) so is left-aligned
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
  mutate(flag_extreme = (abs(zdiff) > 3 & abs(rollz) > 3)) # any z greater than 3
# plot
ggplot(combine_stepdev, aes(clean_datetime, vwc)) +
  geom_line(aes(group = portid)) +
  geom_point(data = subset(combine_stepdev, flag_extreme),aes(clean_datetime, vwc), col = "orchid", size = 2) +
  facet_wrap(~portid)

# try it with two loggers
step_test2 <- check_step(testcase2)
dev_test2 <- check_deviance(testcase2)
stepdev2 <-  left_join(step_test2, dev_test2) %>%
  mutate(flag_extreme = (abs(zdiff) > 3 & abs(rollz) > 3))
# plot
ggplot(stepdev2, aes(clean_datetime, vwc)) +
  geom_line(aes(group = portid), alpha = 0.5) +
  geom_point(data = subset(stepdev2, flag_extreme),aes(clean_datetime, vwc), col = "orchid", size = 2) +
  facet_wrap(~fulltrt)
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
  tempdat$stepchange <- abs(tempdat$vwc / lag(tempdat$vwc))
  tempdat$abschange <- abs(tempdat$vwc - lag(tempdat$vwc))
  # 2. second deriv outside 0.8 to 1.2 (based on calibrated data in Doriga et al. 2013)
  # calculate rel change between lag and lead
  tempdat$step2derv <- abs(lag(tempdat$deriv2) / lead(tempdat$deriv2))
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
facet_wrap(~portid)
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
  tempdat$relchange <- abs((tempdat$vwc - lag(tempdat$vwc)) / tempdat$vwc)
  tempdat$abschange <- abs(tempdat$vwc - lag(tempdat$vwc))
  
  # 2. first deriv at target t not more than 10x the 24hr average 2nd deriv, centered at t
  # calculate 24 hr average of 2nd deriv
  tempdat$mu24_2derv <- rollmean(tempdat$deriv2, 12, align = "center", fill = NA)
  
  # 3. ratios of 2nd derivs (two to crunch)
  tempdat$rat1 <- with(tempdat, round(abs(deriv2 / lead(deriv2)), 0))
  tempdat$rat2 <- with(tempdat, round(abs(lead(deriv2) / lead(deriv2,2)), 0))
  
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
      mutate(diffvwc = round(vwc - lag(vwc),2),
             # increase of at least 15%, or anticipates a wetup of at least 15%, or is still wetting up in timestep after 0.15 increase
             wetup = diffvwc >= 0.15 | (diffvwc > 0.01 & lead(diffvwc) >= 0.15) | (diffvwc >= 0.01 & lag(diffvwc) > 0.15),
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



# -- QA/QC ALL DATA -----
# logic checks:
# plausible range check
smdat_qa <- check_range(smdat)
smdat_qa <- cbind(smdat, flag_range = smdat_qa$flag_range)
# retain original vwc as raw_vwc
smdat_qa$raw_vwc <- smdat_qa$vwc
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
  mutate(flag_extreme = (abs(zdiff) > 4 & abs(rollz) > 4))

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
  mutate(flag_sensor = (flag_extreme & spike) | (flag_extreme & flagbreak) | (spike & flagbreak)) #


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
  facet_wrap(~fulltrt+portid)

# look at outside precip flags only
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
# flag 3 = deviation from near-in-space sensors
# > check on a moving average to smooth out tiny differences between timesteps
# > also screen for absolute spikes


# deviation from same treatment sensors



# -- FINISHING -----
#write.csv(smdat_qa_out, paste0(datpath, "SoilMoisture/SoilMoisture_CleanedData/SoilMoisture_all_clean.csv"), row.names = F)
