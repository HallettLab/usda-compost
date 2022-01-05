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


group_by(smdat, portid) %>%
  subset(is.na(qa_note)) %>%
  summarise(nobs = length(vwc)) %>%
  arrange(nobs) %>%
  data.frame()

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
check_range <- function(dat, id = c("rowid", "portid", "clean_datetime", "cleanorder"), hi = 0.6, low = -0.04){
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
outside_precip

# -- LOGIC CHECKS -----
# flag 1 = value outside of plausible range? ([0-])
# > allow for slight dip below 0 that can be shifted to 0 or up (manual says accurate to 0.03 m3/m3 [3% VWC])
smdat_qa <- subset(smdat, select = c(logger:vwc)) %>%
  mutate(rowid = rownames(.),
         range_flag = ifelse(vwc > 0.6, "high warning",
                             ifelse(vwc < 0 & vwc > -0.04, "low warning", 
                                    ifelse(vwc <= -0.04, "outside range", NA))),
         # start qa_vwc
         qa_vwc = ifelse(grepl("high|outside", range_flag), NA, vwc))
  


# flag 2 = wetup outside of precip event? (and not irrigation)

# track wetup for soilmoisture
rundf <- data.frame()
for(p in unique(smdat$portid)){
  tempdat <- subset(smdat_qa, portid == p) %>%
    mutate(diffvwc = round(vwc - lag(vwc),2),
           wetup = diffvwc >= 0.15 | (diffvwc > 0.01 & lead(diffvwc) > 0.15) | (diffvwc >= 0.01 & lag(diffvwc) > 0.15),
           wetup_event = as.numeric(NA))
  temprun <- list2DF(rle(tempdat$wetup))
  
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
      tempdat$wetup_event[r:tempend] <- wetupevent
      # increment r and rainevent
      r <- tempend+1
      wetupevent <- wetupevent + 1
    }
  }
  rundf <- rbind(rundf, tempdat)
}

# clean up first wetup (installation)
#rundf$wetup_event[which(rundf$wetup_event == 0)] <- NA


accumulated <- group_by(rundf, portid, waterYear) %>%
  mutate(accumulated_sm = cumsum(vwc)) %>% # don't think it sums for portids with missing data
  ungroup() %>%
  left_join(cimis_ppt2hr) %>%
  # screen for wetup events outside of rain events
  mutate(flagwetup = !is.na(wetup_event) & is.na(rainevent)) %>%
  # try more conservative wetupflag
  group_by(date, portid) %>%
  mutate(dryday = sum(ppt_mm, na.rm = T) < 0.02) %>%
  ungroup() %>%
  mutate(flagwetup2 = !is.na(wetup_event) & dryday)

# plot accumulate
ggplot(accumulated, aes(wY_accumulated_ppt, accumulated_sm)) +
  geom_line(aes(col = logger, lty = nut_trt, group = portid)) +
  facet_grid(waterYear~ppt_trt, scales = "free")

# plot soil moisture with flagged values in red (this is not yet taking into account irrigation for wet plots)
ggplot() +
  geom_line(data = subset(accumulated, !flagwetup2), aes(clean_datetime, qa_vwc, group = portid)) +
  geom_point(data = subset(accumulated, flagwetup), aes(clean_datetime, vwc), col = "red", alpha = 0.5) +
  geom_point(data = subset(accumulated, flagwetup2), aes(clean_datetime, vwc), col = "yellow", alpha = 0.5) +
  facet_grid(nut_trt~ppt_trt)

# zoom in on suspect loggers
ggplot() +
  geom_line(data = subset(accumulated, portid %in% unique(accumulated$portid[accumulated$flagwetup])), aes(clean_datetime, vwc, group = portid)) +
  geom_point(data = subset(accumulated, flagwetup), aes(clean_datetime, vwc), col = "red", alpha = 0.5) +
  geom_point(data = subset(accumulated, flagwetup2), aes(clean_datetime, vwc), col = "yellow", alpha = 0.5) +
  facet_grid(nut_trt~ppt_trt)


# -- CONGRUENCY CHECKS ----
# flag 3 = deviation from near-in-space sensors
# > check on a moving average to smooth out tiny differences between timesteps
# > also screen for absolute spikes

drought_congruency <- subset(smdat, ppt_trt == "D")
ggplot(drought_congruency, aes(clean_datetime, vwc))+
  geom_line(aes(group = portid, col = comp_trt)) +
  geom_hline(aes(yintercept = -0.03)) +
  facet_grid(block~nut_trt)
wet_congruency <- subset(smdat, ppt_trt == "W")

control_congruency <- subset(smdat, ppt_trt == "XC")
# deviation from same treatment sensors



# -- FINISHING -----
write.csv(smdat_qa_out, paste0(datpath, "SoilMoisture/SoilMoisture_CleanedData/SoilMoisture_all_clean.csv"), row.names = F)
