# correct timestamps in usda-compost ongoing soil moisture dataset
# author(s): CTW
# questions?: caitlin.t.white@colorado.edu

# SCRIPT PURPOSE:
## in yrs 2 and 3 of 3-yr usda-compost project, timestamps for several loggers got off synch (e.g. years jumped to unrealistic year) so are unreliable
## not all loggers got off synch in the same way at the same times, and not all ports within a given logger got off synch in the same way
## thus, corrections to timestamps are largely based on visual review and manual coding

## code for 2020 corrections were previously in the same script as the raw data compilation loop, but since there are so many corrections splitting workflow into two scripts
## to preserve timestamp corrections made in 2020, this script splits raw data into:
## > 1) all data in by july 2020 (first iteration of timestamp coding corrections)
## > 2) all data appended july 2020 and to end of project (second iteration of timestamp coding corrections)
## there were no major timestamp issues in yr 1 of project (just gaps in data, but timestamps still reliable)


# SCRIPT WORKFLOW:
## i. read in datasets from compilation loop:
## > 1) raw compiled data
## > 2) dataframe noting all timestamp issues and data gaps
## ii. split data into project start -> data in data download june 2020 (period 1); and data downloaded after june 2020 -> project end (period 2)
## > match soilmoisture data to correct timestamps by logger and port:
## iii. treat period 1 as done in summer 2020 (re-run code), compile good-timestamp data for period 1
## iv. treat period 2 (new code), compile good-timestamp data for period 2
## > *note: as of 10/2021, still missing 2021 data for 2 loggers
## v. compile period 1 and period 2 data into one df, run logic checks to make sure all data present as expected
# > compile clean master soil moisture dataset with raw timestamp, correct timestamp, qa notes on affected rows
## >> contains: correct timestamp, data collection sequence order, logger, port, plot and treatment info, vwc, raw timestamp, qa note
## vi. write out .csv to dropbox: USDA-compost/Data/SoilMoisture/SoilMoisture_CleanedData/

# NOTES:
## this script does *not* address true missing data (i.e., infilling) or bad data values (e.g., unrealistic values)
## additional cleaning script(s) should be written for those (as time/desire allows). CTW's suggested workflow is:
## > soilmoisture3_flag_soilmoisture.R (script to flag suspect values: exceed logger sensor range, strongly deviates from trends in other loggers, other suspect values with definable criteria)
## > soilmoisture4_infill_soilmoisture.R (script to infill bad and missing soil moisture values)



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

# read in datasets
## compiled raw soil moisture data
soilmoisture_all <- read.csv(paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/SoilMoisture_compiled_raw.csv"), na.strings = na_vals, strip.white = T)
## table noting what's off overall
adjustdates_all <- read.csv(paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/SoilMoisture_raw_timestampbreaks.csv"), na.strings = na_vals, strip.white = T)

# > auxiliary datasets
## treatment key
trtkey <- read.csv(paste0(datpath, "Compost_TreatmentKey.csv"), na.strings = na_vals, strip.white = T)
## soil data logger lookup table
loggerkey <- read.csv(paste0(datpath, "SoilMoisture/SoilMoisture_RawData/decagon_logger_key.csv"), na.strings = na_vals, strip.white = T)
## prepped hourly CIMIS met data for QA checks after compilation
cimis_hrly <- list.files(paste0(datpath, "CIMIS"), full.names = T) %>%
  subset(grepl(".csv",.)) %>%
  read.csv(na.strings = na_vals) 

# check how data read in
str(soilmoisture_all)
str(adjustdates_all)
# convert date_time to POSIX -- use UTC timezone to avoid NAs with Daylight Savings Time switch at 2am
soilmoisture_all$date_time <- as.POSIXct(soilmoisture_all$date_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
soilmoisture_all$date <- as.Date(soilmoisture_all$date, format = "%Y-%m-%d")
# check date conversion as expected
summary(date(soilmoisture_all$date_time) == soilmoisture_all$date) # yes
summary(soilmoisture_all)

adjustdates_all$date_time <- as.POSIXct(adjustdates_all$date_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
adjustdates_all$filemo <- as.Date(adjustdates_all$filemo, format = "%Y-%m-%d")
summary(adjustdates_all)



# -- SPLIT CORRECTION PERIODS -----
unique(adjustdates_all$filemo)
unique(soilmoisture_all$filename[grepl("B2L4", soilmoisture_all$logger)])
# > all files have nomenclature "[last two digits of year]-" before time, subset by that
# period 1 = project beginning to june 2020 download file
soilmoisture_p1 <- subset(soilmoisture_all, grepl("19-|20-", filename))
# period 2 = data downloaded after june 2020 to project end
soilmoisture_p2 <- subset(soilmoisture_all, grepl("21-", filename))
# check all there (should be TRUE)
(nrow(soilmoisture_p1) + nrow(soilmoisture_p2)) == nrow(soilmoisture_all)
# be sure no NAs in converted date-times
summary(is.na(as.POSIXct(unique(soilmoisture_p1$date_time))))
summary(is.na(as.POSIXct(unique(soilmoisture_all$date_time))))


# -- CORRECT PERIOD 1 -----
## timestamp corrections CTW made in July 2020 for all data downloaded up through June 11 2020 files
# > covers project start (Nov 2018) - June 11 2020
datecheck_p1 <- soilmoisture_p1 %>%
  mutate(timeinterval = as.difftime(time, format = "%H:%M:%S", units = "hours"),
         filemo = substr(filename,6,12),
         filemo = as.Date(filemo, format = "%d%b%y"),
         # create timeid to plot with ppt
         timeid = as.numeric(paste(doy, substr(time,1,2), sep = "."))) %>%
  # crunch interval
  arrange(logger, portid, datorder) %>%
  group_by(portid) %>%
  mutate(timediff = as.numeric(date_time - lag(date_time,1)),
         timediff = ifelse(is.na(timediff) & datorder == 1, 2, timediff),
         diffvwc = round(vwc - lag(vwc, 1),2),
         wetup = ifelse(diffvwc > 0.1, 1, 0),
         increase = ifelse(diffvwc > lag(diffvwc), 1, 0)) %>%
  # check for wetup events) %>%
  ungroup()

str(datecheck_p1)

# add interval event cols
rundf_p1 <- select(datecheck_p1, logger, portid, date_time, datorder, timediff) %>%
  arrange(portid, datorder) %>%
  mutate(intevent = NA, 
         qa_note = NA)

for(i in unique(rundf_p1$portid)){
  if(all(rundf_p1$timediff[rundf_p1$portid == i] == 2)){
    rundf_p1$intevent[rundf_p1$portid == i] <- 1
  }else{
    event <- 1
    tempstart <- which(rundf_p1$portid == i & rundf_p1$datorder == 1)
    tempend <- max(which(rundf_p1$portid == i))
    tempbreaks <- c(which(rundf_p1$portid ==i & rundf_p1$timediff != 2), tempend)
    for(t in tempbreaks){
      if(t == tempbreaks[length(tempbreaks)]){
        rundf_p1$intevent[tempstart:tempend] <- event
      }else{
        rundf_p1$intevent[tempstart:(t-1)] <- event
      }
      # add qa note
      if(event != 1){
        tempnote <- ifelse(rundf_p1$timediff[tempstart] >2 & rundf_p1$timediff[tempstart] <= 8, "needs NA infill", 
                           ifelse(t == tempbreaks[length(tempbreaks)], "begin last run", "needs correct timestamp"))
        rundf_p1$qa_note[tempstart] <- tempnote
      }
      event <- event+1
      tempstart <- t
    }
  }
}

# attach interval event to datecheck
datecheck_p1 <- left_join(datecheck_p1, rundf_p1)

# pull dates that need an adjustment
adjustdates_p1 <- subset(rundf_p1, !is.na(qa_note)) %>%
  # attach filename info
  left_join(distinct(datecheck_p1[c("portid", "fulltrt", "filename", "filemo", "date_time", "datorder")]))


# 2.a. Infill missing intervals with NA ----
# set expected mintime (project starts after 2017)
clean_mintime_p1 <- min(soilmoisture_p1$date_time[year(soilmoisture_p1$date_time) > 2017])
# set expected maxtime (spring 2021 slated as last field season)--can adjust this if needed
clean_maxtime_p1 <- max(soilmoisture_p1$date_time[year(soilmoisture_p1$date_time) > 2017 & year(soilmoisture_p1$date_time) < 2022])   
soilmoisture_clean_p1 <- data.frame(date_time = rep(seq.POSIXt(clean_mintime_p1, clean_maxtime_p1, by = "2 hours"), times = length(unique(soilmoisture_p1$portid)))) %>%
  mutate(portid = rep(unique(soilmoisture_p1$portid), each = length(unique(date_time)))) %>%
  group_by(portid) %>%
  mutate(cleanorder = seq(1, length(date_time), 1)) %>%
  ungroup() %>%
  left_join(soilmoisture_p1[c("portid", "logger", "port", "plotid", "block", "fulltrt", "date_time", "filename", "vwc")]) %>%
  group_by(portid) %>%
  fill(filename, .direction = "downup") %>%
  ungroup()


# 2.b. Triage B3L4 24Apr20 -----
b3l4trt <- gsub("[0-9]", "", unique(datecheck_p1$fulltrt[datecheck_p1$logger == "B3L4"]))

# compare vwc across similar treatments by data order rather than date-time
## actual vwc values
subset(datecheck_p1, grepl(b3l4trt[1], fulltrt)) %>% # trackseq > 2250
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_line(alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates_p1, grepl(b3l4trt[1], fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  ggtitle(paste("troubleshoot B3L4", b3l4trt[1], "date jump (dotted vert line = break)")) +
  geom_smooth(aes(fill = portid)) +
  facet_grid(logger~.)

## t2-t1 (diffed) values
subset(datecheck_p1, grepl(b3l4trt[1], fulltrt)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(datorder, diffvwc, col = portid)) +
  geom_vline(data = mutate(subset(adjustdates_p1, grepl(b3l4trt[1], fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  geom_line(alpha = 0.4) +
  ggtitle(paste("troubleshoot B3L4", b3l4trt[1], "date jump (dotted vert line = break)")) +
  #geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  facet_grid(logger~., scales = "free_x", space = "free_x")


subset(datecheck_p1, grepl(b3l4trt[2], fulltrt)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_vline(data = mutate(subset(adjustdates_p1, grepl(b3l4trt[2], fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  geom_line(alpha = 0.4) +
  geom_smooth(aes(fill = portid)) +
  ggtitle(paste("troubleshoot B3L4", b3l4trt[2], "date jump (dotted vert line = break)")) +
  facet_grid(logger~.)
# > data gap after break, shift series so end lines up with end timestamp collected (and check spikes align with B2L1)

subset(datecheck_p1, grepl(b3l4trt[2], fulltrt)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(datorder, diffvwc, col = portid)) +
  geom_vline(data = mutate(subset(adjustdates_p1, grepl(b3l4trt[2], fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  geom_line(alpha = 0.4) +
  ggtitle(paste("troubleshoot B3L4", b3l4trt[2], "date jump (dotted vert line = break)")) +
  #geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  facet_grid(logger~., scales = "free_x", space = "free_x")



# rule 1: if time of day from end of bad file 2 hrs before time of day of following file, assign new times backwards end to last bad date break
# rule 2: look for moisture spike alignment with comparison

b3l4dates <- subset(datecheck_p1, logger == "B3L4") %>%
  group_by(portid, intevent, filename) %>%
  mutate(minorder = min(datorder),
         maxorder = max(datorder)) %>%
  subset(datorder > 1 & (datorder == minorder | datorder == maxorder)) %>%
  ungroup() %>%
  arrange(portid, datorder) %>%
  mutate(timeinterval = as.numeric(timeinterval),
         hrdiff = ifelse(timeinterval==0, 24-lag(timeinterval), 
                         ifelse(lag(timeinterval) > timeinterval, 24+ (timeinterval  - lag(timeinterval)), timeinterval- lag(timeinterval))),
         hrdiff = ifelse(is.na(qa_note), NA, hrdiff)) %>%
  distinct(logger, date_time, filename, datorder, intevent, qa_note, hrdiff) %>%
  # add clean dates
  left_join(distinct(soilmoisture_clean_p1[c("date_time", "cleanorder")])) %>%
  mutate(difforder = cleanorder - datorder)
unique(b3l4dates$difforder)
# difference was 1 before date break, then 235 off from clean order after (so gap of 234 2-hr intervals.. about 19.5 days)


# all ports have same time breaks
# > 1 is NA to infill for missing timestep
# > other is a date sequence that needs adjustment -- end of sequence seems to match?

# compare differences against ports with good data to get a sense of typical difference between logger-ports 
good_b3l4 <- subset(datecheck_p1, grepl(b3l4trt[1], fulltrt)) %>%
  subset(date_time > b3l4dates$date_time[3]-as.difftime(30, units = "days") & date_time <=b3l4dates$date_time[3]) %>%
  dplyr::select(portid, date_time, vwc) %>%
  distinct() %>%
  spread(portid, vwc) %>%
  gather(badport, vwc, B3L4_1,B3L4_2) %>%
  gather(goodport, goodvwc, B1L3_4:B2L1_5) %>%
  group_by(goodport, badport) %>%
  mutate(vwcrun = cumsum(vwc),
         goodrun = cumsum(goodvwc)) %>%
  ungroup() %>%
  mutate(diffvwc = vwc - goodvwc,
         diffrun = vwcrun - goodrun)

ggplot(subset(good_b3l4, goodport != "B2L1_5"), aes(date_time, diffvwc, col = goodport)) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_line() +
  facet_wrap(~badport)

ggplot(subset(good_b3l4, goodport != "B2L1_5"), aes(date_time, diffrun, col = goodport)) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_line() +
  facet_wrap(~badport)
# clean up
rm(good_b3l4)


searchdata <- subset(datecheck_p1, logger == "B3L4" & datorder >= with(b3l4dates, datorder[grepl("correct timestamp", qa_note)]-1) & filename == with(b3l4dates, filename[grepl("timestamp", qa_note)])) %>%
  group_by(datorder) %>%
  mutate(nobs = sum(wetup, na.rm = T)) %>%
  ungroup() %>%
  subset(datorder >= (min(datorder[wetup == 1 & nobs > 2], na.rm = T)-20) & datorder <= (min(datorder[wetup == 1 & nobs > 2], na.rm = T)+80))

refdata <- subset(datecheck_p1, grepl(str_flatten(b3l4trt, collapse = "|"), fulltrt)  & logger != "B3L4" & datorder >= with(b3l4dates, datorder[grepl("correct timestamp", qa_note)]-1) & grepl("Apr20", filename)) %>% # loggers have different april timestamps as b3l4
  group_by(date_time) %>%
  mutate(nobs = sum(wetup, na.rm = T)) %>%
  ungroup() %>%
  subset(date_time >= (min(date_time[wetup == 1 & nobs > 1], na.rm = T)-as.difftime(20*2, units = "hours")) & date_time <= (min(date_time[wetup == 1 & nobs > 1], na.rm = T)+as.difftime(80*2, units = "hours"))) %>%
  # join clean order
  left_join(distinct(soilmoisture_clean_p1[c("date_time", "cleanorder")]))
unique(refdata$timeid)
unique(searchdata$timeid)
searchdata$timeid2 <- searchdata$timeid+259
plot_grid(ggplot(refdata, aes(timeid, vwc, col = portid)) +
            geom_line() +
            facet_grid(~gsub("[0-9]+", "", fulltrt)),
          ggplot(searchdata, aes(timeid2, vwc, col = portid)) +
            geom_line() +
            facet_grid(~gsub("[0-9]+", "", fulltrt)),
          nrow = 2)

ggplot(refdata, aes(timeid, vwc, col = portid)) +
  geom_line(data = searchdata, aes(timeid2, vwc, group = portid), col = "grey50", lwd = 1.5, alpha = 0.6) +
  geom_line() +
  facet_grid(~gsub("[0-9]+", "", fulltrt))

# does difforder support visual suggestions of new date_time lineup?
min(refdata$cleanorder) - min(searchdata$datorder) # same as difforder after correct timestamps resume in B3L4 loggers
unique(b3l4dates$difforder)
# test plotting by adjusting datorder
plot_grid(ggplot(refdata, aes(cleanorder, vwc, col = portid)) +
            geom_line() +
            facet_grid(~gsub("[0-9]+", "", fulltrt)),
          ggplot(searchdata, aes(datorder+235, vwc, col = portid)) +
            geom_line() +
            facet_grid(~gsub("[0-9]+", "", fulltrt)),
          nrow = 2)
ggplot(refdata, aes(cleanorder, vwc, col = portid)) +
  geom_line(data = searchdata, aes(datorder+235, vwc, group = portid), col = "grey50", lwd = 1.5, alpha = 0.6) +
  geom_line() +
  facet_grid(~gsub("[0-9]+", "", fulltrt))

# what happens if subtract 1 backwards in cleanorder from first good timestamp after time break?
startorder <- with(b3l4dates, cleanorder[grepl("last run", qa_note)])
b3l4test <- left_join(datecheck_p1, distinct(soilmoisture_clean_p1[c("date_time", "cleanorder")])) %>%
  subset(logger == "B3L4") %>%
  arrange(portid, desc(cleanorder)) %>%
  group_by(portid, intevent) %>%
  mutate(cleanorder2 = ifelse(intevent == 3, datorder + with(b3l4dates, difforder[grepl("last run", qa_note)]), cleanorder), 
         cleanorder3 = ifelse(intevent == 3, seq(startorder -length(vwc), startorder-1, 1),cleanorder),
         diffcheck = cleanorder2 - datorder) %>%
  ungroup() %>%
  rename(raw_datetime = date_time)


ggplot(b3l4test, aes(cleanorder2, vwc, col = port, group = portid)) +
  geom_point(data = subset(soilmoisture_clean_p1, grepl(str_flatten(b3l4trt, collapse = "|"), fulltrt) & logger != "B3L4"), aes(cleanorder, vwc, group = portid, col = port)) +
  geom_point() +
  geom_vline(data = mutate(subset(adjustdates_p1, grepl(str_flatten(b3l4trt, collapse = "|"), fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  facet_grid(logger~fulltrt)

ggplot(subset(b3l4test, intevent == 3), aes(cleanorder2, vwc, col = portid, group = portid)) +
  geom_line(data = subset(soilmoisture_clean_p1, grepl(str_flatten(b3l4trt, collapse = "|"), fulltrt) & logger != "B3L4" & cleanorder %in% b3l4test$cleanorder2[b3l4test$intevent == 3]), aes(cleanorder, vwc, group = portid), col = "grey50") +
  geom_line() +
  facet_grid(.~fulltrt) # looks okay now


# correct B3L4 in soilmoisture_clean_p1
clean_b3l4 <- subset(soilmoisture_clean_p1, grepl("B3L4", portid)) %>%
  select(date_time:cleanorder, filename) %>%
  left_join(select(b3l4test, portid, cleanorder3, vwc, raw_datetime, timeinterval, timediff, intevent, qa_note, datorder), by = c("portid", "cleanorder" = "cleanorder3")) %>%
  left_join(distinct(soilmoisture_p1, portid, logger, port, block, plotid, fulltrt))

# create temp master with timestamp-corrected data (NAs not yet addressed)..
soilmoisture_master_p1 <- subset(soilmoisture_clean_p1, !grepl("B3L4", portid), c(date_time, portid, cleanorder, filename, vwc)) %>%
  left_join(distinct(select(soilmoisture_p1, portid, logger, port, block, plotid, fulltrt))) %>%
  select(names(soilmoisture_clean_p1)) %>%
  rbind(clean_b3l4[names(.)])

# check out NAs by logger
with(soilmoisture_master_p1, sapply(split(vwc, portid), function(x) summary(is.na(x)))) #small numbers are just missing data for 2-hr intervals
# clean up
rm(searchdata, refdata, b3l4dates, b3l4trt, b3l4test, startorder)



# 2.c. Triage B2L4 24Apr20 -----
b2l4trt <- gsub("[0-9]", "", unique(datecheck_p1$fulltrt[datecheck_p1$logger == "B2L4"]))

subset(datecheck_p1, grepl(b2l4trt[1], fulltrt) & grepl("Apr20", filename)) %>% # & grepl("Apr20", filename) & trackseq > 2000 trackseq > 2250
  ggplot(aes(datorder, vwc, col = paste(logger, port, sep = "_"))) +
  geom_line(aes(group = paste(logger, port, sep = "_")), alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates_p1, grepl("B2L4", portid)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  #geom_vline(data = subset(datecheck_p1, logger == "B2L4" & timediff != 2 & trackseq > 2000), aes(xintercept = trackseq), lty = 2) +
  geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  ggtitle(paste("troubleshoot B2L4 date jump", b2l4trt[1], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l4trt[1])
#facet_grid(.~filemo, scales = "free_x", space = "free_x")

subset(datecheck_p1, grepl(b2l4trt[1], fulltrt) & grepl("Apr20", filename)) %>% # & grepl("Apr20", filename) & datorder > 2000 datorder > 2250
  ggplot(aes(datorder, diffvwc, col = paste(logger, port, sep = "_"))) +
  geom_line(aes(group = paste(logger, port, sep = "_")), alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates_p1, grepl("B2L4", portid)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  #geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  ggtitle(paste("troubleshoot B2L4 date jump", b2l4trt[1], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l4trt[1])

subset(datecheck_p1, grepl(b2l4trt[2], fulltrt) & grepl("Apr20", filename)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(datorder, vwc, col = paste(logger, port, sep = "_"), group = paste(logger, port, sep = "_"))) +
  geom_line(aes(group = paste(logger, port, sep = "_")), alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates_p1, grepl("B2L4", portid)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  ggtitle(paste("troubleshoot B2L4 date jump", b2l4trt[2], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l4trt[2])
#facet_grid(.~filemo, scales = "free_x", space = "free_x")

subset(datecheck_p1, grepl(b2l4trt[2], fulltrt) & grepl("Apr20", filename)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(datorder, diffvwc, col = paste(logger, port, sep = "_"), group = paste(logger, port, sep = "_"))) +
  geom_line(aes(group = paste(logger, port, sep = "_")), alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates_p1, grepl("B2L4", portid)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  #geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  ggtitle(paste("troubleshoot B2L4 date jump", b2l4trt[2], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l4trt[1])

# > slight data gap (spikes later in sequence should line up on same day [rain event])
b2l4dates <- subset(datecheck_p1, grepl("B2L4", portid)) %>%
  subset(datorder %in% unique(c(datorder[!is.na(qa_note)], datorder[!is.na(qa_note)]-1))) %>%
  arrange(portid, datorder) %>%
  mutate(port = as.numeric(substr(portid, 6,7)),
         timeinterval = as.numeric(timeinterval),
         hrdiff = ifelse(timeinterval==0, 24-lag(timeinterval), 
                         ifelse(lag(timeinterval) > timeinterval, 24+ (timeinterval  - lag(timeinterval)), timeinterval- lag(timeinterval))),
         hrdiff = ifelse(is.na(qa_note), NA, hrdiff)) %>%
  distinct(logger, portid, date_time, filename, datorder, intevent, qa_note, hrdiff) %>%
  # add clean dates
  left_join(distinct(soilmoisture_clean_p1[c("date_time", "cleanorder")])) %>%
  mutate(difforder = cleanorder - datorder)

# visually assess
searchdata_b2l4 <- subset(datecheck_p1, logger == "B2L4" & datorder >= 4106 & grepl("Apr20", filename)) %>%
  #subset(wetup == 1) %>%
  group_by(datorder) %>%
  mutate(nobs = sum(wetup, na.rm = T)) %>%
  ungroup() %>%
  subset(datorder >= (min(datorder[wetup == 1 & nobs >= 2], na.rm = T)-20) & datorder <= (min(datorder[wetup == 1 & nobs >= 2], na.rm = T)+80))

refdata_b2l4 <- subset(datecheck_p1, grepl(str_flatten(b2l4trt, collapse = "|"), fulltrt)  & logger != "B2L4" & datorder >= 4106 & grepl("Apr20", filename)) %>%
  group_by(date_time) %>%
  mutate(nobs = sum(wetup, na.rm = T)) %>%
  ungroup() %>%
  subset(date_time >= (min(date_time[wetup == 1 & nobs > 1], na.rm = T)-as.difftime(20*2, units = "hours")) & date_time <= (min(date_time[wetup == 1 & nobs > 1], na.rm = T)+as.difftime(80*2, units = "hours")))
unique(refdata_b2l4$timeid)
unique(searchdata_b2l4$timeid)
searchdata_b2l4$datorder2 <- searchdata_b2l4$datorder+16
plot_grid(ggplot(refdata_b2l4, aes(datorder, vwc, col = portid)) +
            geom_line() +
            facet_grid(~gsub("[0-9]+", "", fulltrt)),
          ggplot(searchdata_b2l4, aes(datorder2, vwc, col = portid)) +
            geom_line() +
            facet_grid(~gsub("[0-9]+", "", fulltrt)),
          nrow = 2)

ggplot(refdata_b2l4, aes(datorder, vwc, col = portid)) +
  geom_line(data = searchdata_b2l4, aes(datorder2, vwc, group = portid), col = "grey50", lwd = 1, alpha = 0.6) +
  geom_line() +
  facet_grid(~gsub("[0-9]+", "", fulltrt))

max(datecheck_p1$datorder[datecheck_p1$logger == "B2L4" & grepl("Apr20", datecheck_p1$filename)])
max(datecheck_p1$datorder[datecheck_p1$logger == "B1L1" & grepl("Apr20", datecheck_p1$filename)]) # if add 16 to B2L2 logger, get max points in download period

# pull first period that needs timestamp adjustment and doesn't have any spikes to crosscheck -- would +16 datorder match there too?
searchdata_b2l4_p1 <-  subset(datecheck_p1, logger == "B2L4" & datorder %in% c(4100:4384) & grepl("Apr20", filename)) %>%
  mutate(datorder2 = ifelse(datorder < 4107, datorder, datorder + 12))
refdata_b2l4_p1 <- subset(datecheck_p1, grepl(str_flatten(b2l4trt,collapse = "|"), fulltrt) & logger != "B2L4" & grepl("Apr20", filename)) %>% # #block == 2 & nut_trt == "C" &
  subset(datorder %in% c(4100:4396))
#subset(date_time > b2l2dates$date_time[1] & date_time < (b2l2dates$date_time[1] + as.difftime(2*350, units = "hours")))

plot_grid(ggplot(refdata_b2l4_p1, aes(datorder, vwc, col = portid)) +
            geom_line() +
            facet_grid(~gsub("[0-9]+", "", fulltrt)),
          ggplot(searchdata_b2l4_p1, aes(datorder2, vwc, col = portid)) +
            geom_line() +
            facet_grid(~gsub("[0-9]+", "", fulltrt)),
          nrow = 2)

ggplot(refdata_b2l4_p1, aes(datorder, vwc, col = portid)) +
  geom_line(data = subset(searchdata_b2l4_p1, portid != "B2L4_3"), aes(datorder2, vwc, group = portid), col = "grey50", alpha = 0.6) +
  geom_line() +
  geom_vline(aes(xintercept = 4107)) +
  facet_grid(~gsub("[0-9]+", "", fulltrt))


# pull first period that needs timestamp adjustment and doesn't have any spikes to crosscheck -- would +16 datorder match there too?
searchdata_b2l4_p2 <-  subset(datecheck_p1, logger == "B2L4" & datorder %in% c(4130:4600) & grepl("Apr20", filename)) %>%
  mutate(datorder2 = ifelse(datorder < 4385, datorder+12, datorder + 16)) #14

refdata_b2l4_p2 <- datecheck_p1 %>%
  #subset(grepl(str_flatten(b2l4trt,collapse = "|"), fulltrt) & logger != "B2L4" & grepl("Apr20", filename)) %>% 
  subset(block == 2 & nut_trt == "C" & logger != "B2L4" & grepl("Apr20", filename)) %>% # #block == 2 & nut_trt == "C" &
  subset(datorder %in% c(4142:4616))

ggplot(refdata_b2l4_p2, aes(datorder, vwc, col = portid)) +
  geom_line(data = subset(searchdata_b2l4_p2, portid != "B2L4_3"), aes(datorder2, vwc, group = portid), col = "grey50", alpha = 0.6) +
  geom_line() +
  geom_vline(aes(xintercept = 4358))


# try infilling backwards and see how that looks, go by interval event
# what happens if subtract 1 backwards in cleanorder from first good timestamp after time break?
startorder_b2l4 <- unique(with(b2l4dates, cleanorder[grepl("last run", qa_note)]))
b2l4test <- left_join(datecheck_p1, distinct(soilmoisture_clean_p1[c("date_time", "cleanorder")])) %>%
  subset(logger == "B2L4") %>%
  arrange(portid, desc(cleanorder)) %>%
  group_by(portid, intevent) %>%
  mutate(cleanorder2 = ifelse(intevent == 4, datorder + with(b2l4dates, difforder[grepl("last run", qa_note)]), cleanorder), 
         cleanorder3 = ifelse(intevent == 4, seq(startorder_b2l4 -length(vwc), startorder_b2l4-1, 1),cleanorder),
         diffcheck = cleanorder2 - datorder)

ggplot(b2l4test, aes(cleanorder2, vwc, col = port, group = portid)) +
  geom_point(data = subset(soilmoisture_clean_p1, grepl(str_flatten(b2l4trt, collapse = "|"), fulltrt) & logger != "B2L4"), aes(cleanorder, vwc, group = portid, col = port)) +
  geom_point() +
  geom_vline(data = mutate(subset(adjustdates_p1, grepl(str_flatten(b2l4trt, collapse = "|"), fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  facet_grid(logger~gsub("[0-9]+", "", fulltrt))

ggplot(subset(b2l4test, intevent == 4), aes(cleanorder2, vwc, col = portid, group = portid)) +
  geom_line(data = subset(soilmoisture_clean_p1, grepl(str_flatten(b2l4trt, collapse = "|"), fulltrt) & logger != "B2L4" & cleanorder %in% b2l4test$cleanorder2[b2l4test$intevent == 4]), aes(cleanorder, vwc, group = portid), col = "grey50") +
  geom_line() +
  facet_grid(.~gsub("[0-9]+", "", fulltrt)) # looks good


# try to work interval 3 backwards from interval 4
startorder3_b2l4 <- unique(with(b2l4test, max(cleanorder2[intevent == 2]))) + 12
b2l4test <- b2l4test %>%
  mutate(cleanorder3 = ifelse(intevent == 3, seq(startorder3_b2l4 +1, startorder3_b2l4+ length(vwc), 1),cleanorder2),
         diffcheck = cleanorder3 - datorder) %>%
  rename(raw_datetime = date_time)
# look at differences between original datorder and final cleanorder
unique(b2l4test$diffcheck)

ggplot(subset(b2l4test, cleanorder3 %in% c(3500:4500)), aes(cleanorder2, vwc, col = port, group = portid)) +
  geom_point(data = subset(soilmoisture_clean_p1, grepl(str_flatten(b2l4trt, collapse = "|"), fulltrt) & logger != "B2L4" & cleanorder %in% c(3500:4500)), aes(cleanorder, vwc, group = portid, col = port)) +
  geom_point() +
  #geom_vline(data = mutate(subset(adjustdates_p1, grepl(str_flatten(b2l4trt, collapse = "|"), fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  facet_grid(logger~gsub("[0-9]+", "", fulltrt))

ggplot(subset(b2l4test, cleanorder3 %in% c(3500:4500) & portid != "B2L4_3"), aes(cleanorder3, vwc, col = is.na(cleanorder2), group = portid)) +
  geom_line(data = subset(soilmoisture_clean_p1, grepl(str_flatten(b2l4trt, collapse = "|"), fulltrt) & logger != "B2L4" & cleanorder %in% c(3500:4500) ), aes(cleanorder, vwc, group = portid), col = "grey50") + #cleanorder %in% b2l4test$cleanorder3[b2l4test$intevent == 3]
  geom_line() +
  facet_grid(.~gsub("[0-9]+", "", fulltrt)) # looks good -- +12 is correct (best fit)

# focus on interval 3 to be sure
ggplot(subset(b2l4test, intevent == 3 & portid != "B2L4_3"), aes(cleanorder3, vwc, col = is.na(cleanorder2), group = portid)) +
  geom_line(data = subset(soilmoisture_clean_p1, grepl(str_flatten(b2l4trt, collapse = "|"), fulltrt) & logger != "B2L4" & cleanorder %in% b2l4test$cleanorder3[b2l4test$intevent == 3]), aes(cleanorder, vwc, group = portid), col = "grey50") + #
  geom_line() +
  facet_grid(.~gsub("[0-9]+", "", fulltrt)) # looks ok


# correct B2L4 in soilmoisture_clean_p1
clean_b2l4 <- subset(soilmoisture_clean_p1, grepl("B2L4", portid)) %>%
  select(date_time:cleanorder, filename) %>%
  left_join(select(b2l4test, portid, cleanorder3, vwc, raw_datetime, timeinterval, timediff, intevent, qa_note, datorder), by = c("portid", "cleanorder" = "cleanorder3")) %>%
  left_join(distinct(soilmoisture_p1, portid, logger, port, block, plotid, fulltrt))

# create temp master with timestamp-corrected data (NAs not yet addressed)..
soilmoisture_master_p1 <- subset(soilmoisture_master_p1, !grepl("B2L4", portid)) %>%
  rbind(clean_b2l4[names(.)])
# check out NAs by logger
with(soilmoisture_master_p1, sapply(split(vwc, portid), function(x) summary(is.na(x))))
# check original for NA count
with(soilmoisture_p1[soilmoisture_p1$logger == "B2L4",], sapply(split(vwc, portid), function(x) summary(is.na(x))))
#> high count for B2L4_3 was already there, not a matter of poor join 

# clean up
rm(searchdata_b2l4, refdata_b2l4, searchdata_b2l4_p1, searchdata_b2l4_p2, b2l4dates, b2l4trt, b2l4test)




# 2.d. Triage B2L2 24Apr20 -----
b2l2trt <- gsub("[0-9]", "", unique(datecheck_p1$fulltrt[datecheck_p1$logger == "B2L2"]))

subset(datecheck_p1, grepl(b2l2trt[1], fulltrt)) %>% # datorder > 2250
  ggplot(aes(datorder, vwc, col = port, group = portid)) +
  geom_line(alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates_p1, grepl("B2L2", portid)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note, col = as.numeric(substr(portid,6,6)))) +
  ggtitle(paste("troubleshoot B2L2 date jump", b2l2trt[1], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l2trt[1])

subset(datecheck_p1, grepl(b2l2trt[1], fulltrt)) %>% # datorder > 2250
  ggplot(aes(datorder, diffvwc, col = port, group = paste(logger, port, sep = "_"))) +
  geom_line(alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates_p1, grepl("B2L2", portid)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note, col = as.numeric(substr(portid,6,6)))) +
  ggtitle(paste("troubleshoot B2L2 date jump", b2l2trt[1], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l2trt[1])

subset(datecheck_p1, grepl(b2l2trt[2], fulltrt)) %>%  #&  & datorder > 2250
  ggplot(aes(datorder, vwc, col = port, group = portid)) +
  geom_vline(data = mutate(subset(adjustdates_p1, grepl(b2l2trt[2], fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note, col = as.numeric(substr(portid,6,6)))) +
  #geom_vline(aes(xintercept = max(unique(with(datecheck_p1, datorder[logger == "B2L2" & as.numeric(timediff) != 2])), na.rm = T)), lty = 2) +
  geom_line(alpha = 0.4) +
  ggtitle(paste("troubleshoot B2L2 date jump", b2l2trt[2], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l2trt[2])

subset(datecheck_p1, grepl(b2l2trt[2], fulltrt) & grepl("Apr20", filename)) %>%  #&  & datorder > 2250
  ggplot(aes(datorder, diffvwc, col = paste(logger, port, sep = "_"), group = port)) +
  geom_vline(data = mutate(subset(adjustdates_p1, grepl("B2L2", portid)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  geom_line(alpha = 0.4) +
  ggtitle(paste("troubleshoot B2L2 date jump", b2l2trt[2], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l2trt[2])

# > think gap happened after break, and end of series for B2L2 should align with end date -- can gauge by spike in soil moisture in mid-3000s
b2l2dates <- subset(datecheck_p1, grepl("B2L2", portid)) %>%
  subset(datorder %in% unique(c(datorder[!is.na(qa_note)], datorder[!is.na(qa_note)]-1))) %>%
  arrange(portid, datorder) %>%
  mutate(port = as.numeric(substr(portid, 6,7)),
         timeinterval = as.numeric(timeinterval),
         hrdiff = ifelse(timeinterval==0, 24-lag(timeinterval), 
                         ifelse(lag(timeinterval) > timeinterval, 24+ (timeinterval  - lag(timeinterval)), timeinterval- lag(timeinterval))),
         hrdiff = ifelse(is.na(qa_note), NA, hrdiff)) %>%
  distinct(logger, portid, date_time, filename, datorder, intevent, qa_note, hrdiff) %>%
  # add clean dates
  left_join(distinct(soilmoisture_clean_p1[c("date_time", "cleanorder")])) %>%
  mutate(difforder = cleanorder - datorder) %>%
  # note # of loggers date_time - datorder combo is applicable to for subsetting
  group_by(date_time, datorder) %>%
  mutate(nports = length(portid)) %>%
  ungroup()
# check out differces between datorder and cleanorder 
unique(b2l2dates$difforder)

# split out by ports -all doing slightly different things (different datorder, cleanorder numbers)
b2l2dates_1 <- subset(b2l2dates, portid == "B2L2_1") %>% subset(datorder %in% unique(c(datorder[!is.na(qa_note)], datorder[!is.na(qa_note)]-1)) | datorder == max(datorder))
b2l2dates_2 <- subset(b2l2dates, portid == "B2L2_2") %>% subset(datorder %in% unique(c(datorder[!is.na(qa_note)], datorder[!is.na(qa_note)]-1)) | datorder == max(datorder))
b2l2dates_3 <- subset(b2l2dates, portid == "B2L2_3" )%>% subset(datorder %in% unique(c(datorder[!is.na(qa_note)], datorder[!is.na(qa_note)]-1)) | datorder == max(datorder))
b2l2dates_4 <- subset(b2l2dates, portid == "B2L2_4") %>% subset(datorder %in% unique(c(datorder[!is.na(qa_note)], datorder[!is.na(qa_note)]-1)) | datorder == max(datorder))
b2l2dates_5 <- subset(b2l2dates, portid == "B2L2_5") %>% subset(datorder %in% unique(c(datorder[!is.na(qa_note)], datorder[!is.na(qa_note)]-1)) | datorder == max(datorder))

# try infilling backwards and see how that looks, go by interval event
# what happens if subtract 1 backwards in cleanorder from first good timestamp after time break?
startorder_b2l2.1 <- unique(with(b2l2dates_1, cleanorder[grepl("last run", qa_note)]))

# try b2l2.1 first 
searchdata_b2l2.1 <- subset(datecheck_p1, portid == "B2L2_1" & datorder >= with(b2l2dates_1, datorder[intevent == 5 & grepl("timestamp", qa_note)])-1) %>% # filename == with(b3l4dates, filename[grepl("timestamp", qa_note)])
  group_by(datorder) %>%
  mutate(nobs = sum(wetup, na.rm = T)) %>%
  ungroup()
#subset(datorder >= (min(datorder[wetup == 1], na.rm = T)-20) & datorder <= (min(datorder[wetup == 1], na.rm = T)+80))

refdata_b2l2 <- subset(datecheck_p1, grepl(gsub("[0-9]", "", searchdata_b2l2.1$fulltrt[1]), fulltrt)  & logger != "B2L2" & datorder >= with(b2l2dates_1, datorder[intevent == 5 & grepl("timestamp", qa_note)])-1) %>% # & grepl("Apr20", filename)) %>% # loggers have different april timestamps as b3l4
  group_by(date_time) %>%
  mutate(nobs = sum(wetup, na.rm = T)) %>%
  ungroup()
#subset(date_time >= (min(date_time[wetup == 1 & nobs > 1], na.rm = T)-as.difftime(20*2, units = "hours")) & date_time <= (min(date_time[wetup == 1 & nobs > 1], na.rm = T)+as.difftime(80*2, units = "hours"))) %>%
# join clean order
#left_join(distinct(soilmoisture_clean_p1[c("date_time", "cleanorder")]))
unique(refdata_b2l2$timeid)
unique(searchdata_b2l2.1$timeid)
plot_grid(ggplot(refdata_b2l2, aes(datorder, vwc, col = portid)) +
            geom_line() +
            scale_x_continuous(limits = c(4500, 6500)) +
            facet_grid(~gsub("[0-9]+", "", fulltrt)),
          ggplot(searchdata_b2l2.1, aes(datorder, vwc, col = portid)) +
            geom_line() +
            scale_x_continuous(limits = c(4500, 6500)) +
            facet_grid(~gsub("[0-9]+", "", fulltrt)),
          nrow = 2)

b2l2.1_test <- left_join(datecheck_p1, distinct(soilmoisture_clean_p1[c("date_time", "cleanorder")])) %>%
  subset(portid == "B2L2_1") %>%
  group_by(portid, intevent) %>%
  mutate(cleanorder2 = ifelse(intevent == 6, datorder + with(b2l2dates_1, difforder[grepl("last run", qa_note)]), cleanorder), 
         diffcheck = cleanorder2 - datorder)
# > doesn't work to infill cleanorder backwards.. there is a gap between final interval run and penultimate..

ggplot(b2l2.1_test, aes(cleanorder2, vwc, group = intevent)) +
  geom_line(data = subset(soilmoisture_clean_p1, grepl(gsub("[0-9]", "", b2l2.1_test$fulltrt[1]), fulltrt) & logger != "B2L2"), aes(cleanorder, vwc, group = portid), alpha = 0.3, col = "grey50") + #col = port
  #geom_point() +
  geom_line(col = "orchid") +
  geom_vline(data = mutate(subset(adjustdates_p1, grepl(gsub("[0-9]", "", b2l2.1_test$fulltrt[1]), fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  facet_grid(logger~gsub("[0-9]+", "", fulltrt))

# .. this is going to be very manual/visual matching..


# look at last period for matching since middle period doesn't overlap with any spikes
searchdata_b2l2_end <- subset(datecheck_p1, portid == "B2L2_1" & datorder > 4392 & grepl("Apr20", filename)) %>%
  group_by(datorder) %>%
  mutate(nobs = sum(wetup, na.rm = T)) %>%
  ungroup()
#subset(datorder >= (min(datorder[wetup == 1 & nobs > 1], na.rm = T)-20) & datorder <=  (min(datorder[wetup == 1 & nobs > 1], na.rm = T)+80))

nrows <- max(searchdata_b2l2_end$datorder) - min(searchdata_b2l2_end$datorder)
tempend <- max(datecheck_p1$datorder[grepl("Apr20", datecheck_p1$filename)])
refdata_b2l2_end <- #subset(datecheck_p1, grepl(str_flatten(b2l2trt, collapse = "|"), fulltrt)  & logger != "B2L2" & datorder >= (tempend-nrows) & grepl("Apr20", filename)) %>%
  subset(datecheck_p1, grepl("FW", fulltrt)  & logger != "B2L2" & datorder >= (tempend-nrows) & grepl("Apr20", filename)) %>%
  group_by(date_time) %>%
  mutate(nobs = sum(wetup, na.rm = T)) %>%
  ungroup()
#subset(date_time >= (min(date_time[wetup == 1 & nobs > 1], na.rm = T)-as.difftime(20*2, units = "hours")) & date_time <= (min(date_time[wetup == 1 & nobs > 1], na.rm = T)+as.difftime(80*2, units = "hours")))
unique(refdata_b2l2_end$timeid)
unique(searchdata_b2l2_end$timeid)
searchdata_b2l2_end$datorder2 <- searchdata_b2l2_end$datorder + 396
plot_grid(ggplot(refdata_b2l2_end, aes(datorder, vwc, col = portid)) +
            geom_line() +
            scale_x_continuous(limits = c(5000, 6500)) +
            facet_grid(~gsub("[0-9]+", "", fulltrt)),
          ggplot(searchdata_b2l2_end, aes(datorder2, vwc, col = portid)) +
            geom_line() +
            scale_x_continuous(limits = c(5000, 6500)) +
            facet_grid(~gsub("[0-9]+", "", fulltrt)),
          align = "vh",
          nrow = 2)

ggplot(subset(refdata_b2l2_end), aes(datorder, vwc, col = portid), alph = 0.3) +
  geom_line(data = searchdata_b2l2_end, aes(datorder2, vwc, group = portid), col = "grey50", alpha = 0.6, lwd = 1.2) +
  geom_line() +
  facet_grid(.~gsub("[0-9]+", "", fulltrt)) # seems correct..


# try to match first three intervals.. all loggers doing the same thing..
searchdata_b2l2_early <- subset(datecheck_p1, logger == "B2L2" & intevent %in% unique(with(b2l2dates, intevent[nports == 5])) & datorder >= unique(min(b2l2dates$datorder))-50) %>%
  #subset(datecheck_p1, logger == "B2L2" & datorder < 4394 & datorder >= unique(min(b2l2dates$datorder))-150) %>%
  mutate(datorder2 = ifelse(intevent %in% c(2,3), datorder+3, datorder))

refdata_b2l2_early <- subset(datecheck_p1, logger != "B2L2" & grepl(str_flatten(b2l2trt, collapse = "|"), fulltrt) &  datorder %in% searchdata_b2l2_early$datorder)

ggplot(data = refdata_b2l2_early, aes(datorder, vwc, group = portid), col = "grey50", alpha = 0.5) +
  geom_line() +
  geom_line(data = searchdata_b2l2_early, aes(datorder2, vwc, group = paste(portid, intevent), col = as.factor(port), lty = as.factor(intevent))) +
  facet_grid(. ~ gsub("[0-9]", "", fulltrt))
#theme(legend.position = "none")

# join correct cleanorder and date time to searchdata by datorder2
b2l2early_test <- subset(datecheck_p1, logger == "B2L2" & intevent %in% unique(with(b2l2dates, intevent[nports == 5]))) %>%
  #subset(datecheck_p1, logger == "B2L2" & datorder < 4394 & datorder >= unique(min(b2l2dates$datorder))-150) %>%
  mutate(datorder2 = ifelse(intevent %in% c(2,3), datorder+3, datorder)) %>%
  mutate(cleanorder = datorder2+1) %>%
  rename(raw_datetime = date_time) %>%
  left_join(subset(soilmoisture_clean_p1, grepl("B2L2", portid), c(date_time, portid, cleanorder)), by = c("cleanorder", "portid"))

# correct B3L4 in soilmoisture_clean_p1
clean_b2l2_early <- subset(soilmoisture_clean_p1, grepl("B2L2", portid)) %>%
  select(date_time:cleanorder, filename) %>%
  left_join(select(b2l2early_test, portid, cleanorder, vwc, raw_datetime, timeinterval, timediff, intevent, qa_note), by = c("portid", "cleanorder")) %>%
  left_join(distinct(soilmoisture_p1, portid, logger, port, block, fulltrt))

ggplot(b2l2early_test, aes(cleanorder, vwc, group = portid, col= as.factor(port))) +
  geom_line(data = subset(soilmoisture_master_p1, grepl(str_flatten(b2l2trt, collapse = "|"), fulltrt) & logger != "B2L2"), aes(cleanorder, vwc, group = portid), alpha = 0.3, lwd = 1.5) +
  geom_line() +
  facet_grid(block ~ gsub("[0-9]", "", fulltrt)) # looks okay

# clean up environment and correct by port
rm(refdata_b2l2, refdata_b2l2_early, refdata_b2l2_end, searchdata_b2l2_early, searchdata_b2l2_end, 
   nrows, tempend, startorder_b2l2.1, b2l2.1_test)



# 2.e. Triage B2L2 by port -----------
# create refdata
refdata_b2l2 <- subset(soilmoisture_master_p1, grepl(str_flatten(b2l2trt, "|"), fulltrt) & logger != "B2L2")

## 2.e.1. B2L2_1 -------
# press on with corrections for individual ports after interval 3 ..
searchdata_b2l2.1 <- subset(datecheck_p1, portid == "B2L2_1") %>%
  mutate(datorder2 = ifelse(intevent == 1, datorder, 
                            ifelse(intevent == 2, datorder+3, #3 
                                   ifelse(intevent == 3, datorder + 5, #5
                                          # interval 4 is an NA -- adjust by amount of intervals off
                                          ifelse(intevent == 4, datorder + 5 + ((b2l2dates_1$hrdiff[grepl("NA", b2l2dates_1$qa_note)]-2)/2),
                                                 ifelse(intevent == max(intevent)-1,datorder + 396,
                                                        ifelse(intevent == max(intevent), datorder + with(b2l2dates_1, difforder[grepl("last run", qa_note)]),NA)))))),
         cleanorder = ifelse(intevent == max(intevent), datorder2, datorder2+1)) %>%
  rename(raw_datetime = date_time) %>%
  left_join(distinct(soilmoisture_clean_p1[c("date_time", "cleanorder")])) %>%
  mutate(difforder = cleanorder - datorder)
unique(searchdata_b2l2.1$difforder)

# plot to see how it looks
ggplot(subset(searchdata_b2l2.1, intevent %in% 2:4), aes(cleanorder, vwc)) +
  geom_line(data = subset(refdata_b2l2, grepl(gsub("[0-9]", "", searchdata_b2l2.1$fulltrt[1]), fulltrt) & cleanorder %in% searchdata_b2l2.1$cleanorder[searchdata_b2l2.1$intevent %in% 2:4]), aes(cleanorder, vwc, group = portid), alpha = 0.3) +
  geom_line(aes(col = as.factor(intevent))) +
  scale_x_continuous(breaks = seq(4100,4400, 25))

ggplot(searchdata_b2l2.1, aes(datorder, vwc, col = as.factor(intevent))) +
  geom_line()
ggplot(subset(searchdata_b2l2.1, intevent %in% 2:5), aes(cleanorder, vwc)) +
  geom_line(data = subset(refdata_b2l2, portid %in% c("B1L4_3", "B3L1_4") & cleanorder >= min(searchdata_b2l2.1$cleanorder[searchdata_b2l2.1$intevent == 2])), aes(cleanorder, vwc, group = portid), alpha = 0.3) + # grepl(gsub("[0-9]", "", searchdata_b2l2.1$fulltrt[1]), fulltrt) & logger != "B2L2"
  geom_line(aes(col = as.factor(intevent)))



## 2.e.2. B2L2_2 ----- 
# port 2 has 5 runs -- 5th is final, 2, 3, 4 are timestamp breaks... 
searchdata_b2l2.2 <- subset(datecheck_p1, portid == "B2L2_2") %>%
  mutate(datorder2 = ifelse(intevent == 1, datorder, 
                            ifelse(intevent == 2, datorder+3, 
                                   ifelse(intevent == 3, datorder + 5,
                                          ifelse(intevent == max(intevent)-1,datorder + 396,
                                                 ifelse(intevent == max(intevent), datorder + with(b2l2dates_2, difforder[grepl("last run", qa_note)]),NA))))),
         cleanorder = ifelse(intevent == max(intevent), datorder2, datorder2+min(b2l2dates_2$difforder, na.rm = T))) %>%
  rename(raw_datetime = date_time) %>%
  left_join(distinct(soilmoisture_clean_p1[c("date_time", "cleanorder")])) %>%
  mutate(difforder = cleanorder - datorder)
unique(searchdata_b2l2.2$difforder)

# plot to see how it looks
ggplot(subset(searchdata_b2l2.2, intevent %in% 2:4), aes(cleanorder, vwc)) +
  geom_line(data = subset(refdata_b2l2, grepl(gsub("[0-9]", "", searchdata_b2l2.2$fulltrt[1]), fulltrt) & cleanorder > min(with(searchdata_b2l2.2, cleanorder[intevent==2])) & cleanorder < max(with(searchdata_b2l2.2, cleanorder[intevent==4]))), aes(cleanorder, vwc, group = portid), alpha = 0.3) +
  geom_line(aes(col = as.factor(intevent)))
#scale_x_continuous(breaks = seq(4100,4400, 25))

ggplot(searchdata_b2l2.2, aes(datorder, vwc, col = as.factor(intevent))) +
  geom_line()
ggplot(subset(searchdata_b2l2.2, intevent %in% 2:3), aes(cleanorder, vwc)) +
  geom_line(data = subset(soilmoisture_master_p1, grepl(gsub("[0-9]", "", searchdata_b2l2.2$fulltrt[1]), fulltrt) & logger != "B2L2" & cleanorder >= min(searchdata_b2l2.2$cleanorder[searchdata_b2l2.2$intevent == 2]) & cleanorder < max(searchdata_b2l2.2$cleanorder[searchdata_b2l2.2$intevent == 3])), aes(cleanorder, vwc, group = portid), alpha = 0.3) + # grepl(gsub("[0-9]", "", searchdata_b2l2.2$fulltrt[1]), fulltrt) & logger != "B2L2"
  geom_line(aes(col = as.factor(intevent))) +
  scale_x_continuous(breaks = seq(4100,4400, 20)) +
  facet_grid(block~., scales = "free_y")


## 2.e.3. B2L2_3 ----- 
# port 3 has 5 runs -- 5th is final, 2, 3, 4 are timestamp breaks... similar to port 2
searchdata_b2l2.3 <- subset(datecheck_p1, portid == "B2L2_3") %>%
  mutate(datorder2 = ifelse(intevent == 1, datorder, 
                            ifelse(intevent == 2, datorder+3, 
                                   ifelse(intevent == 3, datorder + 5,
                                          ifelse(intevent == max(intevent)-1,datorder + 396,
                                                 ifelse(intevent == max(intevent), datorder + with(b2l2dates_3, difforder[grepl("last run", qa_note)]),NA))))),
         cleanorder = ifelse(intevent == max(intevent), datorder2, datorder2+min(b2l2dates_3$difforder, na.rm = T))) %>%
  rename(raw_datetime = date_time) %>%
  left_join(distinct(soilmoisture_clean_p1[c("date_time", "cleanorder")])) %>%
  mutate(difforder = cleanorder - datorder)
unique(searchdata_b2l2.3$difforder)

# plot to see how it looks
ggplot(subset(searchdata_b2l2.3, intevent %in% 2:4), aes(cleanorder, vwc)) +
  geom_line(data = subset(refdata_b2l2, grepl(gsub("[0-9]", "", searchdata_b2l2.3$fulltrt[1]), fulltrt) & cleanorder > min(with(searchdata_b2l2.3, cleanorder[intevent==2])) & cleanorder < max(with(searchdata_b2l2.3, cleanorder[intevent==4]))), aes(cleanorder, vwc, group = portid), alpha = 0.3) +
  geom_line(aes(col = as.factor(intevent)))

ggplot(searchdata_b2l2.3, aes(datorder, vwc, col = as.factor(intevent))) +
  geom_line()
ggplot(subset(searchdata_b2l2.3, intevent %in% 2:3), aes(cleanorder, vwc)) +
  geom_line(data = subset(refdata_b2l2, grepl(gsub("[0-9]", "", searchdata_b2l2.3$fulltrt[1]), fulltrt)  & cleanorder >= min(searchdata_b2l2.3$cleanorder[searchdata_b2l2.3$intevent == 2]) & cleanorder < max(searchdata_b2l2.3$cleanorder[searchdata_b2l2.3$intevent == 3])), aes(cleanorder, vwc, group = portid), alpha = 0.3) + # grepl(gsub("[0-9]", "", searchdata_b2l2.3$fulltrt[1]), fulltrt) & logger != "B2L2"
  geom_line(aes(col = as.factor(intevent))) +
  scale_x_continuous(breaks = seq(4100,4400, 20)) +
  facet_grid(block~., scales = "free_y")



## 2.e.4. B2L2_4 ----- 
# port 4 has 6 runs -- 6 is final; 2, 3, and 5 are time breaks; 4 is NA skip
searchdata_b2l2.4 <- subset(datecheck_p1, portid == "B2L2_4") %>%
  mutate(datorder2 = ifelse(intevent == 1, datorder, 
                            ifelse(intevent == 2, datorder+3, #3 
                                   ifelse(intevent == 3, datorder + 5, #5
                                          # interval 4 is an NA -- adjust by amount of intervals off
                                          ifelse(intevent == 4, datorder + 5 + ((b2l2dates_4$hrdiff[grepl("NA", b2l2dates_4$qa_note)]-2)/2),
                                                 ifelse(intevent == max(intevent)-1,datorder + 396,
                                                        ifelse(intevent == max(intevent), datorder + with(b2l2dates_4, difforder[grepl("last run", qa_note)]),NA)))))),
         cleanorder = ifelse(intevent == max(intevent), datorder2, datorder2+ min(b2l2dates_4$difforder, na.rm = T))) %>%
  rename(raw_datetime = date_time) %>%
  left_join(distinct(soilmoisture_clean_p1[c("date_time", "cleanorder")])) %>%
  mutate(difforder = cleanorder - datorder)
unique(searchdata_b2l2.4$difforder)


# plot to see how it looks
ggplot(subset(searchdata_b2l2.4, intevent %in% 2:5), aes(cleanorder, vwc)) +
  geom_line(data = subset(refdata_b2l2, grepl(gsub("[0-9]", "", searchdata_b2l2.4$fulltrt[1]), fulltrt) & cleanorder > min(with(searchdata_b2l2.4, cleanorder[intevent==2])) & cleanorder < max(with(searchdata_b2l2.4, cleanorder[intevent==5]))), aes(cleanorder, vwc, group = portid), alpha = 0.3) +
  geom_line(aes(col = as.factor(intevent)))


ggplot(searchdata_b2l2.4, aes(datorder, vwc, col = as.factor(intevent))) +
  geom_line()
ggplot(subset(searchdata_b2l2.4, intevent %in% 2:4), aes(cleanorder, vwc)) +
  geom_line(data = subset(refdata_b2l2, grepl(gsub("[0-9]", "", searchdata_b2l2.4$fulltrt[1]), fulltrt)  & cleanorder >= min(searchdata_b2l2.4$cleanorder[searchdata_b2l2.4$intevent == 2]) & cleanorder < max(searchdata_b2l2.4$cleanorder[searchdata_b2l2.4$intevent == 4])), aes(cleanorder, vwc, group = portid), alpha = 0.3) + # grepl(gsub("[0-9]", "", searchdata_b2l2.4$fulltrt[1]), fulltrt) & logger != "B2L2"
  geom_line(aes(col = as.factor(intevent))) +
  #scale_x_continuous(breaks = seq(4100,4400, 20)) +
  facet_grid(block~gsub("[0-9]", "", fulltrt), scales = "free_y")

# look at end of int 1 through int 3
ggplot(subset(searchdata_b2l2.4, cleanorder %in% 4000:4400), aes(cleanorder, vwc)) +
  geom_line(data = subset(refdata_b2l2, grepl(gsub("[0-9]", "", searchdata_b2l2.4$fulltrt[1]), fulltrt)  & cleanorder %in% 4000:4400), aes(cleanorder, vwc, group = portid), alpha = 0.3) + # grepl(gsub("[0-9]", "", searchdata_b2l2.4$fulltrt[1]), fulltrt) & logger != "B2L2"
  geom_line(aes(col = as.factor(intevent))) +
  scale_x_continuous(breaks = seq(4000,4400, 40)) +
  facet_grid(block~., scales = "free_y") # looks good

# look at end --5 and 6
ggplot(subset(searchdata_b2l2.4, intevent %in% 5:6), aes(cleanorder, vwc)) +
  geom_line(data = subset(refdata_b2l2, grepl(gsub("[0-9]", "", searchdata_b2l2.4$fulltrt[1]), fulltrt)  & cleanorder >= min(searchdata_b2l2.4$cleanorder[searchdata_b2l2.4$intevent == 5]) & cleanorder < max(searchdata_b2l2.4$cleanorder[searchdata_b2l2.4$intevent == 6])), aes(cleanorder, vwc, group = portid), alpha = 0.3) + # grepl(gsub("[0-9]", "", searchdata_b2l2.4$fulltrt[1]), fulltrt) & logger != "B2L2"
  geom_line(aes(col = as.factor(intevent))) +
  #scale_x_continuous(breaks = seq(4100,4400, 20)) +
  facet_grid(block~., scales = "free_y") # looks good



## 2.e.5. B2L2_5 ----- 
# port 5 similar to port 4: has 6 runs -- 6 is final; 2, 3, and 5 are time breaks; 4 is NA skip
searchdata_b2l2.5 <- subset(datecheck_p1, portid == "B2L2_5") %>%
  mutate(datorder2 = ifelse(intevent == 1, datorder, 
                            ifelse(intevent == 2, datorder+3, #3 
                                   ifelse(intevent == 3, datorder + 5, #5
                                          # interval 4 is an NA -- adjust by amount of intervals off
                                          ifelse(intevent == 4, datorder + 5 + ((b2l2dates_5$hrdiff[grepl("NA", b2l2dates_5$qa_note)]-2)/2),
                                                 ifelse(intevent == max(intevent)-1,datorder + 396,
                                                        ifelse(intevent == max(intevent), datorder + with(b2l2dates_5, difforder[grepl("last run", qa_note)]),NA)))))),
         cleanorder = ifelse(intevent == max(intevent), datorder2, datorder2+ min(b2l2dates_5$difforder, na.rm = T))) %>%
  rename(raw_datetime = date_time) %>%
  left_join(distinct(soilmoisture_clean_p1[c("date_time", "cleanorder")])) %>%
  mutate(difforder = cleanorder - datorder)
unique(searchdata_b2l2.5$difforder)


# plot to see how it looks
ggplot(subset(searchdata_b2l2.5, intevent %in% 2:5), aes(cleanorder, vwc)) +
  geom_line(data = subset(refdata_b2l2, grepl(gsub("[0-9]", "", searchdata_b2l2.5$fulltrt[1]), fulltrt) & cleanorder > min(with(searchdata_b2l2.5, cleanorder[intevent==2])) & cleanorder < max(with(searchdata_b2l2.5, cleanorder[intevent==5]))), aes(cleanorder, vwc, group = portid), alpha = 0.3) +
  geom_line(aes(col = as.factor(intevent)))


ggplot(searchdata_b2l2.5, aes(datorder, vwc, col = as.factor(intevent))) +
  geom_line()
ggplot(subset(searchdata_b2l2.5, intevent %in% 2:4), aes(cleanorder, vwc)) +
  geom_line(data = subset(refdata_b2l2, grepl(gsub("[0-9]", "", searchdata_b2l2.5$fulltrt[1]), fulltrt)  & cleanorder >= min(searchdata_b2l2.5$cleanorder[searchdata_b2l2.5$intevent == 2]) & cleanorder < max(searchdata_b2l2.5$cleanorder[searchdata_b2l2.5$intevent == 4])), aes(cleanorder, vwc, group = portid), alpha = 0.3) + # grepl(gsub("[0-9]", "", searchdata_b2l2.5$fulltrt[1]), fulltrt) & logger != "B2L2"
  geom_line(aes(col = as.factor(intevent))) +
  #scale_x_continuous(breaks = seq(4100,4400, 20)) +
  facet_grid(block~.)

# look at end of int 1 through int 3
ggplot(subset(searchdata_b2l2.5, cleanorder %in% 4000:4400), aes(cleanorder, vwc)) +
  geom_line(data = subset(refdata_b2l2, grepl(gsub("[0-9]", "", searchdata_b2l2.5$fulltrt[1]), fulltrt)  & cleanorder %in% 4000:4400), aes(cleanorder, vwc, group = portid), alpha = 0.3) + # grepl(gsub("[0-9]", "", searchdata_b2l2.5$fulltrt[1]), fulltrt) & logger != "B2L2"
  geom_line(aes(col = as.factor(intevent))) +
  scale_x_continuous(breaks = seq(4000,4400, 40)) +
  facet_grid(block~., scales = "free_y") # looks good

# look at end --5 and 6
ggplot(subset(searchdata_b2l2.5, intevent %in% 5:6), aes(cleanorder, vwc)) +
  geom_line(data = subset(refdata_b2l2, grepl(gsub("[0-9]", "", searchdata_b2l2.5$fulltrt[1]), fulltrt)  & cleanorder >= min(searchdata_b2l2.5$cleanorder[searchdata_b2l2.5$intevent == 5]) & cleanorder < max(searchdata_b2l2.5$cleanorder[searchdata_b2l2.5$intevent == 6])), aes(cleanorder, vwc, group = portid), alpha = 0.3) + # grepl(gsub("[0-9]", "", searchdata_b2l2.5$fulltrt[1]), fulltrt) & logger != "B2L2"
  geom_line(aes(col = as.factor(intevent))) +
  #scale_x_continuous(breaks = seq(4100,4400, 20)) +
  facet_grid(block~., scales = "free_y") # looks good


## 2.e.6 Compile all B2L2 ports -----
masterb2l2_clean <- rbind(searchdata_b2l2.1, searchdata_b2l2.2) %>%
  rbind(rbind(searchdata_b2l2.3, searchdata_b2l2.4)) %>% 
  rbind(searchdata_b2l2.5)
unique(masterb2l2_clean$portid)
# attach all possible timestamps
masterb2l2_clean2 <- subset(soilmoisture_clean_p1, grepl("B2L2", portid), c(date_time, cleanorder, portid)) %>%
  left_join(masterb2l2_clean)


# plot everything
ggplot(masterb2l2_clean, aes(cleanorder, vwc, group = paste(portid, intevent), col = as.factor(port))) +
  geom_line() +
  ggtitle("Timestamp-corrected B2L2 soil moisture data") +
  facet_grid(gsub("[0-9]", "", fulltrt)~., scales = "free_y")



# -- 3. COMPILE PERIOD 1 CLEAN DATA ------
# add qa notes for NAs
# want cleanorder, correct timestamp, treatment info, water year info, vwc, qa notes and raw timestamp, and rawdatorder

names(masterb2l2_clean)
names(clean_b2l4)
names(clean_b3l4)
names(soilmoisture_master_p1) # has clean b3l4 and b2l4, but not clean b2l2
names(soilmoisture_clean_p1)
names(soilmoisture_p1)

# remake soilmoisture_master_p1 with all colnames from soilmoisture_p1 + 2 qa cols (qa note and raw timestamp)
# > start with loggers than aren't b3l4, b2l4 and b2l2
# > also reattach logger key and trt key info

soilmoisture_master_p1 <- select(soilmoisture_clean_p1, date_time:cleanorder, filename, vwc) %>%
  subset(!grepl("B3L4|B2L4|B2L2", portid)) %>%
  # add rawdatorder back in 
  left_join(subset(soilmoisture_p1, !grepl("B3L4|B2L4|B2L2", portid), c(date_time, portid, datorder))) %>%
  #summary(is.na(soilmoisture_master_p1))
  #View(subset(soilmoisture_master_p1, is.na(datorder)))
  rename(clean_datetime = date_time) %>%
  # add original datetime in so data user can see where was recorded
  left_join(subset(soilmoisture_p1, !grepl("B3L4|B2L4|B2L2", portid), c(date_time, portid, datorder))) %>%
  rename(raw_datetime = date_time) %>%
  # add in col for qa_note
  mutate(qa_note = NA) %>%
  rbind(rename(clean_b3l4, clean_datetime = date_time)[names(.)]) %>%
  rbind(rename(clean_b2l4, clean_datetime = date_time)[names(.)]) %>%
  rbind(rename(masterb2l2_clean, clean_datetime = date_time)[names(.)]) %>%
  rename(raw_datorder = datorder) %>%
  mutate(qa_note2 = NA)

# infill qa_notes
for(i in unique(adjustdates_p1$portid)){
  tempadjust <- subset(adjustdates_p1, portid == i)
  for(r in 1:nrow(tempadjust)){
    # append qa note dependent on what done
    ## NA infill
    if(grepl("NA infill", tempadjust$qa_note[r])){
      # specify start and end points for note
      ## end = where noted in qa_note
      end_raworder <- tempadjust$datorder[r]
      tempend <- with(soilmoisture_master_p1, which(portid == i & raw_datorder == tempadjust$datorder[r]))
      #grab hr diff
      temphr <- (as.numeric(tempadjust$timediff[r])-2)/2 #this is the number of timestep adjustments from expected 2hr interval
      ## start = work backwards
      tempstart <- tempend - temphr
      # note only goes to spot before where indicated
      soilmoisture_master_p1$qa_note2[tempstart:(tempend-1)] <- paste("Data missing, NA added, no timestamp recorded by", i, "for", temphr, "2hr-timesteps")
    }
    ## timestamp shift
    if(grepl("timestamp", tempadjust$qa_note[r])){
      start_raworder <- tempadjust$datorder[r]
      if(r == nrow(tempadjust)){
        end_raworder <- with(soilmoisture_master_p1, max(raw_datorder[portid == i], na.rm =T))
      }else{
        end_raworder <- with(tempadjust, datorder[(r+1)])-1
      }
      # id rows in soilmoisture_master_p1
      tempstart <- with(soilmoisture_master_p1, which(portid == i & raw_datorder == start_raworder))
      tempend <- with(soilmoisture_master_p1, which(portid == i & raw_datorder == end_raworder))
      # pull adjustment amount for note
      tempdiff <- unique(with(soilmoisture_master_p1, cleanorder[tempstart:tempend] - raw_datorder[tempstart:tempend])) %>% na.exclude()
      # annotate
      soilmoisture_master_p1$qa_note2[tempstart:tempend] <- paste("Timestamp off, data sequence shifted", tempdiff, "2hr-timesteps to correct timestamp")
    }
    
  }
  # annotate first timestamp data recorded for each logger in qa_note for cleanorder == 1
  # grab first timestamp
  firstdat <- min(with(soilmoisture_master_p1, clean_datetime[portid == i & raw_datorder >= 1 & !is.na(vwc)]))
  stopifnot(length(firstdat)==1 & firstdat != Inf)
  # id row in soilmoisture_master_p1
  firstrow <- with(soilmoisture_master_p1, which(portid == i & cleanorder == 1))
  # annotate
  soilmoisture_master_p1$qa_note2[firstrow] <- paste("First non-NA data recording for", i, "starts", firstdat) 
  
  # annotate last timestep if no vwc recorded (should be present on next download date, unless at project stop)
  lastdat <- max(with(soilmoisture_master_p1, clean_datetime[portid == i & raw_datorder >= 1 & !is.na(vwc)]))
  stopifnot(length(firstdat)==1 & firstdat != Inf)
  # id row in soilmoisture_master_p1
  lastrow <- with(soilmoisture_master_p1, which(portid == i & cleanorder == max(cleanorder)))
  notelast <- with(soilmoisture_master_p1, is.na(vwc[lastrow] & is.na(raw_datorder[lastrow])))
  # if no vwc present and no raw_datorder, annotate
  if(notelast){
    # annotate
    soilmoisture_master_p1$qa_note2[lastrow] <- paste("Last non-NA data recording for", i, "ends", lastdat) 
  }
}
# add note for any breaks in data (not first or last clean order, filename, vwc and raw_datorder/raw_datetime will be NA)
checkNA <- apply(select(soilmoisture_master_p1, filename:qa_note2), 1, function(x) all(is.na(x)))
View(soilmoisture_master_p1[checkNA,]) # looks okay for infilling
soilmoisture_master_p1$qa_note2[checkNA] <- "Missing data, break in data logger recording"

# create lookup table for all dates, full water years
firstWY <- as.Date(paste0(min(year(soilmoisture_master_p1$clean_datetime)), "-10-1"))
lastWY <- as.Date(paste0(max(year(soilmoisture_master_p1$clean_datetime))+1, "-9-30"))
date_lookup <- data.frame(date = seq.Date(firstWY, lastWY, 1)) %>%
  # add in year, doy, join water year info, etc back in..
  mutate(mon = month(date), 
         doy = yday(date),
         waterYear = ifelse(mon %in% 10:12, year(date)+1, year(date))) %>%
  # create day of water year
  arrange(waterYear, date) %>%
  group_by(waterYear) %>%
  mutate(dowy = seq(1, length(date), 1)) %>%
  ungroup()

# infill filename where missing, add in trt info, logger info, water year info, and clean up
soilmoisture_master_p1 <- soilmoisture_master_p1 %>%
  group_by(portid) %>%
  fill(filename, .direction= "downup") %>%
  ungroup() %>%
  # add in trt info, logger info..
  left_join(distinct(select(soilmoisture_p1, logger:comp_trt)), by = "portid") %>%
  # add in fulltrt info, year, doy, join water year info, etc back in..
  mutate(date = as.Date(clean_datetime),
         time = format(clean_datetime, format = "%H:%M:%S")) %>%
  left_join(date_lookup, by = "date") %>%
  select(-qa_note) %>%
  rename(qa_note = qa_note2,
         sourcefile = filename) %>%
  # rearrange by portid, clean_datetime, and reorder cols
  arrange(portid, cleanorder) %>%
  select(logger, port, portid, plot, plotid, fulltrt, block, nut_trt, ppt_trt, comp_trt, 
         cleanorder, clean_datetime, date, time, mon, doy, waterYear, dowy,
         vwc, raw_datorder, raw_datetime, qa_note, sourcefile) %>%
  data.frame()


# check for NAs and plot data to be sure all looks good
summary(soilmoisture_master_p1)
summary(is.na(soilmoisture_master_p1))

# remake plots from initial compilation loop
logger_plot_clean <- soilmoisture_master_p1 %>%
  ggplot(aes(dowy, vwc, col = as.factor(port))) + #date_time, 
  geom_line(alpha = 0.5) +
  ggtitle(paste0("Compost prelim plot (", Sys.Date(),  "): cleaned soil moisture (VWC),\nby logger, by water year, colored by port")) +
  scale_color_discrete(name = "port") +
  facet_grid(logger~waterYear, scales = "free_x")
logger_plot_clean

trt_plot_clean <- soilmoisture_master_p1 %>%
  ggplot(aes(dowy, vwc, group = portid, col = as.factor(block))) +
  geom_line(alpha = 0.5) +
  ggtitle(paste0("Compost prelim plot (", Sys.Date(),") : cleaned soil moisture (VWC),\nby water year, nutrient x drought treatment, colored by block")) +
  scale_color_discrete(name = "block") +
  facet_grid(fulltrt~waterYear, scales = "free_x")
trt_plot_clean


# write out clean data preliminary plots to if desired
# ggsave(filename = paste0(datpath, "SoilMoisture/SoilMoisture_Figures/Compost_cleanVWC_bylogger.pdf"),
#        plot = logger_plot_clean,
#        width = 8, height = 8, units = "in")
# ggsave(filename = paste0(datpath, "SoilMoisture/SoilMoisture_Figures/Compost_cleanVWC_bytreatment.pdf"),
#        plot = trt_plot_clean,
#        width = 8, height = 8, unit = "in")


# clean up environment before moving on to period 2
rm(searchdata_b2l2.1, searchdata_b2l2.2, searchdata_b2l2.3, searchdata_b2l2.4, searchdata_b2l2.5,
   tempadjust, checkNA, b2l2trt, clean_maxtime_p1, clean_mintime_p1, end_raworder, event, firstdat, firstrow,
   firstWY, i, lastdat, lastrow, lastWY, notelast, r, start_raworder, startorder_b2l4, startorder3_b2l4,
   t, tempbreaks, tempdiff, tempend, tempnote, tempstart, b2l2dates, b2l2early_test, clean_b2l2_early, clean_b2l4,
   clean_b3l4, masterb2l2_clean, masterb2l2_clean2, refdata_b2l2, refdata_b2l4_p1, refdata_b2l4_p2, temphr,
   b2l2dates_1, b2l2dates_2, b2l2dates_3, b2l2dates_4, b2l2dates_5)



# -- CORRECT PERIOD 2 -----
## timestamp corrections for all data downloaded after June 11 2020 download (Sep 2021 files - ongoing)
# > covers June 11 2020 - project end
datecheck_p2 <- soilmoisture_p2 %>%
  mutate(timeinterval = as.difftime(time, format = "%H:%M:%S", units = "hours"),
         filemo = substr(filename,6,12),
         filemo = as.Date(filemo, format = "%d%b%y"),
         # create timeid to plot with ppt
         timeid = as.numeric(paste(doy, substr(time,1,2), sep = "."))) %>%
  # crunch interval
  arrange(logger, portid, datorder) %>%
  group_by(portid) %>%
  # create new datorder col that starts at this period
  mutate(raw_datorder = datorder,
         datorder = 1:length(date),
         # carry on prep as done for period 1
         timediff = as.numeric(date_time - lag(date_time,1)),
         timediff = ifelse(is.na(timediff) & datorder == 1, 2, timediff),
         diffvwc = round(vwc - lag(vwc, 1),2),
         wetup = ifelse(diffvwc > 0.1, 1, 0),
         increase = ifelse(diffvwc > lag(diffvwc), 1, 0)) %>%
  # check for wetup events) %>%
  ungroup()

str(datecheck_p2)

# add interval event cols
rundf_p2 <- select(datecheck_p2, logger, portid, date_time, datorder, timediff) %>%
  arrange(portid, datorder) %>%
  mutate(intevent = NA, 
         qa_note = NA)

for(i in unique(rundf_p2$portid)){
  if(all(rundf_p2$timediff[rundf_p2$portid == i] == 2)){
    rundf_p2$intevent[rundf_p2$portid == i] <- 1
  }else{
    event <- 1
    tempstart <- which(rundf_p2$portid == i & rundf_p2$datorder == 1)
    tempend <- max(which(rundf_p2$portid == i))
    tempbreaks <- c(which(rundf_p2$portid ==i & rundf_p2$timediff != 2), tempend)
    for(t in tempbreaks){
      if(t == tempbreaks[length(tempbreaks)]){
        rundf_p2$intevent[tempstart:tempend] <- event
      }else{
        rundf_p2$intevent[tempstart:(t-1)] <- event
      }
      # add qa note
      if(event != 1){
        tempnote <- ifelse(rundf_p2$timediff[tempstart] >2 & rundf_p2$timediff[tempstart] <= 8, "needs NA infill", 
                           ifelse(t == tempbreaks[length(tempbreaks)], "begin last run", "needs correct timestamp"))
        rundf_p2$qa_note[tempstart] <- tempnote
      }
      event <- event+1
      tempstart <- t
    }
  }
}

# attach interval event to datecheck
datecheck_p2 <- left_join(datecheck_p2, rundf_p2)

# pull dates that need an adjustment
adjustdates_p2 <- subset(rundf_p2, !is.na(qa_note)) %>%
  # attach filename info
  left_join(distinct(datecheck_p2[c("portid", "fulltrt", "filename", "filemo", "date_time", "datorder")]))


# 4.a. Infill missing intervals with NA ----
# set expected mintime (period 2 starts at min datorder number as long as year is 2020)
clean_mintime_p2 <- with(datecheck_p2, min(date_time[datorder == 1 & year(date_time) > 2019]))
# set expected maxtime -- last field visit was in sep (NS downloaded logger dat for AS)
clean_maxtime_p2 <- with(datecheck_p2, max(date_time[year(date_time) > 2019 & year(date_time)<2022]))
soilmoisture_clean_p2 <- data.frame(date_time = rep(seq.POSIXt(clean_mintime_p2, clean_maxtime_p2, by = "2 hours"), times = length(unique(soilmoisture_p2$portid)))) %>%
  mutate(portid = rep(unique(soilmoisture_p2$portid), each = length(unique(date_time)))) %>%
  group_by(portid) %>%
  mutate(cleanorder = seq(1, length(date_time), 1)) %>%
  ungroup() %>%
  left_join(soilmoisture_p2[c("portid", "logger", "port", "plotid", "block", "fulltrt", "date_time", "filename", "vwc")]) %>%
  group_by(portid) %>%
  fill(filename, .direction = "downup") %>%
  ungroup()


# who needs triage due to NAs or timestamp jump?
distinct(datecheck_p2[grepl("time", datecheck_p2$qa_note), c("logger", "filename")])
# bad years
distinct(datecheck_p2[!year(datecheck_p2$date) %in% c(2020,2021), c("logger", "filename")])
# > b3l4 got off track at last run, so that's how it got past the qa flagging -- yrs are 2044, 2045

# create clean refdat for comparisons below -- only loggers that had NA infills needed or no issues
with(subset(datecheck_p2), lapply(split(qa_note, logger), unique))
goodrefs <- unique(subset(datecheck_p2, select = c(logger, qa_note))) %>%
  group_by(logger) %>%
  filter(all(is.na(qa_note)))
# subset to good ref loggers
soilmoisture_master_p2 <- subset(soilmoisture_clean_p2, logger %in% goodrefs$logger)
# double check NAs
summary(is.na(soilmoisture_master_p2))



# 4.b. Triage B1L2 15Sep21 -----
# little to no missing data, same # nobs per port (5543)
# time difference of 10hrs at break where need correct timestamp (06-27-2020 2am)
b1l2trt <- gsub("[0-9]", "", unique(datecheck_p2$fulltrt[datecheck_p2$logger == "B1L2"]))

# compare vwc across similar treatments by data order rather than date-time
## actual vwc values
subset(datecheck_p2, grepl(b1l2trt[1], fulltrt)) %>% # trackseq > 2250
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_line(alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates_p2, grepl(b1l2trt[1], fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  ggtitle(paste("troubleshoot B1L2", b1l2trt[1], "date jump (dotted vert line = break)")) +
  geom_smooth(aes(fill = portid)) +
  facet_grid(logger~.) # looks like it's lined up where it should be, just switched to bad timestamps

## t2-t1 (diffed) values
subset(datecheck_p2, grepl(b1l2trt[1], fulltrt)) %>% #
  ggplot(aes(datorder, diffvwc, col = portid)) +
  geom_vline(data = subset(adjustdates_p2, grepl(b1l2trt[1], fulltrt)), aes(xintercept = datorder, lty = qa_note)) +
  geom_line(alpha = 0.4) +
  ggtitle(paste("troubleshoot B1L2", b1l2trt[1], "date jump (dotted vert line = break)")) +
  #geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  facet_grid(logger~., scales = "free_x", space = "free_x")

# why do diffs look weird compared to reference logger?
ggplot(subset(datecheck_p2, logger == "B1L2" & datorder < 100), aes(datorder, diffvwc, col = portid)) +
  geom_line() # don't look weird here

# plot second treatment type logger record
subset(datecheck_p2, grepl(b1l2trt[2], fulltrt)) %>% #
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_vline(data = subset(adjustdates_p2, grepl(b1l2trt[2], fulltrt)), aes(xintercept = datorder, lty = qa_note)) +
  geom_line(alpha = 0.4) +
  geom_smooth(aes(fill = portid)) +
  ggtitle(paste("troubleshoot B1L2", b1l2trt[2], "date jump (vert line = break)")) +
  facet_grid(logger~.)
# > data gap after break, shift series so end lines up with end timestamp collected (and check spikes align with B2L1)

# zoom in on period that's off to see if just needs NAs inserted or if time got off..
searchorder <- with(subset(datecheck_p2, grepl(b1l2trt[2], fulltrt) & logger == "B1L2"), min(datorder[grepl("correct", qa_note)]))
searchrange <- (searchorder-10):(searchorder +10)
subset(datecheck_p2,  grepl(b1l2trt[2], fulltrt) & datorder %in% searchrange) %>% #  & grepl("Apr20", filename)
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_vline(data = mutate(subset(adjustdates_p2, grepl(b1l2trt[2], fulltrt) & datorder %in% searchrange), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  geom_line(alpha = 0.4) +
  #geom_smooth(aes(fill = portid)) +
  ggtitle(paste("troubleshoot B1L2", b1l2trt[2], "date jump (dotted vert line = break)"))
# looks pretty spot on

# what are the missing dates in the clean master?
missingdates_b1l2<- unique(with(soilmoisture_clean_p2, date_time[is.na(vwc) & grepl("B1L2", portid)]))
missingdates_b1l2
# what are the date times around the off times in datecheck?
offtimes_b1l2 <- unique(with(datecheck_p2, date_time[logger == "B1L2" & datorder %in% searchrange]))
offtimes_b1l2 # jumps from 4pm on 6/26 to 2am on 6/27

# what does the clean order look like if timestamps left as is?
subset(soilmoisture_clean_p2, portid %in% unique(datecheck_p2$portid[datecheck_p2$fulltrt %in% b1l2trt])) %>%
  #subset to dates of interest
  subset(date_time %in% min(offtimes_b1l2):max(offtimes_b1l2)) %>%
  # create col to pull out logger to correct
  mutate(triage = grepl("B1L2", portid)) %>%
  ggplot(aes(date_time, vwc, col = portid, size = triage)) +
  #geom_vline(data = mutate(subset(adjustdates_p1, grepl(b3l4trt[2], fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  geom_point(alpha = 0.75) 
# seems like the datetimes should just be shifted back

# if everything is off by 10 hrs, how does the end of the sequence look?
subset(soilmoisture_clean_p2, portid %in% unique(datecheck_p2$portid[datecheck_p2$fulltrt %in% b1l2trt])) %>%
  #subset to dates of interest
  subset(date(date_time) > date(max(missingdates_b1l2))-5) %>%
  # create col to pull out logger to correct
  mutate(triage = grepl("B1L2", portid)) %>%
  ggplot(aes(date_time, vwc, col = portid, size = triage, shape = filename)) +
  #geom_vline(data = mutate(subset(adjustdates_p1, grepl(b3l4trt[2], fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  geom_point(alpha = 0.75) 
# datetimes for 16 Sep 2021 look fine, just need to adjust times in 15Sep21

# shift sep 15 file back 10 hrs starting at break. the datorder is fine, but the timestamp needs to be corrected
# think can just keep datorder and drop timestamp
b1l2test <- subset(datecheck_p2, logger == "B1L2") %>%
  # create cleanorder where 15Sep datorder is as is, 16Sep datorder is shifted to be the datorder corresponding to +10hrs (5 2-hr timesteps)
  left_join(distinct(subset(soilmoisture_clean_p2, select = c(date_time, cleanorder)))) %>%
  mutate(cleanorder2 = ifelse(grepl("15Sep", filename), datorder, cleanorder)) %>%
  rename(raw_datetime = date_time)

ggplot(b1l2test, aes(cleanorder2, vwc, col = port, group = portid)) +
  geom_point(data = subset(soilmoisture_clean_p2, grepl(str_flatten(b1l2trt, collapse = "|"), fulltrt) & logger != "B1L2"), aes(cleanorder, vwc, group = portid, col = port)) +
  geom_point() +
  geom_vline(data = mutate(subset(adjustdates_p2, grepl(str_flatten(b1l2trt, collapse = "|"), fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  facet_grid(logger~fulltrt)


# correct B3L4 in soilmoisture_clean_p1
clean_b1l2 <- subset(soilmoisture_clean_p2, grepl("B1L2", portid)) %>%
  select(date_time:cleanorder, filename) %>%
  left_join(select(b1l2test, portid, cleanorder2, vwc, raw_datetime, timeinterval, timediff, intevent, qa_note, datorder), by = c("portid", "cleanorder" = "cleanorder2")) %>%
  left_join(distinct(soilmoisture_p2, portid, logger, port, block, plotid, fulltrt))

# replot
subset(clean_b1l2, date(date_time) > date(max(missingdates_b1l2))-5) %>%
  ggplot(aes(date_time, vwc, fill = portid, shape = filename)) +
  # plot old points first
  geom_point(data = subset(datecheck_p2, logger == "B1L2" & date(date_time) > date(max(missingdates_b1l2))-5), aes(date_time, vwc, col = portid)) +
  geom_point(col = "black", pch = 21) # looks good

# append clean b1l2 to good loggers
soilmoisture_master_p2 <- rbind(soilmoisture_master_p2, clean_b1l2[names(soilmoisture_master_p2)])

# clean up
rm(b1l2trt, b1l2test, missingdates_b1l2, offtimes_b1l2, searchorder, searchrange, t, tempend, tempnote, tempstart, tempbreaks)



# 4.c. Triage B3L4 16Sep21 -----
# all ports same number obs (5543), get off track in timestamp (year) at same time (last run)
b3l4trt <- gsub("[0-9]", "", unique(datecheck_p2$fulltrt[datecheck_p2$logger == "B3L4"])) # only CW

# compare vwc across similar treatments by data order rather than date-time
## actual vwc values
subset(datecheck_p2, grepl(b3l4trt[1], fulltrt)) %>%
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_line(alpha = 0.4) +
  geom_vline(data = subset(adjustdates_p2, grepl(b3l4trt[1], fulltrt)), aes(xintercept = datorder, lty = qa_note)) +
  ggtitle(paste("troubleshoot B3L4", b3l4trt[1], "date jump (vert line = break)")) +
  geom_smooth(aes(fill = portid)) +
  facet_grid(logger~.) # looks like all lines up on dataorder, just datetime off

# color in each portid's line by year
subset(datecheck_p2, grepl(b3l4trt[1], fulltrt)) %>%
  ggplot(aes(datorder, vwc, group = portid)) +
  geom_vline(data = subset(adjustdates_p2, grepl(b3l4trt[1], fulltrt)), aes(xintercept = datorder, lty = qa_note)) +
  ggtitle(paste("troubleshoot B3L4", b3l4trt[1], "date jump (vert line = break)")) +
  geom_smooth(aes(fill = portid), col ="grey50", alpha = 0.5) +
  geom_line(aes(col = factor(year(date)))) +
  facet_grid(logger~.) # easy fix


# keep datorder, assign date_time as raw_datetime and join the master
b3l4test <- subset(datecheck_p2, logger == "B3L4") %>%
  # create cleanorder where 15Sep datorder is as is, 16Sep datorder is shifted to be the datorder corresponding to +10hrs (5 2-hr timesteps)
  left_join(distinct(subset(soilmoisture_clean_p2, select = c(date_time, cleanorder))))
# double check clean order diffs by one for all with dat order
summary(b3l4test$datorder - b3l4test$cleanorder) # yep
b3l4test <- mutate(b3l4test, cleanorder2 = datorder +1) %>%
  rename(raw_datetime = date_time)

# plot with good refdat to be sure adjustment looks good
ggplot(b3l4test, aes(cleanorder2, vwc, group = portid)) +
  geom_vline(data = subset(adjustdates_p2, fulltrt %in% b3l4trt), aes(xintercept = datorder, lty = qa_note)) +
  # add in refdat
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l4trt), aes(cleanorder, vwc, group = portid), alpha = 0.5) +
  geom_smooth(aes(col = portid, fill = portid), alpha = 0.5) +
  geom_line(aes(col = portid), alpha = 0.5) + 
  facet_grid(fulltrt~.) # looks fine
# zoom in
ggplot(subset(b3l4test, cleanorder2 %in% 750:850), aes(cleanorder2, scale(vwc), group = portid)) +
  geom_vline(data = subset(adjustdates_p2, fulltrt %in% b3l4trt), aes(xintercept = datorder, lty = qa_note)) +
  # add in refdat
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l4trt & cleanorder %in% 750:850), aes(cleanorder, scale(vwc), group = portid), alpha = 0.5) +
  #geom_smooth(aes(col = portid, fill = portid), alpha = 0.5) +
  geom_line(aes(col = portid), alpha = 0.5) + 
  facet_grid(fulltrt~.) # looks a little off after adjustment, look at a perdio with wetup events
# zoom in wetup
ggplot(subset(b3l4test, cleanorder2 %in% 2500:3550), aes(cleanorder2+1, vwc, group = portid)) +
  #geom_vline(data = subset(adjustdates_p2, fulltrt %in% b3l4trt), aes(xintercept = datorder, lty = qa_note)) +
  # add in refdat
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l4trt & cleanorder %in% 2500:3550), aes(cleanorder, vwc, group = portid), alpha = 0.5) +
  #geom_smooth(aes(col = portid, fill = portid), alpha = 0.5) +
  geom_line(aes(col = portid), alpha = 0.5) + 
  facet_grid(fulltrt~.) # adjusting by one more after break matches better to ref sources as it did pre break
# there are only 2 events, first one is fine, second needs an extra timestep
b3l4test$cleanorder3 <- with(b3l4test, ifelse(intevent > 1, cleanorder2+1, cleanorder2))

# correct B3L4 in soilmoisture_clean_p2
clean_b3l4 <- subset(soilmoisture_clean_p2, grepl("B3L4", portid)) %>%
  select(date_time:cleanorder, filename) %>%
  left_join(select(b3l4test, portid, cleanorder3, vwc, raw_datetime, timeinterval, timediff, intevent, qa_note, datorder), by = c("portid", "cleanorder" = "cleanorder3")) %>%
  left_join(distinct(soilmoisture_p2, portid, logger, port, block, plotid, fulltrt))
# check NAs in identifying info
summary(is.na(clean_b3l4))
View(subset(clean_b3l4, is.na(intevent))) # at first and at break, fine

# append
soilmoisture_master_p2 <- subset(soilmoisture_master_p2, !grepl("B3L4", portid)) %>%
  rbind(clean_b3l4[names(.)])

# plot to review
ggplot(soilmoisture_master_p2, aes(date(date_time), vwc, group = portid, col = logger, lty = logger %in% goodrefs$logger)) +
  geom_line(alpha = 0.5) +
  facet_grid(substr(fulltrt, 1,1)~gsub("^(C|F|N)", "", fulltrt)) 

# check NAs per port
with(soilmoisture_master_p2, sapply(split(vwc, portid), summary))
# b3l4_3 has a lot of missing.. check against raw (thought they had same nobs, but maybe I included NAs in those counts -- if row present, counted)
with(subset(soilmoisture_master_p2, grepl("B3L4", logger) & is.na(vwc)), sapply(split(vwc, portid), length))
with(subset(soilmoisture_p2, grepl("B3L4", logger) & is.na(vwc)), sapply(split(vwc, portid), length)) # checks out

# clean up environment
rm(b3l4trt, b3l4test)



# 4.d. Triage B2L5 15Sep21 -----
# data end in march 2021?
b2l5trt <- gsub("[0-9]", "", unique(datecheck_p2$fulltrt[datecheck_p2$logger == "B2L5"])) # only CW

# compare vwc across similar treatments by data order rather than date-time
## actual vwc values
subset(datecheck_p2, grepl(b2l5trt[1], fulltrt)) %>%
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_line(alpha = 0.4) +
  geom_vline(data = subset(adjustdates_p2, grepl(b2l5trt[1], fulltrt)), aes(xintercept = datorder, lty = qa_note)) +
  ggtitle(paste("troubleshoot B2L5", b2l5trt[1], "date jump (dotted vert line = break)")) +
  geom_smooth(aes(fill = portid)) +
  facet_grid(logger~.) # eyeballing, looks like it lines up around 4200 mark

#  zoom in
subset(datecheck_p2, grepl(b2l5trt[1], fulltrt) & datorder > 3400) %>%
  ggplot(aes(datorder, vwc, group = portid)) +
  geom_vline(data = mutate(subset(adjustdates_p2, grepl(b2l5trt[1], fulltrt)), qa_note2 = ifelse(grepl("correct", qa_note), qa_note, NA)), aes(xintercept = datorder, lty = qa_note2, col = factor(timediff))) +
  geom_line(alpha = 0.4) +
  ggtitle(paste("troubleshoot B2L5", b2l5trt[1], "date jump (dotted vert line = break)")) +
  #geom_smooth(aes(fill = portid)) +
  facet_grid(logger~.) # eyeballing, looks like it lines up around 4200 mark

# do all ports have the same issue?
datecheck_p2_b2l5 <- subset(datecheck_p2, logger == "B2L5" & !is.na(qa_note), select = c(logger, port, date_time, date, timediff, intevent, qa_note)) %>%
  unique()
View(datecheck_p2_b2l5) # same up thru int event 5. port 5 has more NAs after that.

# try plotting just 2021 dates to see if jump to may 2021 at intevent 2 matches up
subset(datecheck_p2, grepl(b2l5trt[1], fulltrt) & year(date)==2021) %>%
  ggplot(aes(date_time, vwc, group = portid, col = portid)) +
  geom_line(alpha = 0.4) # maybe? can try adjusting data from 2021 date_time in intevent 2; subsequent data seem continuous, just timestamps off
# b2l5 will need some NA added for additional breaks (perhaps?)

b2l5test_thru5 <- subset(datecheck_p2, logger == "B2L5" & intevent <= 5) %>%
  # join clean order
  left_join(distinct(subset(soilmoisture_clean_p2, select = c(date_time, cleanorder))))
# how off is clean order from the b2l5 datorder (ie., in off period and does it start at same time)
group_by(b2l5test_thru5, intevent) %>%
  subset(!is.na(cleanorder)) %>%
  summarise(diff = datorder - cleanorder) %>%
  distinct() # interval 1 matches, interval 2 behind by 789 timesteps

# adjust clean order accordingly and proceed
b2l5test_thru5 <- subset(b2l5test_thru5, select = -cleanorder) %>%
  mutate(cleanorder2 = ifelse(intevent == 1, datorder, datorder + 789)) %>%
  rename(raw_datetime = date_time) %>%
  # rejoin clean datorder and datetime
  left_join(distinct(subset(soilmoisture_clean_p2, select = c(date_time, cleanorder))), by = c("cleanorder2" = "cleanorder"))

# plot to see how it looks
subset(datecheck_p2, grepl(b2l5trt[1], fulltrt) & year(date)==2021 & logger != "B2L5") %>%
  ggplot(aes(date_time, vwc, group = portid, col = portid)) +
  geom_line(data = b2l5test_thru5, alpha = 0.5) +
  geom_line(alpha = 0.5)

int2date <- unique(datecheck_p2_b2l5$date[datecheck_p2_b2l5$intevent == 2])
# zoom in
subset(datecheck_p2, grepl(b2l5trt[1], fulltrt) & year(date)==2021 & logger != "B2L5" & date >= int2date) %>%
  ggplot(aes(date_time, vwc, group = portid)) +
  geom_line(data = subset(b2l5test_thru5, date(date_time) >= int2date), aes(col = factor(intevent)), alpha = 0.5) +
  geom_line(col = "grey30", alpha = 0.5) # not quite right

subset(datecheck_p2, grepl(b2l5trt[1], fulltrt) & year(date)==2021 & logger != "B2L5" & date >= int2date) %>%
  ggplot(aes(datorder, scale(vwc), group = portid)) +
  geom_line(data = subset(b2l5test_thru5, date(date_time) >= int2date), aes(cleanorder2, scale(vwc), col = factor(intevent)), alpha = 0.5) +
  geom_line(col = "grey30", alpha = 0.5) +
  scale_x_continuous(breaks = seq(4000, 5600, 100)) +
  facet_wrap(~portid, nrow = 3, scale = "free_y")
# match events 2-4 on b2l5_4 (can see visual matches best)
# match event 5 on b2l5_1 (longer record)

b2l5_subtest <- subset(b2l5test_thru5, date(date_time) >= int2date) %>%
  mutate(cleanorder3 = cleanorder2)
# try nudging
b2l5_subtest$cleanorder3[b2l5_subtest$intevent == 2] <- c(b2l5_subtest$cleanorder2[b2l5_subtest$intevent == 2] + 10) #+10 looks good
b2l5_subtest$cleanorder3[b2l5_subtest$intevent == 3] <- c(b2l5_subtest$cleanorder2[b2l5_subtest$intevent == 3] + 36) #+36 looks good
b2l5_subtest$cleanorder3[b2l5_subtest$intevent == 4] <- c(b2l5_subtest$cleanorder2[b2l5_subtest$intevent == 4] + 53) #+53 looks good
b2l5_subtest$cleanorder3[b2l5_subtest$intevent == 5] <- c(b2l5_subtest$cleanorder2[b2l5_subtest$intevent == 5] + 215) #+215 looks good

ggplot(subset(b2l5_subtest, intevent %in% c(2:5) & cleanorder3 < 4400), aes(cleanorder3, scale(vwc))) + #portid == "B2L5_2"
  # reference
  geom_point(data = subset(soilmoisture_master_p2, fulltrt %in% b2l5trt & date(date_time) >= int2date & cleanorder < 4400), #& datorder < 4800 
            aes(cleanorder, scale(vwc), group = portid), col = "grey30", alpha = 0.85) +
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b2l5trt & date(date_time) >= int2date & cleanorder < 4400), #& datorder < 4800 
             aes(cleanorder, scale(vwc), group = portid), col = "grey30", alpha = 0.85) +
  geom_point(aes(col = factor(intevent), group = paste0(portid,intevent))) +
  geom_line(aes(col = factor(intevent), group = paste0(portid,intevent))) +
  scale_x_continuous(breaks = seq(4000, 5600, 100))
  #facet_grid(portid~., scale = "free_y")

# how did they line up before break?
ggplot(subset(datecheck_p2, logger == "B2L5" & intevent == 1 & datorder %in% 750:1000), aes(datorder, scale(vwc))) +
  # reference
  geom_point(data = subset(soilmoisture_master_p2, fulltrt %in% b2l5trt & cleanorder %in% 750:1000), #& datorder < 4800 
             aes(cleanorder, scale(vwc), group = portid), col = "grey30", alpha = 0.85) +
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b2l5trt & cleanorder %in% 750:1000), #& datorder < 4800 
            aes(cleanorder, scale(vwc), group = portid), col = "grey30", alpha = 0.85) +
  geom_point(aes(col = factor(portid), group = paste0(portid,intevent))) +
  geom_line(aes(col = factor(portid), group = paste0(portid,intevent))) # lines up with b2l5_2 best

# create different dataframe for ports 4 and 5 for intevents 5+
b2l5test_5up <- subset(datecheck_p2, logger == "B2L5" & port >= 4 & intevent >= 5) %>%
  mutate(cleanorder2 = ifelse(intevent == 1, datorder, datorder + 789)) %>%
  rename(raw_datetime = date_time) %>%
  # rejoin clean datorder and datetime
  left_join(distinct(subset(soilmoisture_clean_p2, select = c(date_time, cleanorder))), by = c("cleanorder2" = "cleanorder"))


# subset port 4 for intervals 6+ (get off at different timestamps)
b2l5_4 <- subset(b2l5test_5up, port == 4 & intevent >= 5) %>%
  mutate(cleanorder3 = cleanorder2)
# try nudging
b2l5_4$cleanorder3[b2l5_4$intevent == 5] <- c(b2l5_4$cleanorder2[b2l5_4$intevent == 5] + 215) # start with adjustment from above
b2l5_4$cleanorder3[b2l5_4$intevent == 6] <- c(b2l5_4$cleanorder2[b2l5_4$intevent == 6] + 216) # timediff = 4 from int 5 (1 extra timestep)
b2l5_4$cleanorder3[b2l5_4$intevent == 7] <- c(b2l5_4$cleanorder2[b2l5_4$intevent == 7] + 217) # timediff = 4 from int 6 (so 2 extra timesteps from 5)
b2l5_4$cleanorder3[b2l5_4$intevent == 8] <- c(b2l5_4$cleanorder2[b2l5_4$intevent == 8] + 218) # 8 is last run (different file)

ggplot(b2l5_4, aes(cleanorder3, vwc)) +
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b2l5trt & cleanorder > 4700), #
            aes(cleanorder, vwc+0.1), col = "grey70", alpha = 0.5, lwd = 1) +
  geom_line(col="grey40") +
  geom_point(aes(col = factor(intevent))) +
  scale_x_continuous(breaks = seq(4000, 5600, 100)) # looks good


# subset port 5 intervals 6+ (get off at different timestamps)
b2l5_5 <- subset(b2l5test_5up, port == 5 & intevent >= 5) %>%
  mutate(cleanorder3 = cleanorder2 + 215)
# try nudging
b2l5_5$cleanorder3[b2l5_5$intevent == 6] <- c(b2l5_5$cleanorder2[b2l5_5$intevent == 6] + 216) # through visual trial and error, diff between each subsequent really is 1 extra time step
b2l5_5$cleanorder3[b2l5_5$intevent == 7] <- c(b2l5_5$cleanorder2[b2l5_5$intevent == 7] + 217) 
b2l5_5$cleanorder3[b2l5_5$intevent == 8] <- c(b2l5_5$cleanorder2[b2l5_5$intevent == 8] + 218) 
b2l5_5$cleanorder3[b2l5_5$intevent == 9] <- c(b2l5_5$cleanorder2[b2l5_5$intevent == 9] + 219)
b2l5_5$cleanorder3[b2l5_5$intevent == 10] <- c(b2l5_5$cleanorder2[b2l5_5$intevent  == 10] + 220)
b2l5_5$cleanorder3[b2l5_5$intevent == 11] <- c(b2l5_5$cleanorder2[b2l5_5$intevent == 11] + 221)
b2l5_5$cleanorder3[b2l5_5$intevent == 12] <- c(b2l5_5$cleanorder2[b2l5_5$intevent == 12] + 222)
b2l5_5$cleanorder3[b2l5_5$intevent == 13] <- c(b2l5_5$cleanorder2[b2l5_5$intevent  == 13] + 223)
b2l5_5$cleanorder3[b2l5_5$intevent == 14] <- c(b2l5_5$cleanorder2[b2l5_5$intevent  == 14] + 224)
b2l5_5$cleanorder3[b2l5_5$intevent >= 15] <- c(b2l5_5$cleanorder2[b2l5_5$intevent >= 15] + 225) # 16 is last run (different file)

ggplot(b2l5_5, aes(cleanorder3, scale(vwc))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b2l5trt & cleanorder > 4700), #
            aes(cleanorder, scale(vwc)), col = "grey70", alpha = 0.5, lwd = 1) +
  geom_line(col="grey40") +
  geom_point(aes(col = factor(intevent))) +
  scale_x_continuous(breaks = seq(4000, 5600, 100)) # looks good


# stack all adjusted dats
b2l5test <- subset(b2l5test_thru5, date(date_time) < int2date) %>% # just everything up thru date of interval 2
  # add a cleanorder 3 col at the end to match other dats
  mutate(cleanorder3 = cleanorder2) %>%
  rbind(b2l5_subtest) %>% # intevents 2-5
  rbind(subset(b2l5_4, intevent > 5)) %>% # for port 4 to end, drop intevent 5
  rbind(subset(b2l5_5, intevent > 5)) %>% # for port 5 to end, drop intevent 5
  distinct()

# check to see if anything missing
with(subset(datecheck_p2, logger == "B2L5"), sapply(split(intevent, portid), unique))
with(subset(datecheck_p2, logger == "B2L5"), sapply(split(intevent, portid), length))

with(b2l5test, sapply(split(intevent, portid), unique)) # intevent 6 missing for ports 1-3
with(b2l5test, sapply(split(intevent, portid), length)) # 10 nobs missing for port 1-3

# adjust intevent 6 for ports 1-3
b2l5_int6 <- subset(datecheck_p2, logger == "B2L5" & port < 4 & intevent >= 5) %>%
  mutate(cleanorder2 = ifelse(intevent == 1, datorder, datorder + 789)) %>%
  rename(raw_datetime = date_time) %>%
  # rejoin clean datorder and datetime
  left_join(distinct(subset(soilmoisture_clean_p2, select = c(date_time, cleanorder))), by = c("cleanorder2" = "cleanorder")) %>%
  # try adjusting int 6 the same as int 5
  mutate(cleanorder3 = cleanorder2 + 215)

# viz to check
ggplot(subset(b2l5_int6, cleanorder3 > 5300), aes(cleanorder3, vwc, group = portid)) +
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b2l5trt & cleanorder > 5300), #
            aes(cleanorder, vwc+0.07), col = "grey70", alpha = 0.5) +
  geom_line(col="grey40") +
  geom_point(aes(col = factor(intevent))) +
  scale_x_continuous(breaks = seq(5400, 5600, 50)) # looks okay
# check that nobs per port are what's expected
with(subset(b2l5_int6, intevent > 5), sapply(split(vwc, portid), length)) # yes
# append
b2l5test <- rbind(b2l5test, subset(b2l5_int6, intevent > 5)) %>%
  distinct() %>%
  arrange(date_time, portid)
# final check on intevents and length
with(b2l5test, sapply(split(intevent, portid), length))
with(b2l5test, sapply(split(intevent, portid), unique)) # looks good

# correct B2L5 in soilmoisture_clean_p2
clean_b2l5 <- subset(soilmoisture_clean_p2, grepl("B2L5", portid)) %>%
  select(date_time:cleanorder, filename) %>%
  left_join(select(b2l5test, portid, cleanorder3, vwc, raw_datetime, timeinterval, timediff, intevent, qa_note, datorder), by = c("portid", "cleanorder" = "cleanorder3")) %>%
  left_join(distinct(soilmoisture_p2, portid, logger, port, block, plotid, fulltrt))
# check again
with(clean_b2l5, sapply(split(intevent, portid), length))
with(clean_b2l5, sapply(split(intevent, portid), unique)) # NAs present for intevent bc NAs are now inserted for blips in data recording
summary(is.na(clean_b2l5))

# append
soilmoisture_master_p2 <- subset(soilmoisture_master_p2, !grepl("B2L5", portid)) %>%
  rbind(clean_b2l5[names(.)])

# review
with(soilmoisture_master_p2, sapply(split(vwc, portid), length))
soilmoisture_master_p2 %>%
  group_by(logger) %>%
  summarise(datrange = str_flatten(as.character(range(date_time)), collapse =  " "))
# some start at 6am and others 8am on 6/11/2020. that's fine. depends on when their period 1 stopped.
# b3l2 stops early, but in june so it's fine. can add on NA vals for it at the end

# clean up environment
rm(list = ls()[grep("^b2l5", ls())])


# 4.e. Triage B3L3 16Sep21 -----
# ports 1, 3, 5 have different number obs (5222,4815,5215 respec.)
# ports 2, 4 have same number (5210)
b3l3trt <- gsub("[0-9]", "", unique(datecheck_p2$fulltrt[datecheck_p2$logger == "B3L3"])) # NW, NXC

# compare vwc across similar treatments by data order rather than date-time
## actual vwc values
subset(datecheck_p2, grepl(b3l3trt[1], fulltrt)) %>%
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_line(alpha = 0.4) +
  geom_vline(data = subset(adjustdates_p2, grepl(b3l3trt[1], fulltrt)), aes(xintercept = datorder, lty = qa_note)) +
  ggtitle(paste("troubleshoot B3L3", b3l3trt[1], "date jump (dotted vert line = break)")) +
  geom_smooth(aes(fill = portid)) +
  facet_grid(logger~.) # seems to match pretty well up until first NA infill needed

# check second treatment
subset(datecheck_p2, grepl(b3l3trt[2], fulltrt)) %>%
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_line(alpha = 0.4) +
  geom_vline(data = subset(adjustdates_p2, grepl(b3l3trt[2], fulltrt)), aes(xintercept = datorder, lty = qa_note)) +
  ggtitle(paste("troubleshoot B3L3", b3l3trt[2], "date jump (dotted vert line = break)")) +
  geom_smooth(aes(fill = portid)) +
  facet_grid(logger~.) # can compare against b2l3 (clean record)
# stinker is no wetup events in off period to match against.. so will be manual again
# gets off after datorder 4000

# for which interval events are ports off in a similar way?
datecheck_p2_b3l3 <- subset(datecheck_p2, logger == "B3L3" & !is.na(qa_note), select = c(logger:port, date_time, date, timediff, intevent, qa_note))
# all get off track similarly for intevent 2, but start to vary after that

# visualize b3l3 data by logger port to see where they get off track
subset(datecheck_p2, logger == "B3L3") %>%
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_line(alpha = 0.4) +
  geom_vline(data = subset(adjustdates_p2, logger == "B3L3"), aes(xintercept = datorder, lty = qa_note)) +
  ggtitle(paste("troubleshoot B3L3", b3l3trt[2], "date jump (dotted vert line = break)")) +
  geom_smooth(aes(fill = portid)) +
  facet_grid(portid~.)
# based on datorder rather than date_time, 5 and 1 seem similar, 2 and 4 seem similar
# even tho port 3 has no data, timestamps generated by logger are still wrong. since all is fine beforehand, could just use cleanorder with NA vwc vals

# treat 1 and 5
# zoom in
subset(datecheck_p2, portid %in% c("B3L3_1", "B3L3_5", "B2L3_4", "B2L3_5", "B1L2") & datorder > 4000) %>%
  ggplot(aes(datorder, vwc, group = portid, col = portid)) + # col = factor(port)
  geom_vline(data = subset(adjustdates_p2, portid %in% c("B3L3_1", "B3L3_5") & datorder > 4000), aes(xintercept = datorder, lty = qa_note, col = portid)) +
  geom_line(alpha = 0.4) +
  ggtitle(paste("troubleshoot B3L3", b3l3trt[1], "date jump (dotted vert line = break)")) +
  #geom_smooth(aes(fill = portid)) +
  facet_grid(logger~.) 
# so they're not quite the same after all.. but adjustments done with probably be similar
# maybe easier to work backwards in adjustments

# 4.e.1. adjust B3L3 port 1 -----
b3l3test_1 <- subset(datecheck_p2, portid == "B3L3_1") %>%
  # join clean order
  left_join(distinct(subset(soilmoisture_clean_p2, dateselect = c(date_time, cleanorder)))) %>%
  mutate(cleanorder2 = ifelse(intevent == 1, cleanorder, datorder)) # cleanorder diffs from datorder by 1, only attached to intevent 1

# how many unique events, and what type?
with(subset(b3l3test_1, !is.na(qa_note)), sapply(split(intevent, qa_note), unique))

# try adjustments
# > adjust port 1 first, then adjust others to port 1, then shift all together
b3l3test_1$cleanorder2[b3l3test_1$intevent == 2] <- b3l3test_1$datorder[b3l3test_1$intevent == 2] + 27
b3l3test_1$cleanorder2[b3l3test_1$intevent == 3] <- b3l3test_1$datorder[b3l3test_1$intevent == 3] + 46
b3l3test_1$cleanorder2[b3l3test_1$intevent == 4] <- b3l3test_1$datorder[b3l3test_1$intevent == 4] + 50
b3l3test_1$cleanorder2[b3l3test_1$intevent == 5] <- b3l3test_1$datorder[b3l3test_1$intevent == 5] + 66
b3l3test_1$cleanorder2[b3l3test_1$intevent == 6] <- b3l3test_1$datorder[b3l3test_1$intevent == 6] + 73
#b3l3test_1$cleanorder2[b3l3test_1$intevent == 7 & b3l3test_1$rowid <= 650067] <- b3l3test_1$datorder[b3l3test_1$intevent == 7 & b3l3test_1$rowid <= 650067] + 81
#b3l3test_1$cleanorder2[b3l3test_1$intevent == 7 & b3l3test_1$rowid %in% 650068:650076] <- b3l3test_1$datorder[b3l3test_1$intevent == 7 & b3l3test_1$rowid %in% 650068:650076] + 83
#b3l3test_1$cleanorder2[b3l3test_1$intevent == 7 & b3l3test_1$rowid > 650076] <- b3l3test_1$datorder[b3l3test_1$intevent == 7 & b3l3test_1$rowid > 650076] + 87 # 7 is a needs NA infill by 1 timestep
b3l3test_1$cleanorder2[b3l3test_1$intevent == 7] <- b3l3test_1$datorder[b3l3test_1$intevent == 7] + 87 # 7 is a needs NA infill by 1 timestep
b3l3test_1$cleanorder2[b3l3test_1$intevent == 8] <- b3l3test_1$datorder[b3l3test_1$intevent == 8] + 103
b3l3test_1$cleanorder2[b3l3test_1$intevent == 9] <- b3l3test_1$datorder[b3l3test_1$intevent == 9] + 106
b3l3test_1$cleanorder2[b3l3test_1$intevent == 10] <- b3l3test_1$datorder[b3l3test_1$intevent == 10] + 120
b3l3test_1$cleanorder2[b3l3test_1$intevent == 11] <- b3l3test_1$datorder[b3l3test_1$intevent == 11] + 123
b3l3test_1$cleanorder2[b3l3test_1$intevent == 12] <- b3l3test_1$datorder[b3l3test_1$intevent == 12] + 127
b3l3test_1$cleanorder2[b3l3test_1$intevent == 13] <- b3l3test_1$datorder[b3l3test_1$intevent == 13] + 130
b3l3test_1$cleanorder2[b3l3test_1$intevent == 14] <- b3l3test_1$datorder[b3l3test_1$intevent == 14] + 133
b3l3test_1$cleanorder2[b3l3test_1$intevent == 15] <- b3l3test_1$datorder[b3l3test_1$intevent == 15] + 136
b3l3test_1$cleanorder2[b3l3test_1$intevent == 16] <- b3l3test_1$datorder[b3l3test_1$intevent == 16] + 139
b3l3test_1$cleanorder2[b3l3test_1$intevent == 17] <- b3l3test_1$datorder[b3l3test_1$intevent == 17] + 141

# just the periods to adjust -- zoom in to datorder of interest as needed to review (i.e., this is manual)
ggplot(subset(b3l3test_1, cleanorder2 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, scale(vwc), col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3trt & cleanorder > 4000), # & cleanorder %in% 4000:4350
            aes(cleanorder, scale(vwc), group = portid), col = "grey50", alpha = 0.4) +
  # plot old as faint line
  #geom_line(data = subset(b3l3test_1, datorder > 4700), aes(datorder, vwc, col = factor(intevent)), alpha = 0.4) +
  # plot adjusted
  #geom_text(alpha = 0.6, aes (label = rowid)) + # use cleanorder %in% 4200:4260 to check adjustments for interval 7
  geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  geom_line(aes(lty = intevent %% 2 == 0), alpha = 0.6)
  #scale_x_continuous(breaks = seq(4000, 4600, 25))

# all dats to see typical relationship
ggplot(subset(b3l3test_1), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, vwc, col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3trt), # & cleanorder %in% 4000:4350
            aes(cleanorder, vwc-0.02, group = portid), col = "grey50", alpha = 0.4) +
  # plot old as faint line
  #geom_line(data = subset(b3l3test_1, datorder > 4700), aes(datorder, vwc, col = factor(intevent)), alpha = 0.4) +
  # plot adjusted
  geom_line(alpha = 0.6)

# create function to try to assign NA infill to avoid doing it manually
# > input should be events that are NA infill
assign_NA <- function(dat, assigncol = "cleanorder3", refcol = "cleanorder2"){
  # ID interval events that are needs NA infill
  NAevents <- with(dat, intevent[grepl("NA infill", qa_note)])
  timestamps_off <- with(dat, intevent[grepl("correct", qa_note)])
  # create column to assign new dat order numbers
  dat[[assigncol]] <- dat[[refcol]]
  for(n in NAevents){
    # what's the timediff from the previous event?
    tempdiff <- abs(with(subset(dat, intevent == n & !is.na(qa_note)), timediff)) #[rowid == min(rowid)]
    # check there's a single value
    stopifnot(length(tempdiff) == 1 & !is.na(tempdiff))
    # where did the count leave off with the last non-NA event?
    lastdatorder <- tail(dat[[assigncol]][dat$intevent == c(n-1)], n = 1)
    # how many rows in NA event?
    n_NA <- nrow(subset(dat, intevent == n))
    # seq to assign
    tempseq <- seq(lastdatorder + (tempdiff/2),length.out = n_NA, by = 1)
    dat[[assigncol]][dat$intevent == n] <- tempseq
  }
  # when done, diff to check
  dat$diffcheck <- dat[[assigncol]] - dat[[refcol]]
  # return
  return(dat)
}

# test it out
copy <- b3l3test_1
b3l3test_1 <- assign_NA(b3l3test_1)

# plot to check
ggplot(subset(b3l3test_1, cleanorder3 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder3, scale(vwc), col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3trt & cleanorder > 4000), # & cleanorder %in% 4000:4350
            aes(cleanorder, scale(vwc), group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  #geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  geom_line(aes(lty = intevent %% 2 == 0), alpha = 0.6)
# worked as intended
# proceed with next


# 4.e.2. adjust B3L3 port 5 -----
# adjust 5 next since similar to port 1
b3l3test_5 <- subset(datecheck_p2, portid == "B3L3_5") %>%
  # join clean order
  left_join(distinct(subset(soilmoisture_clean_p2, dateselect = c(date_time, cleanorder)))) %>%
  mutate(cleanorder2 = ifelse(intevent == 1, cleanorder, datorder)) # cleanorder diffs from datorder by 0, only attached to intevent 1

# how many unique events, and what type?
with(subset(b3l3test_5, !is.na(qa_note)), sapply(split(intevent, qa_note), unique))

# ID the events that need timestamp corrections + last run
b3l3_5_timestamps <- b3l3test_5$intevent[grepl("timestamp|last", b3l3test_5$qa_note)]

# adjust
b3l3test_5$cleanorder2[b3l3test_5$intevent == 2] <- b3l3test_5$datorder[b3l3test_5$intevent == 2] + 27
b3l3test_5$cleanorder2[b3l3test_5$intevent == 3] <- b3l3test_5$datorder[b3l3test_5$intevent == 3] + 45
#b3l3test_5$cleanorder2[b3l3test_5$intevent == 4] <- b3l3test_5$datorder[b3l3test_5$intevent == 4] + 48
b3l3test_5$cleanorder2[b3l3test_5$intevent == 5] <- b3l3test_5$datorder[b3l3test_5$intevent == 5] + 53
#b3l3test_5$cleanorder2[b3l3test_5$intevent == 6] <- b3l3test_5$datorder[b3l3test_5$intevent == 6] + 73
b3l3test_5$cleanorder2[b3l3test_5$intevent == 7] <- b3l3test_5$datorder[b3l3test_5$intevent == 7] + 71
#b3l3test_5$cleanorder2[b3l3test_5$intevent == 8] <- b3l3test_5$datorder[b3l3test_5$intevent == 8] + 86
b3l3test_5$cleanorder2[b3l3test_5$intevent == 9] <- b3l3test_5$datorder[b3l3test_5$intevent == 9] + 81
#b3l3test_5$cleanorder2[b3l3test_5$intevent == 10] <- b3l3test_5$datorder[b3l3test_5$intevent == 10] + 94
b3l3test_5$cleanorder2[b3l3test_5$intevent == 11] <- b3l3test_5$datorder[b3l3test_5$intevent == 11] + 109
b3l3test_5$cleanorder2[b3l3test_5$intevent == 12] <- b3l3test_5$datorder[b3l3test_5$intevent == 12] + 113
b3l3test_5$cleanorder2[b3l3test_5$intevent == 13] <- b3l3test_5$datorder[b3l3test_5$intevent == 13] + 127
b3l3test_5$cleanorder2[b3l3test_5$intevent == 14] <- b3l3test_5$datorder[b3l3test_5$intevent == 14] + 130
b3l3test_5$cleanorder2[b3l3test_5$intevent == 15] <- b3l3test_5$datorder[b3l3test_5$intevent == 15] + 134
b3l3test_5$cleanorder2[b3l3test_5$intevent == 16] <- b3l3test_5$datorder[b3l3test_5$intevent == 16] + 137
b3l3test_5$cleanorder2[b3l3test_5$intevent == 17] <- b3l3test_5$datorder[b3l3test_5$intevent == 17] + 140
b3l3test_5$cleanorder2[b3l3test_5$intevent == 18] <- b3l3test_5$datorder[b3l3test_5$intevent == 18] + 143
b3l3test_5$cleanorder2[b3l3test_5$intevent == 19] <- b3l3test_5$datorder[b3l3test_5$intevent == 19] + 146
b3l3test_5$cleanorder2[b3l3test_5$intevent == 20] <- b3l3test_5$datorder[b3l3test_5$intevent == 20] + 147

# all dats to see typical relationship
ggplot(subset(b3l3test_5, intevent %in% b3l3_5_timestamps & cleanorder2 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350 #%in% 4000:4325
       aes(cleanorder2, vwc, col = factor(intevent), group = intevent)) + #group = factor(intevent)
  # reference
  # geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3trt & cleanorder > 4000), # & cleanorder %in% 4000:4350 #%in% 4000:4325
  #           aes(cleanorder, vwc, group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  geom_line(aes(lty = intevent %% 2 == 0)) + # very similar to port 1
  # add port 1 to help guide
  geom_line(data = subset(b3l3test_1, cleanorder3 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350
            aes(cleanorder2, vwc+0.03, group = intevent), col = "grey50", alpha = 0.6) # seems ok

# run NA infill
b3l3test_5 <- assign_NA(b3l3test_5)

# replot to check
ggplot(subset(b3l3test_5, cleanorder3 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350 #%in% 4000:4325
       aes(cleanorder3, vwc, col = factor(intevent), group = intevent)) + #group = factor(intevent)
  # reference
  # geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3trt & cleanorder > 4000), # & cleanorder %in% 4000:4350 #%in% 4000:4325
  #           aes(cleanorder, vwc, group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  geom_line(aes(lty = intevent %in% b3l3_5_timestamps)) + # very similar to port 1
  # add port 1 to help guide
  geom_line(data = subset(b3l3test_1, cleanorder3 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350
            aes(cleanorder2, vwc, group = intevent), col = "grey50", alpha = 0.6)
# circle back to event 7


# 4.e.2. adjust B3L3 port 2 -----
b3l3test_2 <- subset(datecheck_p2, portid == "B3L3_2") %>%
  # join clean order
  left_join(distinct(subset(soilmoisture_clean_p2, dateselect = c(date_time, cleanorder)))) %>%
  mutate(cleanorder2 = ifelse(intevent == 1, cleanorder, datorder)) # cleanorder diffs from datorder by 1, only attached to intevent 1

# how many unique events, and what type?
with(subset(b3l3test_2, !is.na(qa_note)), sapply(split(intevent, qa_note), unique))

# adjust
b3l3test_2$cleanorder2[b3l3test_2$intevent == 26] <- b3l3test_2$datorder[b3l3test_2$intevent == 26] + 148
b3l3test_2$cleanorder2[b3l3test_2$intevent == 27] <- b3l3test_2$datorder[b3l3test_2$intevent == 27] + 150


ggplot(subset(b3l3test_2, cleanorder2 > 4500), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, vwc+0.03, col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3trt & cleanorder > 4500), # & cleanorder %in% 4000:4350
            aes(cleanorder, vwc, group = portid), col = "grey50", alpha = 0.4) +
  # plot old as faint line
  #geom_line(data = subset(b3l3test_1, datorder > 4700), aes(datorder, vwc, col = factor(intevent)), alpha = 0.4) +
  # plot adjusted
  geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  geom_line(aes(lty = intevent %% 2 == 0), alpha = 0.6) # very similar to port 1


# 4.e.4. adjust B3L3 port 4 -----
b3l3test_4 <- subset(datecheck_p2, portid == "B3L3_4") %>%
  # join clean order
  left_join(distinct(subset(soilmoisture_clean_p2, dateselect = c(date_time, cleanorder)))) %>%
  mutate(cleanorder2 = ifelse(intevent == 1, cleanorder, datorder)) # cleanorder diffs from datorder by 0, only attached to intevent 1

# how many unique events, and what type?
with(subset(b3l3test_4, !is.na(qa_note)), sapply(split(intevent, qa_note), unique))

ggplot(subset(b3l3test_4), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, vwc, col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3trt), # & cleanorder %in% 4000:4350
            aes(cleanorder, vwc, group = portid), col = "grey50", alpha = 0.4) +
  # plot old as faint line
  #geom_line(data = subset(b3l3test_1, datorder > 4700), aes(datorder, vwc, col = factor(intevent)), alpha = 0.4) +
  # plot adjusted
  geom_line(alpha = 0.6) # very similar to port 1


# 4.f. Triage B2L4 15Sep21 -----
# all ports have different number obs
# 1 = 5091, 2 = 5493, 3 = 5497, 4 = 5481, 5 = 5504
# port 1 doesn't have much data from roughly sep 2019-sep 2021
# all others missing data sometime spring 2021 - sep 2021
b2l4trt <- gsub("[0-9]", "", unique(datecheck_p2$fulltrt[datecheck_p2$logger == "B2L4"])) # CXC, CD

# compare vwc across similar treatments by data order rather than date-time
## actual vwc values
subset(datecheck_p2, grepl(b2l4trt[1], fulltrt)) %>%
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_line(alpha = 0.4) +
  geom_vline(data = subset(adjustdates_p2, grepl(b2l4trt[1], fulltrt)), aes(xintercept = datorder, lty = qa_note)) +
  ggtitle(paste("troubleshoot B2L4", b2l4trt[1], "date jump (dotted vert line = break)")) +
  geom_smooth(aes(fill = portid)) +
  facet_grid(logger~.) # seems to match pretty well up until first NA infill needed

subset(datecheck_p2, grepl(b2l4trt[2], fulltrt)) %>%
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_line(alpha = 0.4) +
  geom_vline(data = subset(adjustdates_p2, grepl(b2l4trt[1], fulltrt)), aes(xintercept = datorder, lty = qa_note)) + #b2l4trt[2]
  ggtitle(paste("troubleshoot B2L4", b2l4trt[2], "date jump (dotted vert line = break)")) +
  geom_smooth(aes(fill = portid)) +
  facet_grid(logger~.) # seems to match pretty well up until first NA infill needed



# -- STACK TREATED PERIODS 1 + 2 -----


# remake plots from initial compilation loop
logger_plot_clean <- soilmoisture_master_p1 %>%
  ggplot(aes(dowy, vwc, col = as.factor(port))) + #date_time, 
  geom_line(alpha = 0.5) +
  ggtitle(paste0("Compost prelim plot (", Sys.Date(),  "): cleaned soil moisture (VWC),\nby logger, by water year, colored by port")) +
  scale_color_discrete(name = "port") +
  facet_grid(logger~waterYear, scales = "free_x")
logger_plot_clean

trt_plot_clean <- soilmoisture_master_p1 %>%
  ggplot(aes(dowy, vwc, group = portid, col = as.factor(block))) +
  geom_line(alpha = 0.5) +
  ggtitle(paste0("Compost prelim plot (", Sys.Date(),") : cleaned soil moisture (VWC),\nby water year, nutrient x drought treatment, colored by block")) +
  scale_color_discrete(name = "block") +
  facet_grid(fulltrt~waterYear, scales = "free_x")
trt_plot_clean


# write out clean data preliminary plots to if desired
ggsave(filename = paste0(datpath, "SoilMoisture/SoilMoisture_Figures/Compost_cleanVWC_bylogger.pdf"),
       plot = logger_plot_clean,
       width = 8, height = 8, units = "in")
ggsave(filename = paste0(datpath, "SoilMoisture/SoilMoisture_Figures/Compost_cleanVWC_bytreatment.pdf"),
       plot = trt_plot_clean,
       width = 8, height = 8, unit = "in")


# -- FINISHING -----
write.csv(soilmoisture_master, paste0(datpath, "SoilMoisture/SoilMoisture_CleanedData/SoilMoisture_all_clean.csv"), row.names = F)
