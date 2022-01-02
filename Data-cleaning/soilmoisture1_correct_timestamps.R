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

# search for duplicates in each period (e.g., if two different download files contain the same data)
# rank filename first by date
downloaddates <- data.frame(filename = unique(soilmoisture_all$filename)) %>%
  mutate(download_date = gsub("-", "", substr(filename, 6,12)),
         download_date = as.Date(download_date, "%d%b%y"),
         download_time = as.numeric(str_extract(filename, "[:digit:]{4}(?=[.]xls)"))) %>%
  arrange(download_date, download_time) %>%
  mutate(download_order = 1:nrow(.))

# period 1
soilmoisture_p1 <- left_join(soilmoisture_p1, downloaddates)
soilmoisture_p1$dupcheck <- duplicated(subset(soilmoisture_p1, select = c(logger:vwc)))

# period 2
soilmoisture_p2 <- left_join(soilmoisture_p2, downloaddates)
soilmoisture_p2$dupcheck <- duplicated(subset(soilmoisture_p2, select = c(logger:vwc)))

test <- subset(soilmoisture_p2, month(date_time) == 9 & year(date_time) == 2021 & portid == "B1L2_1", select = c(logger:vwc))
nrow(distinct(test))
test2 <- subset(test, select = c(date_time, vwc)) 
nrow(distinct(test2))
# > the duplicate entries occur at different times.. once they are timestamp matched correctly (below), that's when the timestamps match up. bleh.
# > dates that have duplicates across different download filenames are for b1l2 logger on 9-15-2021

# punt this problem for later in the script
soilmoisture_p1 <- subset(soilmoisture_p1, select = -c(download_date:dupcheck))
soilmoisture_p2 <- subset(soilmoisture_p2, select = -c(download_date:dupcheck))
rm(test, test2)



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



# -- 3. Compile period 1 clean data ------
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

# remake plots from initial compilation loop to be sure looks as expected
logger_plot_clean_p1 <- soilmoisture_master_p1 %>%
  ggplot(aes(dowy, vwc, col = as.factor(port))) + #date_time, 
  geom_line(alpha = 0.5) +
  ggtitle(paste0("Compost prelim plot (", Sys.Date(),  "): cleaned soil moisture (VWC),\nby logger, by water year, colored by port")) +
  scale_color_discrete(name = "port") +
  facet_grid(logger~waterYear, scales = "free_x")
logger_plot_clean_p1

trt_plot_clean_p2 <- soilmoisture_master_p1 %>%
  ggplot(aes(dowy, vwc, group = portid, col = as.factor(block))) +
  geom_line(alpha = 0.5) +
  ggtitle(paste0("Compost prelim plot (", Sys.Date(),") : cleaned soil moisture (VWC),\nby water year, nutrient x drought treatment, colored by block")) +
  scale_color_discrete(name = "block") +
  facet_grid(fulltrt~waterYear, scales = "free_x")
trt_plot_clean_p2




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
# not prepping this clean df as above because there are multiple instances of a logger-port reverting to date-times already recorded. bloats the data.frame
# just making a "clean" df of all expected date-times per logger-port
soilmoisture_clean_p2 <- data.frame(date_time = rep(seq.POSIXt(clean_mintime_p2, clean_maxtime_p2, by = "2 hours"), times = length(unique(soilmoisture_p2$portid)))) %>%
  mutate(portid = rep(unique(soilmoisture_p2$portid), each = length(unique(date_time)))) %>%
  group_by(portid) %>%
  mutate(cleanorder = seq(1, length(date_time), 1)) %>%
  ungroup() %>%
  # join logger/plot trt info
  left_join(distinct(subset(soilmoisture_p2, select = c(logger:comp_trt))))

# review for NAs
summary(soilmoisture_clean_p2) # good to go


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
soilmoisture_master_p2 <- subset(soilmoisture_clean_p2, logger %in% goodrefs$logger) %>%
  left_join(soilmoisture_master_p2[c("portid", "date_time", "filename", "vwc")])
# double check NAs
summary(is.na(soilmoisture_master_p2))
View(subset(soilmoisture_master_p2, is.na(filename))) # just instances where no data recorded, leaving as is for now
with(subset(soilmoisture_master_p2, is.na(filename)), sapply(split(vwc, portid), length))
# 1 could be if no data recorded on first or last date-time of all possible timesteps
with(subset(soilmoisture_master_p2, is.na(filename) & logger != "B3L2"), sapply(split(as.character(date_time), portid), unique))
# yep, all first possible timestep. it's probably in period 1. logger B3L2 stopped recording data in july 2021 it seems..



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

# zoom in on break -- ignore b3l3 since has breaks
subset(datecheck_p2, grepl(b1l2trt[2], fulltrt) & logger != "B3L3" & date(date_time) < as.Date("2020-07-01")) %>% # trackseq > 2250
  ggplot(aes(date_time, vwc, col = portid)) +
  geom_point(alpha = 0.4) +
  geom_line(aes(lty = factor(intevent)), alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates_p2, grepl(b1l2trt[1], fulltrt) & logger != "B3L3" & date(date_time) < as.Date("2020-07-01")), logger = substr(portid, 1,4)), aes(xintercept = date_time, lty = qa_note)) +
  ggtitle(paste("troubleshoot B1L2", b1l2trt[1], "date jump (dotted vert line = break)")) +
  geom_smooth(aes(fill = portid))

# zoom in on period that's off to see if just needs NAs inserted or if time got off..
searchorder <- with(subset(datecheck_p2, grepl(b1l2trt[2], fulltrt) & logger == "B1L2"), min(datorder[grepl("correct", qa_note)]))
searchrange <- (searchorder-10):(searchorder +10)
subset(datecheck_p2,  grepl(b1l2trt[2], fulltrt) & datorder %in% searchrange) %>% #  & grepl("Apr20", filename)
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_vline(data = mutate(subset(adjustdates_p2, grepl(b1l2trt[2], fulltrt) & datorder %in% searchrange), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  geom_point(aes(shape = factor(intevent)), alpha = 0.7) +
  geom_line(aes(lty = factor(intevent)), alpha = 0.7) +
  #geom_smooth(aes(fill = portid)) +
  ggtitle(paste("troubleshoot B1L2", b1l2trt[2], "date jump (dotted vert line = break)"))
# looks pretty spot on

# zoom to last run
lastrunb1l2 <- date(with(adjustdates_p2, date_time[portid == "B1L2_1" & grepl("last", qa_note)]))
# look at all trts
subset(datecheck_p2, fulltrt %in% b1l2trt & logger != "B3L3" & date(date_time) >= lastrunb1l2) %>% # trackseq > 2250
  ggplot(aes(date_time, vwc, col = portid)) +
  geom_point(alpha = 0.4) +
  geom_line(aes(lty = filename), alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates_p2, grepl(b1l2trt[1], fulltrt) & logger != "B3L3" & date(date_time) >= lastrunb1l2), logger = substr(portid, 1,4)), aes(xintercept = date_time, lty = qa_note)) +
  ggtitle(paste("troubleshoot B1L2", b1l2trt[1], "date jump (dotted vert line = break)"))
geom_smooth(aes(fill = portid))
# maybe if shift back interval 2, will connect to interval 3 fine (last interval lines up fine, no need to adjust)


# shift sep 15 file back 10 hrs starting at break. the datorder is fine, but the timestamp needs to be corrected. leave sep 16 as is.
b1l2test <- subset(datecheck_p2, logger == "B1L2") %>%
  rename(raw_datetime = date_time) %>%
  # join cleandat to make date_time
  left_join(distinct(subset(soilmoisture_clean_p2, select = c(date_time, cleanorder))), by = c("datorder" = "cleanorder")) %>%
  mutate(clean_datetime = ifelse(intevent == 2, as.character(date_time), as.character(raw_datetime)),
         clean_datetime = as.POSIXct(clean_datetime, tz = "UTC")) %>%
  # join cleandat for to make clean datorder
  left_join(distinct(subset(soilmoisture_clean_p2, select = c(date_time, cleanorder))), by = c("clean_datetime" = "date_time")) %>%
  rename(cleanorder2 = cleanorder)

# review adjustments
ggplot(b1l2test, aes(cleanorder2, vwc, col = portid)) +
  geom_line(aes(lty = factor(intevent)), alpha = 0.6) +
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b1l2trt & logger == "B2L3"), aes(cleanorder, vwc, group = portid), alpha = 0.6) +
  geom_vline(data = subset(adjustdates_p2, fulltrt %in% b1l2trt & logger == "B1L2"), aes(xintercept = datorder, lty = qa_note))

ggplot(b1l2test, aes(clean_datetime, vwc, col = portid)) +
  geom_line(aes(lty = factor(intevent)), alpha = 0.6) +
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b1l2trt & logger == "B2L3"), aes(date_time, vwc, group = portid), alpha = 0.6) +
  geom_vline(data = subset(adjustdates_p2, fulltrt %in% b1l2trt & logger == "B1L2"), aes(xintercept = date_time, lty = qa_note))

# zoom in by date
ggplot(subset(b1l2test, date(date_time) < as.Date("2020-07-01")), aes(clean_datetime, vwc)) +
  geom_line(aes(lty = factor(intevent), col = portid), alpha = 0.6) +
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b1l2trt & logger == "B2L3" & date(date_time) < as.Date("2020-07-01")), aes(date_time, vwc, group = portid), col = "grey50", alpha = 0.6) +
  geom_vline(data = subset(adjustdates_p2, fulltrt %in% b1l2trt & logger == "B1L2" & date(date_time) < as.Date("2020-07-01")), aes(xintercept = date_time, lty = qa_note))
# looks fine by date
ggplot(subset(b1l2test, date(date_time) >= lastrunb1l2), aes(clean_datetime, vwc)) +
  geom_line(aes(lty = factor(intevent), col = portid), alpha = 0.6) +
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b1l2trt & logger == "B2L3" & date(date_time) >= lastrunb1l2), aes(date_time, vwc, group = portid), col = "grey50", alpha = 0.6) +
  geom_vline(data = subset(adjustdates_p2, fulltrt %in% b1l2trt & logger == "B1L2" & date(date_time) >= lastrunb1l2), aes(xintercept = date_time, lty = qa_note)) 
# looks good

# zoom in by cleanorder
ggplot(subset(b1l2test, cleanorder2 < 500), aes(cleanorder2, vwc, col = portid)) +
  geom_line(aes(lty = factor(intevent)), alpha = 0.6) +
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b1l2trt & logger == "B2L3" & cleanorder < 500), aes(cleanorder, vwc, group = portid), alpha = 0.6) +
  geom_vline(data = subset(adjustdates_p2, fulltrt %in% b1l2trt & logger == "B1L2" & datorder < 500), aes(xintercept = datorder, lty = qa_note))
ggplot(subset(b1l2test, cleanorder2 > 5400), aes(cleanorder2, vwc, col = portid)) +
  geom_line(aes(lty = factor(intevent)), alpha = 0.6) +
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b1l2trt & logger == "B2L3" & cleanorder > 5400), aes(cleanorder, vwc, group = portid), alpha = 0.6) +
  geom_vline(data = subset(adjustdates_p2, fulltrt %in% b1l2trt & logger == "B1L2" & datorder > 5400), aes(xintercept = datorder, lty = qa_note))

# correct B1L2 in soilmoisture_clean_p1
clean_b1l2 <- subset(soilmoisture_clean_p2, grepl("B1L2", portid)) %>%
  left_join(select(b1l2test, portid, filename, clean_datetime, vwc, raw_datetime, timeinterval, timediff, intevent, qa_note, datorder, rowid, raw_datorder), by = c("portid", "date_time" = "clean_datetime"))

# check NAs
summary(clean_b1l2)
data.frame(subset(clean_b1l2, is.na(rowid))) # this is fine. timesteps that have no data between int events 2+3 after timestamp adjustment for period 2

# replot
ggplot(subset(clean_b1l2, date(date_time) <= as.Date("2020-07-01")), aes(date_time, vwc, col = portid)) +
  # plot old points first
  geom_line(data = subset(datecheck_p2, logger == "B1L2" & date(date_time) <= as.Date("2020-07-01")), aes(date_time, vwc, col = portid, group = paste(portid, "old")), lty = 1, alpha = 0.6) +
  geom_line(aes(group = portid), lty = 2) # looks good
# replot
ggplot(subset(clean_b1l2, date(date_time) >= as.Date("2021-09-01")), aes(date_time, vwc, col = portid)) +
  # plot old points first
  geom_line(data = subset(datecheck_p2, logger == "B1L2" & date(date_time) >= as.Date("2021-09-01")), aes(date_time, vwc, col = portid, group = paste(portid, "old")), lty = 1, alpha = 0.6) +
  geom_line(aes(group = portid), lty = 2) # looks good

# append clean b1l2 to good loggers
# > add rowid and global datorder to soilmoisture_master first to keep track of data order
# also add intevent, qa_note, and timediff for easier annotations below
soilmoisture_master_p2 <- left_join(soilmoisture_master_p2, datecheck_p2[c("portid", "datorder", "intevent", "qa_note", "filename", "rowid", "date_time", "raw_datorder")]) %>%
  mutate(raw_datetime = date_time) %>%
  rbind(clean_b1l2[names(.)])

# review
summary(soilmoisture_master_p2)
# who are the NAs?
View(subset(soilmoisture_master_p2, is.na(raw_datetime)))
# is fine -- times when b1l2 stopped recording. can address during annotation at the end.

# clean up
rm(b1l2trt, b1l2test, searchorder, searchrange, lastrunb1l2,
   t, tempend, tempnote, tempstart, tempbreaks, i, event)



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
  left_join(select(b3l4test, portid, cleanorder3, vwc, raw_datetime, timeinterval, timediff, intevent, qa_note, datorder, filename, rowid, raw_datorder), by = c("portid", "cleanorder" = "cleanorder3"))

# check NAs in identifying info
summary(is.na(clean_b3l4))
View(subset(clean_b3l4, is.na(intevent))) # at first and at break, fine

# append
soilmoisture_master_p2 <- rbind(soilmoisture_master_p2, clean_b3l4[names(soilmoisture_master_p2)])

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
  left_join(select(b2l5test, portid, cleanorder3, vwc, raw_datetime, timeinterval, timediff, intevent, qa_note, datorder, filename, rowid, raw_datorder), by = c("portid", "cleanorder" = "cleanorder3"))

# check again
with(clean_b2l5, sapply(split(intevent, portid), length))
with(clean_b2l5, sapply(split(intevent, portid), unique)) # NAs present for intevent bc NAs are now inserted for blips in data recording
summary(is.na(clean_b2l5))
# review
ggplot(clean_b2l5, aes(date_time, vwc)) +
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b2l5trt), aes(date_time, vwc, group = portid), col = "grey50") +
  geom_line(aes(col = portid, lty = intevent %% 2 == 0, group = paste(portid, intevent)))
# zoom in
ggplot(subset(clean_b2l5, date(date_time) > as.Date("2021-03-01")), aes(date_time, vwc)) +
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b2l5trt & date(date_time) > as.Date("2021-03-01")), aes(date_time, vwc+0.05, group = portid), col = "grey50", alpha = 0.7) +
  geom_line(aes(col = portid, lty = intevent %% 2 == 0, group = paste(portid, intevent)), alpha = 0.7)

# append
soilmoisture_master_p2 <- rbind(soilmoisture_master_p2, clean_b2l5[names(soilmoisture_master_p2)])

# review
with(soilmoisture_master_p2, sapply(split(vwc, portid), length))
# all there

# clean up environment
rm(list = ls()[grep("^b2l5", ls())], datecheck_p2_b2l5)


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
b3l3test_1$cleanorder2[b3l3test_1$intevent == 6] <- b3l3test_1$datorder[b3l3test_1$intevent == 6] + 85
#b3l3test_1$cleanorder2[b3l3test_1$intevent == 7 & b3l3test_1$rowid <= 650067] <- b3l3test_1$datorder[b3l3test_1$intevent == 7 & b3l3test_1$rowid <= 650067] + 81
#b3l3test_1$cleanorder2[b3l3test_1$intevent == 7 & b3l3test_1$rowid %in% 650068:650076] <- b3l3test_1$datorder[b3l3test_1$intevent == 7 & b3l3test_1$rowid %in% 650068:650076] + 83
#b3l3test_1$cleanorder2[b3l3test_1$intevent == 7 & b3l3test_1$rowid > 650076] <- b3l3test_1$datorder[b3l3test_1$intevent == 7 & b3l3test_1$rowid > 650076] + 87 # 7 is a needs NA infill by 1 timestep
b3l3test_1$cleanorder2[b3l3test_1$intevent == 7] <- b3l3test_1$datorder[b3l3test_1$intevent == 7] + 87 # 7 is a needs NA infill by 1 timestep
b3l3test_1$cleanorder2[b3l3test_1$intevent == 8] <- b3l3test_1$datorder[b3l3test_1$intevent == 8] + 103
b3l3test_1$cleanorder2[b3l3test_1$intevent == 9] <- b3l3test_1$datorder[b3l3test_1$intevent == 9] + 106
b3l3test_1$cleanorder2[b3l3test_1$intevent == 10] <- b3l3test_1$datorder[b3l3test_1$intevent == 10] + 120
b3l3test_1$cleanorder2[b3l3test_1$intevent == 11] <- b3l3test_1$datorder[b3l3test_1$intevent == 11] + 123
b3l3test_1$cleanorder2[b3l3test_1$intevent == 12] <- b3l3test_1$datorder[b3l3test_1$intevent == 12] + 127
b3l3test_1$cleanorder2[b3l3test_1$intevent == 13] <- b3l3test_1$datorder[b3l3test_1$intevent == 13] + 129
b3l3test_1$cleanorder2[b3l3test_1$intevent == 14] <- b3l3test_1$datorder[b3l3test_1$intevent == 14] + 133
b3l3test_1$cleanorder2[b3l3test_1$intevent == 15] <- b3l3test_1$datorder[b3l3test_1$intevent == 15] + 135
b3l3test_1$cleanorder2[b3l3test_1$intevent == 16] <- b3l3test_1$datorder[b3l3test_1$intevent == 16] + 137
b3l3test_1$cleanorder2[b3l3test_1$intevent == 17] <- b3l3test_1$datorder[b3l3test_1$intevent == 17] + 140

# just the periods to adjust -- zoom in to datorder of interest as needed to review (i.e., this is manual)
ggplot(subset(b3l3test_1, cleanorder2 %in% 4000:5500), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, scale(vwc), col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3trt & cleanorder  %in% 4000:5500), # & cleanorder %in% 4000:4350
            aes(cleanorder, scale(vwc)+0.01, group = portid), col = "grey50", alpha = 0.4) +
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
b3l3test_1 <- assign_NA(b3l3test_1)

# plot to check
ggplot(subset(b3l3test_1, cleanorder3 %in% 4000:5000), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder3, vwc)) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3test_1$fulltrt & cleanorder %in% 4000:5000), # & cleanorder %in% 4000:4350
            aes(cleanorder, vwc, group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  geom_point(aes(shape = intevent %% 2 == 0, col = factor(intevent)), alpha = 0.6) +
  geom_line(alpha = 0.6) #aes(lty = intevent %% 2 == 0)
# worked as intended
# proceed with next

# check lags in cleanorder3 (should never be negative)
sort(unique(b3l3test_1$cleanorder3 - lag(b3l3test_1$cleanorder3)))

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

# adjust periods that need timestamp correct or are the last run
b3l3test_5$cleanorder2[b3l3test_5$intevent == 2] <- b3l3test_5$datorder[b3l3test_5$intevent == 2] + 27
b3l3test_5$cleanorder2[b3l3test_5$intevent == 3] <- b3l3test_5$datorder[b3l3test_5$intevent == 3] + 45
b3l3test_5$cleanorder2[b3l3test_5$intevent == 5] <- b3l3test_5$datorder[b3l3test_5$intevent == 5] + 53
b3l3test_5$cleanorder2[b3l3test_5$intevent == 7] <- b3l3test_5$datorder[b3l3test_5$intevent == 7] + 71
b3l3test_5$cleanorder2[b3l3test_5$intevent == 9] <- b3l3test_5$datorder[b3l3test_5$intevent == 9] + 92
b3l3test_5$cleanorder2[b3l3test_5$intevent == 11] <- b3l3test_5$datorder[b3l3test_5$intevent == 11] + 109
b3l3test_5$cleanorder2[b3l3test_5$intevent == 12] <- b3l3test_5$datorder[b3l3test_5$intevent == 12] + 113
b3l3test_5$cleanorder2[b3l3test_5$intevent == 13] <- b3l3test_5$datorder[b3l3test_5$intevent == 13] + 127
b3l3test_5$cleanorder2[b3l3test_5$intevent == 14] <- b3l3test_5$datorder[b3l3test_5$intevent == 14] + 130
b3l3test_5$cleanorder2[b3l3test_5$intevent == 15] <- b3l3test_5$datorder[b3l3test_5$intevent == 15] + 134
b3l3test_5$cleanorder2[b3l3test_5$intevent == 16] <- b3l3test_5$datorder[b3l3test_5$intevent == 16] + 136
b3l3test_5$cleanorder2[b3l3test_5$intevent == 17] <- b3l3test_5$datorder[b3l3test_5$intevent == 17] + 140
b3l3test_5$cleanorder2[b3l3test_5$intevent == 18] <- b3l3test_5$datorder[b3l3test_5$intevent == 18] + 142
b3l3test_5$cleanorder2[b3l3test_5$intevent == 19] <- b3l3test_5$datorder[b3l3test_5$intevent == 19] + 144
b3l3test_5$cleanorder2[b3l3test_5$intevent == 20] <- b3l3test_5$datorder[b3l3test_5$intevent == 20] + 147

# all dats to see typical relationship
ggplot(subset(b3l3test_5, intevent %in% b3l3_5_timestamps & cleanorder2 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350 #%in% 4000:4325
       aes(cleanorder2, vwc, col = factor(intevent), group = intevent)) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3test_5$fulltrt & cleanorder > 4000), # & cleanorder %in% 4000:4350 #%in% 4000:4325
            aes(cleanorder, vwc, group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  geom_line(aes(lty = intevent %% 2 == 0)) + # very similar to port 1
  # add port 1 to help guide
  geom_line(data = subset(b3l3test_1, cleanorder3 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350
            aes(cleanorder3, vwc+0.03, group = intevent), col = "grey50", alpha = 0.6) # seems ok

# run NA infill
b3l3test_5 <- assign_NA(b3l3test_5)
# check diffs (no negatives should be present)
sort(unique(b3l3test_5$cleanorder3 - lag(b3l3test_5$cleanorder3))) #ok

# replot to check
ggplot(subset(b3l3test_5, cleanorder3 %in% 4000:4500), #, cleanorder2 > 4000 & cleanorder2 < 4350 #%in% 4000:4325
       aes(cleanorder3, vwc, col = factor(intevent), group = intevent)) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3test_5$fulltrt & cleanorder %in% 4000:4500), # & cleanorder %in% 4000:4350 #%in% 4000:4325
            aes(cleanorder, vwc, group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  geom_line(aes(lty = intevent %in% b3l3_5_timestamps)) + # very similar to port 1
  # add port 1 to help guide
  geom_line(data = subset(b3l3test_1, cleanorder3 %in% 4000:4500), #, cleanorder2 > 4000 & cleanorder2 < 4350
            aes(cleanorder3, vwc+0.03, group = intevent), col = "grey50", alpha = 0.6)




# 4.e.3. adjust B3L3 port 2 -----
b3l3test_2 <- subset(datecheck_p2, portid == "B3L3_2") %>%
  # join clean order
  left_join(distinct(subset(soilmoisture_clean_p2, dateselect = c(date_time, cleanorder)))) %>%
  mutate(cleanorder2 = ifelse(intevent == 1, cleanorder, datorder)) # cleanorder diffs from datorder by 1, only attached to intevent 1

# how many unique events, and what type?
with(subset(b3l3test_2, !is.na(qa_note)), sapply(split(intevent, qa_note), unique))
b3l3_2_timestamps <- with(subset(b3l3test_2, !is.na(qa_note)), intevent[grepl("timestamp|last", qa_note)])
b3l3_2_NAs <- with(subset(b3l3test_2, !is.na(qa_note)), intevent[grepl("NA", qa_note)])

# review
ggplot(subset(b3l3test_2, intevent %in% b3l3_2_timestamps & cleanorder2 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, vwc, col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3trt & cleanorder > 4000), # & cleanorder %in% 4000:4350
            aes(cleanorder, vwc, group = portid), col = "grey80", alpha = 0.4) +
  # plot adjusted
  geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  geom_line(aes(lty = intevent %% 2 == 0), alpha = 0.6) +
  # add port 1 to help guide
  geom_line(data = subset(b3l3test_1, cleanorder3 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350
            aes(cleanorder2, vwc, group = intevent), col = "grey30", alpha = 0.6) # seems ok

# try adjusting last pos by difference between port 1 and port 2 (line up the ends)
port2_offset <- max(b3l3test_1$cleanorder3, na.rm =T) - max(b3l3test_2$cleanorder2, na.rm =T)
# adjust
b3l3test_2$cleanorder2[b3l3test_2$intevent == 2] <- b3l3test_2$datorder[b3l3test_2$intevent == 2] + 27
b3l3test_2$cleanorder2[b3l3test_2$intevent == 3] <- b3l3test_2$datorder[b3l3test_2$intevent == 3] + 46
b3l3test_2$cleanorder2[b3l3test_2$intevent == 4] <- b3l3test_2$datorder[b3l3test_2$intevent == 4] + 48
b3l3test_2$cleanorder2[b3l3test_2$intevent == 6] <- b3l3test_2$datorder[b3l3test_2$intevent == 6] + 68
b3l3test_2$cleanorder2[b3l3test_2$intevent == 9] <- b3l3test_2$datorder[b3l3test_2$intevent == 9] + (port2_offset-61)
b3l3test_2$cleanorder2[b3l3test_2$intevent == 10] <- b3l3test_2$datorder[b3l3test_2$intevent == 10] + (port2_offset-44)
b3l3test_2$cleanorder2[b3l3test_2$intevent == 11] <- b3l3test_2$datorder[b3l3test_2$intevent == 11] + (port2_offset-43)
b3l3test_2$cleanorder2[b3l3test_2$intevent == 13] <- b3l3test_2$datorder[b3l3test_2$intevent == 13] + (port2_offset-28)
b3l3test_2$cleanorder2[b3l3test_2$intevent == 14] <- b3l3test_2$datorder[b3l3test_2$intevent == 14] + (port2_offset-24)
b3l3test_2$cleanorder2[b3l3test_2$intevent == 15] <- b3l3test_2$datorder[b3l3test_2$intevent == 15] + (port2_offset-21)
b3l3test_2$cleanorder2[b3l3test_2$intevent == 16] <- b3l3test_2$datorder[b3l3test_2$intevent == 16] + (port2_offset-18)
b3l3test_2$cleanorder2[b3l3test_2$intevent == 17] <- b3l3test_2$datorder[b3l3test_2$intevent == 17] + (port2_offset-15)
b3l3test_2$cleanorder2[b3l3test_2$intevent == 18] <- b3l3test_2$datorder[b3l3test_2$intevent == 18] + (port2_offset-12)
b3l3test_2$cleanorder2[b3l3test_2$intevent == 19] <- b3l3test_2$datorder[b3l3test_2$intevent == 19] + (port2_offset-8)
b3l3test_2$cleanorder2[b3l3test_2$intevent == 27] <- b3l3test_2$datorder[b3l3test_2$intevent == 27] + port2_offset

# review
ggplot(subset(b3l3test_2, intevent %in% c(1,b3l3_2_timestamps) & cleanorder2 %in% 4000:4500), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, scale(vwc), col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3test_2$fulltrt & cleanorder %in% 4000:4500), # & cleanorder %in% 4000:4350
            aes(cleanorder, scale(vwc), group = portid), col = "grey80", alpha = 0.4) +
  # plot adjusted
  #geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  geom_line(aes(lty = intevent %% 2 == 0), alpha = 0.6) +
  # add port 1 to help guide
  geom_line(data = subset(b3l3test_1, cleanorder3 %in% 4000:4500), #, cleanorder2 > 4000 & cleanorder2 < 4350
            aes(cleanorder2, scale(vwc), group = intevent), col = "grey30", alpha = 0.6) # seems ok

# assign clean order for NA jumps
b3l3test_2 <- assign_NA(b3l3test_2)
# check diffs (no negatives should be present)
sort(unique(b3l3test_2$cleanorder3 - lag(b3l3test_2$cleanorder3))) #ok


# replot to check
ggplot(subset(b3l3test_2, cleanorder3 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350 #%in% 4000:4325
       aes(cleanorder3, vwc, col = factor(intevent), group = intevent)) + #group = factor(intevent)
  # reference
  # geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3trt & cleanorder > 4000), # & cleanorder %in% 4000:4350 #%in% 4000:4325
  #           aes(cleanorder, vwc, group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  geom_line(aes(lty = intevent %in% b3l3_2_timestamps)) + # very similar to port 1
  # add port 1 to help guide
  geom_line(data = subset(b3l3test_1, cleanorder3 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350
            aes(cleanorder2, vwc, group = intevent), col = "grey50", alpha = 0.6)
# looks good



# 4.e.4. adjust B3L3 port 4 -----
b3l3test_4 <- subset(datecheck_p2, portid == "B3L3_4") %>%
  # join clean order
  left_join(distinct(subset(soilmoisture_clean_p2, dateselect = c(date_time, cleanorder)))) %>%
  mutate(cleanorder2 = ifelse(intevent == 1, cleanorder, datorder)) # cleanorder diffs from datorder by 0, only attached to intevent 1

# how many unique events, and what type?
with(subset(b3l3test_4, !is.na(qa_note)), sapply(split(intevent, qa_note), unique))
b3l3_4_timestamps <- with(subset(b3l3test_4, !is.na(qa_note)), intevent[grepl("timestamp|last", qa_note)])
b3l3_4_NAs <- with(subset(b3l3test_4, !is.na(qa_note)), intevent[grepl("NA", qa_note)])

# review
ggplot(subset(b3l3test_4, intevent %in% c(1,b3l3_4_timestamps) & cleanorder2 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, vwc, col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3trt & cleanorder > 4000), # & cleanorder %in% 4000:4350
            aes(cleanorder, vwc, group = portid), col = "grey80", alpha = 0.4) +
  # plot adjusted
  geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  geom_line(aes(lty = intevent %% 2 == 0), alpha = 0.6) +
  # add port 1 to help guide
  geom_line(data = subset(b3l3test_1, cleanorder3 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350
            aes(cleanorder2, vwc, group = intevent), col = "grey30", alpha = 0.6) # seems ok

# try adjusting last pos by difference between port 1 and port 2 (line up the ends)
port4_offset <- max(b3l3test_1$cleanorder3, na.rm =T) - max(b3l3test_4$cleanorder2, na.rm =T)

# adjust
b3l3test_4$cleanorder2[b3l3test_4$intevent == 2] <- b3l3test_4$datorder[b3l3test_4$intevent == 2] + 27
b3l3test_4$cleanorder2[b3l3test_4$intevent == 3] <- b3l3test_4$datorder[b3l3test_4$intevent == 3] + 46
b3l3test_4$cleanorder2[b3l3test_4$intevent == 5] <- b3l3test_4$datorder[b3l3test_4$intevent == 5] + 48
b3l3test_4$cleanorder2[b3l3test_4$intevent == 7] <- b3l3test_4$datorder[b3l3test_4$intevent == 7] + 68
b3l3test_4$cleanorder2[b3l3test_4$intevent == 10] <- b3l3test_4$datorder[b3l3test_4$intevent == 10] + (port4_offset-61)
b3l3test_4$cleanorder2[b3l3test_4$intevent == 12] <- b3l3test_4$datorder[b3l3test_4$intevent == 12] + (port4_offset-45)
b3l3test_4$cleanorder2[b3l3test_4$intevent == 13] <- b3l3test_4$datorder[b3l3test_4$intevent == 13] + (port4_offset-42)
b3l3test_4$cleanorder2[b3l3test_4$intevent == 14] <- b3l3test_4$datorder[b3l3test_4$intevent == 14] + (port4_offset-28)
b3l3test_4$cleanorder2[b3l3test_4$intevent == 15] <- b3l3test_4$datorder[b3l3test_4$intevent == 15] + (port4_offset-25)
b3l3test_4$cleanorder2[b3l3test_4$intevent == 16] <- b3l3test_4$datorder[b3l3test_4$intevent == 16] + (port4_offset-21)
b3l3test_4$cleanorder2[b3l3test_4$intevent == 17] <- b3l3test_4$datorder[b3l3test_4$intevent == 17] + (port4_offset-18)
b3l3test_4$cleanorder2[b3l3test_4$intevent == 18] <- b3l3test_4$datorder[b3l3test_4$intevent == 18] + (port4_offset-15)
b3l3test_4$cleanorder2[b3l3test_4$intevent == 19] <- b3l3test_4$datorder[b3l3test_4$intevent == 19] + (port4_offset-12)
b3l3test_4$cleanorder2[b3l3test_4$intevent == 20] <- b3l3test_4$datorder[b3l3test_4$intevent == 20] + (port4_offset-8)
b3l3test_4$cleanorder2[b3l3test_4$intevent == 26] <- b3l3test_4$datorder[b3l3test_4$intevent == 26] + port4_offset

# review
ggplot(subset(b3l3test_4, intevent %in% c(1,b3l3_4_timestamps) & cleanorder2 %in% 4000:4500), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, scale(vwc), col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3trt & cleanorder %in% 4000:4500), # & cleanorder %in% 4000:4350
            aes(cleanorder, scale(vwc), group = portid), col = "grey80", alpha = 0.4) +
  # plot adjusted
  #geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  geom_line(aes(lty = intevent %% 2 == 0)) +
  # add port 1 to help guide
  geom_line(data = subset(b3l3test_2, cleanorder3 %in% 4000:4500), #, cleanorder2 > 4000 & cleanorder2 < 4350
            aes(cleanorder3, scale(vwc), group = intevent), col = "grey30", alpha = 0.6) # is lined up

# assign clean datorder for NA jumps
b3l3test_4 <- assign_NA(b3l3test_4)
# check diffs (no negatives should be present)
sort(unique(b3l3test_4$cleanorder3 - lag(b3l3test_4$cleanorder3))) #ok

# replot to check
ggplot(subset(b3l3test_4, cleanorder3 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350 #%in% 4000:4325
       aes(cleanorder3, vwc, col = factor(intevent), group = intevent)) + #group = factor(intevent)
  # reference
  # geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b3l3trt & cleanorder > 4000), # & cleanorder %in% 4000:4350 #%in% 4000:4325
  #           aes(cleanorder, vwc, group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  geom_line(aes(lty = intevent %in% b3l3_4_timestamps)) + # very similar to port 1
  # add port 1 to help guide
  geom_line(data = subset(b3l3test_1, cleanorder3 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350
            aes(cleanorder2, vwc, group = intevent), col = "grey50", alpha = 0.6)
# looks good



# 4.e.5. prep B3L3 port 3 -----
b3l3test_3 <- subset(datecheck_p2, portid == "B3L3_3") %>%
  # join clean order
  left_join(distinct(subset(soilmoisture_clean_p2, dateselect = c(date_time, cleanorder)))) %>%
  mutate(cleanorder3 = ifelse(intevent == 1, cleanorder, datorder)) %>% # cleanorder diffs from datorder by 0, only attached to intevent 1
  #subset to just intevent 1
  subset(intevent == 1)


# 4.e.6. stack all cleaned B3L3 ports -----
# port 3 is off on timestamps on dates it doesn't have any data anyway
# stack via rbind
clean_b3l3 <- data.frame()
for(i in 1:5){
  tempdat <- subset(soilmoisture_clean_p2, portid == paste0("B3L3_",i)) %>%
    rename(clean_datetime = date_time) %>%
    left_join(select(get(paste0("b3l3test_",i)), portid, cleanorder3, vwc, date_time, timeinterval, timediff, intevent, qa_note, datorder, filename, rowid, raw_datorder), by = c("portid", "cleanorder" = "cleanorder3")) %>%
    rename(raw_datetime = date_time, date_time = clean_datetime)
  clean_b3l3 <- rbind(clean_b3l3, tempdat)
}

# check length by port
with(clean_b3l3, sapply(split(date_time, portid), length)) # all there
summary(clean_b3l3)
# where are the NAs in datorder? (should be mostly port 3)
with(clean_b3l3, sapply(split(datorder, portid), function(x) summary(is.na(x)))) # NAs in others are from inserts

# review
ggplot(clean_b3l3, aes(date_time, vwc, group = portid, col = fulltrt)) +
  geom_line()

# add clean b3l3 to master
soilmoisture_master_p2 <- rbind(soilmoisture_master_p2, clean_b3l3[names(soilmoisture_master_p2)])
# check NAs in VWC
with(soilmoisture_master_p2, sapply(split(vwc, portid), function(x) summary(is.na(x))))
# how many vwc vals present in b3l3_3 event 1? (looking for 438)
summary(is.na(b3l3test_3$vwc)) # matches

# clean up environment
rm(list = ls()[grep("^b3l3|^port", ls())], i, int2date, tempdat, datecheck_p2_b3l3)



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

# follow same approach as last logger. fix one port, then match others to that.
# b2l4_2 looks like first complete

# 4.f.1. adjust B2L4 port 2 -----
b2l4test_2 <- subset(datecheck_p2, portid == "B2L4_2") %>%
  # join clean order
  left_join(distinct(subset(soilmoisture_clean_p2, dateselect = c(date_time, cleanorder)))) %>%
  mutate(cleanorder2 = ifelse(intevent == 1, cleanorder, datorder)) # cleanorder diffs from datorder by 0, only attached to intevent 1

# how many unique events, and what type?
with(subset(b2l4test_2, !is.na(qa_note)), sapply(split(intevent, qa_note), unique))
b2l4_2_timestamps <- with(subset(b2l4test_2, !is.na(qa_note)), intevent[grepl("timestamp|last", qa_note)])
b2l4_2_NAs <- with(subset(b2l4test_2, !is.na(qa_note)), intevent[grepl("NA", qa_note)])

# review
ggplot(subset(b2l4test_2, intevent %in% c(1,b2l4_2_timestamps) & cleanorder2 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, vwc, col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, grepl(unique(b2l4test_2$ppt_trt), fulltrt) & logger %in% goodrefs$logger & cleanorder > 4000), # & cleanorder %in% 4000:4350
            aes(cleanorder, vwc, group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  #geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  geom_line(aes(lty = intevent %% 2 == 0), alpha = 0.7) +
  facet_grid(fulltrt~logger) # try using b2l3

# period 1 gets off track with b2l3, try looking at its actual treatment match
ggplot(subset(b2l4test_2, intevent %in% c(1,b2l4_2_timestamps) & cleanorder2 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, scale(vwc), col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, fulltrt %in% b2l4test_2$fulltrt & cleanorder > 4000), # & cleanorder %in% 4000:4350
            aes(cleanorder, scale(vwc), group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  #geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  geom_line(aes(lty = intevent %% 2 == 0), alpha = 0.7)
# > also gets off track in period 1, but by interval 2 is back on track. match as usual


# try shifting last datorder by last datorder in reference (line up the ends)
b2l4port2_offset <- with(subset(soilmoisture_master_p2, logger == "B2L3" & grepl(unique(b2l4test_2$ppt_trt), fulltrt)), max(cleanorder, na.rm =T)) - max(b2l4test_2$cleanorder2, na.rm =T)

b2l4test_2$cleanorder2[b2l4test_2$intevent == 2] <- b2l4test_2$datorder[b2l4test_2$intevent == 2] + 11
b2l4test_2$cleanorder2[b2l4test_2$intevent == 3] <- b2l4test_2$datorder[b2l4test_2$intevent == 3] + 15
b2l4test_2$cleanorder2[b2l4test_2$intevent == 4] <- b2l4test_2$datorder[b2l4test_2$intevent == 4] + 18
b2l4test_2$cleanorder2[b2l4test_2$intevent == 5] <- b2l4test_2$datorder[b2l4test_2$intevent == 5] + 20
b2l4test_2$cleanorder2[b2l4test_2$intevent == 12] <- b2l4test_2$datorder[b2l4test_2$intevent == 12] + 30
b2l4test_2$cleanorder2[b2l4test_2$intevent == 27] <- b2l4test_2$datorder[b2l4test_2$intevent == 27] + b2l4port2_offset

# review
ggplot(subset(b2l4test_2, intevent %in% c(1,b2l4_2_timestamps) & cleanorder2 %in% 4500:5000), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, scale(vwc), col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, logger == "B2L3" & grepl(unique(b2l4test_2$ppt_trt), fulltrt) & cleanorder %in% 4500:5000), # & cleanorder %in% 4000:4350
            aes(cleanorder, scale(vwc), group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  geom_line(aes(), alpha = 0.7) #lty = intevent %% 2 == 0

# apply NA infill
b2l4test_2 <- assign_NA(b2l4test_2)

# review
ggplot(subset(b2l4test_2, cleanorder3 > 4000), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder3, scale(vwc))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, logger == "B2L3" & grepl(unique(b2l4test_2$ppt_trt), fulltrt) & cleanorder > 4000), # & cleanorder %in% 4000:4350
            aes(cleanorder, scale(vwc), group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  # geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  # geom_line(aes(lty = intevent %in% c(1,b2l4_2_timestamps)), alpha = 0.7) #lty = intevent %% 2 == 
  geom_line(col = "grey30", alpha = 0.7) +
  geom_point(aes(col = factor(intevent), shape = intevent %in% c(1,b2l4_2_timestamps)), alpha = 0.6)
# looks okay

# check diffs between cleanorder to be sure
sort(unique(b2l4test_2$cleanorder3 - lag(b2l4test_2$cleanorder3))) #good



# 4.f.2. adjust B2L4 port 3 -----
b2l4test_3 <- subset(datecheck_p2, portid == "B2L4_3") %>%
  # join clean order
  left_join(distinct(subset(soilmoisture_clean_p2, dateselect = c(date_time, cleanorder)))) %>%
  mutate(cleanorder2 = ifelse(intevent == 1, cleanorder, datorder)) # cleanorder diffs from datorder by 0, only attached to intevent 1

# how many unique events, and what type?
with(subset(b2l4test_3, !is.na(qa_note)), sapply(split(intevent, qa_note), unique))
b2l4_3_timestamps <- with(subset(b2l4test_3, !is.na(qa_note)), intevent[grepl("timestamp|last", qa_note)])
b2l4_3_NAs <- with(subset(b2l4test_3, !is.na(qa_note)), intevent[grepl("NA", qa_note)])

# try shifting last datorder by last datorder in reference (line up the ends)
b2l4port3_offset <- with(subset(soilmoisture_master_p2, logger == "B2L3" & grepl(unique(b2l4test_3$ppt_trt), fulltrt)), max(cleanorder, na.rm =T)) - max(b2l4test_3$cleanorder2, na.rm =T)

b2l4test_3$cleanorder2[b2l4test_3$intevent == 2] <- b2l4test_3$datorder[b2l4test_3$intevent == 2] + 11
b2l4test_3$cleanorder2[b2l4test_3$intevent == 3] <- b2l4test_3$datorder[b2l4test_3$intevent == 3] + 15
b2l4test_3$cleanorder2[b2l4test_3$intevent == 4] <- b2l4test_3$datorder[b2l4test_3$intevent == 4] + 18
b2l4test_3$cleanorder2[b2l4test_3$intevent == 5] <- b2l4test_3$datorder[b2l4test_3$intevent == 5] + 21
b2l4test_3$cleanorder2[b2l4test_3$intevent == 13] <- b2l4test_3$datorder[b2l4test_3$intevent == 13] + 33
b2l4test_3$cleanorder2[b2l4test_3$intevent == 23] <- b2l4test_3$datorder[b2l4test_3$intevent == 23] + 46
b2l4test_3$cleanorder2[b2l4test_3$intevent == 25] <- b2l4test_3$datorder[b2l4test_3$intevent == 25] + (b2l4port3_offset+1)

# review
ggplot(subset(b2l4test_3, intevent %in% c(1,b2l4_3_timestamps) & cleanorder2 %in% 4000:5000), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, scale(vwc), col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, logger == "B2L3" & grepl(unique(b2l4test_3$ppt_trt), fulltrt) & cleanorder %in% 4000:5000), # & cleanorder %in% 4000:4350
            aes(cleanorder, scale(vwc), group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  geom_line(aes(), alpha = 0.7)+  #lty = intevent %% 2 == 0
  # add port 2 to help guide
  geom_line(data = subset(b2l4test_2, cleanorder3 %in% 4000:5000), #, cleanorder2 > 4000 & cleanorder2 < 4350
            aes(cleanorder3, scale(vwc), group = intevent), col = "grey30", alpha = 0.6) # seems ok

# apply NA infill
b2l4test_3 <- assign_NA(b2l4test_3)
# review
ggplot(subset(b2l4test_3, cleanorder3 %in% 4000:5600), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder3, scale(vwc))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, logger == "B2L3" & grepl(unique(b2l4test_2$ppt_trt), fulltrt) & cleanorder %in% 4000:5600), # & cleanorder %in% 4000:4350
            aes(cleanorder, scale(vwc), group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  # geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  # geom_line(aes(lty = intevent %in% c(1,b2l4_2_timestamps)), alpha = 0.7) #lty = intevent %% 2 == 
  geom_line(col = "grey30", alpha = 0.7) +
  geom_point(aes(col = factor(intevent), shape = intevent %in% c(1,b2l4_3_timestamps)), alpha = 0.6)
# looks okay

# check diffs between cleanorder to be sure
sort(unique(b2l4test_3$cleanorder3 - lag(b2l4test_3$cleanorder3))) # no negatives. good.



# 4.f.3. adjust B2L4 port 4 -----
b2l4test_4 <- subset(datecheck_p2, portid == "B2L4_4") %>%
  # join clean order
  left_join(distinct(subset(soilmoisture_clean_p2, dateselect = c(date_time, cleanorder)))) %>%
  mutate(cleanorder2 = ifelse(intevent == 1, cleanorder, datorder)) # cleanorder diffs from datorder by 1

# how many unique events, and what type?
with(subset(b2l4test_4, !is.na(qa_note)), sapply(split(intevent, qa_note), unique))
b2l4_4_timestamps <- with(subset(b2l4test_4, !is.na(qa_note)), intevent[grepl("timestamp|last", qa_note)])
b2l4_4_NAs <- with(subset(b2l4test_4, !is.na(qa_note)), intevent[grepl("NA", qa_note)])

# try shifting last datorder by last datorder in reference (line up the ends)
b2l4port4_offset <- with(subset(soilmoisture_master_p2, logger == "B2L3" & grepl(unique(b2l4test_4$ppt_trt), fulltrt)), max(cleanorder, na.rm =T)) - max(b2l4test_4$cleanorder2, na.rm =T)

# adjust
b2l4test_4$cleanorder2[b2l4test_4$intevent == 2] <- b2l4test_4$datorder[b2l4test_4$intevent == 2] + 11
b2l4test_4$cleanorder2[b2l4test_4$intevent == 3] <- b2l4test_4$datorder[b2l4test_4$intevent == 3] + 15
b2l4test_4$cleanorder2[b2l4test_4$intevent == 4] <- b2l4test_4$datorder[b2l4test_4$intevent == 4] + 18
b2l4test_4$cleanorder2[b2l4test_4$intevent == 6] <- b2l4test_4$datorder[b2l4test_4$intevent == 6] + 21
b2l4test_4$cleanorder2[b2l4test_4$intevent == 16] <- b2l4test_4$datorder[b2l4test_4$intevent == 16] + 33
b2l4test_4$cleanorder2[b2l4test_4$intevent == 21] <- b2l4test_4$datorder[b2l4test_4$intevent == 21] + 45
b2l4test_4$cleanorder2[b2l4test_4$intevent == 38] <- b2l4test_4$datorder[b2l4test_4$intevent == 38] + (b2l4port4_offset)

# review
ggplot(subset(b2l4test_4, intevent %in% c(1,b2l4_4_timestamps) & cleanorder2 %in% 4750:5000), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, scale(vwc), col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, logger == "B2L3" & grepl(unique(b2l4test_4$ppt_trt), fulltrt) & cleanorder %in% 4750:5000), # & cleanorder %in% 4000:4350
            aes(cleanorder, scale(vwc), group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  geom_line(aes(), alpha = 0.7)+  #lty = intevent %% 2 == 0
  # add port 2 to help guide
  geom_line(data = subset(b2l4test_2, cleanorder3 %in% 4750:5000), #, cleanorder2 > 4000 & cleanorder2 < 4350
            aes(cleanorder3, scale(vwc), group = intevent), col = "grey30", alpha = 0.6) # seems ok

# review
ggplot(subset(b2l4test_4, cleanorder2 %in% 4750:5000), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, scale(vwc), col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, logger == "B2L3" & grepl(unique(b2l4test_4$ppt_trt), fulltrt) & cleanorder %in% 4750:5000), # & cleanorder %in% 4000:4350
            aes(cleanorder, scale(vwc), group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  geom_point(aes(shape = intevent %in% c(1,b2l4_4_timestamps)), alpha = 0.6) +
  geom_line(aes(), alpha = 0.7)+  #lty = intevent %% 2 == 0
  # add port 2 to help guide
  geom_line(data = subset(b2l4test_2, cleanorder3 %in% 4750:5000), #, cleanorder2 > 4000 & cleanorder2 < 4350
            aes(cleanorder3, scale(vwc), group = intevent), col = "grey30", alpha = 0.6)


# apply NA infill
b2l4test_4 <- assign_NA(b2l4test_4)
# review
ggplot(subset(b2l4test_4, cleanorder3 %in% 4000:5600), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder3, vwc+0.02)) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, logger == "B2L3" & grepl(unique(b2l4test_4$ppt_trt), fulltrt) & cleanorder %in% 4000:5600), # & cleanorder %in% 4000:4350
            aes(cleanorder, vwc, group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  # geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  # geom_line(aes(lty = intevent %in% c(1,b2l4_2_timestamps)), alpha = 0.7) #lty = intevent %% 2 == 
  geom_line(aes(col = factor(intevent)), alpha = 0.7) + #col = "grey30", 
  geom_point(aes(col = factor(intevent), shape = intevent %in% c(1,b2l4_4_timestamps)), alpha = 0.6)
# looks okay

# check diffs between cleanorder to be sure
sort(unique(b2l4test_4$cleanorder3 - lag(b2l4test_4$cleanorder3)))

# negative initially present at event 21, played with adjustments and plots above until negative diffs gone 
#b2l4test_4$diff <- b2l4test_4$cleanorder3 - lag(b2l4test_4$cleanorder3)
#b2l4test_4 <- subset(b2l4test_4, select = -diff)



# 4.f.4. adjust B2L4 port 5 -----
b2l4test_5 <- subset(datecheck_p2, portid == "B2L4_5") %>%
  # join clean order
  left_join(distinct(subset(soilmoisture_clean_p2, dateselect = c(date_time, cleanorder)))) %>%
  mutate(cleanorder2 = ifelse(intevent == 1, cleanorder, datorder)) # cleanorder diffs from datorder by 1

# how many unique events, and what type?
with(subset(b2l4test_5, !is.na(qa_note)), sapply(split(intevent, qa_note), unique))
b2l4_5_timestamps <- with(subset(b2l4test_5, !is.na(qa_note)), intevent[grepl("timestamp|last", qa_note)])
b2l4_5_NAs <- with(subset(b2l4test_5, !is.na(qa_note)), intevent[grepl("NA", qa_note)])

# try shifting last datorder by last datorder in reference (line up the ends)
b2l4port5_offset <- with(subset(soilmoisture_master_p2, logger == "B2L3" & grepl(unique(b2l4test_5$ppt_trt), fulltrt)), max(cleanorder, na.rm =T)) - max(b2l4test_5$cleanorder2, na.rm =T)

# adjust
b2l4test_5$cleanorder2[b2l4test_5$intevent == 2] <- b2l4test_5$datorder[b2l4test_5$intevent == 2] + 11
b2l4test_5$cleanorder2[b2l4test_5$intevent == 3] <- b2l4test_5$datorder[b2l4test_5$intevent == 3] + 15
b2l4test_5$cleanorder2[b2l4test_5$intevent == 4] <- b2l4test_5$datorder[b2l4test_5$intevent == 4] + 18
b2l4test_5$cleanorder2[b2l4test_5$intevent == 5] <- b2l4test_5$datorder[b2l4test_5$intevent == 5] + 20
b2l4test_5$cleanorder2[b2l4test_5$intevent == 11] <- b2l4test_5$datorder[b2l4test_5$intevent == 11] + 25
b2l4test_5$cleanorder2[b2l4test_5$intevent == 21] <- b2l4test_5$datorder[b2l4test_5$intevent == 21] + 39
b2l4test_5$cleanorder2[b2l4test_5$intevent == 24] <- b2l4test_5$datorder[b2l4test_5$intevent == 24] + (b2l4port5_offset)


# review
ggplot(subset(b2l4test_5, intevent %in% c(1,b2l4_5_timestamps) & cleanorder2 %in% 4000:5000), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, vwc, col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, logger == "B2L3" & grepl(unique(b2l4test_5$ppt_trt), fulltrt) & cleanorder %in% 4000:5000), # & cleanorder %in% 4000:4350
            aes(cleanorder, vwc, group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  #geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  geom_line(aes(), alpha = 0.7)+  #lty = intevent %% 2 == 0
  # add port 2 to help guide
  geom_line(data = subset(b2l4test_2, cleanorder3 %in% 4000:5000), #, cleanorder2 > 4000 & cleanorder2 < 4350
            aes(cleanorder3, vwc, group = intevent), col = "grey30", alpha = 0.6) # seems ok

# look at everything to help with timestamp adjustments
# review
ggplot(subset(b2l4test_5, cleanorder2 %in% 4750:5000), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, scale(vwc), col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, logger == "B2L3" & grepl(unique(b2l4test_5$ppt_trt), fulltrt) & cleanorder %in% 4750:5000), # & cleanorder %in% 4000:4350
            aes(cleanorder, scale(vwc), group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  geom_point(aes(shape = intevent %in% c(1,b2l4_5_timestamps)), alpha = 0.6) +
  geom_line(aes(), alpha = 0.7)+  #lty = intevent %% 2 == 0
  # add port 2 to help guide
  geom_line(data = subset(b2l4test_2, cleanorder3 %in% 4750:5000), #, cleanorder2 > 4000 & cleanorder2 < 4350
            aes(cleanorder3, scale(vwc), group = intevent), col = "grey30", alpha = 0.6)


# apply NA infill
b2l4test_5 <- assign_NA(b2l4test_5)
# review
ggplot(subset(b2l4test_5, cleanorder3 %in% 5000:5600), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder3, vwc)) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, logger == "B2L3" & portid != "B2L3_1" & grepl(unique(b2l4test_5$ppt_trt), fulltrt) & cleanorder %in% 5000:5600), # & cleanorder %in% 4000:4350
            aes(cleanorder, vwc+0.04, group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  # geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  geom_line(aes(col = factor(intevent), lty = intevent %in% c(1,b2l4_2_timestamps)), alpha = 0.7) #lty = intevent %% 2 == 
#geom_line(col = "grey30", alpha = 0.7) + #col = "grey30", 
#geom_point(aes(col = factor(intevent), shape = intevent %in% c(1,b2l4_5_timestamps)), alpha = 0.6) +
#scale_linetype_manual(values = c(2, 1))
# looks okay

# check diffs between cleanorder to be sure
sort(unique(b2l4test_5$cleanorder3 - lag(b2l4test_5$cleanorder3)))
# b2l4test_5$diff <- b2l4test_5$cleanorder3 - lag(b2l4test_5$cleanorder3)
# b2l4test_5 <- subset(b2l4test_5, select = -diff)



# 4.f.5. prep B2L4 port 1 -----
b2l4test_1 <- subset(datecheck_p2, portid == "B2L4_1") %>%
  # join clean order
  left_join(distinct(subset(soilmoisture_clean_p2, dateselect = c(date_time, cleanorder)))) %>%
  mutate(cleanorder2 = ifelse(intevent == 1, cleanorder, datorder)) # cleanorder diffs from datorder by 1

# how many unique events, and what type?
with(subset(b2l4test_1, !is.na(qa_note)), sapply(split(intevent, qa_note), unique))
b2l4_1_timestamps <- with(subset(b2l4test_1, !is.na(qa_note)), intevent[grepl("timestamp|last", qa_note)])
b2l4_1_NAs <- with(subset(b2l4test_1, !is.na(qa_note)), intevent[grepl("NA", qa_note)])

# review
ggplot(subset(b2l4test_1, cleanorder2 %in% 2000:5600), #, cleanorder2 > 4000 & cleanorder2 < 4350 
       aes(cleanorder2, scale(vwc), col = factor(intevent))) + #group = factor(intevent)
  # reference
  geom_line(data = subset(soilmoisture_master_p2, logger == "B2L3" & grepl(unique(b2l4test_1$ppt_trt), fulltrt) & cleanorder %in% 2000:5600), # & cleanorder %in% 4000:4350
            aes(cleanorder, scale(vwc), group = portid), col = "grey50", alpha = 0.4) +
  # plot adjusted
  geom_point(aes(shape = intevent %% 2 == 0), alpha = 0.6) +
  geom_line(aes(), alpha = 0.7)+  #lty = intevent %% 2 == 0
  # add port 2 to help guide
  geom_line(data = subset(b2l4test_2, cleanorder3 %in% 2000:5600), #, cleanorder2 > 4000 & cleanorder2 < 4350
            aes(cleanorder3, scale(vwc), group = intevent), col = "grey30", alpha = 0.6) 

with(b2l4test_1, sapply(split(vwc, intevent), summary))
# no data to adjust after interval event 1, just prep to stack

b2l4test_1 <- rename(b2l4test_1, cleanorder3 = cleanorder2) %>%
  #subset to just intevent 1
  subset(intevent == 1)


# 4.f.6. stack all cleaned B2L4  -----
# stack via rbind
clean_b2l4 <- data.frame()
for(i in 1:5){
  tempdat <- subset(soilmoisture_clean_p2, portid == paste0("B2L4_",i)) %>%
    rename(clean_datetime = date_time) %>%
    left_join(select(get(paste0("b2l4test_",i)), portid, cleanorder3, vwc, date_time, timeinterval, timediff, intevent, qa_note, datorder, filename, rowid, raw_datorder), by = c("portid", "cleanorder" = "cleanorder3")) %>%
    rename(raw_datetime = date_time, date_time = clean_datetime)
  clean_b2l4 <- rbind(clean_b2l4, tempdat)
}

# check length by port
with(clean_b2l4, sapply(split(date_time, portid), length)) # all there
summary(clean_b2l4)
# where are the NAs in datorder? (should be mostly port 3)
with(clean_b2l4, sapply(split(datorder, portid), function(x) summary(is.na(x)))) # NAs in others are from inserts

# review
ggplot(clean_b2l4, aes(date_time, vwc, group = portid, col = portid)) +
  geom_line(alpha = 0.6) +
  facet_grid(fulltrt~., scales = "free_y")

# add clean b2l4 to master
soilmoisture_master_p2 <- rbind(soilmoisture_master_p2, clean_b2l4[names(soilmoisture_master_p2)])

# clean up environment
rm(list = ls()[grep("^b2l4", ls())], tempdat, i)



# -- 5. Compile period 2 clean  data -----
# add qa notes for NAs
# want cleanorder, correct timestamp, treatment info, water year info, vwc, qa notes, raw timestamp, rawdatorder, and sourcefile (filename)
names(soilmoisture_master_p1)
names(soilmoisture_master_p2)
# make copy before proceeding in case need to start over
#copy <- soilmoisture_master_p2
soilmoisture_master_p2 <- copy

# make additional qa_note col for annotations
soilmoisture_master_p2$qa_note2 <- NA
class(soilmoisture_master_p2$qa_note2) <- "character" # reclass from logical


# infill qa_notes
# ignore b3l3 port 3 and b2l4 port 1 after intevent 1 (they stopped measuring vwc even tho data records present -- not possible to correct timestamps with no vwc data)
for(i in unique(adjustdates_p2$portid)){
  # subset from soilmoisture_master since some things have changed
  #tempadjust <- subset(soilmoisture_master_p2, portid == i & !is.na(qa_note))
  tempadjust <- subset(adjustdates_p2, portid == i)
  
  for(r in 1:nrow(tempadjust)){
    # skip b3l3 port 1 and b2l4 port 1 -- will annotate after intevent 1 below
    if(tempadjust$portid[r] %in% c("B3L3_3", "B2L4_1") & tempadjust$intevent[r] > 1){
      next
    }
    # append qa note dependent on what done
    ## NA infill
    if(grepl("NA infill", tempadjust$qa_note[r])){
      # id where interval event with NA flag starts
      tempend <- min(with(soilmoisture_master_p2, which(portid == i & intevent == tempadjust$intevent[r] & !is.na(intevent))))
      # id where previous interval event ends
      tempstart <- max(with(soilmoisture_master_p2, which(portid == i & intevent == (tempadjust$intevent[r]-1) & !is.na(intevent))))
      # id NA rows created between the two that needs annotation
      NArows <- (tempstart+1):(tempend-1)
      # note only goes to spot before where indicated
      soilmoisture_master_p2$qa_note2[NArows] <- paste("Data missing, NA added, no timestamp recorded by", i, "for", (length(NArows)), "2hr-timesteps")
      next
    }
    # pull interval event to infill
    temp_event <- tempadjust$intevent[r]
    # pull where interval begins
    tempstart <- with(soilmoisture_master_p2, which(portid == i & datorder == tempadjust$datorder[r] & raw_datetime == tempadjust$date_time[r]))
    # pull clean order of pertinent interval
    temporder <- soilmoisture_master_p2$cleanorder[tempstart]
    # pull max clean order of previous interval
    temp_prevorder <- max(with(soilmoisture_master_p2, cleanorder[portid == i & intevent == (temp_event-1) & !is.na(intevent)]))
    # pull adjustment amount for note
    tempdiff <- soilmoisture_master_p2$date_time[tempstart] - as.POSIXlt(soilmoisture_master_p2$raw_datetime[tempstart], tz = "UTC")
    temptime <- round(as.numeric(tempdiff),2)
    tempunit <- units(tempdiff)
    if(length(tempdiff) >1){
      print(paste("tempdiff > 1:", i))
    }
    # pull rows to annotate
    temprows <- with(soilmoisture_master_p2, which(portid == i & intevent == temp_event))
    ## timestamp shift or last run
    if(grepl("timestamp|last", tempadjust$qa_note[r])){
      # annotate
      if(tempdiff <= 0){
        soilmoisture_master_p2$qa_note2[temprows] <- paste("Timestamp off, data-time shifted backwards", paste(temptime, tempunit), paste0("(data sequence shifted ", (temporder - temp_prevorder)-1," 2hr-timesteps from last collection run  before break)"), "to correct timestamp")
      }else{
        soilmoisture_master_p2$qa_note2[temprows] <- paste("Timestamp off, data-time shifted forwards", paste(temptime, tempunit), paste0("(data sequence shifted ", (temporder - temp_prevorder)-1," 2hr-timesteps from last collection run before break)"), "to correct timestamp")
      }
    }
  }
  
  # annotate b3l3 port 3 and b2l4 port 1 after intevent
  if(tempadjust$portid[1] %in% c("B3L3_3", "B2L4_1")){
    # pull rows to annotate
    startorder <- with(soilmoisture_master_p2, max(cleanorder[portid == i & intevent == 1], na.rm = T))
    endorder <- with(soilmoisture_master_p2, max(cleanorder[portid == i], na.rm = T))
    tempstart <- with(soilmoisture_master_p2, which(portid == i & cleanorder == (startorder+1)))
    tempend <- with(soilmoisture_master_p2, which(portid == i & cleanorder == endorder))
    # add note
    soilmoisture_master_p2$qa_note2[tempstart:tempend] <- "Timestamp off in raw data and no VWC value present, can't correct timestamp"
  }
  
  # annotate first timestamp data recorded for each logger in qa_note for cleanorder == 1
  # grab first timestamp
  firstdat <- min(with(soilmoisture_master_p2, date_time[portid == i & cleanorder >= 1 & !is.na(vwc)]))
  stopifnot(length(firstdat)==1 & firstdat != Inf)
  # id row in soilmoisture_master_p2
  firstrow <- with(soilmoisture_master_p2, which(portid == i & date_time == firstdat))
  # annotate
  soilmoisture_master_p2$qa_note2[firstrow] <- paste("First non-NA data recording for", i, "starts", firstdat) 
  
  # annotate last timestep where vwc recorded (whether before project end or at project stop)
  lastdat <- max(with(soilmoisture_master_p2, date_time[portid == i & cleanorder >= 1 & !is.na(vwc)]))
  stopifnot(length(lastdat)==1 & lastdat != Inf)
  # id row in soilmoisture_master_p2
  lastrow <- with(soilmoisture_master_p2, which(portid == i & date_time == lastdat))
  # annotate end of logger-port recording for project
  soilmoisture_master_p2$qa_note2[lastrow] <- paste("Last non-NA data recording for", i, "ends", as.character(lastdat)) 
}


# add note for any breaks in data (not first or last clean order, filename, vwc and raw_datorder/raw_datetime will be NA)
checkNA <- apply(select(soilmoisture_master_p2, filename:rowid), 1, function(x) all(is.na(x)))
emptyrows <- soilmoisture_master_p2[checkNA,]
# review
sort(with(emptyrows, sapply(split(vwc, portid), length)))
View(soilmoisture_master_p2[checkNA & is.na(soilmoisture_master_p2$qa_note2),])
# visualize to see if it seems sensical
ggplot() +
  geom_point(data = soilmoisture_master_p2[checkNA,], aes(date_time, 0), alpha = 0.5) +
  geom_point(data = soilmoisture_master_p2[!checkNA,], aes(date_time, 1), alpha = 0.5, col = "lightblue") +
  facet_wrap(~portid) 
# zoom in
ggplot() +
  geom_point(data = soilmoisture_master_p2[checkNA & soilmoisture_master_p2$cleanorder > 3000,], aes(cleanorder, 0), alpha = 0.5, position = position_jitter(width=0, height = 0.1)) +
  geom_point(data = soilmoisture_master_p2[!checkNA& soilmoisture_master_p2$cleanorder > 3000,], aes(cleanorder, 0), alpha = 0.5, col = "lightblue", position = position_jitter(width=0, height = 0.1)) +
  facet_wrap(~portid)
with(soilmoisture_master_p2, sapply(split(cleanorder, portid), function(x) sum(duplicated(x)))) # no cleanorder duplicated, so okay to infill as missing
# infill note as missing/break in data -- will be incorrect for loggers that don't start at cleanorder 1 but can take care of that when rbinding
soilmoisture_master_p2$qa_note2[checkNA] <- "Missing data, break in data logger recording"

# infill missing filenames
# write a loop to infill
soilmoisture_master_p2$sourcefile <- soilmoisture_master_p2$filename
# who has missing filenames?
with(distinct(soilmoisture_master_p2[c("sourcefile", "date_time", "logger")]), sapply(split(sourcefile, logger), function(x) sum(is.na(x))))


# infill filename, by first not-NA date to next not-NA (or last, as applicable)
for(p in unique(soilmoisture_master_p2$portid)){
  tempdat <- subset(soilmoisture_master_p2, portid == p)
  tempfiles <- unique(subset(tempdat, select = c(portid, date_time, cleanorder, filename, vwc, qa_note2))) %>%
    group_by(filename) %>%
    mutate(startorder = min(cleanorder)) %>%
    ungroup() %>%
    subset(!is.na(filename) & cleanorder == startorder)
  if(nrow(tempfiles) == 1){
    needsfill <- with(soilmoisture_master_p2, which(portid == p & is.na(sourcefile)))
    soilmoisture_master_p2$sourcefile[needsfill] <- tempfiles$filename
    next
  }
  # set counter
  i <- 1
  while(i < nrow(tempfiles)){
    tempstart <- tempfiles$startorder[i]
    tempend <- tempfiles$startorder[i+1] -1
    needsfill <-  with(soilmoisture_master_p2, which(portid == p & is.na(sourcefile) & cleanorder %in% tempstart:tempend))
    soilmoisture_master_p2$sourcefile[needsfill] <- tempfiles$filename[i]
    i <- i+1
  }
  # for last period infill last filename
  needsfill <- with(soilmoisture_master_p2, which(portid == p & is.na(sourcefile) & cleanorder >= tempfiles$startorder[nrow(tempfiles)]))
  # if needsfill is empty, this next line will do nothing; otherwise will infill to last entry correctly
  soilmoisture_master_p2$sourcefile[needsfill] <- tempfiles$filename[nrow(tempfiles)]
}

# review (only NAs in sourcename should be if logger didn't start recording at cleanorder 1)
summary(is.na(soilmoisture_master_p2))
View(subset(soilmoisture_master_p2, is.na(sourcefile)))
View(subset(soilmoisture_master_p2, is.na(filename) & !is.na(sourcefile)))

# verify no more than 3 distinct non-NA filenames per logger
View(distinct(soilmoisture_master_p2[c("logger", "filename", "sourcefile")])) # seems fine

# clean up colnames and append water year for rbinding
soilmoisture_master_p2_tobind <- mutate(soilmoisture_master_p2, date = date(date_time),
                                        time = as.character(str_extract(date_time, pattern = "[:digit:]{2}:.*$"))) %>%
  # drop qa_note and keep qa_note2
  subset(select = -qa_note) %>%
  left_join(date_lookup) %>%
  rename(clean_datetime = date_time, qa_note = qa_note2)




# -- STACK TREATED PERIODS 1 + 2 -----
# what's missing from p2 to bind with p1?
names(soilmoisture_master_p2_tobind)[!names(soilmoisture_master_p2_tobind) %in% names(soilmoisture_master_p1)]
names(soilmoisture_master_p1)[!names(soilmoisture_master_p1) %in% names(soilmoisture_master_p2_tobind)]

soilmoisture_master_all <- rbind(cbind(soilmoisture_master_p1, period = 1), 
                                 cbind(soilmoisture_master_p2_tobind[names(soilmoisture_master_p1)], period = 2)) %>%
  arrange(logger, port, clean_datetime, period) %>%
  mutate(stacked_rowid = rownames(.)) %>%
  #look for duplicate date_time entries, if present prioritize the one from period 1
  group_by(portid) %>%
  mutate(checkentry = duplicated(clean_datetime),
         checkentry2 = clean_datetime %in% clean_datetime[checkentry]) %>%
  ungroup()

# review
triage <- subset(soilmoisture_master_all, checkentry2) %>%
  group_by(portid, clean_datetime) %>%
  mutate(keep = ifelse(any(!is.na(vwc)), cleanorder[!is.na(vwc)], NA),
         keep2 = ifelse(is.na(keep) & any(!is.na(raw_datorder)), cleanorder[!is.na(raw_datorder)], keep),
         keep3 = ifelse(is.na(keep) & is.na(keep2) & any(!is.na(qa_note)), cleanorder[grepl("Last non-NA", qa_note)], keep2),
         keep_final = cleanorder == keep3) %>%
  ungroup()
drop <- triage[!triage$keep_final,]
# remove rows to drop
nrow(soilmoisture_master_all)
soilmoisture_master_all <- soilmoisture_master_all[!soilmoisture_master_all$stacked_rowid %in% drop$stacked_rowid,]
nrow(soilmoisture_master_all) # applied correctly

# need to clean up notes for "last non-NA" in period 1 when period 2 data are present, and first non-NA in period 2
soilmoisture_master_all <- soilmoisture_master_all %>%
  mutate(qa_note3 = ifelse(grepl("First|Last", qa_note), qa_note, NA)) %>%
  group_by(portid) %>%
  mutate(NAnote = (grepl("Last", qa_note3) & period == 1 & (any(!is.na(vwc) & period == 2))) | (grepl("First", qa_note3) & period == 2 & (any(!is.na(vwc) & period == 1))),
         # also note who is missing in p2 (just in case data are lost)
         npers = length(unique(period)),
         lastorder = ifelse(npers==1, max(cleanorder), NA),
         maxdate = ifelse(npers == 1 & cleanorder == lastorder, as.character(max(clean_datetime[!is.na(vwc)])), NA)) %>%
  ungroup() %>%
  # clean up notes
  mutate(qa_note = ifelse(NAnote, NA, qa_note),
         qa_note = ifelse(!is.na(maxdate), paste("Last non-NA data recording for", portid, maxdate), qa_note))


# reassign clean dat order and global rowid
firstdatetime <- min(soilmoisture_master_all$clean_datetime)
lastdatetime <- max(soilmoisture_master_all$clean_datetime)
all_datetime <- data.frame(clean_datetime = seq.POSIXt(firstdatetime, lastdatetime, by = "2 hours")) %>%
  mutate(global_cleanorder = 1:nrow(.))

# prep data frame to write out
soilmoisture_master_out <- left_join(soilmoisture_master_all, all_datetime) %>%
  arrange(logger, portid, clean_datetime) %>%
  mutate(cleanorder = global_cleanorder,
         # also realize I don't have a calendar yr col
         yr = year(clean_datetime)) %>%
  subset(select = c(logger:doy, yr, waterYear:sourcefile)) %>%
  data.frame() %>%
  distinct()

# screen for NAs
summary(soilmoisture_master_out)
summary(is.na(soilmoisture_master_out)) # looks good. no NAs present in columns where there shouldn't be NAs

# final pass on date-checking file source
downloadfiles_df <- distinct(subset(soilmoisture_master_out, select = c(logger, portid, port, sourcefile, clean_datetime, qa_note))) %>%
  group_by(portid, sourcefile) %>%
  mutate(firstdate = min(clean_datetime),
         lastdate = max(clean_datetime)) %>%
  subset(clean_datetime == firstdate | clean_datetime == lastdate) %>%
  ungroup() %>%
  group_by(logger, sourcefile) %>%
  mutate(logger_firstdate = min(clean_datetime),
         logger_lastdate = max(clean_datetime)) %>%
  ungroup() %>%
  mutate(firstcheck = firstdate == logger_firstdate,
         lastcheck = lastdate == logger_lastdate )

soilmoisture_master_filecheck <- left_join(soilmoisture_master_out[c("logger", "portid", "clean_datetime", "sourcefile", "cleanorder")], distinct(downloadfiles_df[c("sourcefile", "portid", "logger_firstdate", "logger_lastdate")]))
soilmoisture_master_filecheck$filecheck <- with(soilmoisture_master_filecheck, clean_datetime >= logger_firstdate & clean_datetime <= logger_lastdate)
summary(soilmoisture_master_filecheck$filecheck) # good to write out



# -- FINISHING -----
# remake plots from initial compilation loop
logger_plot_clean <- soilmoisture_master_all %>%
  ggplot(aes(dowy, vwc, col = as.factor(port))) + #date_time, 
  geom_line(alpha = 0.5) +
  ggtitle(paste0("Compost prelim plot (", Sys.Date(),  "): cleaned soil moisture (VWC),\nby logger, by water year, colored by port")) +
  scale_color_discrete(name = "port") +
  facet_grid(logger~waterYear, scales = "free_x")
logger_plot_clean

trt_plot_clean <- soilmoisture_master_all %>%
  ggplot(aes(dowy, vwc, group = portid, col = as.factor(block))) +
  geom_line(alpha = 0.5) +
  ggtitle(paste0("Compost prelim plot (", Sys.Date(),"): cleaned soil moisture (VWC),\nby water year, nutrient x drought treatment, colored by block")) +
  scale_color_discrete(name = "block") +
  facet_grid(fulltrt~waterYear, scales = "free_x")
trt_plot_clean


# write out clean data preliminary plots to if desired
ggsave(filename = paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/PrelimQA_Figures/Compost_timecorrectedVWC_bylogger.pdf"),
       plot = logger_plot_clean,
       width = 8, height = 8, units = "in")
ggsave(filename = paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/PrelimQA_Figures/Compost_timecorrectedVWC_bytreatment.pdf"),
       plot = trt_plot_clean,
       width = 8, height = 8, unit = "in")

# write out dataset
write.csv(soilmoisture_master_out, paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/SoilMoisture_compiled_time-corrected.csv"), row.names = F)
# also write out as cleaned dataset for time being until other data processing scripts made
write.csv(soilmoisture_master_out, paste0(datpath, "SoilMoisture/SoilMoisture_CleanedData/SoilMoisture_all_clean.csv"), row.names = F)
