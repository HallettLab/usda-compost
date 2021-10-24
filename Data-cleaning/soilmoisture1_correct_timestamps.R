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
adjustdates <- read.csv(paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/SoilMoisture_raw_timestampbreaks.csv"), na.strings = na_vals, strip.white = T)

# > auxiliary datasets
## treatment key
trtkey <- read.csv(paste0(datpath, "Compost_TreatmentKey.csv"), na.strings = na_vals, strip.white = T)
## soil data logger lookup table
loggerkey <- read.csv(paste0(datpath, "SoilMoisture/SoilMoisture_RawData/decagon_logger_key.csv"), na.strings = na_vals, strip.white = T)
## prepped hourly CIMIS met data for QA checks after compilation
cimis_hrly <- list.files(paste0(datpath, "CIMIS"), full.names = T) %>%
  subset(grepl(".csv",.)) %>%
  read.csv(na.strings = na_vals) 


# 2.b. Triage B3L4 by download period-----
## 2.b.1 24Apr20 -----
b3l4trt <- gsub("[0-9]", "", unique(datecheck$fulltrt[datecheck$logger == "B3L4"]))

# compare vwc across similar treatments by data order rather than date-time
## actual vwc values
subset(datecheck, grepl(b3l4trt[1], fulltrt)) %>% # trackseq > 2250
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_line(alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates, grepl(b3l4trt[1], fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  ggtitle(paste("troubleshoot B3L4", b3l4trt[1], "date jump (dotted vert line = break)")) +
  geom_smooth(aes(fill = portid)) +
  facet_grid(logger~.)

## t2-t1 (diffed) values
subset(datecheck, grepl(b3l4trt[1], fulltrt)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(datorder, diffvwc, col = portid)) +
  geom_vline(data = mutate(subset(adjustdates, grepl(b3l4trt[1], fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  geom_line(alpha = 0.4) +
  ggtitle(paste("troubleshoot B3L4", b3l4trt[1], "date jump (dotted vert line = break)")) +
  #geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  facet_grid(logger~., scales = "free_x", space = "free_x")


subset(datecheck, grepl(b3l4trt[2], fulltrt)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_vline(data = mutate(subset(adjustdates, grepl(b3l4trt[2], fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  geom_line(alpha = 0.4) +
  geom_smooth(aes(fill = portid)) +
  ggtitle(paste("troubleshoot B3L4", b3l4trt[2], "date jump (dotted vert line = break)")) +
  facet_grid(logger~.)
# > data gap after break, shift series so end lines up with end timestamp collected (and check spikes align with B2L1)

subset(datecheck, grepl(b3l4trt[2], fulltrt)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(datorder, diffvwc, col = portid)) +
  geom_vline(data = mutate(subset(adjustdates, grepl(b3l4trt[2], fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  geom_line(alpha = 0.4) +
  ggtitle(paste("troubleshoot B3L4", b3l4trt[2], "date jump (dotted vert line = break)")) +
  #geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  facet_grid(logger~., scales = "free_x", space = "free_x")


# rule 1: if time of day from end of bad file 2 hrs before time of day of following file, assign new times backwards end to last bad date break
# rule 2: look for moisture spike alignment with comparison


b3l4dates <- subset(datecheck, logger == "B3L4") %>%
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
  left_join(distinct(soilmoisture_clean[c("date_time", "cleanorder")])) %>%
  mutate(difforder = cleanorder - datorder)
unique(b3l4dates$difforder)
# difference was 1 before date break, then 235 off from clean order after (so gap of 234 2-hr intervals.. about 19.5 days)
# then 236 in 2021 data download (another 20d gap)

# all ports have same time breaks
# > 1 is NA to infill for missing timestep
# > other is a date sequence that needs adjustment -- end of sequence seems to match?

# compare differences against ports with good data to get a sense of typical difference between logger-ports 
good_b3l4 <- subset(datecheck, grepl(b3l4trt[1], fulltrt)) %>%
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


## go by download file since there are multiple timestamp corrections needed
# Apr 24 2020
searchdata <- subset(datecheck, logger == "B3L4" & datorder >= with(b3l4dates, datorder[grepl("correct timestamp", qa_note)]-1) & filename == with(b3l4dates, filename[grepl("correct timestamp", qa_note)])) %>%
  group_by(datorder) %>%
  mutate(nobs = sum(wetup, na.rm = T)) %>%
  ungroup() %>%
  subset(datorder >= (min(datorder[wetup == 1 & nobs > 1], na.rm = T)-20) & datorder <= (min(datorder[wetup == 1 & nobs > 1], na.rm = T)+80))

refdata <- subset(datecheck, grepl(str_flatten(b3l4trt, collapse = "|"), fulltrt)  & logger != "B3L4" & datorder >= with(b3l4dates, datorder[grepl("correct timestamp", qa_note)]-1) & grepl("Apr20", filename)) %>% # loggers have different april timestamps as b3l4
  group_by(date_time) %>%
  mutate(nobs = sum(wetup, na.rm = T)) %>%
  ungroup() %>%
  subset(date_time >= (min(date_time[wetup == 1 & nobs > 1], na.rm = T)-as.difftime(20*2, units = "hours")) & date_time <= (min(date_time[wetup == 1 & nobs > 1], na.rm = T)+as.difftime(80*2, units = "hours"))) %>%
  # join clean order
  left_join(distinct(soilmoisture_clean[c("date_time", "cleanorder")]))


# Apr 24 2020
unique(refdata$timeid)
unique(searchdata$timeid)
searchdata$timeid2 <- searchdata$timeid+259.02 # test adjustment
unique(searchdata$timeid2)
plot_grid(ggplot(refdata, aes(timeid, vwc, col = portid)) +
            geom_line() +
            labs(subtitle = "reference logger data") +
            facet_grid(~gsub("[0-9]+", "", fulltrt)),
          ggplot(searchdata, aes(timeid2, vwc, col = portid)) +
            geom_line() +
            labs(subtitle = "targer logger to fix") +
            facet_grid(~gsub("[0-9]+", "", fulltrt)),
          nrow = 2)
# overlay for comparison
ggplot(refdata, aes(timeid, vwc, col = portid)) +
  geom_line(data = searchdata, aes(timeid2, vwc, group = portid), col = "grey50", lwd = 1.5, alpha = 0.6) +
  geom_line() +
  labs(subtitle = "Does adjustment overlap? grey line = target; colored = reference") +
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
# overlay
ggplot(refdata, aes(cleanorder, vwc, col = portid)) +
  geom_line(data = searchdata, aes(datorder+235, vwc, group = portid), col = "grey50", lwd = 1.5, alpha = 0.6) +
  geom_line() +
  facet_grid(~gsub("[0-9]+", "", fulltrt))



# what happens if subtract 1 backwards in clean order from first good timestamp after time break?
startorder <- with(b3l4dates, cleanorder[grepl("last run", qa_note)])
startorder <- 6262
b3l4test <- left_join(datecheck, distinct(soilmoisture_clean[c("date_time", "cleanorder")])) %>%
  subset(logger == "B3L4") %>%
  arrange(portid, desc(cleanorder)) %>%
  group_by(portid, intevent) %>%
  mutate(cleanorder2 = ifelse(intevent == 3, datorder + 235, cleanorder), 
         #cleanorder2 = ifelse(intevent == 3, datorder + with(b3l4dates, difforder[grepl("last run", qa_note)]), cleanorder), 
         cleanorder3 = ifelse(intevent == 3, seq(startorder -length(vwc), startorder-1, 1),cleanorder),
         diffcheck = cleanorder2 - datorder) %>%
  ungroup() %>%
  rename(raw_datetime = date_time)


ggplot(b3l4test, aes(cleanorder2, vwc, col = port, group = portid)) +
  geom_point(data = subset(soilmoisture_clean, grepl(str_flatten(b3l4trt, collapse = "|"), fulltrt) & logger != "B3L4"), aes(cleanorder, vwc, group = portid, col = port)) +
  geom_point() +
  geom_vline(data = mutate(subset(adjustdates, grepl(str_flatten(b3l4trt, collapse = "|"), fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  facet_grid(logger~fulltrt)

ggplot(subset(b3l4test, intevent == 3), aes(cleanorder2, vwc, col = portid, group = portid)) +
  geom_line(data = subset(soilmoisture_clean, grepl(str_flatten(b3l4trt, collapse = "|"), fulltrt) & logger != "B3L4" & cleanorder %in% b3l4test$cleanorder2[b3l4test$intevent == 3]), aes(cleanorder, vwc, group = portid), col = "grey50") +
  geom_line() +
  facet_grid(.~fulltrt) # looks okay now


# correct B3L4 in soilmoisture_clean
clean_b3l4 <- subset(soilmoisture_clean, grepl("B3L4", portid)) %>%
  select(date_time:cleanorder, filename) %>%
  left_join(select(b3l4test, portid, cleanorder3, vwc, raw_datetime, timeinterval, timediff, intevent, qa_note, datorder), by = c("portid", "cleanorder" = "cleanorder3")) %>%
  left_join(distinct(soilmoisture_all, portid, logger, port, block, plotid, fulltrt))

# create temp master with timestamp-corrected data (NAs not yet addressed)..
soilmoisture_master <- subset(soilmoisture_clean, !grepl("B3L4", portid), c(date_time, portid, cleanorder, filename, vwc)) %>%
  left_join(distinct(select(soilmoisture_all, portid, logger, port, block, plotid, fulltrt))) %>%
  select(names(soilmoisture_clean)) %>%
  rbind(clean_b3l4[names(.)])

# check out NAs by logger
with(soilmoisture_master, sapply(split(vwc, portid), function(x) summary(is.na(x)))) #small numbers are just missing data for 2-hr intervals
# clean up
rm(searchdata, refdata, b3l4dates, b3l4trt, b3l4test, startorder)



# 2.c. Triage B2L4 24Apr20 -----
b2l4trt <- gsub("[0-9]", "", unique(datecheck$fulltrt[datecheck$logger == "B2L4"]))

subset(datecheck, grepl(b2l4trt[1], fulltrt) & grepl("Apr20", filename)) %>% # & grepl("Apr20", filename) & trackseq > 2000 trackseq > 2250
  ggplot(aes(datorder, vwc, col = paste(logger, port, sep = "_"))) +
  geom_line(aes(group = paste(logger, port, sep = "_")), alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates, grepl("B2L4", portid)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  #geom_vline(data = subset(datecheck, logger == "B2L4" & timediff != 2 & trackseq > 2000), aes(xintercept = trackseq), lty = 2) +
  geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  ggtitle(paste("troubleshoot B2L4 date jump", b2l4trt[1], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l4trt[1])
#facet_grid(.~filemo, scales = "free_x", space = "free_x")

subset(datecheck, grepl(b2l4trt[1], fulltrt) & grepl("Apr20", filename)) %>% # & grepl("Apr20", filename) & datorder > 2000 datorder > 2250
  ggplot(aes(datorder, diffvwc, col = paste(logger, port, sep = "_"))) +
  geom_line(aes(group = paste(logger, port, sep = "_")), alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates, grepl("B2L4", portid)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  #geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  ggtitle(paste("troubleshoot B2L4 date jump", b2l4trt[1], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l4trt[1])

subset(datecheck, grepl(b2l4trt[2], fulltrt) & grepl("Apr20", filename)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(datorder, vwc, col = paste(logger, port, sep = "_"), group = paste(logger, port, sep = "_"))) +
  geom_line(aes(group = paste(logger, port, sep = "_")), alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates, grepl("B2L4", portid)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  ggtitle(paste("troubleshoot B2L4 date jump", b2l4trt[2], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l4trt[2])
#facet_grid(.~filemo, scales = "free_x", space = "free_x")

subset(datecheck, grepl(b2l4trt[2], fulltrt) & grepl("Apr20", filename)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(datorder, diffvwc, col = paste(logger, port, sep = "_"), group = paste(logger, port, sep = "_"))) +
  geom_line(aes(group = paste(logger, port, sep = "_")), alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates, grepl("B2L4", portid)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  #geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  ggtitle(paste("troubleshoot B2L4 date jump", b2l4trt[2], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l4trt[1])

# > slight data gap (spikes later in sequence should line up on same day [rain event])
b2l4dates <- subset(datecheck, grepl("B2L4", portid)) %>%
  subset(datorder %in% unique(c(datorder[!is.na(qa_note)], datorder[!is.na(qa_note)]-1))) %>%
  arrange(portid, datorder) %>%
  mutate(port = as.numeric(substr(portid, 6,7)),
         timeinterval = as.numeric(timeinterval),
         hrdiff = ifelse(timeinterval==0, 24-lag(timeinterval), 
                         ifelse(lag(timeinterval) > timeinterval, 24+ (timeinterval  - lag(timeinterval)), timeinterval- lag(timeinterval))),
         hrdiff = ifelse(is.na(qa_note), NA, hrdiff)) %>%
  distinct(logger, portid, date_time, filename, datorder, intevent, qa_note, hrdiff) %>%
  # add clean dates
  left_join(distinct(soilmoisture_clean[c("date_time", "cleanorder")])) %>%
  mutate(difforder = cleanorder - datorder)

# visually assess
searchdata_b2l4 <- subset(datecheck, logger == "B2L4" & datorder >= 4106 & grepl("Apr20", filename)) %>%
  #subset(wetup == 1) %>%
  group_by(datorder) %>%
  mutate(nobs = sum(wetup, na.rm = T)) %>%
  ungroup() %>%
  subset(datorder >= (min(datorder[wetup == 1 & nobs >= 2], na.rm = T)-20) & datorder <= (min(datorder[wetup == 1 & nobs >= 2], na.rm = T)+80))

refdata_b2l4 <- subset(datecheck, grepl(str_flatten(b2l4trt, collapse = "|"), fulltrt)  & logger != "B2L4" & datorder >= 4106 & grepl("Apr20", filename)) %>%
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

max(datecheck$datorder[datecheck$logger == "B2L4" & grepl("Apr20", datecheck$filename)])
max(datecheck$datorder[datecheck$logger == "B1L1" & grepl("Apr20", datecheck$filename)]) # if add 16 to B2L2 logger, get max points in download period

# pull first period that needs timestamp adjustment and doesn't have any spikes to crosscheck -- would +16 datorder match there too?
searchdata_b2l4_p1 <-  subset(datecheck, logger == "B2L4" & datorder %in% c(4100:4384) & grepl("Apr20", filename)) %>%
  mutate(datorder2 = ifelse(datorder < 4107, datorder, datorder + 12))
refdata_b2l4_p1 <- subset(datecheck, grepl(str_flatten(b2l4trt,collapse = "|"), fulltrt) & logger != "B2L4" & grepl("Apr20", filename)) %>% # #block == 2 & nut_trt == "C" &
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
searchdata_b2l4_p2 <-  subset(datecheck, logger == "B2L4" & datorder %in% c(4130:4600) & grepl("Apr20", filename)) %>%
  mutate(datorder2 = ifelse(datorder < 4385, datorder+12, datorder + 16)) #14

refdata_b2l4_p2 <- datecheck %>%
  #subset(grepl(str_flatten(b2l4trt,collapse = "|"), fulltrt) & logger != "B2L4" & grepl("Apr20", filename)) %>% 
  subset(block == 2 & nut_trt == "C" & logger != "B2L4" & grepl("Apr20", filename)) %>% # #block == 2 & nut_trt == "C" &
  subset(datorder %in% c(4142:4616))

ggplot(refdata_b2l4_p2, aes(datorder, vwc, col = portid)) +
  geom_line(data = subset(searchdata_b2l4_p2, portid != "B2L4_3"), aes(datorder2, vwc, group = portid), col = "grey50", alpha = 0.6) +
  geom_line() +
  geom_vline(aes(xintercept = 4358))


# try infilling backwards and see how that looks, go by interval event
# what happens if subtract 1 backwards in cleanorder from first good timestamp after time break?
startorder_b2l4 <- with(b2l4dates, cleanorder[grepl("last run", qa_note)])
b2l4test <- left_join(datecheck, distinct(soilmoisture_clean[c("date_time", "cleanorder")])) %>%
  subset(logger == "B2L4") %>%
  arrange(portid, desc(cleanorder)) %>%
  group_by(portid, intevent) %>%
  mutate(cleanorder2 = ifelse(intevent == 4, datorder + with(b2l4dates, difforder[grepl("last run", qa_note)]), cleanorder), 
         cleanorder3 = ifelse(intevent == 4, seq(startorder_b2l4 -length(vwc), startorder_b2l4-1, 1),cleanorder),
         diffcheck = cleanorder2 - datorder)

ggplot(b2l4test, aes(cleanorder2, vwc, col = port, group = portid)) +
  geom_point(data = subset(soilmoisture_clean, grepl(str_flatten(b2l4trt, collapse = "|"), fulltrt) & logger != "B2L4"), aes(cleanorder, vwc, group = portid, col = port)) +
  geom_point() +
  geom_vline(data = mutate(subset(adjustdates, grepl(str_flatten(b2l4trt, collapse = "|"), fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  facet_grid(logger~gsub("[0-9]+", "", fulltrt))

ggplot(subset(b2l4test, intevent == 4), aes(cleanorder2, vwc, col = portid, group = portid)) +
  geom_line(data = subset(soilmoisture_clean, grepl(str_flatten(b2l4trt, collapse = "|"), fulltrt) & logger != "B2L4" & cleanorder %in% b2l4test$cleanorder2[b2l4test$intevent == 4]), aes(cleanorder, vwc, group = portid), col = "grey50") +
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
  geom_point(data = subset(soilmoisture_clean, grepl(str_flatten(b2l4trt, collapse = "|"), fulltrt) & logger != "B2L4" & cleanorder %in% c(3500:4500)), aes(cleanorder, vwc, group = portid, col = port)) +
  geom_point() +
  #geom_vline(data = mutate(subset(adjustdates, grepl(str_flatten(b2l4trt, collapse = "|"), fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  facet_grid(logger~gsub("[0-9]+", "", fulltrt))

ggplot(subset(b2l4test, cleanorder3 %in% c(3500:4500) & portid != "B2L4_3"), aes(cleanorder3, vwc, col = is.na(cleanorder2), group = portid)) +
  geom_line(data = subset(soilmoisture_clean, grepl(str_flatten(b2l4trt, collapse = "|"), fulltrt) & logger != "B2L4" & cleanorder %in% c(3500:4500) ), aes(cleanorder, vwc, group = portid), col = "grey50") + #cleanorder %in% b2l4test$cleanorder3[b2l4test$intevent == 3]
  geom_line() +
  facet_grid(.~gsub("[0-9]+", "", fulltrt)) # looks good -- +12 is correct (best fit)

# focus on interval 3 to be sure
ggplot(subset(b2l4test, intevent == 3 & portid != "B2L4_3"), aes(cleanorder3, vwc, col = is.na(cleanorder2), group = portid)) +
  geom_line(data = subset(soilmoisture_clean, grepl(str_flatten(b2l4trt, collapse = "|"), fulltrt) & logger != "B2L4" & cleanorder %in% b2l4test$cleanorder3[b2l4test$intevent == 3]), aes(cleanorder, vwc, group = portid), col = "grey50") + #
  geom_line() +
  facet_grid(.~gsub("[0-9]+", "", fulltrt)) # looks ok


# correct B2L4 in soilmoisture_clean
clean_b2l4 <- subset(soilmoisture_clean, grepl("B2L4", portid)) %>%
  select(date_time:cleanorder, filename) %>%
  left_join(select(b2l4test, portid, cleanorder3, vwc, raw_datetime, timeinterval, timediff, intevent, qa_note, datorder), by = c("portid", "cleanorder" = "cleanorder3")) %>%
  left_join(distinct(soilmoisture_all, portid, logger, port, block, plotid, fulltrt))

# create temp master with timestamp-corrected data (NAs not yet addressed)..
soilmoisture_master <- subset(soilmoisture_master, !grepl("B2L4", portid)) %>%
  rbind(clean_b2l4[names(.)])
# check out NAs by logger
with(soilmoisture_master, sapply(split(vwc, portid), function(x) summary(is.na(x))))
# check original for NA count
with(soilmoisture_all[soilmoisture_all$logger == "B2L4",], sapply(split(vwc, portid), function(x) summary(is.na(x))))
#> high count for B2L4_3 was already there, not a matter of poor join 

# clean up
rm(searchdata_b2l4, refdata_b2l4, searchdata_b2l4_p1, searchdata_b2l4_p2, b2l4dates, b2l4trt, b2l4test)




# 2.d. Triage B2L2 24Apr20 -----
b2l2trt <- gsub("[0-9]", "", unique(datecheck$fulltrt[datecheck$logger == "B2L2"]))

subset(datecheck, grepl(b2l2trt[1], fulltrt)) %>% # datorder > 2250
  ggplot(aes(datorder, vwc, col = port, group = portid)) +
  geom_line(alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates, grepl("B2L2", portid)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note, col = as.numeric(substr(portid,6,6)))) +
  ggtitle(paste("troubleshoot B2L2 date jump", b2l2trt[1], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l2trt[1])

subset(datecheck, grepl(b2l2trt[1], fulltrt)) %>% # datorder > 2250
  ggplot(aes(datorder, diffvwc, col = port, group = paste(logger, port, sep = "_"))) +
  geom_line(alpha = 0.4) +
  geom_vline(data = mutate(subset(adjustdates, grepl("B2L2", portid)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note, col = as.numeric(substr(portid,6,6)))) +
  ggtitle(paste("troubleshoot B2L2 date jump", b2l2trt[1], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l2trt[1])

subset(datecheck, grepl(b2l2trt[2], fulltrt)) %>%  #&  & datorder > 2250
  ggplot(aes(datorder, vwc, col = port, group = portid)) +
  geom_vline(data = mutate(subset(adjustdates, grepl(b2l2trt[2], fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note, col = as.numeric(substr(portid,6,6)))) +
  #geom_vline(aes(xintercept = max(unique(with(datecheck, datorder[logger == "B2L2" & as.numeric(timediff) != 2])), na.rm = T)), lty = 2) +
  geom_line(alpha = 0.4) +
  ggtitle(paste("troubleshoot B2L2 date jump", b2l2trt[2], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l2trt[2])

subset(datecheck, grepl(b2l2trt[2], fulltrt) & grepl("Apr20", filename)) %>%  #&  & datorder > 2250
  ggplot(aes(datorder, diffvwc, col = paste(logger, port, sep = "_"), group = port)) +
  geom_vline(data = mutate(subset(adjustdates, grepl("B2L2", portid)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  geom_line(alpha = 0.4) +
  ggtitle(paste("troubleshoot B2L2 date jump", b2l2trt[2], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l2trt[2])

# > think gap happened after break, and end of series for B2L2 should align with end date -- can gauge by spike in soil moisture in mid-3000s
b2l2dates <- subset(datecheck, grepl("B2L2", portid)) %>%
  subset(datorder %in% unique(c(datorder[!is.na(qa_note)], datorder[!is.na(qa_note)]-1))) %>%
  arrange(portid, datorder) %>%
  mutate(port = as.numeric(substr(portid, 6,7)),
         timeinterval = as.numeric(timeinterval),
         hrdiff = ifelse(timeinterval==0, 24-lag(timeinterval), 
                         ifelse(lag(timeinterval) > timeinterval, 24+ (timeinterval  - lag(timeinterval)), timeinterval- lag(timeinterval))),
         hrdiff = ifelse(is.na(qa_note), NA, hrdiff)) %>%
  distinct(logger, portid, date_time, filename, datorder, intevent, qa_note, hrdiff) %>%
  # add clean dates
  left_join(distinct(soilmoisture_clean[c("date_time", "cleanorder")])) %>%
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
searchdata_b2l2.1 <- subset(datecheck, portid == "B2L2_1" & datorder >= with(b2l2dates_1, datorder[intevent == 5 & grepl("timestamp", qa_note)])-1) %>% # filename == with(b3l4dates, filename[grepl("timestamp", qa_note)])
  group_by(datorder) %>%
  mutate(nobs = sum(wetup, na.rm = T)) %>%
  ungroup()
#subset(datorder >= (min(datorder[wetup == 1], na.rm = T)-20) & datorder <= (min(datorder[wetup == 1], na.rm = T)+80))

refdata_b2l2 <- subset(datecheck, grepl(gsub("[0-9]", "", searchdata_b2l2.1$fulltrt[1]), fulltrt)  & logger != "B2L2" & datorder >= with(b2l2dates_1, datorder[intevent == 5 & grepl("timestamp", qa_note)])-1) %>% # & grepl("Apr20", filename)) %>% # loggers have different april timestamps as b3l4
  group_by(date_time) %>%
  mutate(nobs = sum(wetup, na.rm = T)) %>%
  ungroup()
#subset(date_time >= (min(date_time[wetup == 1 & nobs > 1], na.rm = T)-as.difftime(20*2, units = "hours")) & date_time <= (min(date_time[wetup == 1 & nobs > 1], na.rm = T)+as.difftime(80*2, units = "hours"))) %>%
# join clean order
#left_join(distinct(soilmoisture_clean[c("date_time", "cleanorder")]))
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

b2l2.1_test <- left_join(datecheck, distinct(soilmoisture_clean[c("date_time", "cleanorder")])) %>%
  subset(portid == "B2L2_1") %>%
  group_by(portid, intevent) %>%
  mutate(cleanorder2 = ifelse(intevent == 6, datorder + with(b2l2dates_1, difforder[grepl("last run", qa_note)]), cleanorder), 
         diffcheck = cleanorder2 - datorder)
# > doesn't work to infill cleanorder backwards.. there is a gap between final interval run and penultimate..

ggplot(b2l2.1_test, aes(cleanorder2, vwc, group = intevent)) +
  geom_line(data = subset(soilmoisture_clean, grepl(gsub("[0-9]", "", b2l2.1_test$fulltrt[1]), fulltrt) & logger != "B2L2"), aes(cleanorder, vwc, group = portid), alpha = 0.3, col = "grey50") + #col = port
  #geom_point() +
  geom_line(col = "orchid") +
  geom_vline(data = mutate(subset(adjustdates, grepl(gsub("[0-9]", "", b2l2.1_test$fulltrt[1]), fulltrt)), logger = substr(portid, 1,4)), aes(xintercept = datorder, lty = qa_note)) +
  facet_grid(logger~gsub("[0-9]+", "", fulltrt))

# .. this is going to be very manual/visual matching..


# look at last period for matching since middle period doesn't overlap with any spikes
searchdata_b2l2_end <- subset(datecheck, portid == "B2L2_1" & datorder > 4392 & grepl("Apr20", filename)) %>%
  group_by(datorder) %>%
  mutate(nobs = sum(wetup, na.rm = T)) %>%
  ungroup()
#subset(datorder >= (min(datorder[wetup == 1 & nobs > 1], na.rm = T)-20) & datorder <=  (min(datorder[wetup == 1 & nobs > 1], na.rm = T)+80))

nrows <- max(searchdata_b2l2_end$datorder) - min(searchdata_b2l2_end$datorder)
tempend <- max(datecheck$datorder[grepl("Apr20", datecheck$filename)])
refdata_b2l2_end <- #subset(datecheck, grepl(str_flatten(b2l2trt, collapse = "|"), fulltrt)  & logger != "B2L2" & datorder >= (tempend-nrows) & grepl("Apr20", filename)) %>%
  subset(datecheck, grepl("FW", fulltrt)  & logger != "B2L2" & datorder >= (tempend-nrows) & grepl("Apr20", filename)) %>%
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
searchdata_b2l2_early <- subset(datecheck, logger == "B2L2" & intevent %in% unique(with(b2l2dates, intevent[nports == 5])) & datorder >= unique(min(b2l2dates$datorder))-50) %>%
  #subset(datecheck, logger == "B2L2" & datorder < 4394 & datorder >= unique(min(b2l2dates$datorder))-150) %>%
  mutate(datorder2 = ifelse(intevent %in% c(2,3), datorder+3, datorder))

refdata_b2l2_early <- subset(datecheck, logger != "B2L2" & grepl(str_flatten(b2l2trt, collapse = "|"), fulltrt) &  datorder %in% searchdata_b2l2_early$datorder)

ggplot(data = refdata_b2l2_early, aes(datorder, vwc, group = portid), col = "grey50", alpha = 0.5) +
  geom_line() +
  geom_line(data = searchdata_b2l2_early, aes(datorder2, vwc, group = paste(portid, intevent), col = as.factor(port), lty = as.factor(intevent))) +
  facet_grid(. ~ gsub("[0-9]", "", fulltrt))
#theme(legend.position = "none")

# join correct cleanorder and date time to searchdata by datorder2
b2l2early_test <- subset(datecheck, logger == "B2L2" & intevent %in% unique(with(b2l2dates, intevent[nports == 5]))) %>%
  #subset(datecheck, logger == "B2L2" & datorder < 4394 & datorder >= unique(min(b2l2dates$datorder))-150) %>%
  mutate(datorder2 = ifelse(intevent %in% c(2,3), datorder+3, datorder)) %>%
  mutate(cleanorder = datorder2+1) %>%
  rename(raw_datetime = date_time) %>%
  left_join(subset(soilmoisture_clean, grepl("B2L2", portid), c(date_time, portid, cleanorder)), by = c("cleanorder", "portid"))

# correct B3L4 in soilmoisture_clean
clean_b2l2_early <- subset(soilmoisture_clean, grepl("B2L2", portid)) %>%
  select(date_time:cleanorder, filename) %>%
  left_join(select(b2l2early_test, portid, cleanorder, vwc, raw_datetime, timeinterval, timediff, intevent, qa_note), by = c("portid", "cleanorder")) %>%
  left_join(distinct(soilmoisture_all, portid, logger, port, block, fulltrt))

ggplot(b2l2early_test, aes(cleanorder, vwc, group = portid, col= as.factor(port))) +
  geom_line(data = subset(soilmoisture_master, grepl(str_flatten(b2l2trt, collapse = "|"), fulltrt) & logger != "B2L2"), aes(cleanorder, vwc, group = portid), alpha = 0.3, lwd = 1.5) +
  geom_line() +
  facet_grid(block ~ gsub("[0-9]", "", fulltrt)) # looks okay

# clean up environment and correct by port
rm(refdata_b2l2, refdata_b2l2_early, refdata_b2l2_end, searchdata_b2l2_early, searchdata_b2l2_end, 
   nrows, tempend, startorder_b2l2.1, b2l2.1_test)



# 2.e. Triage B2L2 by port -----------
# create refdata
refdata_b2l2 <- subset(soilmoisture_master, grepl(str_flatten(b2l2trt, "|"), fulltrt) & logger != "B2L2")

## 2.e.1. B2L2_1 -------
# press on with corrections for individual ports after interval 3 ..
searchdata_b2l2.1 <- subset(datecheck, portid == "B2L2_1") %>%
  mutate(datorder2 = ifelse(intevent == 1, datorder, 
                            ifelse(intevent == 2, datorder+3, #3 
                                   ifelse(intevent == 3, datorder + 5, #5
                                          # interval 4 is an NA -- adjust by amount of intervals off
                                          ifelse(intevent == 4, datorder + 5 + ((b2l2dates_1$hrdiff[grepl("NA", b2l2dates_1$qa_note)]-2)/2),
                                                 ifelse(intevent == max(intevent)-1,datorder + 396,
                                                        ifelse(intevent == max(intevent), datorder + with(b2l2dates_1, difforder[grepl("last run", qa_note)]),NA)))))),
         cleanorder = ifelse(intevent == max(intevent), datorder2, datorder2+1)) %>%
  rename(raw_datetime = date_time) %>%
  left_join(distinct(soilmoisture_clean[c("date_time", "cleanorder")])) %>%
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
searchdata_b2l2.2 <- subset(datecheck, portid == "B2L2_2") %>%
  mutate(datorder2 = ifelse(intevent == 1, datorder, 
                            ifelse(intevent == 2, datorder+3, 
                                   ifelse(intevent == 3, datorder + 5,
                                          ifelse(intevent == max(intevent)-1,datorder + 396,
                                                 ifelse(intevent == max(intevent), datorder + with(b2l2dates_2, difforder[grepl("last run", qa_note)]),NA))))),
         cleanorder = ifelse(intevent == max(intevent), datorder2, datorder2+min(b2l2dates_2$difforder, na.rm = T))) %>%
  rename(raw_datetime = date_time) %>%
  left_join(distinct(soilmoisture_clean[c("date_time", "cleanorder")])) %>%
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
  geom_line(data = subset(soilmoisture_master, grepl(gsub("[0-9]", "", searchdata_b2l2.2$fulltrt[1]), fulltrt) & logger != "B2L2" & cleanorder >= min(searchdata_b2l2.2$cleanorder[searchdata_b2l2.2$intevent == 2]) & cleanorder < max(searchdata_b2l2.2$cleanorder[searchdata_b2l2.2$intevent == 3])), aes(cleanorder, vwc, group = portid), alpha = 0.3) + # grepl(gsub("[0-9]", "", searchdata_b2l2.2$fulltrt[1]), fulltrt) & logger != "B2L2"
  geom_line(aes(col = as.factor(intevent))) +
  scale_x_continuous(breaks = seq(4100,4400, 20)) +
  facet_grid(block~., scales = "free_y")


## 2.e.3. B2L2_3 ----- 
# port 3 has 5 runs -- 5th is final, 2, 3, 4 are timestamp breaks... similar to port 2
searchdata_b2l2.3 <- subset(datecheck, portid == "B2L2_3") %>%
  mutate(datorder2 = ifelse(intevent == 1, datorder, 
                            ifelse(intevent == 2, datorder+3, 
                                   ifelse(intevent == 3, datorder + 5,
                                          ifelse(intevent == max(intevent)-1,datorder + 396,
                                                 ifelse(intevent == max(intevent), datorder + with(b2l2dates_3, difforder[grepl("last run", qa_note)]),NA))))),
         cleanorder = ifelse(intevent == max(intevent), datorder2, datorder2+min(b2l2dates_3$difforder, na.rm = T))) %>%
  rename(raw_datetime = date_time) %>%
  left_join(distinct(soilmoisture_clean[c("date_time", "cleanorder")])) %>%
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
searchdata_b2l2.4 <- subset(datecheck, portid == "B2L2_4") %>%
  mutate(datorder2 = ifelse(intevent == 1, datorder, 
                            ifelse(intevent == 2, datorder+3, #3 
                                   ifelse(intevent == 3, datorder + 5, #5
                                          # interval 4 is an NA -- adjust by amount of intervals off
                                          ifelse(intevent == 4, datorder + 5 + ((b2l2dates_4$hrdiff[grepl("NA", b2l2dates_4$qa_note)]-2)/2),
                                                 ifelse(intevent == max(intevent)-1,datorder + 396,
                                                        ifelse(intevent == max(intevent), datorder + with(b2l2dates_4, difforder[grepl("last run", qa_note)]),NA)))))),
         cleanorder = ifelse(intevent == max(intevent), datorder2, datorder2+ min(b2l2dates_4$difforder, na.rm = T))) %>%
  rename(raw_datetime = date_time) %>%
  left_join(distinct(soilmoisture_clean[c("date_time", "cleanorder")])) %>%
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
searchdata_b2l2.5 <- subset(datecheck, portid == "B2L2_5") %>%
  mutate(datorder2 = ifelse(intevent == 1, datorder, 
                            ifelse(intevent == 2, datorder+3, #3 
                                   ifelse(intevent == 3, datorder + 5, #5
                                          # interval 4 is an NA -- adjust by amount of intervals off
                                          ifelse(intevent == 4, datorder + 5 + ((b2l2dates_5$hrdiff[grepl("NA", b2l2dates_5$qa_note)]-2)/2),
                                                 ifelse(intevent == max(intevent)-1,datorder + 396,
                                                        ifelse(intevent == max(intevent), datorder + with(b2l2dates_5, difforder[grepl("last run", qa_note)]),NA)))))),
         cleanorder = ifelse(intevent == max(intevent), datorder2, datorder2+ min(b2l2dates_5$difforder, na.rm = T))) %>%
  rename(raw_datetime = date_time) %>%
  left_join(distinct(soilmoisture_clean[c("date_time", "cleanorder")])) %>%
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
masterb2l2_clean <- subset(soilmoisture_clean, grepl("B2L2", portid), c(date_time, cleanorder, portid)) %>%
  left_join(masterb2l2_clean)


# plot everything
ggplot(masterb2l2_clean, aes(cleanorder, vwc, group = paste(portid, intevent), col = as.factor(port))) +
  geom_line() +
  ggtitle("Timestamp-corrected B2L2 soil moisture data") +
  facet_grid(gsub("[0-9]", "", fulltrt)~., scales = "free_y")



# -- 3. COMPILE CLEAN DATA ------
# add qa notes for NAs
# want cleanorder, correct timestamp, treatment info, water year info, vwc, qa notes and raw timestamp, and rawdatorder

names(masterb2l2_clean)
names(clean_b2l4)
names(clean_b3l4)
names(soilmoisture_master) # has clean b3l4 and b2l4, but not clean b2l2
names(soilmoisture_clean)
names(soilmoisture_all)

# remake soilmoisture_master with all colnames from soilmoisture_all + 2 qa cols (qa note and raw timestamp)
# > start with loggers than aren't b3l4, b2l4 and b2l2
# > also reattach logger key and trt key info

soilmoisture_master <- select(soilmoisture_clean, date_time:cleanorder, filename, vwc) %>%
  subset(!grepl("B3L4|B2L4|B2L2", portid)) %>%
  # add rawdatorder back in 
  left_join(subset(soilmoisture_all, !grepl("B3L4|B2L4|B2L2", portid), c(date_time, portid, datorder))) %>%
  #summary(is.na(soilmoisture_master))
  #View(subset(soilmoisture_master, is.na(datorder)))
  rename(clean_datetime = date_time) %>%
  # add original datetime in so data user can see where was recorded
  left_join(subset(soilmoisture_all, !grepl("B3L4|B2L4|B2L2", portid), c(date_time, portid, datorder))) %>%
  rename(raw_datetime = date_time) %>%
  # add in col for qa_note
  mutate(qa_note = NA) %>%
  rbind(rename(clean_b3l4, clean_datetime = date_time)[names(.)]) %>%
  rbind(rename(clean_b2l4, clean_datetime = date_time)[names(.)]) %>%
  rbind(rename(masterb2l2_clean, clean_datetime = date_time)[names(.)]) %>%
  rename(raw_datorder = datorder) %>%
  mutate(qa_note2 = NA)

# infill qa_notes
for(i in unique(adjustdates$portid)){
  tempadjust <- subset(adjustdates, portid == i)
  for(r in 1:nrow(tempadjust)){
    # append qa note dependent on what done
    ## NA infill
    if(grepl("NA infill", tempadjust$qa_note[r])){
      # specify start and end points for note
      ## end = where noted in qa_note
      end_raworder <- tempadjust$datorder[r]
      tempend <- with(soilmoisture_master, which(portid == i & raw_datorder == tempadjust$datorder[r]))
      #grab hr diff
      temphr <- (as.numeric(tempadjust$timediff[r])-2)/2 #this is the number of timestep adjustments from expected 2hr interval
      ## start = work backwards
      tempstart <- tempend - temphr
      # note only goes to spot before where indicated
      soilmoisture_master$qa_note2[tempstart:(tempend-1)] <- paste("Data missing, NA added, no timestamp recorded by", i, "for", temphr, "2hr-timesteps")
    }
    ## timestamp shift
    if(grepl("timestamp", tempadjust$qa_note[r])){
      start_raworder <- tempadjust$datorder[r]
      if(r == nrow(tempadjust)){
        end_raworder <- with(soilmoisture_master, max(raw_datorder[portid == i], na.rm =T))
      }else{
        end_raworder <- with(tempadjust, datorder[(r+1)])-1
      }
      # id rows in soilmoisture_master
      tempstart <- with(soilmoisture_master, which(portid == i & raw_datorder == start_raworder))
      tempend <- with(soilmoisture_master, which(portid == i & raw_datorder == end_raworder))
      # pull adjustment amount for note
      tempdiff <- unique(with(soilmoisture_master, cleanorder[tempstart:tempend] - raw_datorder[tempstart:tempend])) %>% na.exclude()
      # annotate
      soilmoisture_master$qa_note2[tempstart:tempend] <- paste("Timestamp off, data sequence shifted", tempdiff, "2hr-timesteps to correct timestamp")
    }
    
  }
  # annotate first timestamp data recorded for each logger in qa_note for cleanorder == 1
  # grab first timestamp
  firstdat <- min(with(soilmoisture_master, clean_datetime[portid == i & raw_datorder >= 1 & !is.na(vwc)]))
  stopifnot(length(firstdat)==1 & firstdat != Inf)
  # id row in soilmoisture_master
  firstrow <- with(soilmoisture_master, which(portid == i & cleanorder == 1))
  # annotate
  soilmoisture_master$qa_note2[firstrow] <- paste("First non-NA data recording for", i, "starts", firstdat) 
  
  # annotate last timestep if no vwc recorded (should be present on next download date, unless at project stop)
  lastdat <- max(with(soilmoisture_master, clean_datetime[portid == i & raw_datorder >= 1 & !is.na(vwc)]))
  stopifnot(length(firstdat)==1 & firstdat != Inf)
  # id row in soilmoisture_master
  lastrow <- with(soilmoisture_master, which(portid == i & cleanorder == max(cleanorder)))
  notelast <- with(soilmoisture_master, is.na(vwc[lastrow] & is.na(raw_datorder[lastrow])))
  # if no vwc present and no raw_datorder, annotate
  if(notelast){
    # annotate
    soilmoisture_master$qa_note2[lastrow] <- paste("Last non-NA data recording for", i, "ends", lastdat) 
  }
}
# add note for any breaks in data (not first or last clean order, filename, vwc and raw_datorder/raw_datetime will be NA)
checkNA <- apply(select(soilmoisture_master, filename:qa_note2), 1, function(x) all(is.na(x)))
View(soilmoisture_master[checkNA,]) # looks okay for infilling
soilmoisture_master$qa_note2[checkNA] <- "Missing data, break in data logger recording"

# create lookup table for all dates, full water years
firstWY <- as.Date(paste0(min(year(soilmoisture_master$clean_datetime)), "-10-1"))
lastWY <- as.Date(paste0(max(year(soilmoisture_master$clean_datetime))+1, "-9-30"))
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
soilmoisture_master <- soilmoisture_master %>%
  group_by(portid) %>%
  fill(filename, .direction= "downup") %>%
  ungroup() %>%
  # add in trt info, logger info..
  left_join(distinct(select(soilmoisture_all, logger:comp_trt)), by = "portid") %>%
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
summary(soilmoisture_master)
summary(is.na(soilmoisture_master))

# remake plots from initial compilation loop
logger_plot_clean <- soilmoisture_master %>%
  ggplot(aes(dowy, vwc, col = as.factor(port))) + #date_time, 
  geom_line(alpha = 0.5) +
  ggtitle(paste0("Compost prelim plot (", Sys.Date(),  "): cleaned soil moisture (VWC),\nby logger, by water year, colored by port")) +
  scale_color_discrete(name = "port") +
  facet_grid(logger~waterYear, scales = "free_x")
logger_plot_clean

trt_plot_clean <- soilmoisture_master %>%
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
