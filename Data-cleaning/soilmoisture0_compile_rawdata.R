# clean and compile usda compost soil moisture
# author(s): ctw (caitlin.t.white@colorado.edu), ashley shaw (ashleysh@uoregon.edu)
# created: 2019-05-08 

# script purpose:
# Part 1:
# > iterate through each raw files soil logger dataset in the compost dropbox SoilMoisture_RawData folder
# > append dataset to master soil moisture data frame
# > return full master soil moisture data frame with plot of data for QA check
# Part 2:
# > assess missing data and time break jumps
# > save prelim QA figs to Dropbox SoilMoisture_DataQA folder (USDA-compost/Data/SoilMoisture/SoilMoisture_DataQA/)
# Part 3:
# > write out two .csvs to Dropbox data QA folder for processing in soilmoisture1_correct_timestamps.R
# i. compiled raw data with plot treatment info
# ii. table of timebreaks that need adjustment


# notes:
# run soilmoisture1_correct_timestamp.R after running this script
# as of Oct 2021 still missing 2021 data for two loggers (B1L1 and B2L2)



# -- SETUP -----
# clear environment
rm(list=ls())
# load libraries needed
library(readxl) # for reading in soil moisture datasets
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

# list files in entered data folder
datfiles <- list.files(paste0(datpath, "SoilMoisture/SoilMoisture_RawData"), full.names = T, pattern = "B[1-4]L[0-9].*[.]xls", ignore.case = T)

# read in treatment key
trtkey <- read.csv(paste0(datpath, "Compost_TreatmentKey.csv"), na.strings = na_vals, strip.white = T)
# read in soil data logger lookup table
loggerkey <- read.csv(paste0(datpath, "SoilMoisture/SoilMoisture_RawData/decagon_logger_key.csv"), na.strings = na_vals, strip.white = T)
# read in prepped hourly CIMIS met data for QA checks after compilation
cimis_hrly <- list.files(paste0(datpath, "CIMIS"), full.names = T) %>%
  subset(grepl(".csv",.)) %>%
  read.csv(na.strings = na_vals) 



# -- 1. COMPILE SOIL MOISTURE RAW FILES -----
# prep logerkey for joining (lower case colnames)
colnames(loggerkey) <- casefold(colnames(loggerkey))
soilmoisture_all <- data.frame()

# order datfiles chronologically and by logger
datfiles_df <- data.frame(filename = datfiles,
                          logger = str_extract(datfiles, "B[0-9]L[0-9]"),
                          filedate = file.info(datfiles)$mtime) %>%
  arrange(logger, filedate)

# initiate data frame for compiling soil moisture data
for(i in datfiles_df$filename){
  print(paste("Compiling ", 
              gsub("-", "", regmatches(i, regexpr("B[0-9].*[0-9]-", i))),
              "soil moisture data!"))
  # read in i dataset
  tempdat <- read_excel(i, sheet = 1, na = na_vals, trim_ws = T,
                        skip = 2,
                        col_types = c("date", rep("text", 5))) #read moisture as text to avoid warnings
  # grab original colnames
  tempdat_names <- names(read_excel(i, sheet = 1))
  # reset port colnames (col 2: ncol)
  colnames(tempdat)[2:ncol(tempdat)] <- tempdat_names[2:length(tempdat_names)]
  colnames(tempdat) <- casefold(colnames(tempdat))
  # convert soil moisture cols to numeric
  tempdat[,2:ncol(tempdat)] <- sapply(tempdat[,2:ncol(tempdat)], as.double)
  print("Here is the summary of the dataset:")
  print(summary(tempdat))
  
  #grab logger from file name
  tempdat$logger <- tempdat_names[1]
  # rename measurement time (col1)
  colnames(tempdat)[1] <- "date_time"
  # split date and time
  tempdat$time <- format(tempdat$date_time, format = "%H:%M:%S")
  tempdat$date <- format(tempdat$date_time, format = "%Y-%m-%d")
  
  #gather port data
  tempdat2 <- gather(tempdat, port, vwc, `port 1`:`port 5`) %>%
    # strip "port" from port to match with logger key
    mutate(port = as.numeric(gsub("port ", "", port)),
           # compute day of year
           doy = yday(date_time),
           waterYear = ifelse(month(date_time) %in% 10:12,
                              year(date_time)+1, year(date_time))) %>%
    # compute water year day of year
    ungroup() %>%
    group_by(waterYear) %>%
    mutate(dowy = ifelse(month(date_time)<10, doy+92, NA),
           maxdoy = max(doy)) %>% #create temporary col to store max day of year
    # ungroup to calculate dowy for Oct-Dec depending on leap year/non leap year total days
    ungroup() %>%
    mutate(dowy = ifelse(is.na(dowy) & maxdoy == 366, doy-274, 
                         ifelse(is.na(dowy) & maxdoy == 365, doy-273, 
                                # if max doy doesn't go to 365 or 366, calculate months 10-12 more manually
                                ifelse(is.na(dowy) & month(date_time) == 10, day(date_time), # if Oct, day of month stay the same
                                       ifelse(is.na(dowy) & month(date_time) == 11, day(date_time)+31, # if Nov, add days in Oct
                                              ifelse(is.na(dowy) & month(date_time) == 12, day(date_time)+(31+30), # if Dec, add days in Oct + Nov 
                                                     dowy))))),
           filename = str_extract(i, "B[0-9]L[0-9].*$"),
           # make date col Date format
           date = as.Date(date, format = "%Y-%m-%d")) %>%
    # join logger key
    left_join(loggerkey, by = c("logger", "port")) %>%
    # add portid
    unite(portid, logger, port, remove = F) %>%
    # unite plot and subplot so joins with treatment key
    unite(plotid, plot, subplot, sep = "") %>% #plot = paste0(plot, subplot)) %>%
    left_join(trtkey, by = "plotid") %>%
    # clarify trt is the composition plot treatment
    rename(comp_trt = trt) %>%
    # reorder cols, remove maxdoy column
    dplyr::select(logger, port, portid, plot, plotid, fulltrt, block, nut_trt, ppt_trt, comp_trt, date_time, date, time, doy:dowy, vwc, filename)
  
  # compile with master soil moisture data frame
  soilmoisture_all <- rbind(soilmoisture_all, tempdat2) %>%
    distinct()
  
  if(i == datfiles_df$filename[nrow(datfiles_df)]){
    #clean up
    rm(tempdat, tempdat2)
    print("All done! Here are 2 plots of the compiled data. Toggle back arrow to see both.")
    logger_plot <- soilmoisture_all %>%
      filter(!is.na(vwc)) %>%
      mutate(full_trt = paste(nut_trt,ppt_trt, sep ="_")) %>%
      ggplot(aes(dowy, vwc, col = as.factor(port))) + #date_time, 
      geom_line(alpha = 0.5) +
      ggtitle(paste0("Compost QA check (", Sys.Date(),"): raw soil moisture (VWC),\nby water year, by logger, colored by port")) +
      scale_color_discrete(name = "port") +
      facet_grid(waterYear~logger)
    
    trt_plot <- soilmoisture_all %>%
      filter(!is.na(vwc)) %>%
      mutate(full_trt = paste(nut_trt,ppt_trt, sep ="_")) %>%
      ggplot(aes(dowy, vwc, group = portid, col = as.factor(block))) +
      geom_line(alpha = 0.5) +
      ggtitle(paste0("Compost treatment check (", Sys.Date(), "): raw soil moisture (VWC),\nby nutrient x drought treatment, by water year, colored by block")) +
      scale_color_discrete(name = "block") +
      facet_grid(full_trt~waterYear)
    
    print(logger_plot)
    print(trt_plot)
    print("Review compiled dataset and if all looks good, write out to Compost Dropbox SoilMoisture_CleanedData folder.")
  }
}

# check cols classed as expected
glimpse(soilmoisture_all)
summary(soilmoisture_all)
sapply(soilmoisture_all, function(x) summary(is.na(x)))
# check all date col vals == date from date_time (no weirdness during loop, forced NAs)
summary(soilmoisture_all$date == date(soilmoisture_all$date_time)) #yips

# check all unique dates have 1 dowy only
wydf <- subset(soilmoisture_all, select = c(date, doy, waterYear, dowy)) %>%
  distinct() %>%
  group_by(date) %>%
  mutate(nobs = length(unique(dowy))) %>%
  ungroup()
summary(wydf) # looks fine
# are the leap years what is expected?
sort(unique(soilmoisture_all$waterYear[soilmoisture_all$dowy == 366])) #2020 only.. that's what it should be of 2019, 2020, and 2021
# > all looks good, proceed!

# order by logger, port, then rowid (chronologically -- even tho datetimes not all correct [yet])
soilmoisture_all <- data.frame(soilmoisture_all) %>%
  mutate(rowid = as.numeric(row.names(.))) %>%
  # add sequence order
  group_by(logger, port) %>%
  mutate(datorder = 1:length(date_time)) %>%
  ungroup()

# write out preliminary plots to qa folder if desired
ggsave(filename = paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/PrelimQA_Figures/Compost_rawVWC_bylogger.pdf"),
       plot = logger_plot,
       width = 8, height = 8, units = "in")
ggsave(filename = paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/PrelimQA_Figures/Compost_rawVWC_bytreatment.pdf"),
       plot = trt_plot,
       width = 8, height = 8, unit = "in")


# > notes:
## > loggers with wonky water years (should only be 2019, 2020, 2021):
## B2L2: 2000, 2001, 2002, (no 2021?)
## B2L4: 2000, 2048 (block 2 comp control + comp dry?)
## B2L5: 2000
## B3L3: 2003, 2004, 2009, 2010
## B3L4: 2002, 2045, 2046 (block 3 fert wet)


# -- 2. PREVIEW DATA ISSUES -----
# add portid to logger key
loggerkey <- unite(loggerkey, portid, logger, port, remove = F)

# do all loggers have the same #obs?
dev.off() # clear out plotting device if needed
#sapply(split(soilmoisture_all$datorder, paste(soilmoisture_all$logger, soilmoisture_all$port, sep = "_")), max)
group_by(soilmoisture_all, logger, port) %>%
  mutate(maxnobs = max(datorder)) %>%
  ungroup() %>%
  subset(datorder == maxnobs) %>%
  ggplot(aes(as.factor(port), as.character(datorder))) +
  geom_point() +
  labs(y = "count", x= "logger port",
       title = "USDA Compost soil moisture QA: # observations per logger-port",
       subtitle = paste0(Sys.Date(), "; if counts differ by +1, needs review")) +
  # shrink y axis text for legibility
  theme(axis.text.y = element_text(size = 5)) +
  facet_wrap(~logger)
# troubleshoot suspect data (e.g. dates, vwc) post-compilation
# save to qa figs
ggsave(paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/PrelimQA_Figures/nobs_perlogger_alldates.pdf"),
       width = 6, height = 6, units = "in", scale = 1.2)

# > notes:
## > loggers to follow up on (* = also flagged for wonky years)
## B2L2*, B2L4*, B2L5*, B3L3* (B3L4 flagged for bad years but has consistent nobs per port)


# 2a. Assess date breaks ----- 
# check if date timestamp in expected year (1 = yes, 0 = no)
logcheck <- dplyr::select(soilmoisture_all, logger, port, portid, date, filename, datorder, rowid) %>%
  distinct() %>%
  mutate(filemo = substr(filename,6,12),
         filemo = as.Date(filemo, format = "%d%b%y"),
         tracker = ifelse(year(date) %in% c(year(filemo), year(filemo)-1), 1, 0),
         yr = substr(date, 1,4)) %>%
  group_by(logger, port, filename) %>%
  mutate(trackseq = 1:length(date)) %>%
  ungroup()

# visualize discrepancies
ggplot(logcheck, aes(trackseq, as.factor(tracker), col = yr)) +
  geom_point() +
  labs(y = "Correct date? (1 = yes, 0 = no)",
       x = "Collection sequence (per file)",
       title = "QA: Soil moisture date sequence consistency, by download file by logger",
       subtitle = paste0(Sys.Date(), "; USDA Compost project")) +
  scale_color_viridis_d() +
  facet_grid(logger~filemo, scale = "free_x", space = "free_x")
# save to qa figs
ggsave(paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/PrelimQA_Figures/years_perlogger_alldates.pdf"),
       width = 8, height = 6, units = "in", scale = 1.5, dpi = 72)
# for apr 2020 files, B2L2, B2L4, B3L4 need corrections from some point in middle all the way through download date
# > B2L2 has two different years: 2000, 2001
# > B2L4 has three different years: 2047, 2048, 2000 .. looks like L2 and L4 switch to final last year on same day?
# > B3L4 has 1 different year: 2002

# can compare vwc values to similar treatments to see if B2L2 and B3L4 stopped on date collected or earlier (i.e. is the data gap in the middle of the sequence or at the end?) 
# go by logger?
unique(soilmoisture_all$fulltrt[soilmoisture_all$logger %in% c("B2L2", "B2L4", "B3L3")])
datecheck <- soilmoisture_all %>%
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


# add interval event cols
rundf <- select(datecheck, logger, portid, date_time, datorder, timediff) %>%
  arrange(portid, datorder) %>%
  mutate(intevent = NA, 
         qa_note = NA)


for(i in unique(rundf$portid)){
  if(all(rundf$timediff[rundf$portid == i] == 2)){
    rundf$intevent[rundf$portid == i] <- 1
  }else{
    event <- 1
    tempstart <- which(rundf$portid == i & rundf$datorder == 1)
    tempend <- max(which(rundf$portid == i))
    tempbreaks <- c(which(rundf$portid ==i & rundf$timediff != 2), tempend)
    for(t in tempbreaks){
      if(t == tempbreaks[length(tempbreaks)]){
        rundf$intevent[tempstart:tempend] <- event
      }else{
        rundf$intevent[tempstart:(t-1)] <- event
      }
      # add qa note
      if(event != 1){
        tempnote <- ifelse(rundf$timediff[tempstart] >2 & rundf$timediff[tempstart] <= 8, "needs NA infill", 
                           ifelse(t == tempbreaks[length(tempbreaks)], "begin last run", "needs correct timestamp"))
        rundf$qa_note[tempstart] <- tempnote
      }
      event <- event+1
      tempstart <- t
    }
  }
}

# check that sequencing worked as expected
# > only plotting logger-ports that have more than 1 event (i.e. was a break in the 2-hr interval recording)
ggplot(subset(rundf, portid %in% portid[!is.na(qa_note)]), aes(datorder, intevent)) +
  geom_point (alpha = 0.5) +
  labs(y = "Run (continuous 2hr interval data)",
       x = "Order data collected",
       title = "USDA Compost soil moisture QA: # continuous 2hr runs for problematic loggers",
       subtitle = paste0(Sys.Date(), "; paneled by logger_port")) +
  facet_wrap(~portid)
# save to qa figs
ggsave(paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/PrelimQA_Figures/datarun_breaks_problemloggers.pdf"),
       width = 7, height = 6, units = "in", scale = 1.5)
# 2020 data update interp:
# > B2L2_1, _4 and _5 look similar; B2L2_2 and _3 are a pair; B2L4 ports look similar; B3L4 all look similar
# > *BUT* each set is doing something different

# 2021 data update interp:
# > B2L4 ports look like they're all doing something different
# > B2L5 ports all doing something different
# > B3L3 ports all doing something different
# > B3L4 all ports look similar
# > again each logger deviates differently than others.. will need to treat manually and individually like last year

# what are the unique collection intervals?
sort(unique(rundf$timediff))
# > 2 is the norm, 4 and 6 hr intervals are just data gaps that can be infilled with NAs for missing 2-hr steps
# > larger magnitude intervals are the loggers to treat

# attach interval event to datecheck
datecheck <- left_join(datecheck, rundf)


# do all ports have same timestamp (time intervals)?
portcheck <- datecheck %>%
  dplyr::select(logger, timediff, datorder) %>%
  distinct() %>%
  group_by(logger) %>%
  mutate(seqcheck = duplicated(datorder)) %>%
  ungroup()

# visualize
# remove last day for each because will have NA diffed time interval
ggplot(subset(portcheck, !is.na(timediff)), aes(datorder, as.factor(timediff))) +
  geom_point(alpha = 0.3) +
  geom_point(data = subset(portcheck, seqcheck), aes(datorder, as.factor(timediff)), col = "turquoise", pch = 1, size = 2) +
  labs(y = "Timestep from last collection (in hours)",
       x = "Collection sequence",
       title = "Compost data QA: Soil moisture timestep interval, blue = +1 interval within logger per collection sequence",
       subtitle = paste0(Sys.Date(), "; programmed for 2hrs, anything not == 2 on y-axis deviates")) +
  theme(plot.title = element_text(size = 11),
        axis.text.y = element_text(size = 5)) +
  facet_wrap(~logger)
# save to qa figs
ggsave(paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/PrelimQA_Figures/collectionintervals_perlogger_alldates.pdf"),
       width = 7, height = 6, units = "in", scale = 1.5)
# 2020 data update notes:
# > just same three loggers that have weird time jumps
# > all loggers have same two days where +1 port has a different timestamp, but that could be from data download/new recording?

# compare total number of recordings in Apr by logger
subset(datecheck, grepl("24Apr", filename)) %>%
  group_by(logger, port) %>%
  summarise(maxseq = length(datorder)) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(as.factor(port), as.character(maxseq), col = as.factor(maxseq))) +
  labs(title = "Compost soil moisture QA: Total # data recordings for 24Apr20 files by logger, by port",
       subtitle = Sys.Date(),
       y = "Total data recordings",
       x = "Logger port") +
  theme(plot.title = element_text(size = 10)) +
  scale_color_discrete(name = "Count") +
  facet_wrap(~logger, scale = "free_y")
# save to qa figs
ggsave(paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/PrelimQA_Figures/nobs_perlogger-port_24Apr20download.pdf"),
       width = 6, height = 6, units = "in")
# 2020 data update interp:
# > most loggers have same numer of collection points for 24apr20 file, but 4 loggers differ, and b2l2 differ from other loggers and also within by port (annoying..)
# > non-problematic ports that have fewer points for 24apr20 file do start with next expected collection point in 11jun20 file

# repeat for Sep 2021 downloads
# compare total number of recordings in Apr by logger
subset(datecheck, grepl("Sep21", filename)) %>%
  group_by(logger, port) %>%
  summarise(maxseq = length(datorder)) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(as.factor(port), as.character(maxseq), col = as.factor(maxseq))) +
  labs(title = "Compost soil moisture QA: Total # data recordings for Sep 2021 files by logger, by port",
       subtitle = Sys.Date(),
       y = "Total data recordings",
       x = "Logger port") +
  theme(plot.title = element_text(size = 10)) +
  scale_color_discrete(name = "Count") +
  facet_wrap(~logger, scale = "free_y")
# save to qa figs
ggsave(paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/PrelimQA_Figures/nobs_perlogger-port_Sep2021download.pdf"),
       width = 6, height = 6, units = "in")
# 2021 interp:
# > yikes. B2L4, B2L5, B3L3 got wacky. Looks like a normal amount of data points per logger should be around 5544-5545
# .. try to assign correct date based on soil moisture patterns in similar treatments?
# > and verify with CIMIS ppt data
# > start with B3L4 logger because that is the simplest correction, then B2L4, then B2L2



# pull date-times that need an adjustment
adjustdates <- subset(rundf, !is.na(qa_note)) %>%
  # attach filename info
  left_join(distinct(datecheck[c("portid", "fulltrt", "filename", "filemo", "date_time", "datorder")]))

# clean up environment
rm(portcheck, logcheck, event, i, t, tempbreaks, tempdat_names, tempend, tempnote, tempstart, rundf)



# 2b. Infill missing intervals with NA ----
# set expected mintime (project starts after 2017)
clean_mintime <- min(soilmoisture_all$date_time[year(soilmoisture_all$date_time) > 2017])
# set expected maxtime (spring 2021 slated as last field season)--can adjust this if needed
clean_maxtime <- max(soilmoisture_all$date_time[year(soilmoisture_all$date_time) > 2017 & year(soilmoisture_all$date_time) < 2022])   
soilmoisture_clean <- data.frame(date_time = rep(seq.POSIXt(clean_mintime, clean_maxtime, by = "2 hours"), times = length(unique(soilmoisture_all$portid)))) %>%
  mutate(portid = rep(unique(soilmoisture_all$portid), each = length(unique(date_time)))) %>%
  group_by(portid) %>%
  mutate(cleanorder = seq(1, length(date_time), 1)) %>%
  ungroup() %>%
  left_join(soilmoisture_all[c("portid", "logger", "port", "plotid", "block", "fulltrt", "date_time", "filename", "vwc")]) %>%
  group_by(portid) %>%
  fill(filename, .direction = "downup") %>%
  ungroup()


# see what's missing
soilmoisture_clean %>%
  replace_na(list(vwc = -999)) %>%
  #subset(grepl("Apr20", filename) & vwc == -999) %>%
  subset(vwc == -999) %>%
  mutate(logger = substr(portid, 1, 4),
         port = as.numeric(substr(portid, 6,6))) %>%
  ggplot(aes(date_time, (vwc + (0.1*port)))) +
  geom_point(aes(col = port, group = port), alpha = 0.4) + #position = position_jitter(height = 0.2)
  labs(title = "Compost soil moisture QA: dates with missing values, by logger and port",
       subtitle = Sys.Date()) +
  scale_x_datetime(date_labels = "%b '%y", breaks = "6 months") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~logger)
# save to qa figs
ggsave(paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/PrelimQA_Figures/missingdata_alllogger-port_allyrs.pdf"),
       width = 7, height = 6, units = "in")


# clean up
rm(clean_maxtime, clean_mintime)



# -- WRITE OUT -----
# quick check dats look good before writing out (e.g. no NAs in timestamp present)
summary(soilmoisture_all)
summary(adjustdates) # yep

# compiled raw data (from loop, not with NAs added in -- will repeat code in timestamp correction script)
write.csv(soilmoisture_all, paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/SoilMoisture_compiled_raw.csv"), row.names = F)
# table of date-time breaks that need adjustment or data missing
write.csv(adjustdates, paste0(datpath, "SoilMoisture/SoilMoisture_DataQA/SoilMoisture_raw_timestampbreaks.csv"), row.names = F)


