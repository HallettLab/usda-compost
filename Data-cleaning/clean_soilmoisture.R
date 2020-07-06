# clean and compile usda compost soil moisture
# author(s): ctw (caitlin.t.white@colorado.edu), ashley shaw
# created: 2019-05-08 

# script purpose
# iterate through each raw files soil logger dataset in the compost dropbox SoilMoisture_RawData folder
# append dataset to master soil moisture data frame
# return full master soil moisture data frame with plot of data for QA check


# notes
# script needs more clean up, basics just taken care of
# e.g. there are duplicate rows in the master dataset (AS thinks might happen if dates overlap in download events)



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
# read in CIMIS met data for QA checks after compilation
cimis <- read.csv(list.files(paste0(datpath, "CIMIS"), full.names = T), na.strings = na_vals) 


# -- COMPILE SOIL MOISTURE RAW FILES -----
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
                                dowy)),
           filename = str_extract(i, "B[0-9]L[0-9].*$")) %>%
    # join logger key
    left_join(loggerkey, by = c("logger", "port")) %>%
    # unite plot and subplot so joins with treatment key
    unite(fulltrt, plot, subplot, sep = "") %>% #plot = paste0(plot, subplot)) %>%
    left_join(trtkey, by = "fulltrt") %>%
    # clarify trt is the composition plot treatment
    rename(comp_trt = trt) %>%
    # reorder cols, remove maxdoy column
    dplyr::select(logger, port, plot, fulltrt, block, nut_trt, ppt_trt, comp_trt, date_time, date, time, doy:dowy, vwc, filename)
  
  # compile with master soil moisture data frame
  soilmoisture_all <- rbind(soilmoisture_all, tempdat2) %>%
    distinct()
  
  if(i == datfiles_df$filename[nrow(datfiles_df)]){
    #clean up
    rm(tempdat, tempdat2)
    print("All done! Here are 2 plots of the compiled data. Toggle back arrow to see both.")
    logger_plot <- soilmoisture_all %>%
      filter(!is.na(vwc)) %>%
      mutate(full_trt = paste(nut_trt,ppt_trt, sep ="_"),
             source = paste(logger, port, sep = "_")) %>%
      ggplot(aes(dowy, vwc, col = as.factor(port))) + #date_time, 
      geom_line(alpha = 0.5) +
      ggtitle("QA check: soil moisture (VWC), by logger, colored by port") +
      scale_color_discrete(name = "port") +
      facet_grid(waterYear~logger)
    
    trt_plot <- soilmoisture_all %>%
      filter(!is.na(vwc)) %>%
      mutate(full_trt = paste(nut_trt,ppt_trt, sep ="_"),
             source = paste(logger, port, sep = "_")) %>%
      ggplot(aes(dowy, vwc, group = source, col = as.factor(block))) +
      geom_line(alpha = 0.5) +
      ggtitle("Treatment check: soil moisture (VWC), by nutrient x drought treatment, colored by block") +
      scale_color_discrete(name = "block") +
      facet_grid(full_trt~waterYear)
    
    print(logger_plot)
    print(trt_plot)
    print("Review compiled dataset and if all looks good, write out to Compost Dropbox SoilMoisture_CleanedData folder.")
  }
}

# order by logger, port, then rowid (chronologically -- even tho datetimes not all correct [yet])
soilmoisture_all <- data.frame(soilmoisture_all) %>%
  mutate(rowid = as.numeric(row.names(.))) %>%
  # add sequence order
  group_by(logger, port) %>%
  mutate(datorder = 1:length(date_time)) %>%
  ungroup()





# -- DATA QA -----
# do all loggers have the same #obs?
#sapply(split(soilmoisture_all$datorder, paste(soilmoisture_all$logger, soilmoisture_all$port, sep = "_")), max)
group_by(soilmoisture_all, logger, port) %>%
  mutate(maxnobs = max(datorder)) %>%
  ungroup() %>%
  subset(datorder == maxnobs) %>%
  ggplot(aes(as.factor(port), as.character(datorder))) +
  geom_point() +
  facet_wrap(~logger)
# troubleshoot suspect data (e.g. dates, vwc) post-compilation


# 1. Prep CIMIS data for comparison -----
# prep cimis
# remove row with all NAs
cimis <- cimis[!apply(cimis,1,function(x) all(is.na(x))),] %>%
  # clean up colnames
  rename_all(function(x) gsub("[.]{2}", ".",x)) %>%
  rename_all(function(x) gsub("[.]$", "", x)) %>%
  # clean up cols to pair with soilmoisture dat
  mutate(date = as.Date(Date, format = "%m/%d/%Y"),
         # adjust 24:00 date to next day (to match SFREC format)
         date = ifelse(Hour.PST == 2400, as.character(date +1), as.character(date)),
         # adjust Julian day
         doy = ifelse(Hour.PST == 2400 & grepl("^12/31", Date), 1, 
                      ifelse(Hour.PST == 2400, Jul+1, Jul)),
         time = ifelse(nchar(Hour.PST)==3, paste0(0, substr(Hour.PST,1,1), ":", substr(Hour.PST, 2,3), ":00"),
                       paste(substr(Hour.PST,1,2), substr(Hour.PST, 3, 4), "00", sep =":")),
         # change 24:00 to 00:00 to match SFRECdat
         time = gsub("24", "00", time),
         #date_time = as.POSIXct(paste(date, time), format = "%Y-%m-%d %H:%M:%S"),
         # create unique continuous field for doy.24hr time to plot both soil moisture and ppt on same hourly x axis
         timeid = as.numeric(paste0(doy, ".", substr(time, 1,2))))

# test <- subset(datecheck, block == 2 & nut_trt == "N" & ppt_trt == "XC") %>%
#   mutate(timeid = as.numeric(paste(doy, substr(time, 1, 2), sep = ".")),
#          yr = year(date_time))
# ggplot(subset(cimis, !is.na(timeid)), aes(timeid, Precip..mm.)) +
#   geom_line(alpha = 0.4, col = "dodgerblue") +
#   geom_line(data = test, aes(timeid, vwc)) +
#   facet_grid(.~yr, scales = "free_x", space = "free_x")
# 
# plot_grid(ggplot(subset(cimis, Date >= min(test$date_time) & Date <= max(test$date_time)), aes(timeid, Air.Temp..C.)) +
#             geom_line(alpha = 0.4, col = "tomato") +
#             #geom_line(data = test, aes(timeid, vwc)) +
#             geom_smooth(col = "tomato", fill = "tomato") +
#             facet_grid(.~yr, scales = "free_x", space = "free_x"),
#           ggplot(subset(cimis, Date >= min(test$date_time) & Date <= max(test$date_time)), aes(timeid, Precip..mm.)) +
#             geom_line(alpha = 0.4, col = "dodgerblue") +
#             #geom_line(data = test, aes(timeid, vwc)) +
#             facet_grid(.~yr, scales = "free_x", space = "free_x"),
#           ggplot(test, aes(timeid, vwc)) +
#             geom_line(alpha = 0.4) +
#             facet_grid(.~yr, scales = "free_x", space = "free_x"),
#           nrow = 3,
#           align = "v")



# 2. Assess date breaks ----- 
# check if date timestamp in expected year (1 = yes, 0 = no)
logcheck <- dplyr::select(soilmoisture_all, logger, port, date, filename, datorder, rowid) %>%
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
       title = "Soil moisture date sequence consistency, by download file by logger") +
  scale_color_viridis_d() +
  facet_grid(logger~filemo, scale = "free_x", space = "free_x")
# save to qa figs
ggsave(paste0(datpath, "SoilMoisture/SoilMoisture_Figures/PrelimQA_Figures/years_perlogger_alldates.pdf"),
       width = 8, height = 6, units = "in")
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
         timeid = as.numeric(paste(doy, substr(time,1,2), sep = ".")),
         # create logger-port id
         portid = paste(logger, port, sep = "_")) %>%
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
rundf <- select(datecheck, portid, date_time, datorder, timediff) %>%
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
    tempbreaks <- which(rundf$portid ==i & rundf$timediff != 2)
    for(t in tempbreaks){
      if(t == tempbreaks[length(tempbreaks)]){
        rundf$intevent[tempstart:tempend] <- event
      }
      rundf$intevent[tempstart:(t-1)] <- event
      # add qa note
      if(event != 1){
        tempnote <- ifelse(rundf$timediff[tempstart] >2 & rundf$timediff[tempstart] <= 8, "needs NA infill", "correct timestamp")
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
  geom_point ()+
  facet_wrap(~portid)

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
       title = "Soil moisture timestep interval, turquoise = +1 interval within logger per collection sequence",
       subtitle = "Programmed for 2hrs, anything not 2 on y-axis deviates") +
  facet_wrap(~logger)
# save to qa figs
ggsave(paste0(datpath, "SoilMoisture/SoilMoisture_Figures/PrelimQA_Figures/collectionintervals_perlogger_alldates.pdf"),
       width = 8, height = 6, units = "in", scale = 1.1)
# > just same three loggers that have weird time jumps
# > all loggers have same two days where +1 port has a different timestamp, but that could be from data download/new recording?

# compare total number of recordings in Apr by logger
subset(datecheck, grepl("24Apr", filename)) %>%
  group_by(logger, port) %>%
  summarise(maxseq = length(datorder)) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(as.factor(port), as.character(maxseq), col = as.factor(maxseq))) +
  labs(title = "Total number data recordings for 24Apr20 files by logger, by port",
       y = "Total data recordings",
       x = "Logger port") +
  scale_color_discrete(name = "Count") +
  facet_wrap(~logger, scale = "free_y")
# save to qa figs
ggsave(paste0(datpath, "SoilMoisture/SoilMoisture_Figures/PrelimQA_Figures/nobs_perlogger-port_24Apr20download.pdf"),
       width = 6, height = 6, units = "in")
# > most loggers have same numer of collection points for 24apr20 file, but 4 loggers differ, and b2l2 differ from other loggers and also within by port (annoying..)
# > non-problematic ports that have fewer points for 24apr20 file do start with next expected collection point in 11jun20 file


# example.. annoying:
View(subset(datecheck, logger == "B2L2" & trackseq %in% 2360:2373)) 
# starts with port 1 having different time than other 4 (starts at trackseq = 2367), then slowly other ports have similar date, BUT they are all wrong timestamp regardless (time from yrs 2000, 2001 and not correct month..)

# .. try to assign correct date based on soil moisture patterns in similar treatments?
# > and verify with CIMIS ppt data
# > start with B3L4 logger because that is the simplest correction, then B2L4, then B2L2

# pull dates that need an adjustment
adjustdates <- subset(rundf, !is.na(qa_note))

# clean up environment
rm(portcheck, logcheck, event, i, t, tempbreaks, tempdat_names, tempend, tempnote, tempstart, rundf)


# 2.a. Triage B3L4 24Apr20 -----
b3l4trt <- gsub("[0-9]", "", unique(datecheck$fulltrt[datecheck$logger == "B3L4"]))


subset(datecheck, grepl(b3l4trt[1], fulltrt)) %>% # trackseq > 2250
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_line(alpha = 0.4) +
  geom_vline(data = subset(datecheck, logger == "B3L4" & timediff != 2), aes(xintercept = datorder), lty = 2) +
  ggtitle(paste("troubleshoot B3L4", b3l4trt[1], "date jump (dotted vert line = break)")) +
  geom_smooth(aes(fill = portid)) +
  facet_grid(logger~.)
# yikes... I guess B2L1 is most reliable for reference?

subset(datecheck, grepl(b3l4trt[1], fulltrt)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(datorder, diffvwc, col = portid)) +
  geom_vline(data = subset(datecheck, logger == "B3L4" & timediff != 2), aes(xintercept = datorder), lty = 2) +
  geom_line(alpha = 0.4) +
  ggtitle(paste("troubleshoot B3L4", b3l4trt[1], "date jump (dotted vert line = break)")) +
  #geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  facet_grid(logger~., scales = "free_x", space = "free_x")


subset(datecheck, grepl(b3l4trt[2], fulltrt)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(datorder, vwc, col = portid)) +
  geom_vline(data = subset(datecheck, logger == "B3L4" & timediff != 2), aes(xintercept = datorder), lty = 2) +
  geom_line(alpha = 0.4) +
  geom_smooth(aes(fill = portid)) +
  ggtitle(paste("troubleshoot B3L4", b3l4trt[2], "date jump (dotted vert line = break)")) +
  facet_grid(logger~.)
# > data gap after break, shift series so end lines up with end timestamp collected (and check spikes align with B2L1)

subset(datecheck, grepl(b3l4trt[2], fulltrt)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(datorder, diffvwc, col = portid)) +
  geom_vline(data = subset(datecheck, logger == "B3L4" & timediff != 2), aes(xintercept = datorder), lty = 2) +
  geom_line(alpha = 0.4) +
  ggtitle(paste("troubleshoot B3L4", b3l4trt[2], "date jump (dotted vert line = break)")) +
  #geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  facet_grid(logger~., scales = "free_x", space = "free_x")



# rule 1: if time of day from end of bad file 2 hrs before time of day of following file, assign new times backwards end to last bad date break
# rule 2: look for moisture spike alignment with comparison

b3l4dates <- subset(datecheck, grepl("B3L4", portid)) %>%
  subset(datorder %in% unique(c(datorder[!is.na(qa_note)], datorder[!is.na(qa_note)]-1))) %>%
  arrange(portid, datorder) %>%
  mutate(logger = substr(portid, 1,4),
         port = as.numeric(substr(portid, 6,7)),
         timeinterval = as.numeric(timeinterval),
         hrdiff = ifelse(timeinterval==0, 24-lag(timeinterval), 
                         ifelse(lag(timeinterval) > timeinterval, 24+ (timeinterval  - lag(timeinterval)), timeinterval- lag(timeinterval))),
         hrdiff = ifelse(is.na(qa_note), NA, hrdiff)) %>%
  dplyr::select(logger, date_time, filename, datorder, intevent, qa_note, hrdiff) %>%
  distinct()

# all ports have some time breaks
# > 1 is NA to infill for missing timestep
# > other is a date sequence that needs adjustment -- end of sequence seems to match?

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
  geom_line() +
  facet_wrap(~badport)

ggplot(subset(good_b3l4, goodport != "B2L1_5"), aes(date_time, diffrun, col = goodport)) +
  geom_line() +
  facet_wrap(~badport)

good_b3l4_trend <- group_by(good_b3l4, goodport, badport) %>%
  summarise(mean_dailydiff = mean(diffvwc)) %>%
  ungroup()

searchdata <- subset(datecheck, logger == "B3L4" & datorder >= 3863 & grepl("Apr20", filename)) %>%
  #subset(wetup == 1) %>%
  group_by(datorder) %>%
  mutate(nobs = sum(wetup, na.rm = T)) %>%
  ungroup() %>%
  subset(datorder >= (min(datorder[wetup == 1 & nobs > 2], na.rm = T)-20) & datorder <= (min(datorder[wetup == 1 & nobs > 2], na.rm = T)+80))

refdata <- subset(datecheck, grepl(str_flatten(b3l4trt, collapse = "|"), fulltrt)  & logger != "B3L4" & datorder >= 3863 & grepl("Apr20", filename)) %>%
  group_by(date_time) %>%
  mutate(nobs = sum(wetup, na.rm = T)) %>%
  ungroup() %>%
  subset(date_time >= (min(date_time[wetup == 1 & nobs > 1], na.rm = T)-as.difftime(20*2, units = "hours")) & date_time <= (min(date_time[wetup == 1 & nobs > 1], na.rm = T)+as.difftime(80*2, units = "hours")))
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




# 2.b. Triage B2L4 24Apr20 -----
b2l4trt <- gsub("[0-9]", "", unique(datecheck$fulltrt[datecheck$logger == "B2L4"]))

subset(datecheck, grepl(b2l4trt[1], fulltrt) & grepl("Apr20", filename)) %>% # & grepl("Apr20", filename) & trackseq > 2000 trackseq > 2250
  ggplot(aes(trackseq, vwc, col = paste(logger, port, sep = "_"))) +
  geom_line(aes(group = paste(logger, port, sep = "_")), alpha = 0.4) +
  geom_vline(data = distinct(subset(datecheck, timediff != 2 & grepl(b2l4trt[1], fulltrt)), logger, plot, timediff, trackseq), aes(xintercept = trackseq), lty = 2) +
  #geom_vline(data = subset(datecheck, logger == "B2L4" & timediff != 2 & trackseq > 2000), aes(xintercept = trackseq), lty = 2) +
  geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  ggtitle(paste("troubleshoot B2L4 date jump", b2l4trt[1], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l4trt[1])
#facet_grid(.~filemo, scales = "free_x", space = "free_x")

subset(datecheck, grepl(b2l4trt[1], fulltrt) & grepl("Apr20", filename)) %>% # & grepl("Apr20", filename) & trackseq > 2000 trackseq > 2250
  ggplot(aes(trackseq, diffvwc, col = paste(logger, port, sep = "_"))) +
  geom_line(aes(group = paste(logger, port, sep = "_")), alpha = 0.4) +
  geom_vline(data = distinct(subset(datecheck, timediff != 2 & grepl(b2l4trt[1], fulltrt)), logger, plot, timediff, trackseq), aes(xintercept = trackseq), lty = 2) +
  #geom_vline(data = subset(datecheck, logger == "B2L4" & timediff != 2 & trackseq > 2000), aes(xintercept = trackseq), lty = 2) +
  #geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  ggtitle(paste("troubleshoot B2L4 date jump", b2l4trt[1], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l4trt[1])

subset(datecheck, grepl(b2l4trt[2], fulltrt) & grepl("Apr20", filename)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(trackseq, vwc, col = paste(logger, port, sep = "_"), group = paste(logger, port, sep = "_"))) +
  geom_line(aes(group = paste(logger, port, sep = "_")), alpha = 0.4) +
  geom_vline(data = distinct(subset(datecheck, timediff != 2 & grepl(b2l4trt[2], fulltrt)), logger, plot, timediff, trackseq), aes(xintercept = trackseq), lty = 2) +
  #geom_vline(data = subset(datecheck, logger == "B2L4" & timediff != 2 & trackseq > 2000), aes(xintercept = trackseq), lty = 2) +
  geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  ggtitle(paste("troubleshoot B2L4 date jump", b2l4trt[2], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l4trt[1])
#facet_grid(.~filemo, scales = "free_x", space = "free_x")

subset(datecheck, grepl(b2l4trt[2], fulltrt) & grepl("Apr20", filename)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(trackseq, diffvwc, col = paste(logger, port, sep = "_"), group = paste(logger, port, sep = "_"))) +
  geom_line(aes(group = paste(logger, port, sep = "_")), alpha = 0.4) +
  geom_vline(data = distinct(subset(datecheck, timediff != 2 & grepl(b2l4trt[2], fulltrt)), logger, plot, timediff, trackseq), aes(xintercept = trackseq), lty = 2) +
  #geom_vline(data = subset(datecheck, logger == "B2L4" & timediff != 2 & trackseq > 2000), aes(xintercept = trackseq), lty = 2) +
  #geom_smooth(aes(fill = paste(logger, port, sep = "_"))) +
  ggtitle(paste("troubleshoot B2L4 date jump", b2l4trt[2], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l4trt[1])

# > slight data gap (spikes later in sequence should line up on same day [rain event])





# 2.c. Triage B2L2 24Apr20 -----
b2l2trt <- gsub("[0-9]", "", unique(datecheck$fulltrt[datecheck$logger == "B2L2"]))

subset(datecheck, grepl(b2l2trt[1], fulltrt) & grepl("Apr20", filename)) %>% # trackseq > 2250
  ggplot(aes(trackseq, vwc, col = paste(logger, port, sep = "_"), group = port)) +
  geom_line(alpha = 0.4) +
  geom_vline(data = distinct(subset(datecheck, timediff != 2 & grepl(b2l2trt[1], fulltrt)), logger, plot, timediff, trackseq), aes(xintercept = trackseq), lty = 2) +
  #geom_vline(aes(xintercept = max(unique(with(datecheck, trackseq[logger == "B2L2" & as.numeric(timediff) != 2])), na.rm = T)), lty = 2) +
  ggtitle(paste("troubleshoot B2L2 date jump", b2l2trt[1], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l2trt[1])

subset(datecheck, grepl(b2l2trt[1], fulltrt) & grepl("Apr20", filename)) %>% # trackseq > 2250
  ggplot(aes(trackseq, diffvwc, col = paste(logger, port, sep = "_"), group = port)) +
  geom_line(alpha = 0.4) +
  geom_vline(data = distinct(subset(datecheck, timediff != 2 & grepl(b2l2trt[1], fulltrt)), logger, plot, timediff, trackseq), aes(xintercept = trackseq), lty = 2) +
  #geom_vline(aes(xintercept = max(unique(with(datecheck, trackseq[logger == "B2L2" & as.numeric(timediff) != 2])), na.rm = T)), lty = 2) +
  ggtitle(paste("troubleshoot B2L2 date jump", b2l2trt[1], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l2trt[1])

subset(datecheck, grepl(b2l2trt[2], fulltrt) & grepl("Apr20", filename)) %>%  #&  & trackseq > 2250
  ggplot(aes(trackseq, vwc, col = paste(logger, port, sep = "_"), group = port)) +
  geom_vline(data = distinct(subset(datecheck, timediff != 2 & grepl(b2l2trt[2], fulltrt)), logger, plot, timediff, trackseq), aes(xintercept = trackseq), lty = 2) +
  #geom_vline(aes(xintercept = max(unique(with(datecheck, trackseq[logger == "B2L2" & as.numeric(timediff) != 2])), na.rm = T)), lty = 2) +
  geom_line(alpha = 0.4) +
  ggtitle(paste("troubleshoot B2L2 date jump", b2l2trt[2], "(dotted vert line = date break)")) +
  facet_grid(logger~b2l2trt[2])

subset(datecheck, grepl(b2l2trt[2], fulltrt) & grepl("Apr20", filename)) %>%  #&  & trackseq > 2250
  ggplot(aes(trackseq, diffvwc, col = paste(logger, port, sep = "_"), group = port)) +
  geom_vline(data = distinct(subset(datecheck, timediff != 2 & grepl(b2l2trt[2], fulltrt)), logger, plot, timediff, trackseq), aes(xintercept = trackseq), lty = 2) +
  #geom_vline(aes(xintercept = max(unique(with(datecheck, trackseq[logger == "B2L2" & as.numeric(timediff) != 2])), na.rm = T)), lty = 2) +
  geom_line(alpha = 0.4) +
  ggtitle(paste("troubleshoot B2L2 date jump", b2l2trt[2], "(dotted vert line = date break)")) +
  scale_x_continuous(breaks = seq(0,4250, 250)) +
  facet_grid(logger~b2l2trt[2])

# > think gap happened after break, and end of series for B2L2 should align with end date -- can gauge by spike in soil moisture in mid-3000s

# first assess relationship of soil moisture loggers to ppt during a reliable period (e.g. first 100 days?)
b2l2refdat <- subset(datecheck, grepl(str_flatten(b2l2trt, collapse = "|"), fulltrt) & grepl("Apr20", filename)) %>%
  filter(trackseq < 250)
# bring in ppt (this may not matter if spikes are due to irrigation)
b2l2refdat <- full_join(b2l2refdat, subset(cimis, date %in% unique(b2l2refdat$date)))

plot_grid(ggplot(b2l2refdat, aes(timeid, Precip.mm)) +
            geom_point()+
            facet_grid("CIMIS"~.),
          ggplot(b2l2refdat, aes(timeid, vwc, group = paste(logger, port, sep = "_"))) +
            geom_line() +
            facet_grid(logger~.),
          nrow = 2)
# > they are all pretty responsive to rainfall, even being irrigated plots


b2l2_fixdates <- distinct(subset(datecheck, timediff != 2 & grepl(b2l2trt[2], fulltrt)), logger, port, plot, fulltrt, comp_trt, date_time, date, time, timeinterval, timediff, trackseq, rowid)

tempdat <- data.frame()
test <- datecheck[(b2l2_fixdates$rowid[1]-20):(b2l2_fixdates$rowid[2]+20),]

# 2.d. Infill missing intervals with NA ----



# -- FINISHING -----
write_csv(soilmoisture_all, paste0(datpath, "SoilMoisture/SoilMoisture_CleanedData/SoilMosture_all_clean.csv"))

# write out preliminary plots if desired
ggsave(filename = paste0(datpath, "SoilMoisture/SoilMoisture_Figures/PrelimQA_Figures/Compost_VWCraw_bylogger.pdf"),
       plot = logger_plot,
       width = 8, height = 8, units = "in")
ggsave(filename = paste0(datpath, "SoilMoisture/SoilMoisture_Figures/PrelimQA_Figures/Compost_VWCraw_bytreatment.pdf"),
       plot = trt_plot,
       width = 8, height = 8, unit = "in")
