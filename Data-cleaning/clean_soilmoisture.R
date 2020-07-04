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

# initiate data frame for compiling soil moisture data
for(i in datfiles){
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
    
    if(i == datfiles[length(datfiles)]){
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


# -- DATA QA -----
# troubleshoot suspect data (e.g. dates, vwc) post-compilation

# check if date timestamp in expected year (1 = yes, 0 = no)
logcheck <- dplyr::select(soilmoisture_all, logger, port, date, filename) %>%
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
  scale_color_viridis_d() +
  facet_grid(logger~filemo, scale = "free_x", space = "free_x")

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
         filemo = as.Date(filemo, format = "%d%b%y")) %>%
  # crunch interval
  group_by(logger, port, filename) %>%
  mutate(timediff = lead(date_time,1)-date_time,
         trackseq = 1:length(date_time)) %>%
  ungroup()

# do all ports have same timestamp (time intervals?
portcheck <- datecheck %>%
  dplyr::select(logger, timediff, trackseq) %>%
  distinct() %>%
  group_by(logger) %>%
  mutate(seqcheck = duplicated(trackseq)) %>%
  ungroup()

# visualize
# remove last day for each because will have NA diffed time interval

ggplot(subset(portcheck, !is.na(timediff)), aes(trackseq, as.factor(timediff))) +
  geom_point(alpha = 0.3) +
  geom_point(data = subset(portcheck, seqcheck), aes(trackseq, as.factor(timediff)), col = "turquoise", pch = 1, size = 2) +
  labs(y = "Deviation from expected 2-hr interval (in days)",
       x = "Collection sequence") +
  facet_wrap(~logger)
# > just same three loggers that have weird time jumps
# > all loggers have same two days where +1 port has a different timestamp, but that could be from data download/new recording?
# > !! HOWEVER, looks like three loggers with weird timestamps had at least 1 port that had correct timestamp? if that is the case, can try to match that?

# example.. annoying:
View(subset(datecheck, logger == "B2L2" & trackseq %in% 2360:2373)) 
# starts with port 1 having different time than other 4 (starts at trackseq = 2367), then slowly other ports have similar date, BUT they are all wrong timestamp regardless (time from yrs 2000, 2001 and not correct month..)

# .. try to assign correct date based on soil moisture patterns in similar treatments?
## B2L2
subset(datecheck, grepl("FW", fulltrt) & grepl("Apr20", filename)) %>% # trackseq > 2250
  ggplot(aes(trackseq, vwc, col = block)) +
  geom_line(alpha = 0.4) +
  geom_vline(aes(xintercept = max(unique(with(datecheck, trackseq[logger == "B2L2" & as.numeric(timediff) != 2])), na.rm = T)), lty = 2) +
  ggtitle("troubleshoot B2L2 date jump (dotted vert line = break)") +
  facet_grid(logger~nut_trt)

subset(datecheck, grepl("NW", fulltrt) & grepl("Apr20", filename) & trackseq > 2250) %>%
  ggplot(aes(trackseq, vwc, col = block)) +
  geom_vline(aes(xintercept = max(unique(with(datecheck, trackseq[logger == "B2L2" & as.numeric(timediff) != 2])), na.rm = T)), lty = 2) +
  geom_point(alpha = 0.4) +
  facet_grid(logger~.)
# think gap happened after break, and end of series for B2L2 should align with end date -- can gauge by spike in soil moisture in mid-3000s

## B2L4
subset(datecheck, grepl("CD", fulltrt)) %>% # & grepl("Apr20", filename) & trackseq > 2000 trackseq > 2250
  ggplot(aes(trackseq, vwc, col = as.factor(block))) +
  geom_line(alpha = 0.4) +
  geom_vline(data = subset(datecheck, logger == "B2L4" & timediff != 2 & trackseq > 2000), aes(xintercept = trackseq), lty = 2) +
  geom_smooth(aes(fill = as.factor(block))) +
  ggtitle("troubleshoot B2L4 (block 2) date jump (dotted vert line = break)") +
  facet_grid(.~filemo, scales = "free_x", space = "free_x")

subset(datecheck, grepl("CXC", fulltrt)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(trackseq, vwc, col = logger)) +
  geom_vline(data = subset(datecheck, logger == "B2L4" & timediff != 2), aes(xintercept = trackseq), lty = 2) +
  geom_line(alpha = 0.2) +
  #geom_vline(data = subset(datecheck, logger == "B2L4" & timediff != 2 & trackseq > 2000), aes(xintercept = trackseq), lty = 2) +
  geom_smooth(aes(fill = logger)) +
  ggtitle("troubleshoot B2L4 (block 2) date jump (dotted vert line = break)") +
  facet_grid(.~filemo, scales = "free_x", space = "free_x")
# > slight data gap (spikes later in sequence should line up on same day [rain event])

# B3L4
subset(datecheck, grepl("FD", fulltrt)) %>% # trackseq > 2250
  ggplot(aes(trackseq, vwc, col = logger)) +
  geom_line(alpha = 0.4) +
  geom_vline(data = subset(datecheck, logger == "B3L4" & timediff != 2), aes(xintercept = trackseq), lty = 2) +
  ggtitle("troubleshoot B3L4 date jump (dotted vert line = break)") +
  geom_smooth(aes(fill = logger)) +
  facet_grid(.~filemo, scales = "free_x", space = "free_x")
# yikes... I guess B2L1 is most reliable for reference?

subset(datecheck, grepl("FXC", fulltrt)) %>% #  & grepl("Apr20", filename)
  ggplot(aes(trackseq, vwc, col = logger)) +
  geom_vline(data = subset(datecheck, logger == "B3L4" & timediff != 2), aes(xintercept = trackseq), lty = 2) +
  geom_line(alpha = 0.4) +
  geom_smooth(aes(fill = logger)) +
  facet_grid(.~filemo, scales = "free_x", space = "free_x")
# > data gap after break, shift series so end lines up with end timestamp collected (and check spikes align with B2L1)


# prep cimis
# create unique continuous field for doy.24hr time to plot both soil moisture and ppt on same hourly x axis
cimis <- subset(cimis, !is.na(Stn.Id))
cimis$timeid <- ifelse(nchar(cimis$Hour..PST.) == 3, paste0(cimis$Jul, ".0", cimis$Hour..PST.), paste0(cimis$Jul, ".", cimis$Hour..PST.))
cimis$timeid <- as.numeric(cimis$timeid)
cimis$Date <- as.Date(cimis$Date, format = "%m/%d/%Y")
cimis$yr <- year(cimis$Date)
test <- subset(datecheck, block == 2 & nut_trt == "N" & ppt_trt == "XC") %>%
  mutate(timeid = as.numeric(paste(doy, substr(time, 1, 2), sep = ".")),
         yr = year(date_time))
ggplot(subset(cimis, !is.na(timeid)), aes(timeid, Precip..mm.)) +
  geom_line(alpha = 0.4, col = "dodgerblue") +
  geom_line(data = test, aes(timeid, vwc)) +
  facet_grid(.~yr, scales = "free_x", space = "free_x")

plot_grid(ggplot(subset(cimis, Date >= min(test$date_time) & Date <= max(test$date_time)), aes(timeid, Air.Temp..C.)) +
            geom_line(alpha = 0.4, col = "tomato") +
            #geom_line(data = test, aes(timeid, vwc)) +
            geom_smooth(col = "tomato", fill = "tomato") +
            facet_grid(.~yr, scales = "free_x", space = "free_x"),
          ggplot(subset(cimis, Date >= min(test$date_time) & Date <= max(test$date_time)), aes(timeid, Precip..mm.)) +
            geom_line(alpha = 0.4, col = "dodgerblue") +
            #geom_line(data = test, aes(timeid, vwc)) +
            facet_grid(.~yr, scales = "free_x", space = "free_x"),
          ggplot(test, aes(timeid, vwc)) +
            geom_line(alpha = 0.4) +
            facet_grid(.~yr, scales = "free_x", space = "free_x"),
          nrow = 3,
          align = "v")

# -- FINISHING -----
write_csv(soilmoisture_all, paste0(datpath, "SoilMoisture/SoilMoisture_CleanedData/SoilMosture_all_clean.csv"))

# write out preliminary plots if desired
ggsave(filename = paste0(datpath, "SoilMoisture/SoilMoisture_Figures/Compost_VWCraw_bylogger.pdf"),
       plot = logger_plot,
       width = 8, height = 8, units = "in")
ggsave(filename = paste0(datpath, "SoilMoisture/SoilMoisture_Figures/Compost_VWCraw_bytreatment.pdf"),
       plot = trt_plot,
       width = 8, height = 8, unit = "in")
