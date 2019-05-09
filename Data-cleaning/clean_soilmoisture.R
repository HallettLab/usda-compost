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
library(dplyr)
library(ggplot2)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals <- c(" ", "", NA, "NA")

# set path to compost data (main level)
# dependent on user, comment out path(s) that aren't pertinent to you
datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/" #ctw's path
#datpath <- "~/Dropbox/USDA-compost/Data/" # should work for LMH and AS
# list files in entered data folder
datfiles <- list.files(paste0(datpath, "SoilMoisture/SoilMoisture_RawData"), full.names = T, pattern = "B[1-4]L[0-9].*[.]xls", ignore.case = T)

# read in treatment key
trtkey <- read.csv(paste0(datpath, "Compost_TreatmentKey.csv"), na.strings = na_vals, strip.white = T)
# read in soil data logger lookup table
loggerkey <- read.csv(paste0(datpath, "SoilMoisture/SoilMoisture_RawData/decagon_logger_key.csv"), na.strings = na_vals, strip.white = T)



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
                        col_types = c("date", rep("numeric", 5)))
  # grab original colnames
  tempdat_names <- names(read_excel(i, sheet = 1))
  # reset port colnames (col 2: ncol)
  colnames(tempdat)[2:ncol(tempdat)] <- tempdat_names[2:length(tempdat_names)]
  colnames(tempdat) <- casefold(colnames(tempdat))
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
    mutate(port = as.numeric(gsub("port ", "", port))) %>%
    # join logger key
    left_join(loggerkey, by = c("logger", "port")) %>%
    # unite plot and subplot so joins with treatment key
    mutate(plot = paste0(plot, subplot)) %>%
    left_join(trtkey, by = "plot") %>%
    # clarify trt is the composition plot treatment
    rename(comp_trt = trt) %>%
    dplyr::select(logger, plot, block, nut_trt, ppt_trt, comp_trt, date_time, date, time, vwc)
  
    # compile with master soil moisture data frame
    soilmoisture_all <- rbind(soilmoisture_all, tempdat2)
    
    if(i == datfiles[length(datfiles)]){
      #clean up
      rm(tempdat, tempdat2)
      print("All done! Here is a plot of the compiled data:")
      prelim_plot <- soilmoisture_all %>%
        mutate(full_trt = paste(nut_trt,ppt_trt, sep ="_")) %>%
        ggplot(aes(date_time, vwc, col = full_trt)) +
           geom_point(alpha = 0.5) +
           facet_wrap(~logger)
      print(prelim_plot)
      print("Review compiled dataset and if all looks good, write out to Compost Dropbox SoilMoisture_CleanedData folder.")
    }
}


# -- FINISHING -----
write.csv(soilmoisture_all, paste0(datpath, "SoilMoisture/SoilMoisture_CleanedData/SoilMosture_all_clean.csv"),
          row.names = F)
