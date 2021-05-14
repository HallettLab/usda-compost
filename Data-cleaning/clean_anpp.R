# read in, compile, and clean all anpp dats for analysis
# author(s): ctw
# date create: may 2021


# -- SETUP --
library(tidyverse)
library(lubridate)
# modify default settings
options(stringsAsFactors = F)
na_vals <- c("" , " ", "NA", NA)

# specify dropbox pathway (varies by user -- ctw can tweak this later)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/"
}else{
  ## probs LMH and AS? (fix if not -- rproj is two levels down in github repo)
  datpath <- "~/Dropbox/USDA-compost/Data/"
}

# list files in entered data folder
anpp_files <- list.files(paste0(datpath, "ANPP/ANPP_EnteredData"), full.names = T, pattern = "20([0-9]{2})", ignore.case = T)
# read in treatment key
trtkey <- read.csv(paste0(datpath, "Compost_TreatmentKey.csv"), na.strings = na_vals, strip.white = T)

# -- COMPILE ANPP ----
#initiate df for all data
anpp_master <- data.frame()

for(i in anpp_files){
  temp_anpp <- read.csv(i, na.strings = na_vals)
  # join trt key
  temp_anpp <- 
}
