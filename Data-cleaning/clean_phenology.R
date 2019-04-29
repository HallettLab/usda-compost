# clean and compile compost phenology dataset
# author(s): ctw (caitlin.t.white@colorado.edu)
# created: april 2019 (will be modified over project lifetime)

# script purpose:
# iterate through each phenology dataset chronologically
# read in phenology entered data and phenology photo key entered data
# merge both and append to master phenology dataset
# write out to compost dropbox: Phenology/Phenology_EnteredData



# -- SETUP -----
rm(list = ls()) # clear environment
# load needed libraries
library(tidyverse)
library(readxl)
# modify default settings
options(stringsAsFactors = F)
na_vals <- c("" , " ", "NA", NA)

# set path to compost data (main level)
# dependent on user, comment out path(s) that aren't pertinent to you
datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/" #ctw's path
#datpath <- "~/Dropbox/USDA-compost/Data/" # should work for LMH and AS

# list files in entered data folder
pheno_files <- list.files(paste0(datpath, "Phenology/Phenology_EnteredData/"), full.names = T, pattern = "_Phenology_.*[.]xl", ignore.case = T)
# read in treatment key
trtkey <- read.csv(paste0(datpath, "Compost_TreatmentKey.csv"), na.strings = na_vals, strip.white = T)


# -- COMPILE ----
#initiate df for all phenology data
pheno_master <- data.frame()

for(i in pheno_files){
  # phenology dat
  temp_pheno <- read_excel(i, sheet = 1, na = na_vals, trim_ws = T)  
  # photo key
  temp_photokey <- read_excel(i, sheet = 2, na = na_vals, trim_ws = T)  
  print(paste("Compiling", unique(temp_pheno$date)[1], "phenology data and photo key"))
  # join both datasets
  temp_dat <- full_join(temp_pheno, temp_photokey, by = c("recorder", "date", "plot", "subplot"))
  # concatenate plot and subplot so plot values match treatment key plot values
  temp_dat$plot <- with(temp_dat, ifelse(!is.na(subplot), paste0(plot, subplot), plot))
  # join treatment data
  temp_dat <- left_join(temp_dat, trtkey, by = "plot")
  # add year
  temp_dat$yr <- substr(as.character(temp_dat$date),1,4) %>% as.numeric()
  # reorder cols
  temp_dat <- dplyr::select(temp_dat, page:date, yr, plot, block:ppt_trt,pct_green:photo_notes)
  
  # rbind to master phenology df
  pheno_master <- rbind(pheno_master, temp_dat)
  
  # if end of loop, print done
  if(i == pheno_files[length(pheno_files)]){
    print("Master phenology dataset compiled! Inspect and if all looks good write out and proceed to analysis! (w00t w00t!)")
  }
}


# -- FINISHING -----
# write out
write.csv(pheno_master, paste0(datpath, "Phenology/Phenology_CleanedData/Compost_Phenology_Clean.csv"), row.names = F)

