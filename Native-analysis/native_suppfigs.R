# supplemental figures
# ctw

# notes:
# just making rain summary for SLM presentation 2022 Apr 5 for now
# more later...


# -- SETUP ----
# modify default settings
options(stringsAsFactors = F)
theme_set(theme_test())
na_vals <- c("" , " ","NA", NA)

# specify dropbox pathway (varies by user)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/"
}else{
  ## probs LMH and AS? (fix if not -- rproj is two levels down in github repo)
  datpath <- "~/Dropbox/USDA-compost/Data/"
}

# climvar datpath for native experiment ther
climvar_path <- "~/Dropbox/ClimVar/DATA/"

# raw cimis data for daily ppt and temp summaries over time
# > there is some temporal overlap here so can stitch together
compost_cimis <- read.csv(paste0(datpath, "CIMIS/Download_RawData/CIMIS_084BrownsValley_hourly_20211021.csv"), na.strings = na_vals)
climvar_cimis <- read.csv(paste0(climvar_path, "Met\ Station\ CIMIS/cimis_brownsvalley_daily_20190320.csv"), na.strings = na_vals)

# read in cleaned soil moisture and rainfall prepped for soil moisture
compost_sm <- read.csv(paste0(datpath, "SoilMoisture/SoilMoisture_CleanedData/SoilMoisture_all_clean.csv"), na.strings = na_vals)
compost_sm_ppt <- read.csv(paste0(datpath, "/CIMIS/CIMIS_084BrownsValley_hourly_prepped.csv"), na.strings = na_vals)


# -- DATA PREP ----
