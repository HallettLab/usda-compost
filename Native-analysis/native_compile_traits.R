# assess coverage of traits for nat recruitment experiment

# -- SETUP -----
# load needed libraries
library(tidyverse)
library(cowplot)

# modify default settings
options(stringsAsFactors = F)
theme_set(theme_test())
na_vals <- c("" , " ","NA", NA)

# path to google drive (Julie's folder) [can make this for more users w access via googledrive package later if wish.. but this less annoying for now]
googlepath <- "/Users/scarlet/Library/CloudStorage/GoogleDrive-cawh3971@colorado.edu/Mi\ unidad/CA_Traits_Collaboration" 
# list all csvs in traits collab
gfiles <- list.files(googlepath, full.names = T, pattern = "csv")
# read all csvs to list to inspect
climvar_traits_list <- lapply(gfiles, read_csv)
# set file names as list element names
names(climvar_traits_list) <- gsub(paste0(googlepath, "/"), "", gfiles)

# path to dropbox
dbpath <- "~/Dropbox/"


# specify dropbox pathway
datpath <- "~/Dropbox/USDA-compost/Data/"
 

# read in various trait data
natlong <- read.csv(paste0(datpath, "Native/Native_CleanedData/Compost_Native_LongClean.csv"), strip.white = T, na.strings = na_vals)
natwide <- read.csv(paste0(datpath, "Native/Native_CleanedData/Compost_Native_WideClean.csv"), strip.white = T, na.strings = na_vals)
# read in treatment key
trtkey <- read.csv(paste0(datpath, "Compost_TreatmentKey.csv"),na.strings = na_vals, strip.white = T)
# read in master spp list (lookup table)
spplist <- read.csv(paste0(datpath, "Compost_SppList.csv"), na.strings = na_vals, strip.white = T)
# read in cleaned cover data for background comparison
coverlong <- read.csv(paste0(datpath, "Cover/Cover_CleanedData/Compost_Cover_LongClean.csv"), strip.white = T, na.strings = na_vals)


# -- COMPILE TRAIT DATA AVAILABLE! -----


# Julie's greenhouse-screened traits for ClimVar competition and recruitment subexperiments in GS 2017
climvar_mature_traits <- read_csv(gfiles[grep("traits_mature", gfiles)]) # has height, biomass, and rootlength per day
climvar_mean_traits_allphase <- read_csv(gfiles[grep("means", gfiles)]) # for this analysis most interested in phase 3 [mature]

