# clean compost cover datasets
# author(s): ctw (caitlin.t.white@colorado.edu)
# created: 2019-04-28 (script will be modified over project lifetime)

# script purpose:
# read in master spp list and treatment lookup table
# iterate through all compost cover datasets to:
# change species names to codes, change "T" abundance values to 0.01
# build 2 types of cover data tables:
## 1) long form (tidy)
## 2) wide form (rows = site x date, cols = species)
## > for wide: transpose species to columns, sites to rows, fill in blank abundance values with 0s
## > for long: append all species descriptive columns in species key
## > for long and wide: add date sampled, recorder, notes, and plot treatment data
# write out clean cover datasets

# notes:
# re-run whenever unknowns are identified and changed in the entered cover dataset
# 2020 adjustment: process and append spring 2020 cover separately (Nikolai collected for us bc of COVID)
# > more unknowns to deal with and potentially bin 


# *** IMPORTANT TO DO BEFORE RUNNING THIS SCRIPT*** ##
# -- >> Run build_spplist.R before running the clean_cover script! << ---
# If new species have been entered, or to be sure, re-create the Compost species list via build_spplist.R
# Species list data appends to the compiled cover data, so important the species list dataset is complete before running this script




# -- SETUP -----
# clear environment
rm(list=ls())
# load libraries needed
library(tidyverse) # for dplyr, tidyr, and stringr
options(stringsAsFactors = F)
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
vegfiles <- list.files(paste0(datpath, "Cover/Cover_EnteredData"), full.names = T, pattern = "_Cover_", ignore.case = T)

# read in master spp list (lookup table)
spplist <- read.csv(paste0(datpath, "Compost_SppList.csv"), na.strings = na_vals, strip.white = T)
# read in treatment key
trtkey <- read.csv(paste0(datpath, "Compost_TreatmentKey.csv"), na.strings = na_vals, strip.white = T)

dat2020 <- read.csv(vegfiles[3], header = F, na.strings = na_vals, strip.white = T)
test <- read.csv(vegfiles[1], header = F, na.strings = na_vals, strip.white = T)

# -- TIDY AND TRANSPOSE ENTERED COVER DATA -----
# loop through cover data to read in, transpose, and append to master cover dataset
# initiate master data frames
cover_master_long <- data.frame()
cover_master_wide <- data.frame()

# loop iterates through each cover dataset and adds new spp to master spp set
# > process 2020 separately
for(i in vegfiles){
  # read in dataset
  vegdat <- read.csv(i, na.strings = na_vals, header = F, blank.lines.skip = T, strip.white = T)
  print(paste("Transposing and tidying", str_extract(i, "(?<=_)[:alnum:]+(?=.csv)"), "composition dataset"))
  # remove blank rows
  vegdat <- vegdat[!apply(vegdat, 1, function(x) all(is.na(x))),]
  # remove common name column if there
  vegdat <- vegdat[,!sapply(vegdat, function(x) any(grep("^Common name", x, ignore.case = T)))]
  # remove Grasses and Forbs section headers if present
  vegdat <- vegdat[!apply(vegdat,1, function(x) all(x %in% c("Grasses", "Forbs", NA))), ]
  # id where plot and cover data start
  plotpos <- grep("plot", vegdat[,1], ignore.case = T)
  
  
  # -- CREATE PLOT-DATE DATA TABLE -----
  # pull out recorder, sample date, and notes into separate data frame
  notes <- vegdat[1:plotpos,]
  notes <- data.frame(t(notes)) #transpose, change matrix class to data frame
  colnames(notes) <- casefold(gsub(":", "", notes[1,])) # set row 1 as column names, make lower case and remove colon
  notes <- notes[-1,] # remove row 1
  rownames(notes) <- seq(1,nrow(notes),1) # rename rownames in numeric sequence
  
  # change date from character to date format and add year column
  notes$date <- as.Date(notes$date, format = "%m/%d/%y")
  notes$yr <- substr(as.character(notes$date), 1,4)
  
  # prep plot col depending on which ID used
  # if plot col is number id, convert to numeric
  if(all(grepl("[[:digit:]]+", notes$plot))){
    notes$plot <- as.numeric(notes$plot)
  }
  # if plot ID is block-nut-ppt ID, rename col to fulltrt to join correctly with trtkey
  # > 2019 we use block-nut-ppt, 2020 use plot #
  if(any(grepl("[:alpha:]", notes$plot))){
    notes <- rename(notes, fulltrt = plot)
  } 
  
  # join plot treatment info
  notes <- left_join(notes, trtkey) #left_join preserves order of sampling, merge alphabetizes plots
  #reorder cols
  notes <- notes[c("plot", "fulltrt", "block", "nut_trt", "ppt_trt", "yr", "date", "recorder", "notes")]
  
  
  # -- PREP ABUNDANCE DATA FOR TIDYING AND TRANSPOSING -----
  # remove notes df rows from vegdat
  vegdat <- vegdat[plotpos:nrow(vegdat),] 
  # id row where species abundance data starts
  litpos <- grep("depth", vegdat[,1])
  # remove any species rows that don't have any entries for abundance value
  allNAs <- apply(vegdat[, 2:ncol(vegdat)], 1, function(x) all(is.na(x)))
  # set rows up to litter depth as FALSE to preserve them (in case 0s not entered and rock or bare never encountered)
  allNAs[1:litpos] <- FALSE 
  # allNAs <- allNAs[!grepl("NA",names(allNAs))] # remove NAs created by ignoring rows 1 through litter depth
  vegdat <- vegdat[!allNAs,] 
  
  # join species codes and descriptive info
  vegdat <- right_join(spplist[c("species", "code4")], vegdat, by = c("species" = "V1"))
  # add codes for non-species cells (e.g. percent grass, percent bare/other)
  vegdat$code4[grepl("%", vegdat$species)] <- with(vegdat[grepl("%", vegdat$species),], 
                                                   paste0("pct_", casefold(str_extract(species, "[A-Z][a-z]+"))))
  
  # ID pct cover column headers as plot or fulltrt depending on whether letter present
  vegdat$code4[1] <- ifelse(grepl("[[:alpha:]]", vegdat$V3[1]), "fulltrt", "plot")
  vegdat$code4[litpos] <- "litter_depth_cm"
  # check to be sure all rows have a code4 value before proceeding
  stopifnot(all(!is.na(vegdat$code4)))
  
  vegdat <- vegdat[,2:ncol(vegdat)] # remove species col (using codes for headers)
  vegdat <- data.frame(t(vegdat)) # transpose, change matrix to data frame
  colnames(vegdat) <- vegdat[1,]
  vegdat <- vegdat[-1,]
  rownames(vegdat) <- seq(1,nrow(vegdat), 1)
  
  #clean up trace values and empty cells for non-species cells (e.g. percent green, percent bare, litter depth)
  vegdat[vegdat == "T"] <- 0.01
  vegdat[colnames(vegdat)[1:7]] <- sapply(vegdat[colnames(vegdat)[1:7]], function(x) ifelse(is.na(x),0,x))
  # make all cols except plot numeric
  vegdat[,2:ncol(vegdat)] <- sapply(vegdat[,2:ncol(vegdat)], as.numeric)
  # reorder species cols alphabetically
  vegdat <- vegdat[c(colnames(vegdat)[1:litpos], sort(colnames(vegdat)[(litpos+1):ncol(vegdat)]))]
  
  # join notes and cover data
  clean_vegdat <- full_join(notes,vegdat, by = "fulltrt")
  # check to make sure all plots accounted for
  stopifnot(all(unique(clean_vegdat$fulltrt) %in% unique(vegdat$fulltrt)))
  # NOTE > can introduce more logic checks here as needed...
  
  
  # -- CREATE TIDY LONG FORM DATASET ----
  vegdat_long <- clean_vegdat %>%
    # gather species cover only (keep pct green through litter depth in their own cols in case want to drop or break out in analysis)
    gather(code4, pct_cover, (grep("depth",colnames(.))+1):ncol(.), na.rm = T) %>%
    # append species descriptive info
    left_join(spplist, by = "code4") %>%
    # order by plot, date, code
    arrange(plot, date, code4)
  
  
  # -- APPEND LONG FORM TO MASTER LONG FORM -----
  cover_master_long <- rbind(cover_master_long, vegdat_long)
  
  
  # -- CREATE MASTER WIDE FORM FROM MASTER LONG FORM -----
  cover_master_wide <- cover_master_long[,1:grep("pct_cover",colnames(cover_master_long), fixed = T)] %>%
    spread(code4, pct_cover, fill = 0) %>%
    arrange(plot, date)
  
  # if end of loop, print done
  if(i == vegfiles[length(vegfiles)]){
    print("Master cover dataset compiled! Inspect and if all looks good write out and proceed to analysis! (w00t w00t!)")
  }
}


# -- FINISHING -----
write.csv(cover_master_long, paste0(datpath, "Cover/Cover_CleanedData/Compost_Cover_LongClean.csv"), row.names = F)
write.csv(cover_master_wide, paste0(datpath, "Cover/Cover_CleanedData/Compost_Cover_WideClean.csv"), row.names = F)
