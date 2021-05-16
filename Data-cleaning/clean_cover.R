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
library(lubridate)
options(stringsAsFactors = F)
na_vals <- c(" ", "", NA, "NA")
theme_set(theme_bw())

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
  rownames(notes) <- NULL # rename rownames in numeric sequence
  
  # change date from character to date format and add year column
  notes$date <- as.Date(notes$date, format = "%m/%d/%y")
  notes$yr <- as.numeric(substr(as.character(notes$date), 1,4))
  
  # prep plot col depending on which ID used
  # if plot col is number id, convert to numeric
  if(all(!grepl("[[:alpha:]]+", notes$plot))){
    notes$plot <- as.numeric(notes$plot)
  }
  # if plot ID is block-nut-ppt ID, rename col to plotid to join correctly with trtkey
  # > 2019 we use block-nut-ppt, 2020 use plot #
  if(any(grepl("[[:alpha:]]", notes$plot))){
    notes <- rename(notes, plotid = plot)
  } 
  
  # join plot treatment info
  notes <- left_join(notes, trtkey) #left_join preserves order of sampling, merge alphabetizes plots
  #reorder cols
  notes <- notes[c("plot", "plotid", "fulltrt", "block", "nut_trt", "ppt_trt", "yr", "date", "recorder", "notes")]
  
  
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
  # rename rownames in sequential order to preserve correct order as in datasheets
  rownames(vegdat) <- 1:nrow(vegdat)
  
  
  # replace entered species name with synonym if present
  synonyms <- vegdat$V1[vegdat$V1 %in% spplist$species[!is.na(spplist$compost_synonym)]]
  if(length(synonyms)>0){
    for(s in synonyms){
      vegdat$V1[vegdat$V1 == s] <- spplist$compost_synonym[spplist$species == s]
    }
  }
  # join species codes and descriptive info
  vegdat <- right_join(spplist[c("species", "code4")], vegdat, by = c("species" = "V1"))
  vegdat <- distinct(vegdat)
  # add codes for non-species cells (e.g. percent grass, percent bare/other)
  vegdat$code4[grepl("%", vegdat$species)] <- with(vegdat[grepl("%", vegdat$species),], 
                                                   paste0("pct_", casefold(str_extract(species, "[A-Z][a-z]+"))))
  
  # ID pct cover column headers as plot or plotid depending on whether letter present
  ## track new litpos and plotpos because something apparently has changed in 2021 (CTW hates new package versions)
  plotpos <- grep("plot", vegdat[,1], ignore.case = T)
  litpos <- grep("depth", vegdat[,1])
  vegdat$code4[plotpos] <- ifelse(grepl("[[:alpha:]]", vegdat$V3[plotpos]), "plotid", "plot")
  vegdat$code4[litpos] <- "litter_depth_cm"
  # check to be sure all rows have a code4 value before proceeding
  stopifnot(all(!is.na(vegdat$code4)))
  
  vegdat <- vegdat[,2:ncol(vegdat)] # remove species col (using codes for headers)
  vegdat <- data.frame(t(vegdat)) # transpose, change matrix to data frame
  colnames(vegdat) <- vegdat[1,]
  vegdat <- vegdat[-1,]
  rownames(vegdat) <- NULL
  
  #clean up trace values and empty cells for non-species cells (e.g. percent green, percent bare, litter depth)
  vegdat[vegdat == "T"] <- 0.01
  # make sure plot and pct cols are first (got reorganized in 2021 with some type of R/package update)
  plotnames <- names(vegdat)[!grepl("[A-Z0-9]{4}", names(vegdat))]
  sppnames <- names(vegdat)[!names(vegdat) %in% plotnames]
  # reorganize
  vegdat <- cbind(vegdat[plotnames],vegdat[sppnames])
  vegdat[plotnames] <- sapply(vegdat[plotnames], function(x) ifelse(is.na(x)|x == "n/a",0,x))
  # make all cols except plot numeric [except if plot is the numeric ID and not fulltrt]
  vegdat[,!grepl("plot", names(vegdat), ignore.case = T)] <- sapply(vegdat[,!grepl("plot", names(vegdat), ignore.case = T)], as.numeric)
  # if plot numeric ID, also make numeric so it joins with trtkey
  if(!all(grepl("[[:alpha:]]", vegdat[,grepl("plot", names(vegdat), ignore.case = T)]))){
    vegdat[,grepl("plot", names(vegdat), ignore.case = T)] <- as.numeric(vegdat[,grepl("plot", names(vegdat), ignore.case = T)])
  }
  # reorder species cols alphabetically
  vegdat <- vegdat[c(plotnames, sort(sppnames))]
  
  # join notes and cover data
  clean_vegdat <- full_join(notes,vegdat, by = names(vegdat)[1])
  
  # NOTE > can introduce more logic checks here as needed...
  
  
  # -- CREATE TIDY LONG FORM DATASET ----
  vegdat_long <- clean_vegdat %>%
    # gather species cover only (keep pct green through litter depth in their own cols in case want to drop or break out in analysis)
    gather(code4, pct_cover, all_of(sppnames), na.rm = T) %>%
    # clean up period and number additions from cols with same code
    mutate(code4 = gsub("[.][0-9]", "", code4)) %>%
    # group by code and sum pct cover [2020 has nikolai's unqiue species names, but correspond to duplicated code4]
    grouped_df(names(.)[names(.) != "pct_cover"]) %>%
    summarise(pct_cover = sum(pct_cover)) %>%
    ungroup() %>%
    # append species descriptive info
    left_join(distinct(subset(spplist, !grepl("red california|different .Navar.|tall branchy", species) & is.na(compost_synonym), select = c(species:nativity))), by = "code4") %>%
    # order by plot, date, code
    arrange(plot, date, code4) %>%
    # rearrange cols
    select(plot, plotid, fulltrt, block:ncol(.))
  
  
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


# -- PRELIM PLOT -----
# quick check that cover dat looks okay
# spp richness
data.frame(cover_master_long) %>%
  mutate(mo = month(date, label = T, abbr = T),
         date_group = paste(mo, yr, sep = " "),
         date_group = factor(date_group, levels = c("Apr 2019", "May 2019", "Apr 2020", "Apr 2021")),
         nut_trt = gsub("N", "XC", nut_trt)) %>%
  grouped_df(c("plot", "nut_trt", "ppt_trt", "date_group")) %>%
  summarise(S = length(unique(code4))) %>%
  ungroup() %>%
  ggplot(aes(plot, S, group = plot)) +
  geom_line(aes(group = plot)) +
  geom_point(aes(fill = date_group), pch = 21, alpha = 0.6) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_viridis_d(name = "Sample date") +
  ggtitle(paste0("USDA Compost data QA (", Sys.Date() ,"):\nprelim spp comp, richness by plot, treatments, and sample date")) +
  facet_grid(nut_trt~ppt_trt) # AS and CTW sampling in 2019 vs NS sampling in 2020.. influence in spp richness
# save to QA figs
ggsave(filename = paste0(datpath, "Cover/Cover_Figures/Prelim_QAFigures/SppComp_S_PlotxYr.pdf"),
       width = 5.25, height = 4, units = "in", scale = 1.1) 

# what about cover by family?
cover_master_long %>%
  mutate(mo = month(date, label = T, abbr = T),
         date_group = paste(mo, yr, sep = " "),
         date_group = factor(date_group, levels = c("Apr 2019", "May 2019", "Apr 2020", "Apr 2021"))) %>%
  grouped_df(c("plot","date_group", "Family")) %>%
  summarise(totcov = sum(pct_cover)) %>%
  ungroup() %>%
  ggplot(aes(plot, totcov)) +
  geom_line(aes(group = plot)) +
  geom_point(aes(fill = date_group), pch = 21, alpha = 0.6) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_viridis_d(name = "Sample date") +
  ggtitle(paste0("USDA Compost data QA (", Sys.Date() ,"): prelim spp comp, Family total cover by plot and sample date")) +
  facet_wrap(~Family, scales = "free") # interesting switch with poa vs forbs in 2020 for lower blocks..
# save to QA figs
ggsave(filename = paste0(datpath, "Cover/Cover_Figures/Prelim_QAFigures/SppComp_Familycover_PlotxYr.pdf"),
       width = 9, height = 6, units = "in", scale = 1.1) 



# -- FINISHING -----
write.csv(cover_master_long, paste0(datpath, "Cover/Cover_CleanedData/Compost_Cover_LongClean.csv"), row.names = F)
write.csv(cover_master_wide, paste0(datpath, "Cover/Cover_CleanedData/Compost_Cover_WideClean.csv"), row.names = F)
