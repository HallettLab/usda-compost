# clean 2021 native diversity sub-experiment data

# script purpose:
# read in apr 2021, may 2021 cover dat + cross walk
# stack apr and may data with spplist
# apply crosswalk to update apr unknowns 
# verify abundances summed by code (no duplicate codes per plot)
# > also add col for gopher disturbance (screen for gopher in notes)
# treatments will be: herbicided, unherbided (both apr + may) + background herbicide (surveyed may 2021 only)

# notes:
# borrowing code from clean_cover script to build
# major difference is will need to add col for herbicided vs. non-herbicided



# -- SETUP -----
rm(list = ls()) # clear environment
# load needed libraries
library(tidyverse)
library(readxl)
# modify default settings
options(stringsAsFactors = F)
na_vals <- c("" , " ",".", "NA", NA)
theme_set(theme_test())

# specify dropbox pathway (varies by user)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/"
}else{
  ## probs LMH and AS? (fix if not -- rproj is two levels down in github repo)
  datpath <- "~/Dropbox/USDA-compost/Data/"
}


# list files in entered data folder
vegfiles <- list.files(paste0(datpath, "Native/Native_EnteredData"), full.names = T)
# read in treatment key
trtkey <- read.csv(paste0(datpath, "Compost_TreatmentKey.csv"),na.strings = na_vals, strip.white = T)
# read in master spp list (lookup table)
spplist <- read.csv(paste0(datpath, "Compost_SppList.csv"), na.strings = na_vals, strip.white = T)
# add generic dichelostemma to spplist

# pull out crosswalk
crosswalk <- read.csv(vegfiles[grepl("crosswalk", vegfiles)], na.strings = na_vals) %>%
  # add herbicide info
  mutate(herbicide = ifelse(grepl("NH", `Plot.`), "Non-herbicided", "Herbicided")) %>%
  # add code for may
  left_join(spplist[c("species", "code4")], by = c("May21_correction" = "species")) %>%
  rename(may_code4 = code4) %>%
  # add code for apr species
  left_join(spplist[c("species", "code4")], by = c("Apr21_species" = "species")) %>%
  rename(apr_code4 = code4) %>%
  dplyr::select(plot, herbicide, Apr21_species:ncol(.)) %>%
  rename_all(casefold) %>%
  # drop any unknown that doesn't have a correction
  filter(!is.na(may_code4)) %>%
  # change T(race) to 0.01 cover
  mutate(cover = gsub("T", "0.01", cover),
         cover = as.numeric(cover))



# -- TIDY AND TRANSPOSE ENTERED COVER DATA -----
# loop through cover data to read in, transpose, and append to master cover dataset
# initiate master data frames
native_master_long <- data.frame()
native_master_wide <- data.frame()

# loop iterates through 2021 native entered data
# > also updates entries for april that were ID'd in May
for(i in vegfiles[!grepl("crosswalk", vegfiles)]){
  # read in dataset
  vegdat <- read.csv(i, na.strings = na_vals, header = F, blank.lines.skip = T, strip.white = T)
  print(paste("Transposing and tidying", str_extract(i, "(?<=_)[:alnum:]+(?=.csv)"), "native dataset"))
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
  notes$mon <- lubridate::month(notes$date)
  
  # split herbicide info from plot
  notes$herbicide <- ifelse(grepl("NH", notes$plot), "Non-herbicided", "Herbicided")
  # note whether treatment native seeded or ambient
  notes$seedtrt <- ifelse(grepl("background", i), "Unseeded", "Native seeded")
  # clean up plot -- data collection only occurred in 2021, so only plot# used
  notes$plot <- parse_number(notes$plot)
  
  # join plot treatment info
  notes <- left_join(notes, trtkey) #left_join preserves order of sampling, merge alphabetizes plots
  #reorder cols
  notes <- notes[c("file", "pages", "plot", "plotid", "fulltrt", "block", "nut_trt", "ppt_trt", "herbicide", "seedtrt", "yr", "mon", "date", "recorder", "notes")]
  
  
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
  rownames(vegdat) <- NULL
  
  
  # replace entered species name with synonym if present
  synonyms <- vegdat$V1[vegdat$V1 %in% spplist$species[!is.na(spplist$compost_synonym)]]
  if(length(synonyms)>0){
    for(s in synonyms){
      vegdat$V1[vegdat$V1 == s] <- spplist$compost_synonym[spplist$species == s]
    }
  }
  # combine rows with same name after synonym swap
  dups <- vegdat$V1[duplicated(vegdat$V1)]
  for(d in dups){
    temprows <- which(vegdat$V1 == d)
    sumcov <- t(vegdat[temprows, 2:ncol(vegdat)])
    # these are still character vals, T value hasn't been assigned
    sumcov_combo <- ifelse(is.na(sumcov[,1]) & !is.na(sumcov[,2]), sumcov[,2], sumcov[,1])
    # assign to both one row and drop the other
    vegdat[temprows[1], 2:ncol(vegdat)] <- sumcov_combo
    vegdat <- vegdat[-(temprows[2]),]
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
  vegdat$code4[plotpos] <- "plot"
  vegdat$code4[litpos] <- "litter_depth_cm"
  # check to be sure all rows have a code4 value before proceeding
  stopifnot(all(!is.na(vegdat$code4)))
  
  vegdat <- vegdat[,2:ncol(vegdat)] # remove species col (using codes for headers)
  vegdat <- data.frame(t(vegdat)) # transpose, change matrix to data frame
  colnames(vegdat) <- vegdat[1,]
  vegdat <- vegdat[-1,]
  rownames(vegdat) <- NULL
  # create column for herbicided
  vegdat$herbicide <- ifelse(grepl("NH", vegdat$plot), "Non-herbicided", "Herbicided")
  # create column for seed trt
  vegdat$seedtrt <-  ifelse(grepl("background", i), "Unseeded", "Native seeded")
  # clean up plot
  vegdat$plot <- parse_number(vegdat$plot)
  
  #clean up trace values and empty cells for non-species cells (e.g. percent green, percent bare, litter depth)
  vegdat[vegdat == "T"] <- 0.01
  # make sure plot and pct cols are first (got reorganized in 2021 with some type of R/package update)
  plotnames <- names(vegdat)[!grepl("[A-Z0-9]{4}", names(vegdat))]
  sppnames <- names(vegdat)[!names(vegdat) %in% plotnames]
  # reorganize
  vegdat <- cbind(vegdat[plotnames],vegdat[sppnames])
  vegdat[plotnames] <- sapply(vegdat[plotnames], function(x) ifelse(is.na(x)|x == "n/a",0,x))
  # make all cols except herbicide numeric
  vegdat[,!grepl("herbic|seedtr", names(vegdat), ignore.case = T)] <- sapply(vegdat[,!grepl("herbic|seedtr", names(vegdat), ignore.case = T)], as.numeric)
  
  # reorder species cols alphabetically
  vegdat <- vegdat[c("plot", "herbicide", "seedtrt", plotnames[grepl("pct|lit", plotnames)], sort(sppnames))]
  
  # join notes and cover data
  clean_vegdat <- full_join(notes,vegdat, by = c("plot", "herbicide", "seedtrt"))
  
  # NOTE > can introduce more logic checks here as needed...
  
  
  # -- CREATE TIDY LONG FORM DATASET ----
  vegdat_long <- clean_vegdat %>%
    # gather species cover only (keep pct green through litter depth in their own cols in case want to drop or break out in analysis)
    gather(code4, pct_cover, all_of(sppnames), na.rm = T) %>%
    # clean up period and number additions from cols with same code
    mutate(code4 = gsub("[.][0-9]", "", code4)) %>%
    # group by code and sum pct cover [2020 has nikolai's unique species names, but correspond to duplicated code4]
    grouped_df(names(.)[names(.) != "pct_cover"]) %>%
    summarise(pct_cover = sum(pct_cover)) %>%
    ungroup()
  
  # if April data, update unknown species
  if(grepl("April", i)){
    vegdat_long$rowid <- as.numeric(rownames(vegdat_long)) 
    for(p in unique(clean_vegdat$plot)){
      crosssub <- subset(crosswalk, plot == p)
      for(h in unique(crosssub$herbicide)){
        crosssub2 <- subset(crosssub, herbicide == h)
        tempnote <- c("Data QA note -- Apr 2021 unknown species updated with May 2021 ID's: ")
        for(r in 1:nrow(crosssub2)){
          tempid <- with(vegdat_long, rowid[plot == p & herbicide == h & code4 == crosssub2$apr_code4[r] & pct_cover == crosssub2$cover[r]])
          # check only pulled 1 row, and that it pulled a row
          stopifnot(length(tempid)==1 & !is.na(tempid))
          # update apr spp with may update
          vegdat_long$code4[tempid] <- crosssub2$may_code4[r]
          # append species updated
          tempnote <- paste0(tempnote, paste0(crosssub2$apr_code4[r], ", ", crosssub2$apr21_species[r], ", to ", crosssub2$may_code4[r], ", ", crosssub2$may21_correction[r], "; "))
        }
        # clean up note before appending to plot notes
        tempnote <- gsub("; $", "", tempnote)
        # add to plot notes
        currentnote <- unique(vegdat_long$notes[vegdat_long$plot == p & vegdat_long$herbicide == h])
        if(is.na(currentnote)){
          # if no note present, assign tempnote
          vegdat_long$notes[vegdat_long$plot == p & vegdat_long$herbicide == h] <- tempnote
        }else{
          # if note present, append tempnote
          vegdat_long$notes[vegdat_long$plot == p & vegdat_long$herbicide == h] <- paste(currentnote, tempnote, sep = "; ")
        }
      }
    }
    # sum cover for any duplicates codes from updated spp, drop rowid
    vegdat_long <- vegdat_long %>%
      grouped_df(names(vegdat_long)[!names(vegdat_long) %in% c("pct_cover", "rowid")]) %>%
      summarise(pct_cover = sum(pct_cover)) %>%
      ungroup() %>%
      data.frame()
  }
  
  
  # continue with data prep
  vegdat_long <- vegdat_long %>%
  # append species descriptive info
  left_join(distinct(subset(spplist, !grepl("red california|different .Navar.|tall branchy", species) & is.na(compost_synonym), select = c(species:nativity))), by = "code4") %>%
    # order by plot, date, code
    arrange(plot, date, code4) %>%
    # rearrange cols
    select(file, pages, plot, plotid, fulltrt, block:ncol(.))
  
  
  
  # -- APPEND LONG FORM TO MASTER LONG FORM -----
  native_master_long <- rbind(native_master_long, vegdat_long)

  
  # -- FINISHING -----
  # if end of loop, make wideform dataset and print done
  if(i == vegfiles[sum(!grepl("crosswalk", vegfiles))]){
    
    # CREATE MASTER WIDE FORM FROM MASTER LONG FORM -----
    native_master_wide <- native_master_long[,1:grep("pct_cover",colnames(native_master_long), fixed = T)] %>%
      spread(code4, pct_cover, fill = 0) %>%
      arrange(plot, date)
    
    print("Master cover dataset compiled! Inspect and if all looks good write out and proceed to analysis! (w00t w00t!)")
  }
}

# -- POST-LOOP CHECKS -----
# check all looks okay
native_master_long <- data.frame(native_master_long)
native_master_wide <- data.frame(native_master_wide)
str(native_master_long)
summary(duplicated(native_master_long))
str(native_master_wide)
summary(duplicated(native_master_wide))
summary(native_master_long)

# check only plots with QA note are those that should have been updated
aprQAcheck <- subset(native_master_wide, grepl("Data QA", notes), select = c(plot, herbicide, seedtrt, yr, mon))
nrow(aprQAcheck) == nrow(distinct(crosswalk, plot, herbicide)) # looks okay
# review notes
with(native_master_wide, unique(notes[grepl("Data QA", notes)])) # okay


ggplot(native_master_long, aes(factor(plot), pct_cover, col = fxnl_grp)) +
  geom_jitter(height = 0) +
  facet_grid(ppt_trt ~ nut_trt)

ggplot(native_master_long, aes(factor(plot), pct_cover, col = nativity)) +
  geom_boxplot() +
  #stat_summary(fun = sum) +
  facet_grid(ppt_trt ~ nut_trt) # native plants seeded maybe like compost w drought??

# add column to long to indicate whether native plant was seeded
# > species added with F. microstachys (FEMI/VUMI), ESCA, BRCA, NEMA, TRCI
native_master_long$native_seeded <- ifelse(native_master_long$code4 %in% c("FEMI", "ESCA", "BRCA", "NEMA", "TRCI"), "Yes", "No")


# -- FINISHING -----
write.csv(native_master_long, paste0(datpath, "Native/Native_CleanedData/Compost_Native_LongClean.csv"), row.names = F)
write.csv(native_master_wide, paste0(datpath, "Native/Native_CleanedData/Compost_Native_WideClean.csv"), row.names = F)

