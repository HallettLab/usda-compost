# clean compost cover datasets
# author(s): ctw (caitlin.t.white@colorado.edu)
# date created: 2019-04-28 (script will be modified over project lifetime)

# script purpose:
# read in master spp list
# iterate through all compost cover datasets
# change species names to codes
# transpose species to columns, sites to rows
# change "T" abundance values to 0.01, fill in blank abundance values with 0s
# add date sampled, recorder, notes and columns
# write out clean cover data

# notes:
# > edit whenever create treatment table (can read that in and join to notes df)
# re-run whenever unknowns are identified and changed in the entered cover dataset


# -- SETUP -----
# clear environment
rm(list=ls())
# load libraries needed
library(dplyr) # for join functions
library(stringr) # for string pattern extraction
options(stringsAsFactors = F)
na_vals <- c(" ", "", NA, "NA")

# set path to compost data (main level)
# dependent on user, comment out path(s) that aren't pertinent to you
datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/" #ctw's path
#datpath <- "~/Dropbox/USDA-compost/Data/" # should work for LMH and AS
# list files in entered data folder
vegfiles <- list.files(paste0(datpath, "Cover/Cover_EnteredData"), full.names = T, pattern = "_Cover_", ignore.case = T)

# read in master spp list (lookup table)
spplist <- read.csv(paste0(datpath, "Compost_SppList.csv"), na.strings = na_vals, strip.white = T)


# -- TRANSPOSE ENTERED COVER DATA -----
# loop through cover data to read in, transpose, and append to master cover dataset
# initiate master data frame
cover_master <- data.frame()
# loop iterates through each cover dataset and adds new spp to master spp set
for(i in vegfiles){
  # read in dataset
  vegdat <- read.csv(i, na.strings = na_vals, header = F, blank.lines.skip = T, strip.white = T)
  print(paste("Transposing and tidying", vegdat$V2[2], "composition dataset"))
  # remove blank rows
  vegdat <- vegdat[!is.na(vegdat$V1),]
  
  # id where plot and cover data start
  plotpos <- grep("plot", vegdat[,1], ignore.case = T)
  
  # pull out recorder, sample date, and notes into separate data frame
  notes <- vegdat[1:plotpos,]  
  notes <- data.frame(t(notes)) #transpose, change matrix class to data frame
  colnames(notes) <- casefold(gsub(":", "", notes[1,])) # set row 1 as column names, make lower case and remove colon
  notes <- notes[-1,] # remove row 1
  rownames(notes) <- seq(1,nrow(notes),1) # rename rownames in numeric sequence
  notes$date <- as.Date(notes$date, format = "%m/%d/%y")
  # break out plot into block, compost treatment, and precip treatment, add year
  notes$block <- substr(notes$plot,1,1)
  notes$nut_trt <- substr(notes$plot,2,2)
  notes$ppt_trt <- substr(notes$plot,3,nchar(notes$plot))
  notes$yr <- substr(as.character(notes$date), 1,4)
  #reorder cols
  notes <- notes[c("plot", "block", "nut_trt", "ppt_trt", "yr", "date", "recorder", "notes")]
  
  # prep spp comp to transpose
  vegdat <- vegdat[plotpos:nrow(vegdat),] # remove notes df rows from vegdat
  # id row where species abundance data starts
  litpos <- grep("depth", vegdat[,1])
  # remove any species rows that don't have any entries for abundance value
  allNAs <- apply(vegdat[, 2:ncol(vegdat)], 1, function(x) all(is.na(x)))
  # set rows up to litter depth as FALSE to preserve them (in case 0s not entered and rock or bare never encountered)
  allNAs[1:litpos] <- FALSE 
  # allNAs <- allNAs[!grepl("NA",names(allNAs))] # remove NAs created by ignoring rows 1 through litter depth
  vegdat <- vegdat[!allNAs,] 
  
  vegdat <- right_join(spplist[c("species", "code4")], vegdat, by = c("species" = "V1"))
  vegdat$code4[grepl("%", vegdat$species)] <- with(vegdat[grepl("%", vegdat$species),], 
                                                     paste0("pct_", casefold(str_extract(species, "[A-Z][a-z]+"))))
  vegdat$code4[1] <- "plot"
  vegdat$code4[litpos] <- "litter_depth_cm"
  vegdat <- vegdat[,2:ncol(vegdat)] # remove species col (using codes for headers)
  vegdat <- data.frame(t(vegdat)) # transpose, change matrix to data frame
  colnames(vegdat) <- vegdat[1,]
  vegdat <- vegdat[-1,]
  rownames(vegdat) <- seq(1,nrow(vegdat), 1)
  
  #clean up trace values and empty cells
  vegdat[vegdat == "T"] <- 0.01
  vegdat[is.na(vegdat)] <- 0
  # reorder species cols alphabetically
  vegdat <- vegdat[c(colnames(vegdat)[1:litpos], sort(colnames(vegdat)[(litpos+1):ncol(vegdat)]))]
  
  # join notes and cover data
  clean_vegdat <- full_join(notes,vegdat, by = "plot")
  # check to make sure all plots accounted for
  stopifnot(all(unique(clean_vegdat$plot) %in% unique(vegdat$plot)))
  # NOTE > can introduce more logic checks here as needed...
  
  # append ith cover dataset to master cover dataset
  if(ncol(cover_master)==0){
    print("Adding first USDA Compost composition survey to master cover dataset")
  cover_master <- rbind(cover_master, clean_vegdat)
  }else{
    # in order to rbind ith cover dataset to master dataset, both sets have to have the same species (i.e. colnames match)
    #id which col is litter_depth_cm (last col before species start -- should be 10 but just in case, make generic)
    # list all spp in both datasets
    all_spp <- sort(unique( # alphabetize and select unique:
      # species in master cover dataset
      c(names(cover_master)[grep("depth", colnames(cover_master))+1:ncol(cover_master)],
        # species in ith cover dataset
        colnames(vegdat)[(litpos+1):ncol(vegdat)])))
    
    # check that all species are in vegdat
    if(!all((all_spp) %in% colnames(clean_vegdat))){
      vegdat_missing <- all_spp[!all_spp %in% colnames(clean_vegdat)]
      print(paste("Past compost species not encountered in", clean_vegdat$date[1], "composition survey:"))
      print(vegdat_missing)
      for(m in vegdat_missing){
        clean_vegdat[vegdat_missing] <- NA
      }
      # re-alphabetize species columns
     clean_vegdat <- clean_vegdat[c(colnames(clean_vegdat)[1:grep("depth_cm", colnames(clean_vegdat))],all_spp)]  
    }
    
    # check that all species are in master
    # check that all species are in vegdat
    if(!all((all_spp) %in% colnames(cover_master))){
      covdat_missing <- all_spp[!all_spp %in% colnames(cover_master)]
      print("New species added to master cover dataset!:")
      print(covdat_missing)
      for(m in covdat_missing){
        cover_master[covdat_missing] <- NA
      }
      # re-alphabetize species columns
      cover_master <- cover_master[c(colnames(cover_master)[1:grep("depth_cm", colnames(cover_master))],all_spp)]  
    }
    # append ith cover dataset to master
    print(paste("Adding", clean_vegdat$date[1], "cover data to master cover dataset!"))
    cover_master <- rbind(cover_master,clean_vegdat)
  }
  if(i == vegfiles[length(vegfiles)]){
    print("Master cover dataset compiled! Inspect and if all looks good write out and proceed to analysis! (w00t w00t!)")
  }
}


# -- FINISHING -----
write.csv(cover_master, paste0(datpath, "Cover/Cover_CleanedData/Compost_Cover_AllClean.csv"), row.names = F)
