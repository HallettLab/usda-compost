# clean and compile compost phenology dataset
# author(s): ctw (caitlin.t.white@colorado.edu)
# created: april 2019 (will be modified over project lifetime)

# script purpose:
# iterate through each phenology dataset chronologically
# read in phenology entered data and phenology photo key entered data
# merge both and append to master phenology dataset
# write out to compost dropbox: Phenology/Phenology_EnteredData

# modification jan 2020:
# below main loop that compiles growing season phenology, code added to read in and clean/compile winter 2020 phenology (format different than main phenology)
# modification for spring 2020:
# prep excel data sent by dustin before appending to 2019 data



# -- SETUP -----
rm(list = ls()) # clear environment
# load needed libraries
library(tidyverse)
library(readxl)
# modify default settings
options(stringsAsFactors = F)
na_vals <- c("" , " ",".", "NA", NA)

# specify dropbox pathway (varies by user -- ctw can tweak this later)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/"
}else{
  ## probs LMH and AS? (fix if not -- rproj is two levels down in github repo)
  datpath <- "~/Dropbox/USDA-compost/Data/"
}

# list files in entered data folder
pheno_files <- list.files(paste0(datpath, "Phenology/Phenology_EnteredData"), full.names = T, pattern = "_Phenology_.*[.]xl", ignore.case = T)
# read in treatment key
trtkey <- read.csv(paste0(datpath, "Compost_TreatmentKey.csv"), na.strings = na_vals, strip.white = T)


# -- COMPILE SPRING PHENOLOGY ----
#initiate df for all phenology data
pheno_master <- data.frame()
pct_cols <- c("pct_green", "pct_brown", "pct_bare")


# 1. spring 2019 + 2021 loop -----
# june 2020: changing this loop to prep 2019 only since winter 2020 and spring 2020 pheno files structured differently
# may 2021: modifying loop to accommodate 2021 data

for(i in pheno_files[grepl("2019|2021", pheno_files)]){
  # except if Dustin's 2021 dat, run in next loop
  if(grepl("20210427", i)){
    next
  }
  # phenology dat
  temp_pheno <- read_excel(i, sheet = 1, na = na_vals, trim_ws = T) %>% data.frame()  
  # photo key
  temp_photokey <- read_excel(i, sheet = 2, na = na_vals, trim_ws = T) %>% data.frame()
  print(paste("Compiling", unique(temp_pheno$date)[1], "phenology data and photo key"))
  
  # clean up pheno dat
  # change Ts to 0.01 and adjust pct cover of green (or if "<" noted)
  temp_pheno[pct_cols] <- sapply(temp_pheno[pct_cols], function(x) ifelse(trimws(x) %in% c("T", "<1"), 0.01,x))
  # check 1: look for < > and adjust based on sum of other cols
  flag1 <- sapply(temp_pheno[pct_cols],function(x) grepl("<|>",x))
  for(f in which(flag1)){
    adjust_pos <- grep("<|>", temp_pheno[f,pct_cols])
    stopifnot(length(adjust_pos)==1) # need to add code if there is more than 1 col that has greater/less than
    temp_pheno[f, pct_cols[adjust_pos]] <- 100-(sum(as.numeric(temp_pheno[f,pct_cols[-adjust_pos]])))
  }
  # make sure all pct cols converted to numeric
  temp_pheno <- mutate_at(temp_pheno, vars(pct_green:pct_bare), as.numeric)
  # check 2: flag rows that don't sum to 100 (or 0)
  flag2 <- apply(temp_pheno[pct_cols],1,function(x) !(sum(as.numeric(x)) %in% c(0, 100) | all(is.na(x))))
  for(f in which(flag2)){
    # logic check sum over 100 due to T value
    if(!0.01 %in% temp_pheno[f, pct_cols]){
      stop(paste0("Pct phenology cover for ", temp_pheno$plot[f], temp_pheno$subplot[f], " doesn't sum to 100 and not due to trace value. Review datasheet and data entry for ", unique(temp_pheno$date)))
    }
    # subtract from whichever column has the greatest number
    tempgreat <- which.max(temp_pheno[f, pct_cols])
    tempmin <- which.min(temp_pheno[f, pct_cols][which(temp_pheno[f, pct_cols]>0)])  # returns named value
    temp_pheno[[names(tempgreat)]][f] <- temp_pheno[[names(tempgreat)]][f] - temp_pheno[[names(tempmin)]][f]
    # check row sums correctly now
    stopifnot(sum(temp_pheno[f, pct_cols])==100)
    # subtract from pct green
    #temp_pheno$pct_green[f] <- with(temp_pheno, 100-(sum(pct_brown[f]+pct_bare[f]))) 
  }
  
  # join both datasets, correct name
  # need to trt 2019 and 2021 a little differently because datasheets structured differently
  if(unique(year(temp_pheno$date))==2019){
    tempdat <- full_join(temp_pheno, temp_photokey, by = c("recorder", "date", "plot", "subplot")) %>%
      #  # concatenate plot and subplot so plot values match treatment key plot values
      mutate(plotid = ifelse(!is.na(subplot), paste0(plot, subplot), plot)) %>%
      # drop plot so pairs correctly with trtkey
      dplyr::select(-c(plot, subplot))
  }else{
    tempdat <- full_join(temp_pheno, temp_photokey) %>%
      # add line to tempdat to match 2019 data. will be same as plot #
      mutate(line = plot) %>%
      #reorder cols
      dplyr::select(page, line, recorder:ncol(.))
  }
  # join treatment data
  tempdat <- left_join(tempdat, trtkey)
  # add year
  tempdat$yr <- substr(as.character(tempdat$date),1,4) %>% as.numeric()
  # reorder cols
  tempdat <- dplyr::select(tempdat, page:date, yr, plot, plotid, fulltrt, block:ppt_trt,pct_green:photo_notes)
  
  # rbind to master phenology df
  pheno_master <- rbind(pheno_master, tempdat)
  
  # if end of loop, print done
  if(i == last(pheno_files[grepl("2019|2021", pheno_files)])){
    # make df data frame
    pheno_master <- data.frame(pheno_master)
    print("Master phenology dataset compiled! Inspect and if all looks good write out and proceed to analysis! (w00t w00t!)")
  }
}

# 2. spring 2020 loop -----
# write loop for spring 2020 data + Dustin's 2021 data
# pull names in spring 2020 + 2021 excel workbooks
sp2020 <- excel_sheets(pheno_files[grep("Spring2020", pheno_files)])
sp2021 <- excel_sheets(pheno_files[grep("20210427", pheno_files)])

for(i in c(sp2020[sp2020 != "Master"], sp2021)){
  if(grepl("-20$", i)){
    tempdat <- read_excel(pheno_files[grep("Spring2020", pheno_files)], sheet = i, na = na_vals)
  }else{
    tempdat <- read_excel(pheno_files[grep("20210427", pheno_files)], sheet = i, na = na_vals)
  }
  tempinfo <- names(tempdat)[grep("[:alpha:]", names(tempdat))]
  # remove any cols that are all NA
  tempdat <- tempdat[,!sapply(tempdat, function(x) all(is.na(x)))]
  tempdat <- data.frame(tempdat)
  
  # id row with correct headers
  startrow <- grep("^Plot", tempdat[[names(tempdat)[1]]])
  # rename cols
  names(tempdat) <- tempdat[startrow, ]
  # start dat where data actually start
  tempdat <- tempdat[(startrow+1):nrow(tempdat),]
  # remove any rows that are all NA
  tempdat <- tempdat[!apply(tempdat, 1, function(x) all(is.na(x))),]
  tempdat$recorder <- trimws(gsub("^R[[:alpha:]]+ *: *", "", tempinfo[grep("Recorder", tempinfo, ignore.case = T)]))
  # clean up i (characters present in 2021)
  i <- str_remove(i, "[:alpha:]+")
  i <- trimws(i)
  tempdat$date <- as.Date(i, format = "%m-%d-%y")
  # standardize colnames
  names(tempdat) <- gsub(" |#", "", casefold(names(tempdat)))
  names(tempdat) <- gsub("%", "pct_", names(tempdat))
  names(tempdat)[grep("^photo", names(tempdat))] <- "photo_subplot" # these are full subplot photos
  
  # if there are notes in col 1, move them to notes col
  noterow <- grep("[[:alpha:]]", tempdat[,1])
  if(length(noterow)>0){
    # search for general notes (in first column) and put them in notes col
    tempdat$notes[noterow] <- paste("General phenology notes:",tempdat$plot[noterow])
    # NA the note in the plot col
    tempdat$plot[noterow] <- NA
  }
  
  # clean up trace or less than numbers (we make it 0.01)
  # check that pct_cols in tempdat
  stopifnot(all(pct_cols %in% names(tempdat)))
  # change Ts to 0.01 and adjust pct cover of green (or if "<" noted)
  tempdat[pct_cols] <- sapply(tempdat[pct_cols], function(x) ifelse(trimws(x) %in% c("T", "<1"), 0.01,x))
  # check 1: look for < > and adjust based on sum of other cols
  flag1 <- sapply(tempdat[pct_cols],function(x) grepl("<|>",x))
  for(f in which(flag1)){
    adjust_pos <- grep("<|>", tempdat[f,pct_cols])
    stopifnot(length(adjust_pos)==1) # need to add code if there is more than 1 col that has greater/less than
    tempdat[f, pct_cols[adjust_pos]] <- 100-(sum(as.numeric(tempdat[f,pct_cols[-adjust_pos]])))
  }
  # make sure all pct cols converted to numeric
  tempdat <- mutate_at(tempdat, vars(plot:pct_bare), as.numeric)
  
  # check 2: flag rows that sum to +100
  flag2 <- apply(tempdat[c("pct_green", "pct_brown", "pct_bare")],1,function(x) sum(as.numeric(x)) > 100)
  for(f in which(flag2)){
    # logic check sum over 100 due to T value
    if(!0.01 %in% tempdat[f, pct_cols]){
      stop(paste0("Pct phenology cover for ", tempdat$plot[f], tempdat$subplot[f], " exceeds 100 and not due to trace value. Review datasheet and data entry."))
    }
    # subtract from whichever column has the greatest number
    tempgreat <- which.max(tempdat[f, pct_cols])
    tempmin <- which.min(tempdat[f, pct_cols][which(tempdat[f, pct_cols]>0)])  # returns named value
    tempdat[[names(tempgreat)]][f] <- tempdat[[names(tempgreat)]][f] - tempdat[[names(tempmin)]][f]
    # check row sums correctly now
    stopifnot(sum(tempdat[f, pct_cols])==100)
  }  
  # finish with columns needed in pheno_master (others finished at end)
  tempdat$page <- 1 # every pheno sampling fits on 1 page
  tempdat$line <- 1:nrow(tempdat)
  tempdat$yr <- as.numeric(substr(unique(tempdat$date), 1, 4)) # just in case this loop gets recycled
  # except any general note should have NA for line #
  tempdat$line[grep("General", tempdat$notes)] <- NA
  # join treatment info
  tempdat <- merge(tempdat, trtkey, all.x = T)
  addnames <- names(pheno_master[!names(pheno_master) %in% names(tempdat)])
  tempdat <- cbind(tempdat, data.frame(matrix(nrow = nrow(tempdat), ncol = length(addnames), dimnames = list(NULL, addnames))))
  # order by plot
  tempdat <- tempdat[order(tempdat$plot),]
  # rbind to master df
  pheno_master <- rbind(pheno_master, tempdat[names(pheno_master)])
}


# 3. final QA check -----
# check range and sum
sapply(split(pheno_master[,pct_cols], pheno_master$date), function(x) range(x, na.rm = T)) # should always be between 0 and 100
sumcheck <- apply(pheno_master[,pct_cols], 1, function(x) unique(sum(x, na.rm = T))) # should be either 0 (NAs in row) or 100
View(pheno_master[sumcheck < 100,]) # if all NA's okay -- just means a general plot note

# visualize data to be sure nothing funky
ggplot(subset(pheno_master, !is.na(pct_green)), aes(date, pct_green, col = nut_trt)) +
  geom_point(alpha = 0.75) +
  geom_smooth(se = F) +
  scale_x_datetime(date_labels = "%m-%d") +
  ggtitle("USDA Compost senesence, all years, QA quickplot") +
  facet_grid(ppt_trt~yr, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# write out
ggsave(paste0(datpath, "Phenology/Phenology_Figures/Compost_senescence_QAprelim.pdf"), 
       width = 6, height = 4)

ggplot(subset(pheno_master, !is.na(pct_bare)), aes(date, pct_bare, col = nut_trt, group = plot)) +
  geom_line() +
  #geom_smooth(se = F) +
  scale_x_datetime(date_labels = "%m-%d") +
  facet_grid(ppt_trt~yr, scales = "free_x") # gopher disturbance causes hi bare

ggplot(subset(pheno_master, !is.na(pct_brown)), aes(date, pct_brown, col = nut_trt)) +
  geom_point() +
  geom_smooth(se = F) +
  ggtitle("USDA Compost brownup, all years, QA quickplot") +
  scale_x_datetime(date_labels = "%m-%d") +
  facet_grid(ppt_trt~yr, scales = "free_x")



# -- FINISHING -----
# write out
write.csv(pheno_master, paste0(datpath, "Phenology/Phenology_CleanedData/Compost_Phenology_Clean.csv"), row.names = F)



# -- WINTER 2020 PHENOLOGY ----
# 1/10/2020: for now, code to clean jan 2020 is in its own section as jan 2020 survey has phenology, veg and litter heights, and time spent on each row
# note: if phenology excel file is open while running code, read_excel line will throw error because pheno_files will contain working copy temp version and actual file (i.e. 2 files)
jan20dat <- read_excel(pheno_files[grepl("202001", pheno_files)], sheet = 1)
jan20foto <- read_excel(pheno_files[grepl("202001", pheno_files)], sheet = 2)

glimpse(jan20dat) # fix time
glimpse(jan20foto)
names(pheno_master) # <-- these are the cols that should be in final + height cols

# clean dat should have block and trts in their own cols, veg and litter heights, notes, and photo info -- keep wide format
# convert T(race) to 0.01
jan20_times <- dplyr::select(jan20dat, start_time, stop_time, date, plot) %>%
  rename(plotid = plot) %>%
  left_join(trtkey) %>%
  group_by(block, nut_trt) %>%
  mutate(start = unique(start_time[!is.na(start_time)]),
         stop = unique(stop_time[!is.na(stop_time)])) %>%
  ungroup() %>%
  dplyr::select(block, nut_trt, start, stop) %>%
  distinct() %>%
  mutate(duration_min = stop - start)

ggplot(jan20_times, aes(nut_trt, duration_min, col = as.factor(block))) +
  geom_jitter(width = 0.1) +
  labs(y = "Minutes", x = "Nutrient treatment", title = "Time to survey phenology per nutrient row (3 plots per row)",
       subtitle = "Generally takes longer as move uphill to grass-dominated blocks, but all compost comparable") +
  scale_color_discrete(name = "Block")

jan20dat_clean <- dplyr::select(jan20dat, date:notes) %>%
  rename(plotid = plot) %>%
  left_join(trtkey) %>%
  left_join(dplyr::select(jan20foto, -c(page_order, order)), by = c("plotid" = "plot", "date", "recorder")) %>%
  mutate(yr = as.numeric(substr(date, 1, 4))) %>%
  # reorder cols
  dplyr::select(page, line, recorder, date, yr, plot, plotid, fulltrt, block:ppt_trt, pct_litter:notes, photo_subplot:photo_notes) %>%
  # convert trace to 0.01
  mutate_at(grep("pct", names(.)), function(x) gsub("T", "0.01", x)) %>%
  # make all pct cols numeric
  mutate_at(grep("pct", names(.)), as.numeric)

# do range and sumcheck on pct cols
sapply(jan20dat_clean[grep("pct", names(jan20dat_clean))], range)
apply(jan20dat_clean[pct_cols], 1, sum) # only green, brown, bare (can be under 100 if litter comprised some)
summary(apply(jan20dat_clean[grep("pct", names(jan20dat_clean))], 1, sum)) # all at least 100 in cover. good.
# range check on heights
sapply(jan20dat_clean[grep("ht_", names(jan20dat_clean))], range) #ok
str(jan20dat_clean) #looks okay

# write out
write_csv(jan20dat_clean, paste0(datpath, "Phenology/Phenology_CleanedData/Compost_PhenologyJan2020_Clean.csv"))

