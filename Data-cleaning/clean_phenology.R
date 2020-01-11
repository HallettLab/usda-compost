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



# -- SETUP -----
rm(list = ls()) # clear environment
# load needed libraries
library(tidyverse)
library(readxl)
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
pheno_files <- list.files(paste0(datpath, "Phenology/Phenology_EnteredData/"), full.names = T, pattern = "_Phenology_.*[.]xl", ignore.case = T)
# read in treatment key
trtkey <- read.csv(paste0(datpath, "Compost_TreatmentKey.csv"), na.strings = na_vals, strip.white = T)


# -- COMPILE ----
#initiate df for all phenology data
pheno_master <- data.frame()

# 1/10/2020: changing main loop to exclude winter 2020 pheno file for now 
for(i in pheno_files[!grepl("202001", pheno_files)]){
  # phenology dat
  temp_pheno <- read_excel(i, sheet = 1, na = na_vals, trim_ws = T)  
  # photo key
  temp_photokey <- read_excel(i, sheet = 2, na = na_vals, trim_ws = T)  
  print(paste("Compiling", unique(temp_pheno$date)[1], "phenology data and photo key"))
  pct_cols <- c("pct_green", "pct_brown", "pct_bare")
  
  # clean up pheno dat
  # change Ts to 0.01 and adjust pct cover of green (or if "<" noted)
  temp_pheno[pct_cols] <- sapply(temp_pheno[pct_cols], function(x) ifelse(trimws(x)=="T", 0.01,x))
  # check 1: look for < > and adjust based on sum of other cols
  flag1 <- sapply(temp_pheno[pct_cols],function(x) grepl("<|>",x))
  for(f in which(flag1)){
    adjust_pos <- grep("<|>", temp_pheno[f,pct_cols])
    stopifnot(length(adjust_pos)==1) # need to add code if there is more than 1 col that has greater/less than
    temp_pheno[f, pct_cols[adjust_pos]] <- 100-(sum(as.numeric(temp_pheno[f,pct_cols[-adjust_pos]])))
  }
  # make sure all pct cols converted to numeric
  temp_pheno <- mutate_at(temp_pheno, vars(pct_green:pct_bare), as.numeric)
  # check 2: flag rows that sum to +100
  flag2 <- apply(temp_pheno[c("pct_green", "pct_brown", "pct_bare")],1,function(x) sum(as.numeric(x)) > 100)
  for(f in which(flag2)){
    # logic check sum over 100 due to T value
    if(!0.01 %in% temp_pheno[f, pct_cols]){
      stop(paste0("Pct phenology cover for ", temp_pheno$plot[f], temp_pheno$subplot[f], " exceeds 100 and not due to trace value. Review datasheet and data entry"))
    }
    # subtract from pct green
    temp_pheno$pct_green[f] <- with(temp_pheno, 100-(sum(pct_brown[f]+pct_bare[f]))) 
  }  
  
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
  left_join(trtkey) %>%
  left_join(dplyr::select(jan20foto, -c(page_order, order))) %>%
  mutate(yr = as.numeric(substr(date, 1, 4))) %>%
  # reorder cols
  dplyr::select(page, line, recorder, date, yr, plot, block:ppt_trt, pct_litter:notes, photo_subplot:photo_notes) %>%
  # convert trace to 0.01
  mutate_at(grep("pct", names(.)), function(x) gsub("T", "0.01", x))

# write out
write_csv(jan20dat_clean, paste0(datpath, "Phenology/Phenology_CleanedData/Compost_PhenologyJan2020_Clean.csv"))

