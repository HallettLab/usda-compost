# clean and compile compost phenology dataset
# author(s): ctw (caitlin.t.white@colorado.edu)
# created: april 2019 (will be modified over project lifetime)

# script purpose:
# iterate through each phenology dataset chronologically
# read in phenology entered data and phenology photo key entered data
# merge both and append to master phenology dataset
# write out to compost dropbox: Phenology/Phenology_EnteredData

# modification jan 2020:
# below main code that compiles growing season pheno, code added to read in and clean/compile winter 2020 phenology (format different than main phenology)



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

# 1/10/2020: changing main code to exclude winter 2020 pheno file for now 
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
# 1/10/2020: for now, code to clean jan 2020 is in its own section as that phenology has veg heights
jan20dat <- read_excel(pheno_files[grepl("202001", pheno_files)], sheet = 1)
jan20foto <- read_excel(pheno_files[grepl("202001", pheno_files)], sheet = 2)

# compile comp plot pheno with mean veg height, mean litter height, and photo info, with notes
jan20dat_clean <- dplyr::select(jan20dat,Plot:Notes) %>%
  gather(met, val, vht_1:lht_4) %>%
  mutate(rep = parse_number(met),
         met = ifelse(grepl("^v", met), "veg", "litter")) %>%
  rename_all(casefold) %>%
  grouped_df(names(.)[grep("pl|pct|not|met", names(.))]) %>%
  summarise(mean = mean(val),
            se = sd(val)/length(sqrt(val))) %>%
  ungroup() %>%
  left_join(trtkey) %>%
  mutate(pct_bare = ifelse(pct_bare == "T", "0.01", pct_bare),
                           pct_bare = parse_number(pct_bare))

ggplot(jan20dat_clean, aes(ppt_trt, mean, col = met)) +
  geom_errorbar(aes(ymax = mean + se, min = mean - se), position = position_dodge(width = 0.3), width = 0.2) +
  geom_point(position = position_dodge(width = 0.3)) +
  labs(y = "Mean height (cm)") +
  scale_color_manual(values = c(veg = "green", litter = "brown")) +
  theme_bw() +
  facet_grid(block ~ nut_trt)

distinct(dplyr::select(jan20dat_clean, plot:pct_bare, block:ppt_trt)) %>%
           gather(met, val, pct_litter:pct_bare) %>%
  ggplot( aes(ppt_trt, val, col = met)) +
  #geom_errorbar(aes(ymax = mean + se, min = mean - se), position = position_dodge(width = 0.3), width = 0.2) +
  geom_point(position = position_dodge(width = 0.3)) +
  labs(y = "Percent") +
  scale_color_manual(values = c(pct_litter = "darkgoldenrod1", pct_green = "green", pct_bare = "brown", pct_brown = "pink")) +
  theme_bw() +
  facet_grid(block ~ nut_trt)




