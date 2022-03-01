## LOAD PACKAGES ##
library(tidyverse)
library(minpack.lm)
library(nlstools)
library(grid)
library(gridExtra)

## READ IN THE DATA ##
#set data pathway 
datpath <- "~/Desktop/USDA-compost/Data/Competition/"  ##NAT TO EDIT 

#datpath <- "~/Research/Nat_Thesis/USDA-compost/Data/Competition/" ## Carmen's file path - please leave :)

# list files in entered data folder
datfiles <- dats <- list.files(paste0(datpath, "Competition_CleanedData"), full.names = T)

# read in raw trifolium data
dat<- read.csv(paste0(datpath, "Competition_CleanedData/competition_trifolium_seeds_2021_updated.csv"), strip.white = T)



################################################################################
## CLEAN UP THE DATA FOR USE IN MODELS ##
    ## keep only the no nutrient treatment (N) and only TRHI (Control) and LOMU (another background species for comparison)
    ## Also filter out wet treatment so we are left with drought and control treatments


## need to filter by the data available
dat_actual <- dat %>%
  filter(!is.na(mean_density_halfm2)) %>% ## filters for stem density data present
  filter(tot_stems > 0) ## filters for plants that grew

## check each background to see how many rows there are for each.
trhi_avba <- dat_actual %>%
  filter(background == "AVBA")
  ## AVBA has the most rows, so it might be best to use this for the purposes of your project. 
  ## Otherwise I think we will be struggling to fit models with not enough data

trhi_taca <- dat_actual %>%
  filter(background == "TACA")

trhi_homu <- dat_actual %>%
  filter(background == "HOMU")

trhi_lomu <- dat_actual %>%
  filter(background == "LOMU")

trhi_erbo <- dat_actual %>%
  filter(background == "ERBO")


## get data ready for modeling
trhi_avba_m <- trhi_avba %>%
  filter(background == "AVBA") %>% ## filter for avba background
  mutate(TRseedin = tot_stems, ## doing this by the # phytometer stems as seeds in (3 in every case, is not enough info for model)
         TRseedout = tot_seeds, ## total seeds out for phytometers
         AVseedin = max_density_halfm2, ## seeds in for avena
         AVstemdens = mean_density_halfm2) ## avena STEM density

## germination and survival fractions
## Trifolium
tg <- 0.2 ## What I think is germination from Andrew's code in the USDA Climvar repository; DOUBLE CHECK this is correct
eg <- .6 ## erodium germination from Lauren's models_no_facilitation script
ag <- .9 ## avena germination from Lauren's models_no_facilitation script


################################################################################
## Model TRHI in AVBA background ##

## NOTE: we have NOT filtered for nutrient treatment. This script is only parsing by precip treatment
## I chose to do this as we do not have enough data for TRHI otherwise
## we are therefore making the assumption that compost application vs. fertilizer application does NOT make a difference

## model formula
m1TinA <- as.formula(log(TRseedout +1) ~  log(tg*(TRseedin+1)*exp(log(lambda)-log((1+aiT*(TRseedin+1)*tg+aiA*(AVstemdens+1)*ag))))) 
  ## NOTE TO NAT: AVstemdens goes into this model near the end


## Run the model in a for loop
treatments <- unique(trhi_avba_m$ppt_trt) ## create vector of unique precipitation treatments

TRoutput <- as.data.frame(matrix(nrow = 0, ncol = 7)) ## create empty data frame to contain output of for-loop
names(TRoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species") ## rename those columns


for (i in 1:length(treatments)){
  #eg <- ifelse(treatments[i] == "consistentDry" | treatments[i] == "fallDry", 0.64, .6)
  m1out <- nlsLM(m1TinA, start=list(lambda=1, aiT = .01, aiA=.01),
                 lower = c(0, 0, 0), upper = c(200, 1, 1),
                 control=nls.lm.control(maxiter=500), trace=T,
                 data = subset(trhi_avba_m, !is.na(TRseedout) & ppt_trt == treatments[i]))
  
  outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
  names(outreport) = c("estimate", "se", "t", "p")
  outreport$params <- row.names(outreport)
  outreport$treatment <- treatments[i]
  outreport$species <- "Trifolium"
  TRoutput <- rbind(TRoutput, outreport)
}






################################################################################
## Nat's code ## unedited by Carmen


## Competitor species

TRIBACKseedin <- 1250 ## THIS IS NOT NEEDED ANYMORE
 
# ai = competition coefficient; effect of that species on the focal phytometer 
  
## Trifolium model - annual plant competition model found in Hallett et al. 2019
m1T <- as.formula(log(TRseedout +1) ~  log(tg*(TRseedin+1)*exp(log(lambda)-log((1+aiT*(TRseedin+1)*tg+aiTB*(TRIBACKseedin+1)*tg)))))
    ## ahh okay the AVseedin part is where we need the background seeds in. 
    ## THIS IS WHERE WE ARE STILL MISSING INFORMATION ## 

## JUST CONTROL TEST RUN
TRIBACKseedin <- 0

m1T <- as.formula(log(TRseedout +1) ~  log(tg*(TRseedin+1)*exp(log(lambda)-log((1+aiT*(TRseedin+1)*tg+aiTB*(TRIBACKseedin+1)*tg)))))
treatments <- unique(dat_ppt$ppt_trt) ## create a vector of unique treatments (control (XC) and dry (D))

TRoutput <- as.data.frame(matrix(nrow = 0, ncol = 7)) ## make an empty data frame to store the outputs of the for loop below
names(TRoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species") ## rename columns fo empty dataframe

for (i in 1:length(treatments)){ ## start the for loop; run once for each treatment
  
  #eg <- ifelse(treatments[i] == "consistentDry" | treatments[i] == "fallDry", 0.64, .6) ## I think this counts for differential germination in wet and dry...
  ## taking this out as we do not have this info, we'll just use one germination estimate
  
  ## this line below is using maximum likely hood and the model we specified above (m1T) to estimate lambda and alphas
  m1out <- nlsLM(m1T, start=list(lambda=1, aiT = .01, aiTB=.01), ## aiA will change depending on which competitive background we use
                 lower = c(0, 0, 0), upper = c(200, 1, 1),
                 control=nls.lm.control(maxiter=500), trace=T,
                 data = subset(tr1, !is.na(TRseedout) & ppt_trt == treatments[i]))
  
  ## lines below are getting the model outputs into a data frame that matches the format & column names of 
  ## TRouput, the empty data frame we created above
  outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
  names(outreport) = c("estimate", "se", "t", "p")
  outreport$params <- row.names(outreport)
  outreport$treatment <- treatments[i] ## make a treatment column in the final data frame
  outreport$species <- "Trifolium" ## make a column for species in the final data frame
  TRoutput <- rbind(TRoutput, outreport) ## appending the info from each loop into the empty dataframe
} ## end the for loop here

## GOAL: different intraspecific comp. in density + precip treatments
