# starter script for analyzing trends in competition recruitment jan 2020
# > can become the main script for analyzing all competition trends (i.e. spring data too), or change file name and code annotation later to reflect competition recruitment only
# author(s): CTW (caitlin.t.white@colorado.edu)
# date create: 21 jan 2020 (this script will be modified over lifetime of project)

# script purpose:
# read in clean competition recruitment dat
# visualize phytometer winter recruitment by covariates (e.g. phyto species , background competitors, treatments, blocks)
# visualize competitor background recruitment by covariates(e.g. comp species, treatments, blocks)
# assess species-species outcomes (who does well/poorly against what competitor(s))?

# notes:




# -- SETUP ----
rm(list = ls()) # clean enviro
#load libraries needed
library(tidyverse)
# change default settings
options(stringsAsFactors = F)
theme_set(theme_bw())

# specify dropbox pathway (varies by user -- ctw can tweak this later)
if(file.exists("~/DropboxCU/Dropbox/USDA-compost/Data/Competition/Competition_EnteredData")){
  ## CTW
  datpath <- "~/DropboxCU/Dropbox/USDA-compost/Data/Competition/Competition_EnteredData"
}else{
  ## probs LMH and AS? (fix if not -- rproj is two levels down in github repo)
  datpath <- "~/Dropbox/USDA-compost/Data/Competition/Competition_EnteredData"
}


# read in clean dats
# clean competitor recruitment density

# clean phytometer recruitment



