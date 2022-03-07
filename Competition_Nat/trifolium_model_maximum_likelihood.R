## LOAD PACKAGES ##
library(tidyverse)
library(minpack.lm)
library(nlstools)
library(grid)
library(gridExtra)

install.packages("dplyr")    # alternative installation of the %>%
library(dplyr) # needs to be run every time you start R and want to use %>%

## READ IN THE DATA ##
#set data pathway 
datpath <- "~/Desktop/USDA-compost/Data/Competition/"  ##NAT TO EDIT 

#datpath <- "~/Research/Nat_Thesis/USDA-compost/Data/Competition/" ## Carmen's file path - please leave :)

# list files in entered data folder
datfiles <- dats <- list.files(paste0(datpath, "Competition_CleanedData"), full.names = T)

# read in raw trifolium data
dat<- read.csv(paste0(datpath, "Competition_CleanedData/competition_trifolium_seeds_2021.csv"), strip.white = T)



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
  filter(background == "AVBA") ##changed AVBA to TACA to test NaN issues 
## AVBA has the most rows, so it might be best to use this for the purposes of your project. 
## Otherwise I think we will be struggling to fit models with not enough data

#trhi_taca <- dat_actual %>%
  #filter(background == "TACA")

#trhi_homu <- dat_actual %>%
  #filter(background == "HOMU")

#trhi_lomu <- dat_actual %>%
  #filter(background == "LOMU")

#trhi_erbo <- dat_actual %>%
  #filter(background == "ERBO")


## get data ready for modeling
trhi_avba_m <- trhi_avba %>%
  filter(background = "AVBA") %>% ## filter for avba background ##I changed AVBA to TACA to see if max_density_halfm2 is the problem...I get an error I think has to do with the fact TACA has 10 rows and AVBA has 9 rows...
  
###REMOVED: used for the trhi_avba_m table 
#mutate(TRseedin = as.numeric(tot_stems), ## doing this by the # phytometer stems as seeds in (3 in every case, is not enough info for model)
         #TRseedout = tot_seeds, ## total seeds out for phytometers
         #AVseedin = max_density_halfm2, ## seeds in for avena; didn't actually use this
         #AVstemdens = mean_density_halfm2) ## avena STEM density

  ###Simulate more data for AVBA and TRHI to fix lambda. 
  ###use means and sd from seed_in (TRHI), seed_out (TRHI), mean_density (AVBA); make sure to find these for the separate treatment
  #XC - TRseed_in (TRHI) = tot_stems, TRseed_out = tot_seeds (TRHI), AVstemdens = mean_density (AVBA)
  #W - seed_in (TRHI), seed_out (TRHI), mean_density (AVBA)
  #D - seed_in (TRHI), seed_out (TRHI), mean_density (AVBA)
  #convert into a new table  

###START SIMULATION
TRseed_in_new = as.numeric(trhi_avba_m$tot_stems) ##numeric 

###JEFF ADVICE FOR SIMULATION

#vector (just one column from the original dataframe)
 TRseed_out_new.mean=mean(trhi_avba_m$tot_seeds)
 TRseed_out_new.sd=sd(trhi_avba_m$tot_seeds)
 TRseed_in_new.mean=mean(TRseed_in_new) 
 TRseed_in_new.sd=sd(trhi_avba_m$tot_stems)
 AVstemdens_new.mean=mean(trhi_avba_m$mean_density_halfm2)
 AVstemdens_new.sd=sd(trhi_avba_m$mean_density_halfm2)
 
 
#N1= TRseed_in
#N2= AVseed_out
 
  # Set arbitrary parameters
 N<-100 # number of data points you want to simulate
 N1.mean<-2.363636;  N1.std<-0.67419;  N1.max<-100  # Trif stats (intraspecific)
 N2.mean<-13.181818; N2.std<-14.972;  N2.max<-100 #Avena stats (interspecific)
 fitness.mean<-80; fitness.std<-60 #fitness = how many seeds produced? estimate avg. Trif seed produced (per plant)
 
 
 ## Simple random 
 dat.simple = tibble(
   N1 = round(runif(N, 0,N1.max)) # seems to me it makes more sense to have a full density gradient rather than just data surrounding the mean.
   ,N2 = round(runif(N, 0,N2.max)) 
   #,N2 = round(rnorm(N, N2.mean,N2.std))# replace 30 with mean and 5 with std
   ,fitness_data = round(rnorm(N, fitness.mean, fitness.std))
   ,treatment = sample(c("XC","W","D"), size=N, replace=T) # randomly assign treatments
 ) #checking data...shows no correlation
 dat.simple
 
 
 # plot relative to N1
 ggplot(dat.simple, aes(x= N1, y=fitness_data)) +
   geom_point() +
   geom_smooth(method = "lm", formula = y ~ x + I(x^2))
 
 # plot relative to N2
 ggplot(dat.simple, aes(x= N2, y=fitness_data)) +
   geom_point() +
   geom_smooth(method = "lm", formula = y ~ x + I(x^2))
 
 
 # Set arbitrary parameters
 N<-100 # number of data points
 lambda<-50 # true lambda (fitness, no competition means exponentially increasing pop.)
      ###in Bev-Holt model = growth rate accounting for effects of comp.  
      ###1vs0 = intraspecific comp. promotes coexistence 
aiT<-.05 # competition coefficient - you could play around with different intra- and inter-sp strengths of competition
 #how was .05 calculated??
aiA<-.01 # competition coefficient; weaker interaction 
N1.mean<-2.363636;  N1.std<-0.67419; N1.max<-100
N2.mean<-13.181818; N2.std<-14.972; N2.max<-100
 
 bev1 <- function(N1,N2, lamdba, aiT,aiA) {lambda / (1+aiT*N1+aiA*N2)} # Beverton-Holt
 
 # simulate data
 datXC = tibble(
   N1 = round(runif(N, 0,N1.max))
   ,N2 = round(runif(N, 0,N2.max))
   ,fitness_sim = bev1(N1,N2, lambda, aiT, aiA) # fitness, using BH
   ,fitness_data = fitness_sim  + rnorm(N,0, 5) # add some noise; adding random amount of error
   ,treatment = "XC")
 
 # simulate data
 datW = tibble(
   N1 = round(runif(N, 0,N1.max))
   ,N2 = round(runif(N, 0,N2.max))
   ,fitness_sim = bev1(N1,N2, 70, aiT, aiA) # fitness, using BH
   ,fitness_data = fitness_sim  + rnorm(N,0, 5) # add some noise; adding random amount of error
   ,treatment = "W")
 
 #lambda = performance in the ABSENCE of comp. 
 
 # simulate data
 datD = tibble(
   N1 = round(runif(N, 0,N1.max))
   ,N2 = round(runif(N, 0,N2.max))
   ,fitness_sim = bev1(N1,N2, 20, aiT, aiA) # fitness, using BH
   ,fitness_data = fitness_sim  + rnorm(N,0, 5) # add some noise; adding random amount of error
   ,treatment = "D")
 
 dat<-bind_rows(datXC,datW,datD)
 
 ggplot(dat, aes(x= N1, y=fitness_data, color=treatment)) +
   geom_point()+
 geom_smooth(aes(col=treatment), method="lm", se=F)
 
 ggplot(dat, aes(x= N2, y=fitness_data, color=treatment)) +
   geom_point()+
   geom_smooth(aes(col=treatment), method="lm", se=F)
###ADD LINES
 dat
 dat$fitness_data[dat$fitness_data<0] <- 0 #lambda values below 0 are just equal to 0
 
#  ,fitness_data=ifelse(treatment %in% c("W"), fitness_data=fitness_data+5,fitness_data=fitness_data)#W add effect of add to fitness and add error; rnorm randomly
#   ,fitness_data=ifelse(treatment %in% c("D"), fitness_data=fitness_data-5,fitness_data=fitness_data#D add effect of subtract from fitness and add error; rnorm assign treatments))
 
 #example ifelse: mydata$y = ifelse(mydata$x3 %in% c("A","D") ,mydata$x1*2,mydata$x1*3)
 #command Z (undo)
 # plot relative to N1
 ggplot(dat, aes(x= N1, y=fitness_data)) +
   geom_point() +
   geom_smooth(method = "lm", formula = y ~ x + I(x^2))+
   stat_smooth(method = "nls",
               formula = y ~ a/(1+b*x), # BH
               method.args = list(start = list(a = lambda, b = .1)),
               se = FALSE, 
               color='orange')
 #randomness about "this" shape (due to Bev-Holt model, as oppossed to just normal.distrib. randomness)
 #lines = show different theoretical models 
 
 # plot relative to N2
 ggplot(dat, aes(x= N2, y=fitness_data)) +
   geom_point() +
   geom_smooth(method = "lm", formula = y ~ x + I(x^2))+
   stat_smooth(method = "nls",
               formula = y ~ a/(1+b*x), # BH
               method.args = list(start = list(a = lambda, b = .1)),
               se = FALSE, 
               color='orange')
 #graphi curved; fitness of N1 related to pop.size N1
 #graph straight line; fitness N1 is not correlated to pop.size to N2
 
# SO, you could try fitting your model to the data in dat

###END SIMULATION 

## germination and survival fractions
## Trifolium
tg <- 0.2 ## What I think is germination from Andrew's code in the USDA Climvar repository; DOUBLE CHECK this is correct
#eg <- .6 ## erodium germination from Lauren's models_no_facilitation script
ag <- .9 ## avena germination from Lauren's models_no_facilitation script

## Model TRHI in AVBA background ##

## NOTE: we have NOT filtered for nutrient treatment. This script is only parsing by precip treatment
## we are therefore making the assumption that compost application vs. fertilizer application does NOT make a difference

## model formula
m1TinA <- as.formula(log(fitness_data +1) ~  log(tg*(N1+1)*exp(log(lambda)-log((1+aiT*(N1+1)*tg+aiA*(N2+1)*ag))))) 
## AVstemdens goes into this model near the end
# ai = competition coefficient; effect of that species on the focal phytometer 

## Run the model in a for loop
#treatments <- unique(trhi_avba_m$ppt_trt) ## create vector of unique precipitation treatments

TRoutput <- as.data.frame(matrix(nrow = 0, ncol = 7)) ## create empty data frame to contain output of for-loop
names(TRoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species") ## rename those columns


for (i in 1:length(treatments)){
  #eg <- ifelse(treatments[i] == "consistentDry" | treatments[i] == "fallDry", 0.64, .6)
  m1out <- nlsLM(m1TinA, start=list(lambda=1, aiT = .01, aiA=.01),
                 lower = c(0, 0, 0), upper = c(200, 1, 1),
                 control=nls.lm.control(maxiter=500), trace=T,
                 data = subset(dat, treatment == treatments[i]))
  
  outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
  names(outreport) = c("estimate", "se", "t", "p")
  outreport$params <- row.names(outreport)
  outreport$treatment <- treatments[i]
  outreport$species <- "Trifolium"
  TRoutput <- rbind(TRoutput, outreport)
}


##GROWTH PLOT

annual_growth_competition <- function(N0, lambda, aiT, aiA, aiAA, aiTA, tg, ag, period){
  # Argument definitions:
N0 = 26 #starting population size
aiAA = 0.001 #affect of Avena vs. Avena 
aiTA = 0.03 #affect of Trifolium vs. Avena 
period = 50
  # lambda = per capita growth rate in the absence of competition
  # alpha = strength of intraspecific competition
  # g = germination rate (between 0 - 1)
  
  t <- c() # create an empty vector to store timesteps in
  Nt_aiT <- c() # create an empty vector to store pop sizes in
  Nt_aiAA <- c() # create an empty vector to store pop#2 sizes in
  fec_aiT <- c() # create an empty vector to store fecundities in
  fec_aiAA <- c() # create an empty vector to store pop#2 fecundities in
  
  Nt_aiT[1] <- N0 # first item in pop size vector is starting pop size
  Nt_aiAA[1] <- N0 # first item in pop size vector is starting pop size
  
  t[1] <- 1 # first item in timestep vector is 1
  fec_aiT[1] <- lambda/(1 + aiT*N0*tg+ aiA*N0*ag) # calculating the fecundity at time 1
  fec_aiAA[1] <- lambda/(1 + aiAA*N0*ag+ aiTA*N0*tg) # calculating the fecundity pop#2 at time 1
  
  for(i in 2:period){
    t[i] <- i # store each time step
    
    fec_aiT[i] <- lambda/(1 + aiT*Nt_aiT[i-1]*tg+aiA*Nt_aiAA[i-1]*ag) # store each new fecundity
    fec_aiAA[i] <- lambda/(1 + aiAA*Nt_aiAA[i-1]*ag+aiTA*Nt_aiAA[i-1]*tg)
    
    Nt_aiT[i] <- (((1-tg))*Nt_aiT[i-1])+tg*Nt_aiT[i-1]*fec_aiT[i-1] # store each pop size
    Nt_aiAA[i] <- (((1-ag))*Nt_aiAA[i-1])+ag*Nt_aiAA[i-1]*fec_aiAA[i-1] # store each pop size
  }
  dat <- cbind.data.frame(t, Nt_aiT, Nt_aiAA) # bind together time and pop size vectors
  
  ggplot(data = dat) + geom_line(aes(x = t, y = Nt_aiT), color = "green") + geom_line(aes(x = t, y = Nt_aiAA), color = "orange") +
    theme_classic() + xlab("Years") + ylab("Population Size")
  
}



#HALLETT PAPER PLOT CODE
#ggplot(TRoutput, aes(x=treatment, y=aiA)) + geom_bar(stat="identity") + facet_wrap(~species)
#ggplot(TRoutput, aes(x=treatment, y=aiT)) + geom_bar(stat="identity") + facet_wrap(~species)


###JEREMY ADVICE FOR CHECKING DATA

#Make historgrams and facet wrap for treatments to see what is wrong
#Make histogram of everything that's a continuous variable
#Make scatterplots (bi-plots)

ggplot(trhi_avba_m, aes(x=background, y=tot_seed_mass)) +
  geom_boxplot() +
  facet_wrap(~ppt_trt)

ggplot(trhi_avba_m, aes(x=background, y=tot_stem_mass)) +
  geom_boxplot() +
  facet_wrap(~ppt_trt)

ggplot(trhi_avba_m, aes(x=background, y=tot_seeds)) +
  geom_boxplot() +
  facet_wrap(~ppt_trt)

ggplot(trhi_avba_m, aes(x=tot_seeds, y=tot_stems, color = ppt_trt))+
  geom_point()
  
ggplot(trhi_avba_m, aes(x=tot_seed_mass, y=tot_stem_mass, color = ppt_trt))+
  geom_point()

ggplot(trhi_avba_m, aes(x=tot_seed_mass, y=mean_density_halfm2, color = ppt_trt))+
  geom_point()

ggplot(trhi_avba_m, aes(x=tot_stem_mass, y=mean_density_halfm2, color = ppt_trt))+
  geom_point()

