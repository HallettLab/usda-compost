#########################################################
###        Module 2 FD, PD, Null models             #####
#########################################################

#### Let's start with loading some of the key packages we will be using today
library(vegan)
library(picante)
# install.packages("FD")
library(FD)
library(phytools)
library(ape)
library(tidyr)
library(pez)
library(matrixStats)
library(NST)
library(data.table)
library(tidyverse)

## Today we will continue working with Harrison CA data.
## We will start be exploring functional and phylogenetic diversity patters
## We will end with some PCM on some other data

## first we will load all of our datasets
## 1. species composition data
dat = read.csv("Community Analysis - Caryn/model examples/CAstatewide_Local_Woody_Comdata_20190918.csv")

## now let's organize and calculate relative abundances
## make a site by species matrix
dat2 = spread(dat, Species, Abundance, fill = 0)
head(dat2)

##now let's check the size
dim(dat2)

## fix is to make the row names the site names and remove the offending columns
row.names(dat2) = dat2$Loc.code
head(dat2)

## and lets remove the Site info columns
abundances = dat2[,6:144]
head(abundances)


## 2. trait data
traits = read.csv("Community Analysis - Caryn/model examples/CAstatewide_Traits_SoilMean20190918.csv")
head(traits)
str(traits)

## important note here, we have trait data measured on different soil types, 
## so we need to do some data indexing to do incorporate this

##subset to relevant trait data (SLA, Wood Density, height (after Anacker and Harrison 2012 Am Nat) )
traits2 = traits[,c("Soil", "Species", "Height.cm", "SLA", "wood_density")]
head(traits2)

## We are going to focus on 3 traits for simplicity
## 1) Height - how tall a species is, influences competition for light
## 2) SLA – specific leaf area, measured as leaf area (cm2) per unit mass (g). This trait is related to a general growth strategy where low SLA species are slow growing and stress tolerant and high SLA species are quick growing resource acquisitive species
## 3) Wood Density - How dense the wood is, higher wood density typically means slower growing more stress tollerant



## 3. environmental data
env = read.csv("Community Analysis - Caryn/model examples/CAstatewide_Enviro_21090918.csv")
head(env)

## let's set the row names
row.names(env) = env$Loc.code

## and look a the data
head(env)
dim(env)

## 4. last we need a phylogentic tree. 
# Load Zanne et al. (2014) phylogeny - DOI:10.1038/nature12872
Ltree <- read.tree("https://github.com/willpearse/pd_app/raw/master/Vascular_Plants_rooted.dated.tre")
str(Ltree)
# plot(Ltree)
# don't plot this unless you have some time to kill

###############    part 1 functional diversity       ##############################################


## First we are going to calculate a few metrics of functional diversity.
## we will be using FD package
help(FD)

## and the main function in the FD package id dbFD
help(dbFD)

##Check to make sure same number of species in both files CRITICAL##
if(dim(abundances)[2]!=dim(traits)[1])stop("error:differentnumberofspeciesin'traits'and'abundances'matrices")

## big nope!

## here are some common issues in running this code
## 1 - names don't match
## 2 - species with no traits
## 3 - traits with no species in comp
## 4 - communities with no species with trait data

## we have a bunch of these common issues to fix that all can be resolved with indexing

## let's look at our traits again real quick
head(traits2)

## so we have traits measured on different soil types, 
## so we need to subset our composition data to only include Serp or Non-Serp soils

##subset woody species to Serp for calculations of FD
WoodSE = subset(dat,Soil =="se" )
head(WoodSE)

# select composition data
WoodSE2 = WoodSE[,c(5:7)]
head(WoodSE2)

## transform to site by species matrix and relativize abundance
WoodSE3 = as.data.frame(spread(WoodSE2, Species,Abundance, fill = 0))
head(WoodSE3)
dim(WoodSE3)

## name rows as regions
row.names(WoodSE3) = WoodSE3$Loc.code

## remove all but species comp data
WoodSE4 = subset(WoodSE3, select = -c(Loc.code) )
head(WoodSE4)
dim(WoodSE4)

## relative abundances
WoodSE5 = decostand(WoodSE4, "total")
head(WoodSE5)

##subset woody species to nonSerp for calculations of FD
WoodNS = subset(dat,Soil =="ns" )
head(WoodNS)

# select composition data
WoodNS2 = WoodNS[,c(5:7)]
head(WoodNS2)

## transform to site by species matrix and relativize abundance
WoodNS3 = as.data.frame(spread(WoodNS2, Species,Abundance, fill = 0))
head(WoodNS3)
dim(WoodNS3)

## name rows as regions
row.names(WoodNS3) = WoodNS3$Loc.code

## remove all but species comp data
WoodNS4 = subset(WoodNS3, select = -c(Loc.code) )
head(WoodNS4)
dim(WoodNS4)

## relativize abundances
WoodNS5 = decostand(WoodNS4, "total")
head(WoodNS5)

#### now that our comp data is subsetted, let's do the same for trait data

##subset trait data to Serp sites for calculations of FD
traitSE = subset(traits2,Soil =="s" )
head(traitSE)
dim(traitSE)

## remove species without complete trait data
traitSE2 = na.omit(traitSE)
head(traitSE2)
dim (traitSE2)

#remove soil column
traitSE3 = traitSE2[,2:5]
head(traitSE3)
dim(traitSE3)

##remove species from Composition data with no trait data
##remove species from trait data that are not present in comp data

##list of species in trait data
SE.species<-sort(unique(traitSE3$Species))

#subsetting columns in comp data frame to only have species with trait data
selcol = (is.element(names(WoodSE5), as.vector(SE.species)))
WoodSE6 = WoodSE5[,selcol]
head(WoodSE6)

###check - row numbers should remain the same, but should be fewer columns
dim(WoodSE5)
dim(WoodSE6)

#remove species with trait data that have no cover data
ind = match(unique(traitSE3$Species),names(WoodSE6))
traitSE4 = traitSE3[-which(is.na(ind==F)),]
head(traitSE4)
dim(traitSE4)

##convert species names into row names for trait data
traitSE5 = traitSE4[,-1]
rownames(traitSE5) = traitSE4[,1]
head(traitSE5)

##Remove communites with poor trait coverage <80% 
##get row sums
R1 = rowSums(WoodSE6)
str(R1)
## check to see what plots are being dropped
drop=as.vector(names(R1[R1<0.80]))
length(drop)

##drop plots
## %in% This operator is used to identify if an element belongs to a vector.
WoodSE7<-WoodSE6[-which(row.names(WoodSE6)%in%names(R1[R1<0.80])),]
dim(WoodSE7)

##Drop species that are not present in at least one community
#remove zero abundance spp from percent cover and trait data
x2<-colSums(WoodSE7)
zero.cols=(!is.element(names(WoodSE7), as.vector(names(x2[x2==0]))))
WoodSE8<-WoodSE7[,zero.cols]
head(WoodSE8)
dim(WoodSE8)

##drop traits for dropped species
traitSE6<-traitSE5[-which(row.names(traitSE5)%in%names(x2[x2==0])),]
dim(traitSE6)

##calculate Functional Diversity - about 200 lines of code to get one function to run
## this is worst case scenario example

#check coherence of number of species in 'traits' and 'abundances'
if(dim(WoodSE8)[2]!=dim(traitSE6)[1])stop("error:differentnumberofspeciesin'traits'and'abundances'matrices")

##RUN FD CODE
resSE0=dbFD(traitSE6,WoodSE8)
str(resSE0)


help(dbFD)
## what I usually run
resSE=dbFD(traitSE6,WoodSE8, corr = c("lingoes"), stand.FRic = TRUE, calc.FDiv=FALSE)

## what did this give us
str(resSE)

## a bunch of metrics and qual.FRic which tells us the quality of the reduced-space representation 
## required to compute FRic. it's like a r-square adn a higher number (up to 1) is better

## let's look at a few things
pairs(resSE)
## notice SP and FRic are correlated
## notice Rao's Q and FDis are VERY correlated

## We can subset out individual items from this list usng the $
cwm=resSE$CWM
head(cwm)

## and then specifically CWM SLA. 
cwmLS = cwm$SLA
head(cwmLS)

## we could also use two dollar signs
cwmLS2=resSE$CWM$SLA
head(cwmLS2)

## One important thing to notice here is that the FD package doesn't calculate FD for
## individual traits so we will need to do some sub-setting

### subset out individual traits
traits_H=subset(traitSE6,select = +Height.cm)
traits_S=subset(traitSE6,select = +SLA)
traits_W=subset(traitSE6,select = +wood_density)

##RUN FD CODE
resHse=dbFD(traits_H,WoodSE8,  stand.FRic = TRUE, calc.FDiv=FALSE)
resSse=dbFD(traits_S,WoodSE8,  stand.FRic = TRUE, calc.FDiv=FALSE)
resWse=dbFD(traits_W,WoodSE8,  stand.FRic = TRUE, calc.FDiv=FALSE)

##write out CSV file with FD calculations
write.csv(resSE, "StatewideFD_MV_SE_3traits.csv")
write.csv(resHse, "StatewideFD_H_SE.csv")
write.csv(resSse, "StatewideFD_SLA_SE.csv")
write.csv(resWse, "StatewideFD_WD_SE.csv")

## Since we have already subsetted our data, let's calculate phylogenetic diversity as well

##################  part 2 phylogenetic diversity ##################################################

##isolate species list
LspeciesSE = row.names(traitSE6)
sppSE <-as.character(LspeciesSE)

##capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x}

sppSE2 = firstup(sppSE)

# Build tree
Ltree2 <- congeneric.merge(Ltree, sppSE2)
Ltree3 <- keep.tip(Ltree2, sppSE2)

## simply typing in the name of your tree (“Ltre3”) 
Ltree3
## tells us some important info about: the number of tip and nodes, 
## whether it is rooted (reflects the most basal ancestor of the tree) 
## or unrooted (does not imply a known ancestral root), 
## and if it has branch length information, 
## or if all branch lengths are the same. ##
## We can use the unroot function to see what an unrooted tree looks like, 
## which we will do after we learn some of the basics of plotting trees.

## First we are going to go through a few ways to visualize phylogenies
plot(Ltree3)

#You can adjust the text size so you can read the species names 
plot(Ltree3, cex=0.7)

## This gives us the simplest view – refered to as a "phylogram" 
##We can also change how we view the phylogeny by specifying one of the 
##following types "cladogram", "fan", "unrooted", "radial"
plot(Ltree3,type="cladogram",cex = 0.9)
plot(Ltree3,type="fan",cex = 0.9)
plot(Ltree3,type="unrooted",cex = 0.9)
plot(Ltree3,type="radial",cex = 0.9)


## We can also change the direction of the phylogeny by specifying one of the 
## following directions “rightwards" (the default), "leftwards", "upwards", and "downwards"
## Let’s stick with an upwards phylogram for now. 
plot(Ltree3, direction = "upwards", cex = 0.9)

# #Let’s quickly compare our tree to an unrooted tree
Ltree4=unroot(Ltree3)
plot(Ltree4, direction = "upwards", cex = 0.9)
## Can you see the difference?

## Now we want to test our key assumption that closely related species have more similar traits. 
## K values of 1 correspond to a Brownian motion process 
## K values closer to zero correspond to a random or convergent pattern of evolution
## K values greater than 1 indicate strong phylogenetic signal and conservatism of traits.


#All of these analyses assume the trait data are in the same order as the phylogeny tip.labels. 
## Let's make sure the trait data are in this order and then measure phylogenetic signal in these data.

## first we need to make our trait names uppercase
head(traitSE6)
dim(traitSE6)

## update regional comp data with uppercase species names
traitSE7 = traitSE6
row.names(traitSE7) = sppSE2
head(traitSE7)


traitSE7 <- traitSE7[Ltree3$tip.label, ]
head(traitSE7)

## Now let's visualize the trait values on the trees by plotting a different color for each trait value. 
## The arguments to the tiplabels function give each trait value a unique color and adjust the size of the trait symbols.
## Let start with column 1 of our trait data and see if it looks like being a canopy or understory species has a phylogenetic signal

## First we plot the tree removing the species names with “show.tip.label=FALSE”
plot(tre, show.tip.label=FALSE, direction = "upwards")

## Now let’s add tip labels corresponding to the traits where “traits[, 1] 
## “ specifies the column of interest in your trait matrix: 1 = Height, 2 = SLA, 3 = WD
tiplabels(pch = 19, col = traitSE7[, 1] )
tiplabels(pch = 19, col = traitSE7[, 2] )
tiplabels(pch = 19, col = traitSE7[, 3] )

## this isn't super obvious with continuous traits, so let's map them onto our phylogeny
## Let's make SLA into a vector with names
SLA = setNames(traitSE7$SLA, rownames(traitSE7))
obj = contMap(Ltree3, SLA )
obj



## Now let's calculate Bloomberg's K for our traits
help("multiPhylosignal")
multiPhylosignal(traitSE7, Ltree3)

##  So let's make our tree fully dichotomous
Ltree4 = multi2di(Ltree3)
plot(Ltree3)
plot(Ltree4)

#Now we can use the mulitphylosignal to test for a significant phylogenetic signal in our four traits
multiPhylosignal(traitSE7, Ltree4)

## K values of 1 correspond to a Brownian motion process 
## K values closer to zero correspond to a random or convergent pattern of evolution
## K values greater than 1 indicate strong phylogenetic signal and conservatism of traits.

## now let's calculate Community phylogenetic patterns

## calculate PD
LSEphydist = cophenetic(Ltree4)

## update comp data with uppercase species names
LsppSE = names(WoodSE8)
LsppSE2 = firstup(LsppSE)

WoodSE9 = WoodSE8
names(WoodSE9) = LsppSE2
head(WoodSE9)

## calculate Mean Pairwise Distance (MPD for each plot)
LSE_ses.mpd.result = ses.mpd(WoodSE9, LSEphydist, 
                             null.model="independentswap", 
                             abundance.weighted=TRUE, runs=999)

## let's see what this gives us. 
str(LSE_ses.mpd.result)
head(LSE_ses.mpd.result)
write.csv(LSE_ses.mpd.result, "LOC_PD_SE.csv")


## this function actually does our null modeling for us! YAY!

## but our FD did not so we need to do it by hand


######################## part 3 null modeling #################################################

## the first step is determining you regional species pool. Realistically we should include serp and non-serp
## species as they are both in the region and we should subset our plots by regions
## this would make things really complicated as we would have to add is traits for species not in our comp data
## we can do this, but it is a more complex model and for simplicity sake today we will not do that

## this is also our first example of a "for loop"

## for loops are a power way in R to do repeated tasks, but they are often computationally intense

## here is a super basic example
for(i in 1:10){
  print("EEOB 230!")
  }
## what this is telling us is for i (which is the number of times), do a thing

## load regional species trait pool to randomize trait values across entire dataset
rsp=traitSE6

##multiply abundances by 100 to get whole numbers for null modeling. Does not affect calculations of FD or CWM

WoodSE8.2 = floor(WoodSE8 * 100 / min(WoodSE8[WoodSE8>0]) )
WoodSE8.2[WoodSE8.2==0] <-1
head(WoodSE8.2)

## set number of reps for null modeling
## scale this down to 3 the first time to make sure it works
## once working scale up to 9999
nreps<-3  

##create matrix to store loop results
## set number of rows to number of plots
dim(WoodSE8.2)
res4 = matrix(data=NA, nrow=114, ncol=nreps)

##run loop
for(i in 1:nreps){
  
  ##Create a null community by randomizing the species matrix, while maintaining sample counts and species richness with each plot
  ##Only works if Relative abundance data is whole numbers##
  
  ##shuffles abundances to different species, but keeps row total and actual abundance numbers
  comm = permatfull(WoodSE8.2, fixedmar = "rows", shuffle = "samp", strata = NULL,
                    mtype = "count", times = 1)
  
  ##Put into usable form for FD package
  mat=data.frame(comm$perm)
  
  ## randomly draw species from RSP and name to match composition dataset
  Ntraits=rsp[sample(nrow(rsp), (dim(mat)[2]), replace=FALSE),]
  ## order traits alphabetically		  
  Ntraits=Ntraits[order(rownames(Ntraits)), ]
  
  ###Calculate Null functional diversity for each plot####
  res2 = dbFD(Ntraits, mat, calc.FRic = FALSE, calc.FDiv = FALSE, print.pco = FALSE, messages = FALSE)
  res3 = as.matrix(res2$FDis)
  
  ## store each run in a matrix, set number of rows to the 
  ## number of your plots
  
  res4[,i]=res3
} 


##double check loop, should be columns equal to nreps
head(res4)
dim(res4)

##add plot names
row.names(res4)=row.names(WoodSE8.2)


##load package for row stats
library(matrixStats)


## create matrix to store final results
resFD = data.frame(cbind(meanFDis=rowMeans(res4), 
                          sdFDis = rowSds(res4),
                         upFDis=rowQuantiles(res4,probs=0.95), 
                         loFDis=rowQuantiles(res4,probs=0.05)))

head(resFD)

## save results to a csv file
write.table(resFD,file="resFD.csv",row.names=TRUE,col.names=TRUE, sep=",")


## this gives us our Null values, but now we want to calculate a SES
SES =(resSE$FDis - resFD$meanFDis)/(resFD$sdFDis)
SESh = (resFD$upFDis - resFD$meanFDis)/(resFD$sdFDis)
SESl = (resFD$loFDis - resFD$meanFDis)/(resFD$sdFDis)


## lets take a look at our results

## now we can plot how FD changes along an environmental gradient

## first we need to subset our env data to match
##list of plots in comp data. Since we have our plots as row names, we can just use the row.names function
env.sites = row.names(res4)
length(env.sites)

#sub-setting columns in ENV data frame to only have plots with comp data
selrow<-(is.element(row.names(env), as.vector(env.sites)))
env2 = env[selrow,]
dim(env2)
summary(env2)


plot(env2$elev,SES, ylab = "SES_FD" )
## add a line for the null
abline (h=0, col="blue")
## add a line for upper and lower CI
abline(h=(mean(SESh)), col="red")
abline(h=(mean(SESl)), col="red")

## Look at the relationship
test = lm(SES~env2$elev)
summary(test)
abline(test)

### what do we find?


## let's take a quick look a the phylogenetic patterns
plot(env2$elev,LSE_ses.mpd.result$mpd.obs, ylab = "obs_PD" )
test2 = lm(LSE_ses.mpd.result$mpd.obs~env2$elev)
summary(test2)
abline(test2)


## now lets quickly calculate some functional-beta
FD.db <- vegdist(cwm, 'euclidean')

envF = rep(x = "F1", times = 114)
Fb = betadisper(FD.db, envF)
Fb
boxplot(Fb)

## this isn't terribly informative, but we will add in the spatial 
## and environmental variables in future lectures


############################  part 4 PCMs    ###################
## Felsenstein (1985) identified (and, to a large extent, solved) a problem that had previously recognized, 
## but that was underappreciated: and that is, that the data for species cannot be treated as 
## independent data points from the point of view of statistical analysis.

## To learn how to use the contrasts method to fit linear models in R, 
## we'll first need to learn some basics in linear model.

##  load some data for gape width & buccal length (different attributes of the mouth) 
## in different species of fish & we will fit a simple linear regression to the data.
dat0<-read.csv("Centrarchidae.csv",row.names=1)
dat0

plot(dat0[,c("buccal.length","gape.width")],
     xlab="relative buccal length",
     ylab="relative gape width",pch=21,bg="grey",
     cex=1.4)

## now lets fit a simple linear model
fit.ols<-lm(gape.width~buccal.length,data=dat0)
fit.ols
summary(fit.ols)
abline(fit.ols,lwd=2,lty="dashed",col="red")

## However, we cannot forget that when our data are phylogenetic, the assumption independent and identical distribution of the residual error does not hold. 
## Consequently, we need to take the phylogeny into account. One way to do this is by using PICs:
cent.tree<-read.tree("Centrarchidae.tre")
plot(cent.tree)
buccal.length<-setNames(dat0[,"buccal.length"],
                        rownames(dat0))
gape.width<-setNames(dat0[,"gape.width"],rownames(dat0))

### 
help(pic)
pic.bl<-pic(buccal.length,cent.tree)
pic.gw<-pic(gape.width,cent.tree)
fit.pic<-lm(pic.gw~pic.bl+0)
fit.pic

## lets see how our model looks
summary(fit.pic)


## and compare it to our original model 
summary(fit.ols)

## and graph the new information
plot(pic.bl,pic.gw,xlab="PICs for buccal length",
     ylab="PICs for gape width",bg="grey",
     cex=1.4,pch=21)
abline(fit.pic,lwd=2,lty="dashed",col="red")



########################## part 5 exercise  ###################################
# Now that you have gone through the guided tasks, apply them to your own dataset.

#For the exercise on Thursday
# 1) Upload and explore your own data
# 2) Calculate alpha functional diversity indices
# 3) see if your FD differs from the null expectation of random assembly
# 3) Calculate Bloomberg's K for your data 
# 4) Calculate phylogenetic diversity for your data 
# 5) if you have environmental data, look at relationships with CWM, FD, PD

## if you don't have data use the Non-serpentine plot data


