#####################################################
###            Module 5 Ordinations II           ####
#####################################################
## load in the required libraries
library(vegan)
library(tidyr)
library(dplyr)
library(ggfortify)

## start with the usual data cleaning

## first we will load all of our datasets
## 1. species composition data
dat = read.csv("model examples/CAstatewide_Local_Woody_Comdata_20190918.csv")

## now let's organize and calculate relative abundances
## make a site by species matrix
dat2 = spread(dat, Species, Abundance, fill = 0)

## fix is to make the row names the site names and remove the offending columns
row.names(dat2) = dat2$Loc.code

## and lets remove the Site info columns
abundances = dat2[,6:144]
head(abundances)

## Because we want to analyze community composition in a relative sense, 

## let’s relativize the data and focus on using Bray-Curtis dissimilarities:
comp <- decostand(abundances, "total")
head(comp)

### Again I am removing some plots here to make life easier, we will come back to this later
comp2=comp
comp2 = comp2[row.names(comp2) != c("ncwmsse"),]
comp2 = comp2[row.names(comp2) != c("kltmnns"),]
comp2 = comp2[row.names(comp2) != c("ncetnse"),]
comp2 = comp2[row.names(comp2) != c("sccisns"),]
comp2 = comp2[row.names(comp2) != c("ncrisns"),]
dim(comp2)

## 2. environmental data
env = read.csv("model examples/CAstatewide_Enviro_21090918.csv")
head(env)

## let's set the row names
row.names(env) = env$Loc.code

## We are also going to subset our data so that all the plots without soil data are excluded
env2 = na.omit(env)
dim(env2)
summary(env2)

##and subset our env data to match the comp data
##list of plots in comp data. Since we have our plots as row names, we can just use the row.names function
env.sites = row.names(comp2)
length(env.sites)

#sub-setting columns in ENV data frame to only have plots with comp data
selrow<-(is.element(row.names(env2), as.vector(env.sites)))
env3 = env2[selrow,]
dim(env3)


## and subset our comp data to match our env data
comp.sites = row.names(env3)
length(comp.sites)

comp3 = selrow<-(is.element(row.names(comp2), as.vector(comp.sites)))
comp3 = comp2[selrow,]
dim(comp3)

## now that everything is matched, we are going to further subset to start some "real analyses" and compare 
## drivers of variation among provinces

KLenv = env3[env3$Province == "kl",]
dim(KLenv)
NCenv = env3[env3$Province == "nc",]
dim(NCenv)
SNenv = env3[env3$Province == "sn",]
dim(SNenv)
SCenv = env3[env3$Province == "sc",]
dim(SCenv)


### now we match the composition data

## Kalamath
KL.sites = row.names(KLenv)
length(KL.sites)
KLcomp = selrow<-(is.element(row.names(comp3), as.vector(KL.sites)))
KLcomp =comp3[selrow,]
dim(KLcomp)
## remove species with no plots
x2<-colSums(KLcomp)
zero.cols=(!is.element(names(KLcomp), as.vector(names(x2[x2==0]))))
KLcomp2<-KLcomp[,zero.cols]
head(KLcomp2)
dim(KLcomp2)

## nor Cal
NC.sites = row.names(NCenv)
length(NC.sites)
NCcomp = selrow<-(is.element(row.names(comp3), as.vector(NC.sites)))
NCcomp =comp3[selrow,]
dim(NCcomp)
## remove species with no plots
x2<-colSums(NCcomp)
zero.cols=(!is.element(names(NCcomp), as.vector(names(x2[x2==0]))))
NCcomp2<-NCcomp[,zero.cols]
head(NCcomp2)
dim(NCcomp2)

## Sierra Nevada
SN.sites = row.names(SNenv)
length(SN.sites)
SNcomp = selrow<-(is.element(row.names(comp3), as.vector(SN.sites)))
SNcomp =comp3[selrow,]
dim(SNcomp)
## remove species with no plots
x2<-colSums(SNcomp)
zero.cols=(!is.element(names(SNcomp), as.vector(names(x2[x2==0]))))
SNcomp2<-SNcomp[,zero.cols]
head(SNcomp2)
dim(SNcomp2)

## SoCal
SC.sites = row.names(SCenv)
length(SC.sites)
SCcomp = selrow<-(is.element(row.names(comp3), as.vector(SC.sites)))
SCcomp =comp3[selrow,]
dim(SCcomp)
## remove species with no plots
x2<-colSums(SCcomp)
zero.cols=(!is.element(names(SCcomp), as.vector(names(x2[x2==0]))))
SCcomp2<-SCcomp[,zero.cols]
head(SCcomp2)
dim(SCcomp2)


######################  Part 1 PCA    ##################################


## There are many PCA methods available for R. 
## Here we will focus on the most basic one (prcomp), 
## a default function from the R base package.

## Let's look at our env data. We don't actually want all the variables, 
## so let's index the ones we want and make a new object
head(KLenv)
soils = KLenv[,c(24,27:32,35)]
head(soils)

## REMEMBER! we can't use a PCA on categorical variables

## now let's run the PCA
Spca = prcomp(soils)
biplot(Spca)

## what is the problem here?
## REMEMBER! scaling your data is SUPER important

Spca2 = prcomp(soils, scale = TRUE)
biplot(Spca2)
## so now that we have a good looking PCA what does it all mean?

##first we want to know how much of our data is explained by each axis
summary(Spca2)

## PC1 explains 40.0%
## PC2 explains 25.4%
## and so on

## if we just look at the object we get a lot of info
str(Spca2)

## rotation tells us our loadings - higher values mean a variable loads more strongly on that axis
## there is one for each variable - 
Spca2$rotation

## x gives us our site scores, useful for 1) analysis, and 2) graphing elsewhere
## you can also use the "scores" function
## so lets extract PC1 for analysis
PC1 = scores(Spca2, choices=c(1), display=c("sites"))
PC1

##changing the choices gives us different axes
PC2 = scores(Spca2, choices=c(2), display=c("sites"))
PC2

## changing the display gives us the "rotation" or loadings
PC1s = scores(Spca2, choices=c(1), display=c("species"))
PC2s = scores(Spca2, choices=c(2), display=c("species"))

PC1s
PC2s

## let's calculate diversity in each plot
div = diversity(KLcomp2)
div

##and ask how diversity changes with our environmental gradient
plot(PC1,div)
lm1 = lm(div~PC1)
abline(lm1)
summary (lm1)

## what about PC2?
plot(PC2,div)
lm2 = lm(div~PC2)
abline(lm2)
summary (lm2)

## this is one of the primary uses of PCAs you will see in community ecology
## often people have many env variables and need to summarize them to examine
## how a variable of interest changes along a env gradient

## lastly we will make a "nicer" plot
## making some nicer PCA plots
#install.packages('ggfortify')
pca=autoplot(Spca2, loadings = TRUE, label = FALSE, shape = TRUE, label.size = 5, loadings.label = TRUE,
         loadings.colour = "black",loadings.label.size=5, loadings.label.colour="black",size=3,
         loadings.label.vjust = -0.6)+
      theme_bw(base_size = 16)+ theme(text = element_text(size = 20)) +
      theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))+
  geom_hline(yintercept=0, color="black",size=1.0, lty="dashed")+
  geom_vline(xintercept=0, color="black",size=1.0, lty="dashed")
pca
## NOTICE THE BRITISH SPELLING!!!!

## you can play around with this more on your own



##################  part 2 CA and DCA   #################################################

## CA and DCA are not really used in ecology anymore as there are more robust analyses
## here we will quickly run some code and give them a look, but we will not spend much time on them

## correspondence analysis is very similar to PCA, but instated we will be using our species data
Sca = cca(KLcomp2)
Sca

## let's see how it looks
plot(Sca)
summary(Sca)



## now lets detrend it 
## here's how to conduct a DCA so that you can see how it works, 
## but DO NOT USE DCA:

Sdca = decorana(KLcomp2)
plot(Sdca)
summary (Sdca)

## NEVER DO A DCA for your analysis ##


################# Part 3 CCA and RDA   #####################################################

## remember, CCA assumes that our species have a unimodal relationship with our variables.
## RDA assumes a linear relationship
## let's start by looking at that assumption
pairs(c(KLcomp2[,5:16],soils))

## so we can;t really do this for all of our species
## so how do we chose between a CCA and RDA?

## we use a DCA!!!

## The length of first DCA axis:
## > 4 indicates heterogeneous dataset on which unimodal methods should be used (CCA), 
##  < 3 indicates homogeneous dataset for which linear methods are suitable (RDA) 
## between 3 and 4 both linear and unimodal methods are OK.
summary (Sdca)

## Looks like CCA is best, but we will run both and compare your results

## CCA does not try to display all variation in the data, 
## but only the part that can be explained by the used constraints. 
## Consequently, the results are strongly dependent on the set of constraints 
## The shotgun method is to use all environmental variables as constraints, 
## however, rely on non-metric multidimensional scaling (metaMDS) and environmental 
## interpretation after analysis (envfit, ordisurf). 
## CCA is a good choice if the user has clear and strong a priori hypotheses on constraints 

head(KLenv)
## let's start with a CCA of several environmental variables
KLcca = cca(KLcomp2 ~  Soil + Aspect + elev + bio1 + bio12 + OM + CA + MG, data = KLenv )

## what do we find?
KLcca
plot(KLcca)

## remember, Scaling
plot(KLcca, scaling = 1)
## this emphasizes relationships among sites
plot(KLcca, scaling = 2)
## this emphasizes relationships among species

## do you see the difference?
par(mfrow = c(1,2))
plot(KLcca, scaling = 1)
## this emphasizes relationships among sites
plot(KLcca, scaling = 2)

## let's go back to trying to interpret this plot
## what do we find?
summary(KLcca)

## 1) partition of scaled chi-squared Constrained = explained by variable,
##    Unconstrained = explained by residuals
## 2) Importance of componenents - which axes contribute
## 3) Scores: LC or WA? LC give best fit if you have no noise in constraining
##    variables, WA scores good with noisy data

## So Is our Constrained inertia our R-sqauare?

## NO, we need to bootstrap to get the real r-square

## how much variation does it explain?
RsquareAdj(KLcca)


## we can then test for significance of the model by permutation
## if it is not significant, 
## it doesn't matter how much of the variation is explained
## Is it significant?
anova(KLcca, permutations = how(nperm=999))

## we can also do a permutation test by axis
anova(KLcca, by = "axis", permutations = how(nperm=999))
## or by terms
anova(KLcca, by = "terms", permutations = how(nperm=999))
## this will help us interpret our CCA and we can see some variable are not significant

## we can use model selection instead of picking variables we think are important

KLcca2 = ordistep(cca(KLcomp2 ~ 1, data = KLenv ),
                 scope=formula(KLcca),
                 direction = "forward",
                 permutations = how(nperm=999))

KLcca2
## this tells us that we really only need a few variables 
## KLcomp2 ~ Soil + bio1 + MG + bio12 

## so lets rerun our model with these and see how it looks

KL.p.cca =  cca(KLcomp2 ~ Soil + bio1 + MG + bio12, data = KLenv )
summary(KL.p.cca)
anova(KL.p.cca, permutations = how(nperm=999))
anova(KL.p.cca, by = "terms", permutations = how(nperm=999))

##What about r-squared

## original model
RsquareAdj(KLcca)

## new model
RsquareAdj(KL.p.cca)

## but we have limited our model selection using "scope" to just what was in our initial model
## let's change that
## subset out all the junk
KLenvM = KLenv[,11:57]
head(KLenvM)

## the "." here mean everything in the dataset
KLcca3 = cca(KLcomp2 ~ ., data = KLenvM )

## now let rerun our model selection
KLcca4 = ordistep(cca(KLcomp2 ~ 1, data = KLenvM ),
                  scope=formula(KLcca3),
                  direction = "forward",
                  permutations = how(nperm=999))

## so our model that was selected is 
# KLcomp2 ~ MG + bio10 + Shrub.Cover + ndvi_aug + bio4 + bio11 + PH + bio6 +
##  Mean.shrub.ht..m + Bare.Soil + Rocks 


KL.p.cca2 =  cca(KLcomp2 ~ MG + bio10 + Shrub.Cover + ndvi_aug + bio4 + bio11 + PH + bio6 +
                   Mean.shrub.ht..m + Bare.Soil + Rocks , data = KLenvM )
plot(KL.p.cca2)
KL.p.cca2
anova(KL.p.cca2, permutations = how(nperm=999))
anova(KL.p.cca2, by = "terms", permutations = how(nperm=999))

## now let's compare this to our original model
RsquareAdj(KL.p.cca)
## new model
RsquareAdj(KL.p.cca2)

### ok, now lets compare our CCA to a RDA and see if it tells us anything different
KLrda = rda(KLcomp2 ~ MG + bio10 + ndvi_aug + Shrub.Cover + elev + PH + bio3 + bio15 + CaMgRatio + K + Mean.shrub.ht..m + Bare.Soil + Scat +Olsen_P , data = KLenvM)
 
## what do we find?
plot(KLrda)

## it's not that different than our CCA, just rotated differently
par(mfrow = c(1,2))
plot(KL.p.cca2)
plot(KLrda)

## Let's see our summary
summary(KLrda)
# anova(KLrda, perm = 999)
## how much variation does it explain?
RsquareAdj(KLrda)


## we can then test for significance of the model by permutation
## if it is not significant, 
## it doesn't matter how much of the variation is explained
## Is it significant?
anova(KLrda, permutations = how(nperm=999))

## we can also do a permutation test by axis as before and use the same code above to do model selection

#Now let's compare the Kalamth with SoCal'
SCenvM = SCenv[,11:57]
head(SCenvM)

## the "." here mean everything in the dataset
SCcca = cca(SCcomp2 ~ ., data = SCenvM)

## now let rerun our model selection
SCcca2 = ordistep(cca(SCcomp2 ~ 1, data = SCenvM ),
                  scope=formula(SCcca),
                  direction = "forward",
                  permutations = how(nperm=999))

SC.p.cca =  cca(SCcomp2 ~ elev + ndvi_aug + Shrub.Cover + MG + Herb.Cover + Moss + bio11 + 
                  bio17 + NA. + Litter.Depth + Aspect + ndvi_mar + K , data = SCenvM )
plot(SC.p.cca)
SC.p.cca
anova(SC.p.cca, permutations = how(nperm=999))
anova(SC.p.cca, by = "terms", permutations = how(nperm=999))
RsquareAdj(SC.p.cca)

## let's comare teh two regions
par(mfrow = c(1,2))
plot(KL.p.cca2)
plot(SC.p.cca)

## so how do we interpret this?
## the env explains 18% of the variation in community comp in the Kalamath
## and it is primarily driven by KLcomp2 ~ MG + bio10 + Shrub.Cover + ndvi_aug + bio4 + bio11 + PH + bio6 +
## Mean.shrub.ht..m + Bare.Soil + Rocks

## the env explains 15% of the variation in community comp in the SoCal
## and it is primarily driven by elev + ndvi_aug + Shrub.Cover + MG + Herb.Cover + Moss + bio11 + 
## bio17 + NA. + Litter.Depth + Aspect + ndvi_mar + K



## that's it for today. Questions?
## in future lectures we will learn about dbRDA and variation partitioning

###################     part 4  exercise     #######################
## 1) use the other two provinces
## 2) run a PCA on the environmental variables and regress against diversity 
## 3) determine if you should run a CCA or RDA
## 4) run a CCA or RDA































 
