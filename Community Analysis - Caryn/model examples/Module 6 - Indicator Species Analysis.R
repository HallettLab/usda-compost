#####################################################
##          Module 6 - Cluster ISA SIMPER          ##
#####################################################

## let's load in our libraries for the day
## please install any packages that are new for you
library(vegan)
library(cluster)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(ggdendro)
library(dendextend)
library(indicspecies)
library(tidyr)

setwd("/Volumes/Memorex USB/eeob230")
## start with the usual data cleaning

## first we will load all of our datasets
## 1. species composition data
dat = read.csv("./model examples/CAstatewide_Local_Woody_Comdata_20190918.csv")

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
env = read.csv("./model examples/CAstatewide_Enviro_21090918.csv")
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

## 3. trait data
traits = read.csv("./model examples/CAstatewide_Traits_SoilMean20190918.csv")
head(traits)
str(traits)

##subset to relevant trait data (SLA, Wood Density, height (after Anacker and Harrison 2012 Am Nat) )
traits2 = traits[,c("Soil", "Species", "Height.cm", "SLA",  "DispersalSyndrome")]
head(traits2)
dim(traits2)
str(traits2)
traits2$DispersalSyndrome = as.factor(traits2$DispersalSyndrome)

## We are going to focus on 3 traits for simplicity
## 1) Height - how tall a species is, influences competition for light
## 2) SLA – specific leaf area, measured as leaf area (cm2) per unit mass (g). This trait is related to a general growth strategy where low SLA species are slow growing and stress tolerant and high SLA species are quick growing resource acquisitive species
## 3) Dispersal Syndrome - How does the species disperse

## and just do a few species for visualization
traits3=traits2[1:35,]
dim(traits3)

###########################  part 1 cluster analysis    #############
## Many times in biology we want to classify things into groups
## In some cases, we may not have data with habitats or groups assigned apriori. 
## If this is the case, then a clustering approach can be used to determine groups. 
## As an example we will ask if there trait groups in the Serpentine data trees


## for hierarchical clustering we will be using the "hclust" function in the cluster package
## let's start by taking a look at it
help(hclust)

## ok, so we need a dissimilarity matrix as our input
## let's start with our first question about habitat types

## instead of using vegdist in VEGAN, we can use the "daisy" function in CLUSTER
help(daisy)
trtD = daisy(traits3[,3:5], metric = "gower", stand = TRUE)
str(trtD)

## now lets run a cluster analysis of the env variables and start by looking at different methods
## We will try "single", "complete", "average" and "ward.D"


#hierarchical clustering
## single linkage 
trtCL1 = hclust(trtD, method="single")
#plot dendogram
plot(trtCL1)
## remember this is prone to chaining 

## Complete linkage
trtCL2 = hclust(trtD, method="complete")
#plot dendogram
plot(trtCL2)

## Average linkage UPGMA
trtCL3 = hclust(trtD, method="average")
#plot dendogram
plot(trtCL3)

## ward's method
trtCL4 = hclust(trtD, method="ward.D")
#plot dendogram
plot(trtCL4)

## as we mentioned different packages do things differently
## cluster also includes ward.D2 which is the same as AGNES ward
## what does AGNES do that hclust doesn't?

## let's find out
trtCL5 = agnes(trtD, diss=TRUE, method = "ward")
## one benefit is that you don't have to run a dissimilariy matrix first
plot(trtCL5)


## the banner plot is a horizontal barplot visualizing the hierarchical clustering
## note labels are only printed in the banner plot when the number of 
## observations is limited less than 35 


## ok we have our cluster analysis done, how many clusters are there?
## we could go through by hand and look or we can use PAM to tell us

## what pam does is Partitioning (clustering) of the data into k clusters 
## around medoids", a more robust version of K-means.

# Calculate silhouette width for many k using PAM
sil_width <- c(NA)

## the number here is for the number of potential clusters

for(i in 2:20){
  
  pam_fit <- pam(trtD,
                 diss = TRUE,
                 k = i)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}

# Plot silhouette width (higher is better, but not always biologically meaningful)

plot(1:20, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:20, sil_width)

### here our highest value arises from 10 clusters 

## we can then set the number of clusters in PAM
pamfit = pam(trtD, diss=TRUE, k=10)
summary (pamfit)
plot(pamfit)

## Range of SC	Interpretation
## 0.71-1.0	A strong structure has been found
## 0.51-0.70	A reasonable structure has been found
## 0.26-0.50	The structure is weak and could be artificial
## < 0.25	No substantial structure has been found

## this plot tells us how good our clusters are


info = pamfit$medoids

## let's visualize our clusters
tsne_obj <- Rtsne(trtD, is_distance = TRUE, perplexity = 10)
##t-SNE is a method for constructing a low dimensional embedding 
## of high-dimensional data, distances or similarities

tsne_data <- as.data.frame(tsne_obj$Y )
names(tsne_data)=(c("X", "Y")) 
cluster = factor(pamfit$clustering)

ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster))

#obtain clusters
clust_trt = pamfit$clustering

## Plot our tree of clusters
dends = as.dendrogram(trtCL4)
dends = color_branches(dends, k=10)
plot(dends, horiz = F)

## add clusters back to trait data
traits3$cluster = clust_trt
head(traits3)
write.csv(traits3, "clusters.csv")

## so clustering can be a useful tool, but it HAS to still be biologically meaningful.
## be very careful when conducting and interpreting cluster analyses
## this is not just something where you can "run the code"


############  Part 2. SIMPER and ISA     ################

## Biologists often ask how species can act as indicators of different habitat types. 
## Two of the most popular analyses available with which to ask this question are SIMPER
## and IndVal (ISA). SIMPER is the older of these methods and is the first demonstrated below,
## but remember, SIMPER is rapidly becoming out of date

## SIMPER stands for 'similarity percentages,' and provides the contribution of 
## each species to Bray-Curtis dissimilarities. 
## The 'simper' function in R performs a pairwise comparison procedure among two or more groups. 
## Here, we'll assess similarity percentages across soil type in the Kalamath region.
sim <- simper(KLcomp2, KLenv$Soil)
summary(sim)

## Of particular interest in the output are the average contribution of each species to 
## overall dissimilarity, their standard deviations, and the cumulative contribution to 
## the dissimilarities. The average abundances per group (A and B) are provided to 
##interpret which group the species indicates.

## Though SIMPER has been used widely, it is difficult to interpret. 
## The difficulty arises in that the similarity percentages primarily 
## reflect the most variable species both within and across groups, 
## instead of distinctive species across groups. 

## instead we can look at the INDICSPECIES package to conduct our indicator species analysis. 
## A newer method, called indicator value (IndVal) analysis, can better separate each species' 
## predictive value across habitat types. 

##For this we can apply the 'multipatt' function from this package:
help("multipatt")
serp.ind = multipatt(KLcomp2, KLenv$Soil, func="IndVal.g", control = how(nperm=999))
summary(serp.ind)

## IndVal.g is based on the product of each species' frequency and relative abundance in each group,
## adjusted for sample size and square-root transformed. These values can be found in the 'stat' 
## column of the summary output. Probability values are also provided, based on a permutational 
## procedure. By default, the function only returns summaries for species with P < 0.05. 
## We can easily see every species' maximum indicator value by adjusting the alpha level:
summary(serp.ind, alpha=1)


###############  part 3  exercise    ###############
## using a different province
## 1) run one cluster analysis topographic variables and one only using edaphic variables (soil pH, soil nutrient content)
## 2) determine what species are indicators of serpentine and non serpentine in that province
colnames(env) # Province, Aspect_NS, elev, Aspect, Slope 

##subset to relevant trait data (SLA, Wood Density, height (after Anacker and Harrison 2012 Am Nat) )
env4 = env[,c("Province", "Aspect_NS", "elev", "Aspect",  "Slope")]
env4 = env4 %>%
  filter(env4$Province =="kl")
env4$Aspect_NS = as.factor(env4$Aspect_NS)
typeof(env4$Aspect_NS)
#instead of using vegdist in VEGAN, we can use the "daisy" function in CLUSTER
help(daisy)
envD = daisy(env4[,2:5], metric = "gower", stand = TRUE) # why isn't it taking in the aspect_NS?
str(envD)

## now lets run a cluster analysis of the env variables and start by looking at different methods
## We will try "single", "complete", "average" and "ward.D"


#hierarchical clustering
## single linkage 
single_env = hclust(envD, method="single")
#plot dendogram
plot(single_env)
## remember this is prone to chaining 

## Complete linkage
complete_env = hclust(envD, method="complete")
#plot dendogram
plot(complete_env)

## Average linkage UPGMA
ave_env = hclust(envD, method="average")
#plot dendogram
plot(ave_env)

## ward's method
ward_env = hclust(envD, method="ward.D")
#plot dendogram
plot(ward_env)

## as we mentioned different packages do things differently
## cluster also includes ward.D2 which is the same as AGNES ward
## what does AGNES do that hclust doesn't?

## let's find out
agnes_env = agnes(envD, diss=TRUE, method = "ward")
## one benefit is that you don't have to run a dissimilariy matrix first
plot(agnes_env)


#----------------------------pam---------------------------#

# Calculate silhouette width for many k using PAM
sil_width <- c(NA)

## the number here is for the number of potential clusters

for(i in 2:20){
  
  pam_fit <- pam(envD,
                 diss = TRUE,
                 k = i)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}

# Plot silhouette width (higher is better, but not always biologically meaningful)

plot(1:20, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:20, sil_width)

### here our highest value arises from 10 clusters 

## we can then set the number of clusters in PAM
pamfit = pam(envD, diss=TRUE, k=10)
summary (pamfit)
plot(pamfit)

## Range of SC	Interpretation
## 0.71-1.0	A strong structure has been found
## 0.51-0.70	A reasonable structure has been found
## 0.26-0.50	The structure is weak and could be artificial
## < 0.25	No substantial structure has been found



## let's visualize our clusters
env_obj <- Rtsne(envD, is_distance = TRUE, perplexity = 10)
##t-SNE is a method for constructing a low dimensional embedding 
## of high-dimensional data, distances or similarities

env_data <- as.data.frame(env_obj$Y )
names(env_data)=(c("X", "Y")) 
cluster = factor(pamfit$clustering)

ggplot(aes(x = X, y = Y), data = env_data) +
  geom_point(aes(color = cluster))

#obtain clusters
clust_trt = pamfit$clustering

## Plot our tree of clusters
dends = as.dendrogram(ward_env)
dends = color_branches(dends, k=10)
plot(dends, horiz = F)

## add clusters back to trait data
env4$cluster = clust_trt
head(env4)
write.csv(traits3, "clusters.csv")

## Indicator Species Analysis ##

## so clustering can be a useful tool, but it HAS to still be biologically meaningful.
## be very careful when conducting and interpreting cluster analyses
## this is not just something where you can "run the code"

site.ind = multipatt(KLcomp2, KLenv$Site, func="IndVal.g", control = how(nperm=999))
summary(site.ind)
unique(KLenv$Site)

