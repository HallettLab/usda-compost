#########################################
##    Module 4 Ordinations 1           ##
#########################################


## Today we are going to start looking at ordinations, focusing in on distance based ordinations 
## Non-metric Multidimensional Scaling and Principal Coordinates Analysis (PCoA).

## Let's start by loading the relevant packages
library(rgl)
library(MASS)
library(ape)
library(vegan)
library(ecodist)
library(scatterplot3d)
library(vegan3d)
library(magick)
library(ggplot2)
library(grid)
library(tidyr)
library(dplyr)
library(pairwiseAdonis)


################## preliminaries #################################################################
## organize our data

## first we will load all of our datasets
## 1. species composition data
dat = read.csv("model examples/CAstatewide_Local_Woody_Comdata_20190918.csv")

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

## Because we want to analyze community composition in a relative sense, 

## let’s relativize the data and focus on using Bray-Curtis dissimilarities:
comp <- decostand(abundances, "total")
head(comp)

### I am removing some plots here to make life easier, we will come back to this later
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
summary(env3)


## and subset our comp data to match our env data
comp.sites = row.names(env3)
length(comp.sites)

comp3 = selrow<-(is.element(row.names(comp2), as.vector(comp.sites)))
comp3 = comp2[selrow,]
dim(comp3)
head(comp3)



##################  part 1 NMDS    #####################################################################

## First let's make sure we understand how NMDS works

## Consider a single axis representing the abundance of a single species. 
## Along this axis, we can plot the communities in which this species appears, 
## based on its abundance within each.

plot(0:10,0:10,type="n",axes=F,xlab="Abundance of Species 1",ylab="") 
axis(1)
points(5,0); text(5.5,0.5,labels="community A")
points(3,0); text(3.2,0.5,labels="community B")
points(0,0); text(0.8,0.5,labels="community C")

## Now consider a second axis of abundance, representing another species. 
##We can now plot each community along the two axes (Species 1 and Species 2).
plot(0:10,0:10,type="n",xlab="Abundance of Species 1",  ylab="Abundance of Species 2")
points(5,5); text(5,4.5,labels="community A")
points(3,3); text(3,3.5,labels="community B")
points(0,5); text(0.8,5.5,labels="community C")


## Now consider a third axis of abundance representing yet another species.
d=scatterplot3d(0:10,0:10,0:10,type="n",xlab="Abundance of Species 1",
                ylab="Abundance of Species 2",zlab="Abundance of Species 3"); d
d$points3d(5,5,0); text(d$xyz.convert(5,5,0.5),labels="community A")
d$points3d(3,3,3); text(d$xyz.convert(3,3,3.5),labels="community B")
d$points3d(0,5,5); text(d$xyz.convert(0,5,5.5),labels="community C")


## we could keep doing this for as many species as we had in our community
## and it would be very hard to visualize

## The goal of NMDS is to represent the original position of communities in 
## multidimensional space as accurately as possible using a reduced number of 
## dimensions that can be easily plotted and visualized (and to spare your thinking).

## Okay, so lets do it.
# NMDS requires a community-by-species matrix, here we will use the same dataset as last time.

## The function `metaMDS` will take care of most of the distance calculations, iterative fitting, etc, We need
## simply to supply our sample unit by species matrix (dune.rel) and number of reduced dimensions (k=#).
serp.nmds0=metaMDS(comp, k=2)

## let's see what it is doing and what we can specify
help("metaMDS")


##   There are several key things we need to know
##   1. distance. Calculate a dissimilarity matrix. The default is Bray-Curtis dissimilarity, and because we did not
##      specify, that is what was used here. Any other dissimilarity index in vegdist can be used using the
##      argument distance. 
##   2. K the number of dimensions
##   3. NMDS with random starts. NMDS gets easily trapped into local optima, so it is important to start
##      NMDS several times from random starts to be confident that you have found the global solution.
##      metaMDS starts NMDS from several random starts (minimum number is given by try and
##      maximum number by trymax). If a solution is better (has lower stress) than the previous best, it is
##      taken as the new standard. Procrustes analysis compares two solutions, and if they are similar
##      (small residuals) they are considered convergent (indicating a global solution). If you want to be
##      more certain of reaching a global solution, you can compare results from several independent runs.
##   4. There is an older alternative NMDS engine (isoMDS in the MASS package), but it does not have
##      several random starts; for this reason, monoMDS (default in metaMDS) is preferable.
##   5. autoTransform. The function will perform a Wisconsin double standardization (wisconsin) if data
##      maximum values are >9, and if the values are very large (>50), a square-root transformation
##      (sqrt). If you want to have full control over data transformations (you should, and particularly
##      should if you have data types other than typical community data) you can set
##      autotransformation=FALSE. In this case we just ran, the maximum data values did not trigger
##      transformations but if would be listed here it they did.
##   6. The last step in metaMDS is that species scores are added to the final solution. Remember the
##      dissimilarity matrix is based on similarity across sample units. Species scores are weighted
##      averages using function wascores. These are helpful to include in visual representations of the NMDS.

## You saw each iteration of the NMDS until a solution is reached (i.e., stress was minimized after
## some number of reconfiguration of the points in 2 dimensions). You can increase the number of default
## runs (each at a random starting point) using the argument "trymax=##"
serp.nmds0=metaMDS(comp, k=2, trymax=25)

## And then we can look at the NMDS object
serp.nmds0

## The results list the stress value 
## Generally, stress < 0.05 provides an excellent representation in reduced dimensions, 
## < 0.1 is great, < 0.2 is good, and stress > 0.3 provides a poor representation. 
## You always should report the stress of a NMDS, as it is important for readers to interpret your analysis.


## Let's examine a Shepard plot, which shows scatter around the regression between the interpoint
## distances in the final configuration (distances between each pair of communities) against their original
## dissimilarities. The stress of this analysis indicates a pretty good representation and the Shepard plot confirms.
stressplot(serp.nmds0)

## Now we can plot the NMDS
plot(serp.nmds0, type="t")

## this tells us we have a few weird plots that are super unique
## for now let's drop those plots and rerun our
serp.nmds=metaMDS(comp3, k=2)
stressplot(serp.nmds)
plot(serp.nmds, type="t")

    
## this is the basic NMDS plot and is not terribly attractive
## vegan has several functions for making "better" ordinations
ordiplot(serp.nmds, type="points")
ordiplot(serp.nmds, type="text")
ordiplot(serp.nmds, type="text", display="sites")
ordiplot(serp.nmds, type="text", display="species")
ordiplot(serp.nmds, type="points", display="sites")

## for full control call no points
ordiplot(serp.nmds, type="n")
with (env3, ordiellipse(serp.nmds, Province, kind="se", conf=0.99, col="blue", lwd=2,
                            label=TRUE))
orditorp (serp.nmds, display="sites", col="black")

## Another visualization is to show convex hulls connecting the vertices of the points made by these
## treatments on the plot.
ordiplot(serp.nmds, type="p")
with (env3, ordihull(serp.nmds, Province, col="blue", lwd=2, label=TRUE))
orditorp (serp.nmds, display="species", col="black", air=3)

## air argument specifies label padding

## Lastly, let’s consider a continuous variable. In this case, map contour lines 
## onto the plot using the function ordisurf:
ordiplot(serp.nmds,type="n")
with (env3, ordisurf(serp.nmds, MG, add = TRUE))
orditorp(serp.nmds,display="sites",col="grey30",air=0.01)

## Using ordisurf, you might notice that this does more than graphically represent contour lines, it also
## includes a gam function (REML is the restricted maximum likelihood estimate). There are also ways to fit
## environmental factors to NMDS results using envfit. These give P-values and are based on permutations,
## but consider them as post-hoc tests because they fit environmental (or other explanatory) vectors onto
## an ordination rather than starting from the data matrix (e.g., PERMANOVA). They also can consider
## categorical groups (factors) such as Management. 
## CONSIDER THIS A confirmatory visualization ONLY!! 
## Rely on hypothesis testing based directly on resemblance matrices (like PERMANOVA).

##SCALE!!!!!
soils = as.data.frame(env3[,24:36])
summary(soils)
head(soils)
soils2 = scale(soils)
head(soils2)
ef<-envfit(serp.nmds, soils, permu=999)
ef

## This function is useful though for adding environmental vectors to you NMDS
ordiplot(serp.nmds, type="points", display="sites")
plot(envfit(serp.nmds,soils))
## it doesn't really work with categorical variables and is best for continuous data

### as you can see with big data sets, its hard to see anything. 
### I like to instead graph mean and SE (or SD) of the treatments

## first we need to get our scores out to plot - note this will work with PCAs as well
NMDS1 = scores(serp.nmds, choices=c(1), display=c("sites"))
NMDS2 = scores(serp.nmds, choices=c(2), display=c("sites"))

##addin them to env object
env3$NMDS1 = NMDS1
env3$NMDS2 = NMDS2
head(env3)

### function for calculating Standard Error
calcSE <- function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

## calculating means and SE for each province
meanNMDS1 = as.data.frame(aggregate(env3$NMDS1~env3$Province, FUN=mean))
names(meanNMDS1) = c("prov", "NMDS1")
meanNMDS2 = as.data.frame(aggregate(env3$NMDS2~env3$Province, FUN=mean))
names(meanNMDS2) = c("prov", "NMDS2")
seNMDS1 = as.data.frame(aggregate(env3$NMDS1~env3$Province, FUN=calcSE))
names(seNMDS1) = c("prov", "NMDS1se")
seNMDS2 = as.data.frame(aggregate(env3$NMDS2~env3$Province, FUN=calcSE))
names(seNMDS2) = c("prov", "NMDS2se")

##merging all the pieces
df = merge(meanNMDS1, meanNMDS2, )
df = merge(df, seNMDS1)
df = merge(df, seNMDS2)

## plot of mean and SE of scores  ##
ggplot(data = df,aes(x = NMDS1,y = NMDS2,shape = prov)) + 
  geom_point(aes(size=5)) + 
  ylim(-0.3, 0.3)+
  theme_bw() +
  theme(axis.line.x.top = element_line(colour = 'black', size=0.5, linetype='solid'),
       axis.line.x.bottom = element_line(colour = 'black', size=0.5, linetype='solid'),
       axis.line.y.left = element_line(colour = 'black', size=0.5, linetype='solid'),
       axis.line.y.right = element_line(colour = 'black', size=0.5, linetype='solid'))+
  geom_errorbar(aes(ymin = (NMDS2-NMDS2se),ymax = (NMDS2+NMDS2se), width = 0.025)) + 
  geom_errorbarh(aes(xmin = (NMDS1-NMDS1se),xmax = (NMDS1+NMDS1se), height = 0.025))+guides(size="none")+
  scale_shape_manual(values=c(15,16,17,18),name="Province", label=c("klamath", "NorCal", "SoCal", "Sierra"))+
  annotate("text", x=-0.25, y=0.30, label = "Stress = 0.16", cex=4)

### does this match our analysis from last week?
##one issue with adonis is that it doesn't do multiple comparisons
adonis2(comp3 ~ as.factor(env3$Province),  permutations = 999, method = "bray")

## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!

pair.mod<-pairwise.adonis(comp3,as.factor(env3$Province))
pair.mod

## why does it look like Sierra and NorCal don't differ in the plot?

## the last thing we will explore with NMDS is 3D plots
ordiplot3d(serp.nmds, display="sites" )

## why doesn't this work? Remember above we only calculated 2 axes k=2
serp.nmds3=metaMDS(comp3, k=3)
serp.nmds3
## notice the stress goes down

ordiplot3d(serp.nmds3, display="sites" )

### A boxed 'pin' version
ordiplot3d(serp.nmds3, type = "h")

### More user control
pl = ordiplot3d(serp.nmds3, scaling = "symmetric", angle=15, type="n")
points(pl, "points", pch=16, col="red", cex = 0.7)

### put labels in better positions
text(pl, "points", col="blue", pos=3)

### Add species using xyz.convert function returned by ordiplot3d
sp <- scores(serp.nmds3, choices=1:3, display="species", scaling="symmetric")
text(pl$xyz.convert(sp), rownames(sp), cex=0.7, xpd=TRUE)

### Two ways of adding fitted variables to ordination plots
ef2 <- envfit(serp.nmds3, soils2, choices = 1:3)
### 1. use argument 'envfit'
ordiplot3d(serp.nmds3, envfit = ef2)

### 2. use returned envfit.convert function for better user control
pl3 <- ordiplot3d(serp.nmds3)
plot(pl3$envfit.convert(ef2), at = pl3$origin)

### envfit.convert() also handles different 'choices' of axes
pl3 <- ordiplot3d(serp.nmds3, choices = c(1,3,2))
plot(pl3$envfit.convert(ef2), at = pl3$origin)

### vegan::ordiXXXX functions can add items to the plot
pl4 <- with(env3, ordiplot3d(serp.nmds3, col = 1:4, pch=16))
with(env3, ordiellipse(pl4, Province, draw = "poly", col = 1:4, alpha = 60))
with(env3, ordispider(pl4, Province, col = 1:4, label = TRUE))

## Bonus code - 3D animated plot

#Create x,y refs 
data.x <- scores(serp.nmds3, choices=c(1), display=c("sites"))
data.y <- scores(serp.nmds3, choices=c(2), display=c("sites"))
data.z <- scores(serp.nmds3, choices=c(3), display=c("sites"))  

#Plot 
plot3d(data.x, data.y, data.z, col=1:4) 

#Animate by spinning on Y & Z axes 
play3d(spin3d(axis=c(0,1,1), rpm=3), duration=10)

## make an animated gif of your ordination
movie3d(spin3d(axis=c(0,1,1), rpm=3), duration=10, 
        dir ="C:/Users/mspas/Dropbox/EEOB 230 Fall2022 Instructor/Module 4"
      , movie = "NMDS")



#############        Part 2 Principal Coordinates Analysis (PCoA).    ######################

## NMDS is best when your goal is to preserve the ordering relationships among objects, 
## particularly into two or three dimensions. However, it does not preserve the exact 
## distance relationships among objects and is computationally intensive; if these considerations are
## important, you might use a PCoA. It should be reserved for cases where no Euclidean measure 
## is appropriate (e.g., lots of ## shared zeros) and you prefer to use a similarity measure.


## Let’s compute a matrix of Bray-Curtis similarities among sites, and subject this matrix to PCoA.
serp.bray<-vegdist(comp3)
serp.pcoa<-cmdscale(serp.bray, k=(nrow(comp3)-1), eig=TRUE)

## Let’s compute a matrix of Bray-Curtis similarities among sites, and subject this matrix to PCoA.
## If the metric is non-euclidean (as in our case), then the PCoA may produce several negative eigenvalues in
## addition to the positive ones. In most applications, this does not affect the representation of the first
## several axes. You will still receive a warning message, though, when this occurs. You also will receive a
## warning about species scores not being available; there is a way to project weighted averages of species
## abundances on a PCoA plot using the function wascores, we will do that too.

## Check to see if negative eigenvalues affect the interpretation of the first several axes
serp.pcoa

## We will discuss eigenvalues more next week as we move further into eigenvalue-based ordination
## methods; for now just know that each ordination axis has an eigenvalue, it is a measure of the strength of
## an axis, the amount of variation explained by the axis. We typically consider the first several axes, as they
## have greater eigenvalues and thus explain much of the variation in a dataset

## how much does axis 1 explain?
serp.pcoa$eig
## 21%

## Now lets graph it.
ordiplot(scores(serp.pcoa)[,c(1,2)], type="t")
abline(h=0, lty=3)
abline(v=0, lty=3)
## this is very different than our NMDS!!


serp.wa <-wascores(serp.pcoa$points[,1:2], comp3)
text(serp.wa, rownames(serp.wa), cex=0.7, col="red")
##  You can also add convex hulls or ellipses, using the same functions as you did with the NMDS.

## but these are sort of ugly. A lot of packages have built in functions

## let's try and make them look better
DPCOA = pcoa(serp.bray)
biplot(DPCOA)

## project species on the PCoA ordination
biplot(DPCOA, comp3)

## project soils on the PCoA ordination
biplot(DPCOA, soils)

## so just CA and MG? or should we scale our variables
biplot(DPCOA, scale(soils))

## these look a little better, but you sill wouldn't want to publish this
## I personally use a different program, but if you want to do it in R I 
## suggest you continue to explore ggplot

## unlike NMDS, we can use PCoA sores for other analyses and
str(DPCOA)
PCoAscores=as.data.frame(DPCOA$vectors[,1:2])

## plot to check
plot(PCoAscores$Axis.2 ~ PCoAscores$Axis.1 )
## these can be used as a predictor variable in other analyses. NMDS scores cannot!!!



###########  part 4  Exercise ##########
## 1) conduct either a NMDS or PCoA
## 2) Graph it
## 3) Overlay species
## 4) Overlay environmental variables









