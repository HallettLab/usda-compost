# native causal inference and other test code

# -- TEST DAGGITY ASSUMPTIONS -----
testImplications <- function( covariance.matrix, sample.size ){
  library(ggm)
  tst <- function(i){ pcor.test( pcor(i,covariance.matrix), length(i)-2, sample.size )$pvalue }
  tos <- function(i){ paste(i,collapse=" ") }
  implications <- list(c("Precip trt","Nutrient addition"),
                       c("Precip trt","Herbicide"),
                       c("Nutrient addition","Herbicide"))
  data.frame( implication=unlist(lapply(implications,tos)),
              pvalue=unlist( lapply( implications, tst ) ) )
  
}

mydag <- daggity('dag {
bb="-3.38,-2.024,2.638,1.911"
FEMI [outcome,pos="-0.541,1.411"]
forb [pos="-0.552,-0.116"]
grass [pos="0.036,-1.524"]
herbicide [exposure,pos="1.749,1.360"]
nut_trt [exposure,pos="-2.879,-0.574"]
ppt_trt [exposure,pos="2.137,-0.609"]
forb -> FEMI
grass -> FEMI
grass -> forb
herbicide -> FEMI
herbicide -> forb
herbicide -> grass
nut_trt -> FEMI
nut_trt -> forb
nut_trt -> grass
ppt_trt -> FEMI
ppt_trt -> forb
ppt_trt -> grass
}
'
)

mydag <- (dagitty('dag {
  bb="-2.968,-3.059,2.203,2.38"
  "FEMI" [outcome,pos="-0.827,-0.860"]
  "forb" [pos="-0.279,0.934"]
  "grass" [pos="1.082,-0.303"]
  "nut_trt" [exposure,pos="0.751,-2.559"]
  "ppt_trt" [exposure,pos="-2.468,0.769"]
  herbicide [exposure,pos="1.703,1.880"]
  "forb" -> "FEMI"
  "grass" -> "FEMI"
  "grass" -> "forb"
  "nut_trt" -> "FEMI"
  "nut_trt" -> "forb"
  "nut_trt" -> "grass"
  "ppt_trt" -> "FEMI"
  "ppt_trt" -> "forb"
  "ppt_trt" -> "grass"
  herbicide -> "FEMI"
  herbicide -> "forb"
  herbicide -> "grass"
}')
)
plot(mydag)
impliedConditionalIndependencies(mydag)
testdat <- neighborhood[c("plot", "herbicide", "ppt_trt", "nut_trt", "forb", "grass", "BRCA")]
testdat$herbicide <- as.numeric(testdat$herbicide)
testdat$ppt_trt <- as.numeric(testdat$ppt_trt)
testdat$nut_trt <- as.numeric(testdat$nut_trt)
str(testdat)
localTests(x = mydag, data = testdat, test = "cis.chisq")

testcorr <- lavaan::lavCor( testdat )
localTests( mydag, sample.cov=testcorr, sample.nobs=nrow( testdat ) )


# i had this in the native_modeling script and not sure why.. may have been typed into wrong script in prep for comps exam?
# notes:
latlong_proj <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"
boulder_coords <- sp::SpatialPoints(coords = data.frame(lat = 40.0101190898142, long = -105.24232781534424))
sp::proj4string(boulder_coords) <- latlong_proj
sfrec_coords <- sp::SpatialPoints(coords = data.frame(lat = 39.253229939065406, long = -121.31334268578593))
sp::proj4string(boulder_coords) <- latlong_proj
sp::spDists(boulder_coords, sfrec_coords, longlat = T)
