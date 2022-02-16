library(SSDM)
library(raster)
library(parallel)
detectCores()

#load environmental rasters
env<-load_var(path="~/Documents/jpv/projet_CEBC/projet_bocage_CEBC/env_variables",format=".asc", verbose=FALSE)
env
projection(env)=CRS("+init=epsg:2154") #assign a projection system, here 2154 for Lambert 93

occ=read.csv("~/Documents/jpv/projet_CEBC/projet_bocage_CEBC/DATA_SPECIES_79/data_species_79_L93_for_analysis_testSSDM.csv")
head(occ)

SSDM<-stack_modelling(c("CTA", "SVM"), occ, env, rep=1, ensemble.thresh=0, Xcol="lon", Ycol="lat", Spcol="species", cv="holdout",cv.param=c(0.7,1),method="pSSDM", cores=6, verbose=F)

plot(SSDM@diversity.map, main="SSDM\nfor Amphibians")



