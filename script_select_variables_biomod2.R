######################################################
#script written by Jean-Pierre Vacher 4 January 2022
#updated 24 February 2022
######################################################

#charge the packages
x=c("here","sf","dismo","ggmap", "rgdal", "rgeos", "maptools", "dplyr", "tidyr", "tmap", "raster","mapdata","sp","spdep","colorRamps","ggplot2","gridExtra","usdm","virtualspecies") #create a list with names of packages
lapply(x, library, character.only=TRUE) #loop that read all the packages from the list

################################
#SELECTION OF THE BEST MODEL####
################################

###################################
#WITH LOG-TRANSFORMED VARIABLES####
###################################

#Read species data####
data.sp=read.csv("DATA_SPECIES_79/data_species_grid100m_79.csv") #read filtered species data

#Read stacked raster files####
#create an object that gather all the raster layers of the environmental variables
#attention pour générer la liste de rasters, il ne faut pas prendre les fichiers .xml qui sont générés au moment de l'export des couches raster en ascii. Donc prévoir un nouveau dossier dans lequel on range les seulement les fichiers .asc et .prj
list.rasters=(list.files("env_variables_log_complete", full.names=T,pattern=".asc")) #provide the path of the directory where the raster layers saved
rasters=stack(list.rasters) #stack the layers
rasters #check what it looks like
projection(rasters)=CRS("+init=epsg:2154") #assign a projection system, here 2154 for Lambert 93
#pdf(file="plot_raster_variables.pdf") #save as a pdf file
jpeg(file="figures/plot_raster_variables_log_transformed.jpg", height=17, width=17, units="cm", res=200) #save as a jpg file
plot(rasters) #plot all the rasters on one panel
dev.off() #save the file in the working directory

#Check collinearity between variables####
#with Pearson's####
#library("virtualspecies")

rasters.reduced=removeCollinearity(rasters, multicollinearity.cutoff=0.7, plot=T, select.variables=T, sample.points=F) #use a Pearson's test 
# - No multicollinearity detected in your data at threshold 0.7
rasters.reduced #check the results

#with VIF####
env.var=read.table("env_variables.txt", h=T) #read the table
env.var.log=env.var
str(env.var.log)
cols=colnames(env.var.log[,c(2:10)])
env.var.log[cols]=lapply(env.var.log[cols]+1, log) #on transforme en log (en ajoutant +1 pour éviter log(0)) les variables de distance
#library(usdm)
vif(env.var.log[,2:10])
#            Variables      VIF
#1       per.pasture  1.400109
#2      per.housing  1.111231
#3       per.forest 10.212438
#4        dens.road  1.125261
#5 dens.forest.edge  3.230023
#6       dens.hedge  1.491303
#7       dens.river  1.123069
#8      dist.forest  7.240955
#9        dist.pond  1.254532

#VIF values from 5 to 10 are considered critical. So we can run vif again without per.forest and dist.forest

vif(env.var.log[,c(2:3,5:8,10)])
#           Variables      VIF
#1        per.pasture 1.455162
#2      per.housing 1.075617
#3        dens.road 1.108439
#4 dens.forest.edge 1.058071
#5       dens.hedge 1.440763
#6       dens.river 1.103888
#7        dist.pond 1.278821

#Now it seems ok

vifcor(env.var.log[,2:10], th=.7) #VIF with a correlation threshold at r = 0.7
#2 variables from the 9 input variables have collinearity problem: 
#
#per.forest dist.forest 
#
#After excluding the collinear variables, the linear correlation coefficients ranges between: 
#min correlation ( dist.pond ~ dens.forest.edge ):  -0.005402574 
#max correlation ( dens.hedge ~ per.pasture ):  0.4561517 
#
#---------- VIFs of the remained variables -------- 
#  Variables      VIF
#1      per.pasture 1.408473
#2      per.housing 1.076855
#3        dens.road 1.122766
#4 dens.forest.edge 1.036860
#5       dens.hedge 1.400765
#6       dens.river 1.075766
#7        dist.pond 1.242366

rasters.selected=subset(rasters, c("per_pastures","per_housing","dens_roads","dens_forest_edges","dens_hedges","dens_river","dist_pond")) #create an object that contains all the variables that were selected
#plot(rasters.selected) #check how it looks like

###########################
#Salamandra salamandra####
##########################
ss=data.sp[data.sp$species=="Salamandra_salamandra",] #select the data for the focal species
ss=ss[!duplicated(ss),] #remove duplicated occurrences to avoid pseudoreplication

#check graphically how it looks like
#plot(rasters.selected$dens_hedges) #plot the raster of hedges density
#points(ss[,2:3], pch=21, bg=alpha("yellow",.5), lwd=.3,cex=.5) #plot the points of occurrence of S. salamandra. The numbers of columns correspond to the longitude and latitude.

#select variables with PCA####
library(ade4)
envdata=data.frame(raster::extract(x=rasters, y=cbind(ss[,2], ss[,3]))) #build a dataframe that contains the values of each variable for the species locations
envdata=na.omit(envdata) #remove the rows with na
str(envdata) #check how it looks like
pca1<-dudi.pca(envdata, scannf = F, nf=2) #run the pca
round(pca1$eig/sum(pca1$eig)*100,2)
plot(pca1$li[,1:2]) #plot the pca to see what it looks like (no outliers)
summary(pca1$li)
s.corcircle(pca1$co) #plot the correlation circle
#per_forest and dist_forest negatively correlated
#per_pastures and dens_hedges correlated
#dist_pond and dens_hedges+per_pastures negatively correlated
#dens_roads and per_housing correlated

#select variables with GLM####
presvals=raster::extract(rasters.selected,as.data.frame(ss)[,2:3]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.11433  -0.32603   0.05096   0.32478   0.80935  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.584947   0.073659   7.941 3.20e-15 ***
#dist_pond           -0.049673   0.011213  -4.430 9.91e-06 ***
#length_river         0.011345   0.006768   1.676   0.0938 .  
#length_forest_edges  0.084966   0.004011  21.185  < 2e-16 ***
#length_hedges       -0.005622   0.005044  -1.115   0.2651    
#length_roads         0.029904   0.004913   6.086 1.37e-09 ***
#per_housing          0.077908   0.012158   6.408 1.81e-10 ***
#per_pastures         0.008251   0.006055   1.363   0.1731    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1915011)
#
#Null deviance: 535.53  on 2152  degrees of freedom
#Residual deviance: 410.77  on 2145  degrees of freedom
#AIC: 2561.3
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.10086  -0.32561   0.06104   0.32612   0.79228  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.605933   0.063620   9.524  < 2e-16 ***
#dist_pond           -0.052674   0.010422  -5.054 4.69e-07 ***
#length_forest_edges  0.086953   0.003887  22.368  < 2e-16 ***
#length_roads         0.028661   0.004805   5.964 2.87e-09 ***
#per_housing          0.074237   0.012020   6.176 7.83e-10 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1916908)
#
#Null deviance: 535.53  on 2152  degrees of freedom
#Residual deviance: 411.75  on 2148  degrees of freedom
#AIC: 2560.5
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          9 2561.323
#reduced.model.present  6 2560.464

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


###########################
#Lissotriton helveticus####
###########################
lh=data.sp[data.sp$species=="Lissotriton_helveticus",] #select the data for the focal species
lh=lh[!duplicated(lh),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(lh)[,2:3]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.2222  -0.4461   0.1810   0.3684   0.8776  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          1.054405   0.056491  18.665  < 2e-16 ***
#dist_pond           -0.110723   0.008601 -12.873  < 2e-16 ***
#length_river        -0.001188   0.006618  -0.179    0.858    
#length_forest_edges  0.050363   0.004007  12.568  < 2e-16 ***
#length_hedges        0.018502   0.004595   4.026 5.83e-05 ***
#length_roads         0.004133   0.005035   0.821    0.412    
#per_housing          0.087567   0.011040   7.932 3.19e-15 ***
#per_pastures         0.006849   0.005378   1.274    0.203    
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.2005061)
#
#Null deviance: 615.24  on 2598  degrees of freedom
#Residual deviance: 519.51  on 2591  degrees of freedom
#AIC: 3209.3
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variable within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.2240  -0.4505   0.1815   0.3709   0.8756  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          1.074093   0.053870  19.939  < 2e-16 ***
#dist_pond           -0.112830   0.008315 -13.570  < 2e-16 ***
#length_forest_edges  0.050166   0.003933  12.754  < 2e-16 ***
#length_hedges        0.021342   0.004092   5.216 1.98e-07 ***
#per_housing          0.087609   0.010351   8.464  < 2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.2004401)
#
#Null deviance: 615.24  on 2598  degrees of freedom
#Residual deviance: 519.94  on 2594  degrees of freedom
#AIC: 3205.4
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         9 3209.270
#reduced.model.present  6 3205.421

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


#########################
#Triturus cristatus####
#########################
tc=data.sp[data.sp$species=="Triturus_cristatus",] #select the data for the focal species
tc=tc[!duplicated(tc),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(tc)[,2:3]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.8588  -0.2939  -0.1402   0.3399   1.1241  
#
#Coefficients:
#    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.899379   0.076182  11.806  < 2e-16 ***
#dist_pond           -0.126380   0.011406 -11.080  < 2e-16 ***
#length_river        -0.032899   0.009474  -3.473 0.000531 ***
#length_forest_edges  0.032136   0.005513   5.830 6.93e-09 ***
#length_hedges        0.009461   0.005756   1.644 0.100486    
#length_roads         0.005996   0.006447   0.930 0.352485    
#per_housing          0.051183   0.015641   3.272 0.001093 ** 
#per_pastures         0.015941   0.006749   2.362 0.018317 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1647093)
#
#Null deviance: 269.54  on 1368  degrees of freedom
#Residual deviance: 224.17  on 1361  degrees of freedom
#AIC: 1425.9
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.8709  -0.2915  -0.1413   0.3434   1.1073  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.944182   0.072472  13.028  < 2e-16 ***
#dist_pond           -0.129832   0.011230 -11.561  < 2e-16 ***
#length_river        -0.031164   0.009435  -3.303 0.000981 ***
#length_forest_edges  0.030999   0.005490   5.647 1.99e-08 ***
#per_housing          0.056703   0.015096   3.756 0.000180 ***
#per_pastures         0.020734   0.006119   3.389 0.000722 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1649911)
#
#Null deviance: 269.54  on 1368  degrees of freedom
#Residual deviance: 224.88  on 1363  degrees of freedom
#AIC: 1426.3
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         9 1425.938
#reduced.model.present  7 1426.289

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "model.present"


#########################
#Triturus marmoratus####
#########################
tm=data.sp[data.sp$species=="Triturus_marmoratus",] #select the data for the focal species
tm=tm[!duplicated(tm),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(tm)[,2:3]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.21129  -0.40479  -0.01392   0.39872   0.99028  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          1.148470   0.065714  17.477  < 2e-16 ***
#dist_pond           -0.145809   0.010049 -14.510  < 2e-16 ***
#length_river        -0.020926   0.008447  -2.477   0.0133 *  
#length_forest_edges  0.044004   0.004813   9.143  < 2e-16 ***
#length_hedges        0.020627   0.005265   3.918 9.23e-05 ***
#length_roads         0.006247   0.005737   1.089   0.2764    
#per_housing          0.108670   0.012198   8.908  < 2e-16 ***
#per_pastures         0.006035   0.006149   0.981   0.3265    
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.2004802)
#
#Null deviance: 500.50  on 2001  degrees of freedom
#Residual deviance: 399.76  on 1994  degrees of freedom
#AIC: 2474.1
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.2166  -0.4058  -0.0156   0.3973   0.9859  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          1.164411   0.062797  18.542  < 2e-16 ***
#dist_pond           -0.147294   0.009704 -15.178  < 2e-16 ***
#length_river        -0.020361   0.008410  -2.421   0.0156 *  
#length_forest_edges  0.043785   0.004810   9.102  < 2e-16 ***
#length_hedges        0.023622   0.004684   5.043 4.99e-07 ***
#per_housing          0.110249   0.011580   9.520  < 2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.2004797)
#
#Null deviance: 500.50  on 2001  degrees of freedom
#Residual deviance: 400.16  on 1996  degrees of freedom
#AIC: 2472.1
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          9 2474.120
#reduced.model.present  7 2472.123

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


#####################
#Triturus blasii####
#####################
tb=data.sp[data.sp$species=="Triturus_blasii",] #select the data for the focal species
tb=tb[!duplicated(tb),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(tb)[,2:3]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.42489  -0.13094  -0.07102  -0.00025   1.06573  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.505969   0.063578   7.958 4.35e-15 ***
#dist_pond           -0.076394   0.009560  -7.991 3.38e-15 ***
#length_river        -0.015632   0.007060  -2.214  0.02702 *  
#length_forest_edges  0.012791   0.004353   2.939  0.00337 ** 
#length_hedges        0.010899   0.004271   2.552  0.01085 *  
#length_roads        -0.002340   0.004961  -0.472  0.63718    
#per_housing          0.030621   0.012261   2.497  0.01266 *  
#per_pastures        -0.002949   0.005095  -0.579  0.56288    
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.07499314)
#
#Null deviance: 90.909  on 1099  degrees of freedom
#Residual deviance: 81.893  on 1092  degrees of freedom
#AIC: 282.24
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.42176  -0.13324  -0.06952  -0.00270   1.06681  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.495511   0.060128   8.241 4.84e-16 ***
#dist_pond           -0.075140   0.009113  -8.245 4.69e-16 ***
#length_river        -0.016030   0.007025  -2.282  0.02270 *  
#length_forest_edges  0.012816   0.004342   2.952  0.00323 ** 
#length_hedges        0.009583   0.003836   2.498  0.01263 *  
#per_housing          0.030110   0.011859   2.539  0.01126 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.07489147)
#
#Null deviance: 90.909  on 1099  degrees of freedom
#Residual deviance: 81.931  on 1094  degrees of freedom
#AIC: 278.76
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         9 282.2410
#reduced.model.present  7 278.7616

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


########################
#Alytes obstetricans####
########################
ao=data.sp[data.sp$species=="Alytes_obstetricans",] #select the data for the focal species
ao=ao[!duplicated(ao),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(ao)[,2:3]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.05079  -0.15459  -0.08315   0.08237   1.03602  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.238913   0.070015   3.412 0.000662 ***
#dist_pond           -0.022356   0.010511  -2.127 0.033606 *  
#length_river         0.011597   0.007443   1.558 0.119434    
#length_forest_edges  0.033549   0.004563   7.352 3.31e-13 ***
#length_hedges        0.005047   0.004662   1.083 0.279175    
#length_roads         0.017884   0.005119   3.494 0.000491 ***
#per_housing          0.218361   0.009252  23.601  < 2e-16 ***
#per_pastures        -0.029285   0.006028  -4.858 1.32e-06 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1212735)
#
#Null deviance: 291.28  on 1410  degrees of freedom
#Residual deviance: 170.15  on 1403  degrees of freedom
#AIC: 1037.4
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.04472  -0.15290  -0.08675   0.07516   1.00808  

Coefficients:
  #Estimate Std. Error t value Pr(>|t|)    
  #(Intercept)          0.267080   0.066901   3.992 6.88e-05 ***
  #dist_pond           -0.025216   0.010336  -2.440 0.014825 *  
  #length_forest_edges  0.034544   0.004456   7.752 1.73e-14 ***
  #length_roads         0.019091   0.005068   3.767 0.000172 ***
  #per_housing          0.217263   0.009214  23.580  < 2e-16 ***
  #per_pastures        -0.025428   0.005477  -4.643 3.76e-06 ***
  #---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  #
  #(Dispersion parameter for gaussian family taken to be 0.1214343)
#
#Null deviance: 291.28  on 1410  degrees of freedom
#Residual deviance: 170.62  on 1405  degrees of freedom
#AIC: 1037.3
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          9 1037.425
#reduced.model.present  7 1037.304

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


########################
#Pelodytes punctatus####
########################
pp=data.sp[data.sp$species=="Pelodytes_punctatus",] #select the data for the focal species
pp=pp[!duplicated(pp),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(pp)[,2:3]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.5965  -0.2351  -0.1466  -0.0449   1.0509  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.487915   0.085200   5.727 1.28e-08 ***
#dist_pond           -0.058593   0.012899  -4.542 6.09e-06 ***
#length_river         0.024313   0.008442   2.880  0.00405 ** 
#length_forest_edges  0.034261   0.005484   6.248 5.68e-10 ***
#length_hedges        0.001724   0.005653   0.305  0.76046    
#length_roads        -0.008674   0.006817  -1.272  0.20347    
#per_housing          0.039344   0.016633   2.365  0.01816 *  
#per_pastures         0.010212   0.006793   1.503  0.13299    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1566135)
#
#Null deviance: 215.07  on 1273  degrees of freedom
#Residual deviance: 198.27  on 1266  degrees of freedom
#AIC: 1263.5
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.59003  -0.22848  -0.15413  -0.05038   1.05884  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.564619   0.071950   7.847 8.97e-15 ***
#dist_pond           -0.068754   0.011758  -5.847 6.35e-09 ***
#length_river         0.026708   0.008307   3.215  0.00134 ** 
#length_forest_edges  0.034407   0.005466   6.295 4.22e-10 ***
#per_housing          0.030444   0.015880   1.917  0.05545 .  
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1568388)
#
#Null deviance: 215.07  on 1273  degrees of freedom
#Residual deviance: 199.03  on 1269  degrees of freedom
#AIC: 1262.3
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          9 1263.467
#reduced.model.present  6 1262.314

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


#####################
#Bufo spinosus####
#####################
bs=data.sp[data.sp$species=="Bufo_spinosus",] #select the data for the focal species
bs=bs[!duplicated(bs),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(bs)[,2:3]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.2823  -0.4207   0.1407   0.3383   0.7139  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.745967   0.053080  14.054  < 2e-16 ***
#dist_pond           -0.057554   0.008027  -7.170 9.37e-13 ***
#length_river         0.020939   0.005263   3.979 7.10e-05 ***
#length_forest_edges  0.041706   0.003539  11.786  < 2e-16 ***
#length_hedges        0.019627   0.003957   4.961 7.42e-07 ***
#length_roads         0.036195   0.003829   9.452  < 2e-16 ***
#per_housing          0.108655   0.008051  13.495  < 2e-16 ***
#per_pastures        -0.002179   0.004888  -0.446    0.656    
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1829041)
#
#Null deviance: 670.29  on 3032  degrees of freedom
#Residual deviance: 553.28  on 3025  degrees of freedom
#AIC: 3464.8
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.2841  -0.4185   0.1410   0.3385   0.7143  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.739824   0.051253  14.435  < 2e-16 ***
#dist_pond           -0.056833   0.007861  -7.229 6.11e-13 ***
#length_river         0.020788   0.005251   3.959 7.71e-05 ***
#length_forest_edges  0.041739   0.003537  11.799  < 2e-16 ***
#length_hedges        0.018925   0.003629   5.215 1.96e-07 ***
#length_roads         0.036378   0.003807   9.556  < 2e-16 ***
#per_housing          0.109226   0.007948  13.743  < 2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1828557)
#
#Null deviance: 670.29  on 3032  degrees of freedom
#Residual deviance: 553.32  on 3026  degrees of freedom
#AIC: 3463
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present           9 3464.831
#reduced.model.present  8 3463.030

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


######################
#Epidalea calamita####
######################
ec=data.sp[data.sp$species=="Epidalea_calamita",] #select the data for the focal species
ec=ec[!duplicated(ec),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(ec)[,2:3]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.26848  -0.10071  -0.06972  -0.03786   1.02496  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.226515   0.064823   3.494 0.000495 ***
#dist_pond           -0.027523   0.009814  -2.804 0.005132 ** 
#length_river        -0.019597   0.007052  -2.779 0.005545 ** 
#length_forest_edges  0.013459   0.004265   3.156 0.001645 ** 
#length_hedges       -0.001114   0.004165  -0.268 0.789067    
#length_roads         0.002737   0.004886   0.560 0.575471    
#per_housing          0.030876   0.012000   2.573 0.010216 *  
#per_pastures         0.003483   0.005006   0.696 0.486646    
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.07191946)
#
#Null deviance: 80.037  on 1086  degrees of freedom
#Residual deviance: 77.601  on 1079  degrees of freedom
#AIC: 233.53
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.27330  -0.09925  -0.07009  -0.03973   1.01750  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.241833   0.054848   4.409 1.14e-05 ***
#dist_pond           -0.029308   0.008936  -3.280  0.00107 ** 
#length_river        -0.019040   0.006946  -2.741  0.00623 ** 
#length_forest_edges  0.013512   0.004237   3.189  0.00147 ** 
#per_housing          0.031800   0.011502   2.765  0.00580 ** 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.07177032)
#
#Null deviance: 80.037  on 1086  degrees of freedom
#Residual deviance: 77.655  on 1082  degrees of freedom
#AIC: 228.29
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          9 233.5322
#reduced.model.present  6 228.2938

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


#####################
#Hyla arborea####
#####################
ha=data.sp[data.sp$species=="Hyla_arborea",] #select the data for the focal species
ha=ha[!duplicated(ha),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(ha)[,2:3]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.2087  -0.4199   0.1288   0.3745   1.0075  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          1.046783   0.060738  17.234  < 2e-16 ***
#dist_pond           -0.125255   0.009282 -13.495  < 2e-16 ***
#length_river         0.005976   0.006822   0.876 0.381118    
#length_forest_edges  0.040640   0.004357   9.328  < 2e-16 ***
#length_hedges        0.022580   0.004917   4.592 4.62e-06 ***
#length_roads         0.008541   0.005230   1.633 0.102590    
#per_housing          0.107913   0.011225   9.614  < 2e-16 ***
#per_pastures         0.020884   0.005559   3.757 0.000177 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1977555)
#
#Null deviance: 571.92  on 2335  degrees of freedom
#Residual deviance: 460.37  on 2328  degrees of freedom
#AIC: 2853.3
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.2249  -0.4197   0.1318   0.3772   1.0000  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          1.043018   0.060719  17.178  < 2e-16 ***
#dist_pond           -0.123910   0.009252 -13.392  < 2e-16 ***
#length_forest_edges  0.040910   0.004307   9.499  < 2e-16 ***
#length_hedges        0.024474   0.004811   5.087 3.93e-07 ***
#per_housing          0.113326   0.010675  10.616  < 2e-16 ***
#per_pastures         0.020497   0.005528   3.708 0.000214 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1978783)
#
#Null deviance: 571.92  on 2335  degrees of freedom
#Residual deviance: 461.06  on 2330  degrees of freedom
#AIC: 2852.7
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          9 2853.257
#reduced.model.present  7 2852.712

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


###################
#Rana dalmatina####
###################
rd=data.sp[data.sp$species=="Rana_dalmatina",] #select the data for the focal species
rd=rd[!duplicated(rd),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(rd)[,2:3]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.1986  -0.4028   0.1600   0.3403   0.7774  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.922457   0.051546  17.896  < 2e-16 ***
#dist_pond           -0.089616   0.007874 -11.381  < 2e-16 ***
#length_river         0.021283   0.005106   4.168 3.16e-05 ***
#length_forest_edges  0.056222   0.003523  15.959  < 2e-16 ***
#length_hedges        0.021546   0.004275   5.039 4.95e-07 ***
#length_roads        -0.005962   0.004860  -1.227     0.22    
#per_housing          0.056679   0.011317   5.008 5.82e-07 ***
#per_pastures         0.020256   0.004911   4.124 3.82e-05 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1871213)
#
#Null deviance: 656.59  on 2911  degrees of freedom
#Residual deviance: 543.40  on 2904  degrees of freedom
#AIC: 3393.4
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.1909  -0.4084   0.1601   0.3394   0.7812  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.922986   0.051549  17.905  < 2e-16 ***
#dist_pond           -0.090172   0.007862 -11.470  < 2e-16 ***
#length_river         0.021348   0.005106   4.181 2.99e-05 ***
#length_forest_edges  0.056311   0.003523  15.986  < 2e-16 ***
#length_hedges        0.020734   0.004224   4.908 9.69e-07 ***
#per_housing          0.052991   0.010911   4.857 1.26e-06 ***
#per_pastures         0.020752   0.004895   4.239 2.31e-05 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1871538)
#
#Null deviance: 656.59  on 2911  degrees of freedom
#Residual deviance: 543.68  on 2905  degrees of freedom
#AIC: 3392.9
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          9 3393.381
#reduced.model.present  8 3392.889

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model
#[1] "reduced.model.present"


####################
#Rana temporaria####
####################
rt=data.sp[data.sp$species=="Rana_temporaria",] #select the data for the focal species
rt=rt[!duplicated(rt),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(rt)[,2:3]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.9889  -0.2271  -0.1628   0.3621   0.9003  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.307409   0.086855   3.539 0.000412 ***
#dist_pond           -0.015891   0.013208  -1.203 0.229093    
#length_river         0.048527   0.006931   7.001  3.7e-12 ***
#length_forest_edges  0.083639   0.004660  17.950  < 2e-16 ***
#length_hedges       -0.014477   0.005541  -2.613 0.009068 ** 
#length_roads        -0.012532   0.006600  -1.899 0.057772 .  
#per_housing          0.029078   0.016799   1.731 0.083654 .  
#per_pastures         0.018669   0.006538   2.855 0.004354 ** 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1815888)
#
#Null deviance: 384.24  on 1623  degrees of freedom
#Residual deviance: 293.45  on 1616  degrees of freedom
#AIC: 1848.1
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.9708  -0.2191  -0.1742   0.3574   0.8758  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.204565   0.018808  10.877  < 2e-16 ***
#length_river         0.048764   0.006928   7.039 2.86e-12 ***
#length_forest_edges  0.083079   0.004660  17.829  < 2e-16 ***
#length_hedges       -0.014296   0.005306  -2.694 0.007127 ** 
#per_pastures         0.020849   0.006305   3.307 0.000964 ***
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1820755)
#
#Null deviance: 384.24  on 1623  degrees of freedom
#Residual deviance: 294.78  on 1619  degrees of freedom
#AIC: 1849.5
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          9 1848.132
#reduced.model.present  6 1849.490

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model
#[1] "model.present"


##########################
#Pelophylax esculentus####
##########################
pe=data.sp[data.sp$species=="Pelophylax_esculentus",] #select the data for the focal species
pe=pe[!duplicated(pe),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(pe)[,2:3]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.1259  -0.3911  -0.1005   0.4136   0.9970  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          1.013342   0.069815  14.515  < 2e-16 ***
#dist_pond           -0.133646   0.010537 -12.683  < 2e-16 ***
#length_river         0.008566   0.007747   1.106    0.269    
#length_forest_edges  0.043698   0.005080   8.602  < 2e-16 ***
#length_hedges        0.026552   0.005530   4.801 1.71e-06 ***
#length_roads        -0.001926   0.006153  -0.313    0.754    
#per_housing          0.103261   0.013700   7.537 7.58e-14 ***
#per_pastures         0.005302   0.006354   0.835    0.404    
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1986126)
#
#Null deviance: 444.14  on 1798  degrees of freedom
#Residual deviance: 355.72  on 1791  degrees of freedom
#AIC: 2207.4
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.1337  -0.3901  -0.1025   0.4145   0.9985  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          1.035555   0.066585  15.552  < 2e-16 ***
#dist_pond           -0.136779   0.010157 -13.467  < 2e-16 ***
#length_forest_edges  0.045161   0.004934   9.154  < 2e-16 ***
#length_hedges        0.028935   0.004928   5.872 5.12e-09 ***
#per_housing          0.099338   0.012887   7.708 2.10e-14 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.198513)
#
#Null deviance: 444.14  on 1798  degrees of freedom
#Residual deviance: 356.13  on 1794  degrees of freedom
#AIC: 2203.5
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          9 2207.421
#reduced.model.present  6 2203.529

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model
#[1] "reduced.model.present"


##########################
#Pelophylax ridibundus####
##########################
pr=data.sp[data.sp$species=="Pelophylax_ridibundus",] #select the data for the focal species
pr=pr[!duplicated(pr),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(pr)[,2:3]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.17920  -0.39832   0.02506   0.37418   0.94962  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.9109251  0.0607151  15.003  < 2e-16 ***
#dist_pond           -0.1190560  0.0091554 -13.004  < 2e-16 ***
#length_river         0.0360347  0.0061364   5.872 5.00e-09 ***
#length_forest_edges  0.0502532  0.0044076  11.402  < 2e-16 ***
#length_hedges        0.0368508  0.0051433   7.165 1.08e-12 ***
#length_roads        -0.0004309  0.0055474  -0.078    0.938    
#per_housing          0.1093033  0.0126070   8.670  < 2e-16 ***
#per_pastures         0.0083107  0.0057306   1.450    0.147    
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1904065)
#
#Null deviance: 518.77  on 2077  degrees of freedom
#Residual deviance: 394.14  on 2070  degrees of freedom
#AIC: 2460.5
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.17460  -0.39377   0.02284   0.37622   0.96395  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.933757   0.058657  15.919  < 2e-16 ***
#dist_pond           -0.121905   0.008946 -13.626  < 2e-16 ***
#length_river         0.036171   0.006135   5.896 4.34e-09 ***
#length_forest_edges  0.050530   0.004399  11.487  < 2e-16 ***
#length_hedges        0.039817   0.004623   8.612  < 2e-16 ***
#per_housing          0.105872   0.011872   8.918  < 2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1904207)
#
#Null deviance: 518.77  on 2077  degrees of freedom
#Residual deviance: 394.55  on 2072  degrees of freedom
#AIC: 2458.7
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          9 2460.535
#reduced.model.present  7 2458.697

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model
#[1] "reduced.model.present"


#########################
#Anguis fragilis####
#########################
af=data.sp[data.sp$species=="Anguis_fragilis",] #select the data for the focal species
af=af[!duplicated(af),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(af)[,2:3]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-0.46532  -0.08333  -0.03522  -0.00226   1.00745  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.2411561  0.0603177   3.998 6.82e-05 ***
#dist_pond           -0.0344597  0.0091131  -3.781 0.000165 ***
#length_river        -0.0101547  0.0064155  -1.583 0.113755    
#length_forest_edges  0.0201919  0.0039139   5.159 2.95e-07 ***
#length_hedges       -0.0050962  0.0038733  -1.316 0.188555    
#length_roads         0.0092194  0.0044166   2.087 0.037081 *  
#per_housing          0.0994277  0.0102757   9.676  < 2e-16 ***
#per_pastures         0.0008659  0.0046686   0.185 0.852892    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.06105928)
#
#Null deviance: 75.786  on 1081  degrees of freedom
#Residual deviance: 65.578  on 1074  degrees of freedom
#AIC: 55.379
#
#Number of Fisher Scoring iterations: 2

sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-0.46692  -0.08104  -0.03539  -0.00548   0.99431  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.196150   0.051197   3.831 0.000135 ***
#dist_pond           -0.029422   0.008323  -3.535 0.000425 ***
#length_forest_edges  0.019480   0.003844   5.068 4.74e-07 ***
#length_roads         0.007785   0.004323   1.801 0.071987 .  
#per_housing          0.100324   0.010235   9.802  < 2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.0611804)
#
#Null deviance: 75.786  on 1081  degrees of freedom
#Residual deviance: 65.891  on 1077  degrees of freedom
#AIC: 54.541
#
#Number of Fisher Scoring iterations: 2

aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         9 55.37851
#reduced.model.present  6 54.54090

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


#########################
#Podarcis muralis####
#########################
pm=data.sp[data.sp$species=="Podarcis_muralis",] #select the data for the focal species
pm=pm[!duplicated(pm),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(pm)[,2:3]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.2701  -0.3934   0.1299   0.3391   0.6722  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.5047416  0.0601365   8.393  < 2e-16 ***
#dist_pond           -0.0223926  0.0090342  -2.479  0.01324 *  
#length_river         0.0165121  0.0052812   3.127  0.00179 ** 
#length_forest_edges  0.0511732  0.0034757  14.723  < 2e-16 ***
#length_hedges        0.0261584  0.0039966   6.545 6.97e-11 ***
#length_roads         0.0129518  0.0040731   3.180  0.00149 ** 
#per_housing          0.1418510  0.0078182  18.144  < 2e-16 ***
#per_pastures        -0.0004443  0.0049272  -0.090  0.92816    
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1803317)
#
#Null deviance: 666.44  on 2997  degrees of freedom
#Residual deviance: 539.19  on 2990  degrees of freedom
#AIC: 3382.5
#
#Number of Fisher Scoring iterations: 2

sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.2705  -0.3931   0.1299   0.3394   0.6722  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.503399   0.058255   8.641  < 2e-16 ***
#dist_pond           -0.022230   0.008851  -2.512  0.01207 *  
#length_river         0.016504   0.005280   3.126  0.00179 ** 
#length_forest_edges  0.051179   0.003474  14.730  < 2e-16 ***
#length_hedges        0.026012   0.003652   7.122 1.32e-12 ***
#length_roads         0.012979   0.004061   3.196  0.00141 ** 
#per_housing          0.141959   0.007724  18.378  < 2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1802719)
#
#Null deviance: 666.44  on 2997  degrees of freedom
#Residual deviance: 539.19  on 2991  degrees of freedom
#AIC: 3380.5
#
#Number of Fisher Scoring iterations: 2

aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         9 3382.499
#reduced.model.present  8 3380.507

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


#########################
#Lacerta bilineata####
#########################
lb=data.sp[data.sp$species=="Lacerta_bilineata",] #select the data for the focal species
lb=lb[!duplicated(lb),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(lb)[,2:3]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.0655  -0.4359   0.1307   0.4211   0.8221  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.520182   0.077098   6.747 1.94e-11 ***
#dist_pond           -0.043325   0.011731  -3.693 0.000227 ***
#length_river        -0.014218   0.007779  -1.828 0.067723 .  
#length_forest_edges  0.068937   0.004424  15.581  < 2e-16 ***
#length_hedges        0.027407   0.005318   5.154 2.79e-07 ***
#length_roads         0.015603   0.005494   2.840 0.004553 ** 
#per_housing          0.079100   0.013180   6.002 2.29e-09 ***
#per_pastures         0.017341   0.006125   2.831 0.004680 ** 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.2122716)
#
#Null deviance: 534.88  on 2149  degrees of freedom
#Residual deviance: 454.69  on 2142  degrees of freedom
#AIC: 2779.2
#
#Number of Fisher Scoring iterations: 2

sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.0783  -0.4374   0.1386   0.4218   0.8211  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.517122   0.077121   6.705 2.56e-11 ***
#dist_pond           -0.042817   0.011734  -3.649  0.00027 ***
#length_forest_edges  0.067096   0.004311  15.566  < 2e-16 ***
#length_hedges        0.026481   0.005297   5.000 6.21e-07 ***
#length_roads         0.016015   0.005492   2.916  0.00358 ** 
#per_housing          0.079201   0.013187   6.006 2.23e-09 ***
#per_pastures         0.016937   0.006124   2.766  0.00573 ** 
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.2125035)
#
#Null deviance: 534.88  on 2149  degrees of freedom
#Residual deviance: 455.39  on 2143  degrees of freedom
#AIC: 2780.5
#
#Number of Fisher Scoring iterations: 2

aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         9 2779.160
#reduced.model.present  8 2780.511

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "model.present"


#########################
#Natrix helvetica####
#########################
nh=data.sp[data.sp$species=="Natrix_helvetica",] #select the data for the focal species
nh=nh[!duplicated(nh),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(nh)[,2:3]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.2074  -0.3879  -0.1032   0.4159   0.9287  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.547306   0.076328   7.170 1.06e-12 ***
#dist_pond           -0.060256   0.011563  -5.211 2.08e-07 ***
#length_river         0.027707   0.006886   4.023 5.96e-05 ***
#length_forest_edges  0.051278   0.004815  10.651  < 2e-16 ***
#length_hedges        0.031784   0.005505   5.774 9.02e-09 ***
#length_roads         0.035643   0.005320   6.699 2.74e-11 ***
#per_housing          0.100181   0.012291   8.151 6.44e-16 ***
#per_pastures         0.003784   0.006358   0.595    0.552    
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.2024449)
#
#Null deviance: 484.80  on 1940  degrees of freedom
#Residual deviance: 391.33  on 1933  degrees of freedom
#AIC: 2418
#
#Number of Fisher Scoring iterations: 2

sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.2038  -0.3862  -0.1033   0.4152   0.9287  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.558311   0.074043   7.540 7.16e-14 ***
#dist_pond           -0.061647   0.011322  -5.445 5.84e-08 ***
#length_river         0.027709   0.006885   4.024 5.93e-05 ***
#length_forest_edges  0.051273   0.004814  10.651  < 2e-16 ***
#length_hedges        0.033255   0.004919   6.761 1.81e-11 ***
#length_roads         0.035362   0.005298   6.674 3.24e-11 ***
#per_housing          0.098982   0.012123   8.165 5.74e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.2023773)
#
#Null deviance: 484.8  on 1940  degrees of freedom
#Residual deviance: 391.4  on 1934  degrees of freedom
#AIC: 2416.3
#
#Number of Fisher Scoring iterations: 2

aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         9 2417.968
#reduced.model.present  8 2416.324

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


##################
#Natrix maura####
##################
nm=data.sp[data.sp$species=="Natrix_maura",] #select the data for the focal species
nm=nm[!duplicated(nm),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(nm)[,2:3]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-0.96363  -0.15241  -0.02214   0.02218   1.04520  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.195648   0.069110   2.831  0.00472 ** 
#dist_pond           -0.029514   0.010506  -2.809  0.00504 ** 
#length_river         0.077442   0.005900  13.125  < 2e-16 ***
#length_forest_edges  0.046414   0.004339  10.697  < 2e-16 ***
#length_hedges        0.005901   0.004582   1.288  0.19805    
#length_roads         0.014197   0.005027   2.824  0.00482 ** 
#per_housing          0.108803   0.011584   9.392  < 2e-16 ***
#per_pastures        -0.015448   0.005489  -2.814  0.00497 ** 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.09687071)
#
#Null deviance: 185.00  on 1226  degrees of freedom
#Residual deviance: 118.09  on 1219  degrees of freedom
#AIC: 627.77
#
#Number of Fisher Scoring iterations: 2

sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-0.96159  -0.15318  -0.02483   0.02077   1.04614  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.218499   0.066812   3.270  0.00110 ** 
#dist_pond           -0.031522   0.010392  -3.033  0.00247 ** 
#length_river         0.078802   0.005806  13.571  < 2e-16 ***
#length_forest_edges  0.045942   0.004325  10.623  < 2e-16 ***
#length_roads         0.015308   0.004953   3.090  0.00204 ** 
#per_housing          0.108789   0.011588   9.388  < 2e-16 ***
#per_pastures        -0.012416   0.004960  -2.503  0.01244 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.09692299)
#
#Null deviance: 185.00  on 1226  degrees of freedom
#Residual deviance: 118.25  on 1220  degrees of freedom
#AIC: 627.44
#
#Number of Fisher Scoring iterations: 2

aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         9 627.7670
#reduced.model.present  8 627.4353

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


###########################
#Hierophis viridiflavus####
###########################
hv=data.sp[data.sp$species=="Hierophis_viridiflavus",] #select the data for the focal species
hv=hv[!duplicated(hv),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(hv)[,2:3]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.2295  -0.4354   0.1205   0.3330   0.7103  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.269372   0.064134   4.200 2.75e-05 ***
#dist_pond            0.009335   0.009580   0.974   0.3299    
#length_river         0.010744   0.005817   1.847   0.0648 .  
#length_forest_edges  0.040167   0.003804  10.560  < 2e-16 ***
#length_hedges        0.032080   0.004256   7.538 6.39e-14 ***
#length_roads         0.052483   0.003909  13.425  < 2e-16 ***
#per_housing          0.096619   0.008735  11.061  < 2e-16 ***
#per_pastures        -0.001233   0.005115  -0.241   0.8095    
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1860754)
#
#Null deviance: 645.26  on 2818  degrees of freedom
#Residual deviance: 523.06  on 2811  degrees of freedom
#AIC: 3269.5
#
#Number of Fisher Scoring iterations: 2

sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.2300  -0.4347   0.1203   0.3327   0.6718  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)         0.328181   0.016961  19.349   <2e-16 ***
#length_forest_edges 0.042010   0.003673  11.438   <2e-16 ***
#length_hedges       0.031530   0.003704   8.513   <2e-16 ***
#length_roads        0.052724   0.003867  13.636   <2e-16 ***
#per_housing         0.096157   0.008637  11.133   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1861646)
#
#Null deviance: 645.26  on 2818  degrees of freedom
#Residual deviance: 523.87  on 2814  degrees of freedom
#AIC: 3267.9
#
#Number of Fisher Scoring iterations: 2

aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         9 3269.525
#reduced.model.present  6 3267.882

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


########################
#Zamenis longissimus####
########################
zl=data.sp[data.sp$species=="Zamenis_longissimus",] #select the data for the focal species
zl=zl[!duplicated(zl),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(zl)[,2:3]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.0998  -0.3104  -0.1630   0.4146   0.9970  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)         -0.1597702  0.0864638  -1.848 0.064805 .  
#dist_pond            0.0486245  0.0129857   3.744 0.000187 ***
#length_river        -0.0009719  0.0085665  -0.113 0.909682    
#length_forest_edges  0.0488443  0.0050986   9.580  < 2e-16 ***
#length_hedges        0.0285082  0.0055661   5.122 3.38e-07 ***
#length_roads         0.0471397  0.0056010   8.416  < 2e-16 ***
#per_housing          0.1251322  0.0123520  10.131  < 2e-16 ***
#per_pastures        -0.0054286  0.0068150  -0.797 0.425812    
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1922397)
#
#Null deviance: 400.84  on 1668  degrees of freedom
#Residual deviance: 319.31  on 1661  degrees of freedom
#AIC: 1994.2
#
#Number of Fisher Scoring iterations: 2

sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.0903  -0.3084  -0.1610   0.4166   0.9821  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)         -0.180479   0.082529  -2.187   0.0289 *  
#dist_pond            0.051378   0.012527   4.101 4.31e-05 ***
#length_forest_edges  0.048630   0.004964   9.797  < 2e-16 ***
#length_hedges        0.026515   0.004984   5.320 1.18e-07 ***
#length_roads         0.047798   0.005539   8.629  < 2e-16 ***
#per_housing          0.126304   0.012262  10.300  < 2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1920852)
#
#Null deviance: 400.84  on 1668  degrees of freedom
#Residual deviance: 319.44  on 1663  degrees of freedom
#AIC: 1990.9
#
#Number of Fisher Scoring iterations: 2

aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         9 1994.197
#reduced.model.present  7 1990.862

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


########################
#Vipera aspis####
########################
va=data.sp[data.sp$species=="Vipera_aspis",] #select the data for the focal species
va=va[!duplicated(va),] #remove duplicated occurrences to avoid pseudoreplication

presvals=raster::extract(rasters.selected,as.data.frame(va)[,2:3]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-0.67131  -0.22868  -0.12079   0.01602   1.07821  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.536009   0.079838   6.714 2.89e-11 ***
#dist_pond           -0.077822   0.012094  -6.435 1.76e-10 ***
#length_river        -0.016427   0.008846  -1.857   0.0635 .  
#length_forest_edges  0.034851   0.005276   6.606 5.87e-11 ***
#length_hedges        0.006652   0.005534   1.202   0.2296    
#length_roads         0.024629   0.005932   4.152 3.52e-05 ***
#per_housing          0.080954   0.014377   5.631 2.22e-08 ***
#per_pastures         0.010652   0.006513   1.636   0.1022    
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1368318)
#
#Null deviance: 194.85  on 1241  degrees of freedom
#Residual deviance: 168.85  on 1234  degrees of freedom
#AIC: 1064.3
#
#Number of Fisher Scoring iterations: 2

sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-0.68185  -0.22357  -0.12098   0.00504   1.06817  
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.638980   0.067165   9.514  < 2e-16 ***
#dist_pond           -0.089597   0.011058  -8.102 1.28e-15 ***
#length_forest_edges  0.033041   0.005189   6.368 2.70e-10 ***
#length_roads         0.026240   0.005800   4.524 6.65e-06 ***
#per_housing          0.078424   0.014338   5.470 5.46e-08 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1375474)
#
#Null deviance: 194.85  on 1241  degrees of freedom
#Residual deviance: 170.15  on 1237  degrees of freedom
#AIC: 1067.8
#
#Number of Fisher Scoring iterations: 2

aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         9 1064.276
#reduced.model.present  6 1067.770

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "model.present"


###################################
#WITH NON-TRANSFORMED VARIABLES####
###################################

############################
#Open species data table####
############################
data.sp=read.table("DATA_SPECIES_79/DONNEES_AMPHREP_79.txt", h=T)
str(data.sp)
coordinates(data.sp)=~lon+lat
proj4string(data.sp)=CRS("+init=epsg:4326")
data.sp=spTransform(data.sp, CRS("+init=epsg:2154"))

#List of species
#sort(unique(data.sp@data$species))
#[1] "Alytes_obstetricans"    "Anguis_fragilis"        "Bombina_variegata"      "Bufo_spinosus"          "Coronella_austriaca"   
# [6] "Emys_orbicularis"       "Epidalea_calamita"      "Hierophis_viridiflavus" "Hyla_arborea"           "Hyla_meridionalis"     
#[11] "Ichthyosaura_alpestris" "Lacerta_bilineata"      "Lissotriton_helveticus" "Lissotriton_vulgaris"   "Natrix_helvetica"      
#[16] "Natrix_maura"           "Pelodytes_punctatus"    "Pelophylax_esculentus"  "Pelophylax_lessonae"    "Pelophylax_ridibundus" 
#[21] "Podarcis_muralis"       "Rana_dalmatina"         "Rana_temporaria"        "Salamandra_salamandra"  "Trachemys_scripta"     
#[26] "Triturus_blasii"        "Triturus_cristatus"     "Triturus_marmoratus"    "Vipera_aspis"           "Xenopus_laevis"        
#[31] "Zamenis_longissimus" 


#Read stacked raster files####
#create an object that gather all the raster layers of the environmental variables
#attention pour générer la liste de rasters, il ne faut pas prendre les fichiers .xml qui sont générés au moment de l'export des couches raster en ascii. Donc prévoir un nouveau dossier dans lequel on range les seulement les fichiers .asc et .prj
list.rasters=(list.files("env_variables", full.names=T,pattern=".asc")) #provide the path of the directory where the raster layers saved
rasters=stack(list.rasters) #stack the layers
rasters #check what it looks like
projection(rasters)=CRS("+init=epsg:2154") #assign a projection system, here 2154 for Lambert 93
#pdf(file="plot_raster_variables.pdf") #save as a pdf file
jpeg(file="plot_raster_variables.jpg", height=17, width=17, units="cm", res=200) #save as a jpg file
plot(rasters) #plot all the rasters on one panel
dev.off() #save the file in the working directory

#Check collinearity between variables####
#with Pearson's####
#library("virtualspecies")
rasters.reduced=removeCollinearity(rasters, multicollinearity.cutoff=0.7, plot=T, select.variables=T, sample.points=F) #use a Pearson's test 
# - No multicollinearity detected in your data at threshold 0.7
rasters.reduced #check the results
rasters.selected=subset(rasters, c("dist_forest","dist_pond","dist_river","length_forest_edges","length_hedges","length_roads","per_housing","per_pastures")) #create an object that contains all the variables that were selected

#with VIF####
vif(env.var[,2:11])
#            Variables      VIF
#1         per.pasture 2.249548
#2            per.agri 2.627156
#3          per.forest 2.356400
#4         per.housing 1.321686
#5         length.road 1.158513
#6  length.forest.edge 1.612014
#7        length.hedge 1.342562
#8          dist.river 1.314031
#9         dist.forest 1.394635
#10          dist.pond 1.334876

#VIF values from 5 to 10 are considered critical.

vifcor(env.var[,2:11], th=.7) #VIF with a correlation threshold at r = 0.7
#No variable from the 10 input variables has collinearity problem. 

#The linear correlation coefficients ranges between: 
#min correlation ( dist.forest ~ length.road ):  0.0003669296 
#max correlation ( length.forest.edge ~ per.forest ):  0.5725145 
#
#---------- VIFs of the remained variables -------- 
#            Variables      VIF
#1         per.pasture 2.214779
#2            per.agri 2.619966
#3          per.forest 2.295791
#4         per.housing 1.338639
#5         length.road 1.143286
#6  length.forest.edge 1.586353
#7        length.hedge 1.335140
#8          dist.river 1.325311
#9         dist.forest 1.390679
#10          dist.pond 1.349130

#All variables are kept at a threshold of .7


##########################
#Salamandra salamandra####
##########################
ss=data.sp[data.sp@data$species=="Salamandra_salamandra",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(ss)[,7:8]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.0607  -0.2918   0.1213   0.2777   1.1498  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          7.172e-01  1.837e-02  39.033  < 2e-16 ***
#dist_forest         -1.109e-03  4.642e-05 -23.890  < 2e-16 ***
#dist_pond           -8.008e-05  1.709e-05  -4.685 2.92e-06 ***
#dist_river           8.767e-06  8.896e-06   0.986  0.32442    
#length_forest_edges  8.696e-04  7.900e-05  11.007  < 2e-16 ***
#length_hedges        2.262e-04  8.793e-05   2.572  0.01014 *  
#length_roads         5.899e-04  1.567e-04   3.764  0.00017 ***
#per_housing          8.524e-03  1.641e-03   5.195 2.17e-07 ***
#per_pastures         1.732e-04  2.371e-04   0.731  0.46494    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1499779)
#
#    Null deviance: 688.47  on 3209  degrees of freedom
#Residual deviance: 480.08  on 3201  degrees of freedom
#AIC: 3030.3
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.0620  -0.2931   0.1233   0.2782   1.1400  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          7.242e-01  1.711e-02  42.335  < 2e-16 ***
#dist_forest         -1.106e-03  4.631e-05 -23.875  < 2e-16 ***
#dist_pond           -7.547e-05  1.554e-05  -4.856 1.26e-06 ***
#length_forest_edges  8.661e-04  7.892e-05  10.975  < 2e-16 ***
#length_hedges        2.369e-04  8.124e-05   2.916 0.003572 ** 
#length_roads         5.780e-04  1.539e-04   3.755 0.000176 ***
#per_housing          8.316e-03  1.624e-03   5.119 3.25e-07 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1499478)
#
#    Null deviance: 688.47  on 3209  degrees of freedom
#Residual deviance: 480.28  on 3203  degrees of freedom
#AIC: 3027.7
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          10 3030.346
#reduced.model.present  8 3027.705

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


############################
#Lissotriton helveticus####
###########################
lh=data.sp[data.sp@data$species=="Lissotriton_helveticus",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(lh)[,7:8]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.11749  -0.00352   0.11957   0.21509   0.95259  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          7.401e-01  1.430e-02  51.751  < 2e-16 ***
#dist_forest         -2.916e-04  3.120e-05  -9.345  < 2e-16 ***
#dist_pond           -1.784e-04  1.271e-05 -14.034  < 2e-16 ***
#dist_river           3.443e-05  7.345e-06   4.687 2.84e-06 ***
#length_forest_edges  6.745e-04  6.936e-05   9.725  < 2e-16 ***
#length_hedges        6.385e-04  5.518e-05  11.572  < 2e-16 ***
#length_roads        -2.438e-05  1.274e-04  -0.191    0.848    
#per_housing          1.038e-02  1.109e-03   9.359  < 2e-16 ***
#per_pastures         1.011e-03  1.499e-04   6.742 1.73e-11 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1320464)
#
#    Null deviance: 809.16  on 5239  degrees of freedom
#Residual deviance: 690.73  on 5231  degrees of freedom
#AIC: 4272.6
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variable within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.11672  -0.00511   0.11962   0.21461   0.95338  
#
#Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          7.400e-01  1.429e-02  51.796  < 2e-16 ***
#dist_forest         -2.918e-04  3.117e-05  -9.362  < 2e-16 ***
#dist_pond           -1.785e-04  1.270e-05 -14.058  < 2e-16 ***
#dist_river           3.440e-05  7.343e-06   4.684 2.88e-06 ***
#length_forest_edges  6.737e-04  6.922e-05   9.733  < 2e-16 ***
#length_hedges        6.366e-04  5.424e-05  11.737  < 2e-16 ***
#per_housing          1.031e-02  1.045e-03   9.870  < 2e-16 ***
#per_pastures         1.014e-03  1.487e-04   6.823 9.90e-12 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1320221)
#
#    Null deviance: 809.16  on 5239  degrees of freedom
#Residual deviance: 690.74  on 5232  degrees of freedom
#AIC: 4270.6
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         10 4272.555
#reduced.model.present  9 4270.592

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


#########################
#Triturus cristatus####
#########################
tc=data.sp[data.sp@data$species=="Triturus_cristatus",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(tc)[,7:8]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.0030  -0.3945  -0.1321   0.4140   1.0703  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          2.495e-01  2.888e-02   8.640  < 2e-16 ***
#dist_forest         -1.257e-04  5.952e-05  -2.112   0.0348 *  
#dist_pond           -1.375e-04  2.110e-05  -6.513 9.55e-11 ***
#dist_river           6.858e-05  1.450e-05   4.730 2.43e-06 ***
#length_forest_edges  1.282e-03  1.590e-04   8.062 1.37e-15 ***
#length_hedges        1.073e-03  1.195e-04   8.984  < 2e-16 ***
#length_roads        -2.981e-04  2.826e-04  -1.055   0.2916    
#per_housing          1.887e-02  2.489e-03   7.581 5.49e-14 ***
#per_pastures         1.728e-03  3.152e-04   5.480 4.85e-08 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.2032213)
#
#    Null deviance: 440.40  on 1786  degrees of freedom
#Residual deviance: 361.33  on 1778  degrees of freedom
#AIC: 2234.8
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.9979  -0.3925  -0.1334   0.4246   1.0751  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          2.473e-01  2.880e-02   8.586  < 2e-16 ***
#dist_forest         -1.268e-04  5.951e-05  -2.131   0.0333 *  
#dist_pond           -1.376e-04  2.110e-05  -6.522 9.03e-11 ***
#dist_river           6.774e-05  1.448e-05   4.678 3.11e-06 ***
#length_forest_edges  1.294e-03  1.587e-04   8.154 6.56e-16 ***
#length_hedges        1.046e-03  1.167e-04   8.966  < 2e-16 ***
#per_housing          1.832e-02  2.434e-03   7.526 8.23e-14 ***
#per_pastures         1.752e-03  3.144e-04   5.572 2.90e-08 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.2032342)
#
#    Null deviance: 440.40  on 1786  degrees of freedom
#Residual deviance: 361.55  on 1779  degrees of freedom
#AIC: 2233.9
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         10 2234.751
#reduced.model.present  9 2233.869

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


#########################
#Triturus marmoratus####
#########################
tm=data.sp[data.sp@data$species=="Triturus_marmoratus",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(tm)[,7:8]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.1157  -0.3684   0.1585   0.2853   0.9458  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          6.257e-01  2.013e-02  31.076  < 2e-16 ***
#dist_forest         -2.569e-04  4.331e-05  -5.931 3.31e-09 ***
#dist_pond           -2.721e-04  1.793e-05 -15.171  < 2e-16 ***
#dist_river           4.025e-05  1.124e-05   3.582 0.000346 ***
#length_forest_edges  8.468e-04  1.086e-04   7.795 8.57e-15 ***
#length_hedges        8.733e-04  7.915e-05  11.034  < 2e-16 ***
#length_roads        -7.010e-05  1.800e-04  -0.389 0.697046    
#per_housing          1.612e-02  1.517e-03  10.631  < 2e-16 ***
#per_pastures         1.018e-03  2.122e-04   4.798 1.67e-06 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1681287)
#
#    Null deviance: 700.06  on 3333  degrees of freedom
#Residual deviance: 559.03  on 3325  degrees of freedom
#AIC: 3527.9
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.1169  -0.3685   0.1589   0.2865   0.9469  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          6.252e-01  2.010e-02  31.113  < 2e-16 ***
#dist_forest         -2.570e-04  4.330e-05  -5.937 3.21e-09 ***
#dist_pond           -2.722e-04  1.793e-05 -15.185  < 2e-16 ***
#dist_river           4.013e-05  1.123e-05   3.573 0.000358 ***
#length_forest_edges  8.452e-04  1.085e-04   7.786 9.13e-15 ***
#length_hedges        8.678e-04  7.783e-05  11.149  < 2e-16 ***
#per_housing          1.591e-02  1.414e-03  11.247  < 2e-16 ***
#per_pastures         1.027e-03  2.111e-04   4.863 1.21e-06 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1680858)
#
#    Null deviance: 700.06  on 3333  degrees of freedom
#Residual deviance: 559.05  on 3326  degrees of freedom
#AIC: 3526
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          10 3527.863
#reduced.model.present  9 3526.015

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


#####################
#Triturus blasii####
#####################
tb=data.sp[data.sp@data$species=="Triturus_blasii",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(tb)[,7:8]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.4133  -0.1738  -0.1052  -0.0341   0.9751  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          7.521e-02  2.512e-02   2.994 0.002815 ** 
#dist_forest         -8.240e-05  5.087e-05  -1.620 0.105513    
#dist_pond           -6.007e-05  1.856e-05  -3.237 0.001243 ** 
#dist_river           4.678e-05  1.270e-05   3.683 0.000241 ***
#length_forest_edges  5.246e-04  1.718e-04   3.053 0.002316 ** 
#length_hedges        6.026e-04  1.147e-04   5.253 1.79e-07 ***
#length_roads        -4.253e-04  2.661e-04  -1.598 0.110223    
#per_housing          7.989e-03  2.722e-03   2.935 0.003406 ** 
#per_pastures         4.128e-04  2.945e-04   1.402 0.161270    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1099154)
#
#    Null deviance: 135.70  on 1156  degrees of freedom
#Residual deviance: 126.18  on 1148  degrees of freedom
#AIC: 739.68
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.41805  -0.16681  -0.09995  -0.04577   0.97364  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          6.612e-02  1.974e-02   3.350 0.000834 ***
#dist_pond           -6.855e-05  1.823e-05  -3.761 0.000178 ***
#dist_river           4.105e-05  1.253e-05   3.275 0.001087 ** 
#length_forest_edges  6.584e-04  1.532e-04   4.297 1.88e-05 ***
#length_hedges        6.088e-04  1.078e-04   5.647 2.05e-08 ***
#per_housing          6.918e-03  2.662e-03   2.599 0.009474 ** 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1103714)
#
#    Null deviance: 135.70  on 1156  degrees of freedom
#Residual deviance: 127.04  on 1151  degrees of freedom
#AIC: 741.49
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         10 739.6811
#reduced.model.present  7 741.4907

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "model.present"


########################
#Alytes obstetricans####
########################
ao=data.sp[data.sp@data$species=="Alytes_obstetricans",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(ao)[,7:8]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.0584  -0.2757  -0.1180   0.3129   1.0357  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          2.902e-01  2.412e-02  12.033  < 2e-16 ***
#dist_forest         -2.926e-04  5.173e-05  -5.657 1.80e-08 ***
#dist_pond           -1.420e-05  1.704e-05  -0.834   0.4046    
#dist_river           2.463e-05  1.052e-05   2.341   0.0194 *  
#length_forest_edges  9.553e-04  1.218e-04   7.846 7.42e-15 ***
#length_hedges        2.547e-04  1.182e-04   2.155   0.0313 *  
#length_roads         9.899e-04  1.936e-04   5.114 3.51e-07 ***
#per_housing          2.392e-02  1.182e-03  20.239  < 2e-16 ***
#per_pastures        -2.380e-03  3.394e-04  -7.014 3.31e-12 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1480019)
#
#    Null deviance: 431.82  on 1759  degrees of freedom
#Residual deviance: 259.15  on 1751  degrees of freedom
#AIC: 1643.1
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.0495  -0.2776  -0.1193   0.3131   1.0380  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          2.837e-01  2.282e-02  12.434  < 2e-16 ***
#dist_forest         -2.986e-04  5.122e-05  -5.830 6.58e-09 ***
#dist_river           2.191e-05  1.000e-05   2.191   0.0286 *  
#length_forest_edges  9.418e-04  1.207e-04   7.805 1.01e-14 ***
#length_hedges        2.638e-04  1.176e-04   2.243   0.0250 *  
#length_roads         9.905e-04  1.936e-04   5.117 3.44e-07 ***
#per_housing          2.382e-02  1.176e-03  20.259  < 2e-16 ***
#per_pastures        -2.341e-03  3.360e-04  -6.966 4.60e-12 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1479761)
#
#    Null deviance: 431.82  on 1759  degrees of freedom
#Residual deviance: 259.25  on 1752  degrees of freedom
#AIC: 1641.8
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          10 1643.108
#reduced.model.present  9 1641.806

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


########################
#Pelodytes punctatus####
########################
pp=data.sp[data.sp@data$species=="Pelodytes_punctatus",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(pp)[,7:8]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.9798  -0.3428  -0.1855   0.4799   1.2893  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          2.518e-01  2.953e-02   8.527  < 2e-16 ***
#dist_forest         -1.708e-04  6.400e-05  -2.670  0.00767 ** 
#dist_pond           -3.249e-05  2.315e-05  -1.403  0.16070    
#dist_river          -4.721e-07  1.596e-05  -0.030  0.97641    
#length_forest_edges  1.458e-03  1.615e-04   9.026  < 2e-16 ***
#length_hedges        7.423e-04  1.256e-04   5.909 4.22e-09 ***
#length_roads        -1.334e-03  3.287e-04  -4.057 5.21e-05 ***
#per_housing          8.644e-03  3.083e-03   2.804  0.00511 ** 
#per_pastures         1.326e-03  3.322e-04   3.992 6.84e-05 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1983029)
#
#    Null deviance: 358.56  on 1558  degrees of freedom
#Residual deviance: 307.37  on 1550  degrees of freedom
#AIC: 1912.8
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.9903  -0.3427  -0.1872   0.4842   1.0450  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          2.294e-01  2.561e-02   8.955  < 2e-16 ***
#dist_forest         -1.854e-04  6.265e-05  -2.960  0.00313 ** 
#length_forest_edges  1.458e-03  1.615e-04   9.028  < 2e-16 ***
#length_hedges        7.794e-04  1.233e-04   6.322 3.37e-10 ***
#length_roads        -1.376e-03  3.270e-04  -4.207 2.73e-05 ***
#per_housing          9.166e-03  3.066e-03   2.990  0.00283 ** 
#per_pastures         1.455e-03  3.207e-04   4.537 6.15e-06 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1983675)
#
#    Null deviance: 358.56  on 1558  degrees of freedom
#Residual deviance: 307.87  on 1552  degrees of freedom
#AIC: 1911.3
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          10 1912.825
#reduced.model.present  8 1911.343

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


#####################
#Bufo spinosus####
#####################
bs=data.sp[data.sp@data$species=="Bufo_spinosus",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(bs)[,7:8]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.21418  -0.07188   0.13291   0.25514   0.86479  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          7.424e-01  1.635e-02  45.418  < 2e-16 ***
#dist_forest         -4.306e-04  3.880e-05 -11.099  < 2e-16 ***
#dist_pond           -6.830e-05  1.231e-05  -5.547 3.07e-08 ***
#dist_river          -5.323e-05  8.876e-06  -5.997 2.17e-09 ***
#length_forest_edges  4.797e-04  7.485e-05   6.409 1.62e-10 ***
#length_hedges        5.904e-04  6.122e-05   9.644  < 2e-16 ***
#length_roads         8.849e-04  1.149e-04   7.700 1.67e-14 ***
#per_housing          1.374e-02  9.818e-04  13.997  < 2e-16 ***
#per_pastures        -1.043e-04  1.946e-04  -0.536    0.592    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1475479)
#
#    Null deviance: 772.05  on 4386  degrees of freedom
#Residual deviance: 645.96  on 4378  degrees of freedom
#AIC: 4065.8
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.21381  -0.07193   0.13348   0.25222   0.86423  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          7.391e-01  1.515e-02  48.773  < 2e-16 ***
#dist_forest         -4.311e-04  3.879e-05 -11.116  < 2e-16 ***
#dist_pond           -6.728e-05  1.216e-05  -5.531 3.37e-08 ***
#dist_river          -5.281e-05  8.840e-06  -5.973 2.51e-09 ***
#length_forest_edges  4.840e-04  7.441e-05   6.505 8.67e-11 ***
#length_hedges        5.854e-04  6.051e-05   9.675  < 2e-16 ***
#length_roads         8.930e-04  1.139e-04   7.839 5.65e-15 ***
#per_housing          1.383e-02  9.674e-04  14.297  < 2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1475239)
#
#    Null deviance: 772.05  on 4386  degrees of freedom
#Residual deviance: 646.01  on 4379  degrees of freedom
#AIC: 4064.1
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present           10 4065.783
#reduced.model.present  9 4064.071

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


######################
#Epidalea calamita####
######################
ec=data.sp[data.sp@data$species=="Epidalea_calamita",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(ec)[,7:8]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.37452  -0.10444  -0.07617  -0.05975   1.05902  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          7.326e-02  2.270e-02   3.227  0.00129 ** 
#dist_forest          3.097e-05  4.646e-05   0.667  0.50514    
#dist_pond           -2.072e-05  1.637e-05  -1.266  0.20573    
#dist_river          -3.372e-06  1.216e-05  -0.277  0.78160    
#length_forest_edges  9.037e-04  1.488e-04   6.075  1.7e-09 ***
#length_hedges       -6.543e-05  1.088e-04  -0.601  0.54789    
#length_roads         3.132e-04  2.403e-04   1.303  0.19277    
#per_housing          6.403e-03  2.389e-03   2.680  0.00746 ** 
#per_pastures         2.633e-04  2.747e-04   0.959  0.33801    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.09185061)
#
#    Null deviance: 107.14  on 1119  degrees of freedom
#Residual deviance: 102.05  on 1111  degrees of freedom
#AIC: 515.28
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-0.36982  -0.08054  -0.07303  -0.07303   0.92697  
#
#Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)    
#(Intercept)         0.0730273  0.0102886   7.098 2.25e-12 ***
#length_forest_edges 0.0008607  0.0001317   6.533 9.75e-11 ***
#per_housing         0.0068748  0.0023178   2.966  0.00308 ** 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.09183092)
#
#    Null deviance: 107.14  on 1119  degrees of freedom
#Residual deviance: 102.58  on 1117  degrees of freedom
#AIC: 509.08
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          10 515.2831
#reduced.model.present  4 509.0753

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


#####################
#Hyla arborea####
#####################
ha=data.sp[data.sp@data$species=="Hyla_arborea",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(ha)[,7:8]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.0911  -0.1626   0.1491   0.2503   1.5329  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          7.106e-01  1.817e-02  39.121  < 2e-16 ***
#dist_forest         -2.452e-04  3.938e-05  -6.225 5.33e-10 ***
#dist_pond           -2.431e-04  1.680e-05 -14.474  < 2e-16 ***
#dist_river          -2.115e-05  1.128e-05  -1.875   0.0609 .  
#length_forest_edges  5.698e-04  9.978e-05   5.710 1.21e-08 ***
#length_hedges        7.312e-04  6.736e-05  10.856  < 2e-16 ***
#length_roads         4.408e-05  1.626e-04   0.271   0.7863    
#per_housing          1.361e-02  1.487e-03   9.153  < 2e-16 ***
#per_pastures         1.203e-03  1.852e-04   6.497 9.25e-11 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1588223)
#
#    Null deviance: 737.74  on 3812  degrees of freedom
#Residual deviance: 604.16  on 3804  degrees of freedom
#AIC: 3816
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.0913  -0.1839   0.1505   0.2554   1.4866  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          7.038e-01  1.776e-02  39.619  < 2e-16 ***
#dist_forest         -2.546e-04  3.900e-05  -6.529 7.50e-11 ***
#dist_pond           -2.534e-04  1.582e-05 -16.013  < 2e-16 ***
#length_forest_edges  5.785e-04  9.965e-05   5.805 6.96e-09 ***
#length_hedges        7.368e-04  6.617e-05  11.135  < 2e-16 ***
#per_housing          1.390e-02  1.404e-03   9.899  < 2e-16 ***
#per_pastures         1.203e-03  1.842e-04   6.530 7.44e-11 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1588871)
#
#    Null deviance: 737.74  on 3812  degrees of freedom
#Residual deviance: 604.72  on 3806  degrees of freedom
#AIC: 3815.6
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          10 3816.012
#reduced.model.present  8 3815.570

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


###################
#Rana dalmatina####
###################
rd=data.sp[data.sp@data$species=="Rana_dalmatina",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(rd)[,7:8]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.14393   0.00205   0.09651   0.19719   0.70004  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          8.095e-01  1.315e-02  61.582  < 2e-16 ***
#dist_forest         -4.472e-04  3.118e-05 -14.343  < 2e-16 ***
#dist_pond           -1.556e-04  1.196e-05 -13.011  < 2e-16 ***
#dist_river          -1.430e-07  6.761e-06  -0.021  0.98313    
#length_forest_edges  5.718e-04  6.078e-05   9.407  < 2e-16 ***
#length_hedges        5.720e-04  4.898e-05  11.678  < 2e-16 ***
#length_roads        -3.667e-04  1.301e-04  -2.818  0.00485 ** 
#per_housing          7.211e-03  1.274e-03   5.661 1.58e-08 ***
#per_pastures         1.014e-03  1.326e-04   7.650 2.34e-14 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1206039)
#
#    Null deviance: 826.84  on 5774  degrees of freedom
#Residual deviance: 695.40  on 5766  degrees of freedom
#AIC: 4184.2
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.14391   0.00207   0.09655   0.19725   0.69951  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          8.095e-01  1.293e-02  62.582  < 2e-16 ***
#dist_forest         -4.472e-04  3.093e-05 -14.459  < 2e-16 ***
#dist_pond           -1.556e-04  1.143e-05 -13.618  < 2e-16 ***
#length_forest_edges  5.718e-04  6.066e-05   9.428  < 2e-16 ***
#length_hedges        5.720e-04  4.890e-05  11.698  < 2e-16 ***
#length_roads        -3.668e-04  1.299e-04  -2.823  0.00477 ** 
#per_housing          7.211e-03  1.273e-03   5.663 1.56e-08 ***
#per_pastures         1.014e-03  1.324e-04   7.658 2.20e-14 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.120583)
#
#    Null deviance: 826.84  on 5774  degrees of freedom
#Residual deviance: 695.40  on 5767  degrees of freedom
#AIC: 4182.2
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          10 4184.201
#reduced.model.present  9 4182.201

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model
#[1] "reduced.model.present"


####################
#Rana temporaria####
####################
rt=data.sp[data.sp@data$species=="Rana_temporaria",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(rt)[,7:8]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.1297  -0.3253   0.1166   0.2728   1.1800  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          6.192e-01  2.108e-02  29.374  < 2e-16 ***
#dist_forest         -9.984e-04  5.083e-05 -19.640  < 2e-16 ***
#dist_pond            2.326e-05  1.867e-05   1.246  0.21295    
#dist_river          -1.128e-04  1.356e-05  -8.317  < 2e-16 ***
#length_forest_edges  1.026e-03  9.867e-05  10.398  < 2e-16 ***
#length_hedges        2.691e-04  9.540e-05   2.820  0.00483 ** 
#length_roads        -6.801e-05  2.196e-04  -0.310  0.75679    
#per_housing         -1.034e-03  2.580e-03  -0.401  0.68859    
#per_pastures         1.927e-03  2.346e-04   8.216 3.26e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1558075)
#
#    Null deviance: 622.36  on 2647  degrees of freedom
#Residual deviance: 411.18  on 2639  degrees of freedom
#AIC: 2602.7
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.1207  -0.3304   0.1202   0.2717   1.1588  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          6.293e-01  1.850e-02  34.020  < 2e-16 ***
#dist_forest         -9.960e-04  5.066e-05 -19.662  < 2e-16 ***
#dist_river          -1.070e-04  1.269e-05  -8.429  < 2e-16 ***
#length_forest_edges  1.024e-03  9.659e-05  10.605  < 2e-16 ***
#length_hedges        2.397e-04  8.894e-05   2.695  0.00709 ** 
#per_pastures         1.925e-03  2.259e-04   8.520  < 2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1557425)
#
#    Null deviance: 622.36  on 2647  degrees of freedom
#Residual deviance: 411.47  on 2642  degrees of freedom
#AIC: 2598.6
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          10 2602.696
#reduced.model.present  7 2598.600

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model
#[1] "reduced.model.present"


##########################
#Pelophylax esculentus####
#########################
pe=data.sp[data.sp@data$species=="Pelophylax_esculentus",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(pe)[,7:8]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.0531  -0.4243   0.1731   0.3300   1.0665  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          5.006e-01  2.393e-02  20.916  < 2e-16 ***
#dist_forest         -2.027e-04  5.138e-05  -3.946 8.16e-05 ***
#dist_pond           -2.204e-04  1.992e-05 -11.065  < 2e-16 ***
#dist_river           1.925e-05  1.314e-05   1.465    0.143    
#length_forest_edges  8.916e-04  1.312e-04   6.798 1.32e-11 ***
#length_hedges        9.626e-04  9.607e-05  10.020  < 2e-16 ***
#length_roads         4.857e-05  2.200e-04   0.221    0.825    
#per_housing          1.631e-02  2.119e-03   7.700 1.93e-14 ***
#per_pastures         1.570e-03  2.562e-04   6.126 1.04e-09 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1927976)
#
#    Null deviance: 605.37  on 2533  degrees of freedom
#Residual deviance: 486.81  on 2525  degrees of freedom
#AIC: 3030.9
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.0615  -0.4271   0.1781   0.3280   1.0283  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          5.078e-01  2.343e-02  21.676  < 2e-16 ***
#dist_forest         -1.931e-04  5.097e-05  -3.789 0.000155 ***
#dist_pond           -2.103e-04  1.872e-05 -11.234  < 2e-16 ***
#length_forest_edges  8.819e-04  1.310e-04   6.734 2.03e-11 ***
#length_hedges        9.689e-04  9.364e-05  10.347  < 2e-16 ***
#per_housing          1.631e-02  2.040e-03   7.993 1.98e-15 ***
#per_pastures         1.556e-03  2.539e-04   6.127 1.03e-09 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1928139)
#
#    Null deviance: 605.37  on 2533  degrees of freedom
#Residual deviance: 487.24  on 2527  degrees of freedom
#AIC: 3029.1
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          10 3030.911
reduced.model.present  8 3029.132

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model
#[1] "reduced.model.present"


##########################
#Pelophylax ridibundus####
##########################
pr=data.sp[data.sp@data$species=="Pelophylax_ridibundus",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(pr)[,7:8]) #extract the values from the rasters
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#create a dataframe with presence/absence values
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.1604  -0.3847   0.1623   0.3043   0.8796  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          6.433e-01  2.133e-02  30.160  < 2e-16 ***
#dist_forest         -3.745e-04  4.831e-05  -7.752 1.22e-14 ***
#dist_pond           -1.691e-04  1.763e-05  -9.595  < 2e-16 ***
#dist_river          -6.607e-05  1.279e-05  -5.166 2.54e-07 ***
#length_forest_edges  9.002e-04  1.022e-04   8.810  < 2e-16 ***
#length_hedges        9.176e-04  8.452e-05  10.856  < 2e-16 ***
#length_roads        -4.412e-05  2.008e-04  -0.220    0.826    
#per_housing          1.402e-02  1.848e-03   7.587 4.32e-14 ***
#per_pastures         8.783e-04  2.212e-04   3.970 7.35e-05 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1744785)
#
#    Null deviance: 673.95  on 3066  degrees of freedom
#Residual deviance: 533.56  on 3058  degrees of freedom
#AIC: 3359.9
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the significant variables within the model
sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.1595  -0.3843   0.1609   0.3056   0.8808  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          6.430e-01  2.128e-02  30.218  < 2e-16 ***
#dist_forest         -3.749e-04  4.827e-05  -7.768 1.08e-14 ***
#dist_pond           -1.691e-04  1.763e-05  -9.596  < 2e-16 ***
#dist_river          -6.628e-05  1.275e-05  -5.199 2.14e-07 ***
#length_forest_edges  9.002e-04  1.022e-04   8.812  < 2e-16 ***
#length_hedges        9.141e-04  8.295e-05  11.020  < 2e-16 ***
#per_housing          1.394e-02  1.810e-03   7.701 1.81e-14 ***
#per_pastures         8.848e-04  2.192e-04   4.037 5.54e-05 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1744242)
#
#    Null deviance: 673.95  on 3066  degrees of freedom
#Residual deviance: 533.56  on 3059  degrees of freedom
#AIC: 3358
#
#Number of Fisher Scoring iterations: 2

#A significant intercept in a model only means that there is a constant in the model, in addition to all other explanatory variables

#select the model with AIC value, the model with the smaller value of AIC will be selected
aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present          10 3359.915
#reduced.model.present  9 3357.963

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model
#[1] "reduced.model.present"


#########################
#Anguis fragilis####
#########################
af=data.sp[data.sp@data$species=="Anguis_fragilis",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(af)[,7:8]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-0.91220  -0.11701  -0.06921  -0.00705   1.02848  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          1.460e-01  2.133e-02   6.847 1.24e-11 ***
#dist_forest         -1.668e-04  4.368e-05  -3.819 0.000142 ***
#dist_pond           -5.009e-05  1.622e-05  -3.087 0.002069 ** 
#dist_river          -1.458e-05  1.129e-05  -1.292 0.196726    
#length_forest_edges  1.987e-04  1.475e-04   1.347 0.178244    
#length_hedges        1.969e-04  1.006e-04   1.957 0.050554 .  
#length_roads         2.786e-05  2.186e-04   0.127 0.898592    
#per_housing          2.329e-02  1.997e-03  11.659  < 2e-16 ***
#per_pastures        -8.999e-04  2.621e-04  -3.434 0.000618 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.07838321)
#
#    Null deviance: 106.345  on 1118  degrees of freedom
#Residual deviance:  87.005  on 1110  degrees of freedom
#AIC: 337.41
#
#Number of Fisher Scoring iterations: 2

sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-0.92762  -0.11700  -0.07144  -0.01284   1.00866  
#
#Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   1.694e-01  1.661e-02  10.199  < 2e-16 ***
#dist_forest  -2.024e-04  3.888e-05  -5.206 2.29e-07 ***
#dist_pond    -6.329e-05  1.455e-05  -4.352 1.48e-05 ***
#per_housing   2.338e-02  1.958e-03  11.940  < 2e-16 ***
#per_pastures -7.291e-04  2.499e-04  -2.918  0.00359 ** 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.0786175)
#
#    Null deviance: 106.34  on 1118  degrees of freedom
#Residual deviance:  87.58  on 1114  degrees of freedom
#AIC: 336.78
#
#Number of Fisher Scoring iterations: 2

aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         10 337.4113
#reduced.model.present  6 336.7761

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


#########################
#Podarcis muralis####
#########################
pm=data.sp[data.sp@data$species=="Podarcis_muralis",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(pm)[,7:8]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.19195  -0.01998   0.11554   0.20360   0.84993  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          7.790e-01  1.380e-02  56.469  < 2e-16 ***
#dist_forest         -3.510e-04  3.369e-05 -10.418  < 2e-16 ***
#dist_pond           -5.269e-05  1.022e-05  -5.157 2.60e-07 ***
#dist_river          -5.162e-05  8.283e-06  -6.232 4.92e-10 ***
#length_forest_edges  4.889e-04  5.663e-05   8.633  < 2e-16 ***
#length_hedges        5.958e-04  4.474e-05  13.319  < 2e-16 ***
#length_roads         2.070e-04  9.918e-05   2.087   0.0369 *  
#per_housing          1.376e-02  7.452e-04  18.470  < 2e-16 ***
#per_pastures        -2.405e-04  1.625e-04  -1.481   0.1388    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1209866)
#
#    Null deviance: 826.84  on 5774  degrees of freedom
#Residual deviance: 697.61  on 5766  degrees of freedom
#AIC: 4202.5
#
#Number of Fisher Scoring iterations: 2

sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
=======
  >>>>>>> 97c60f6e7d19899df0e5a0202bcd0a1e6ad26ed3
#-1.19205  -0.01786   0.11469   0.20431   0.83368  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          7.708e-01  1.265e-02  60.934  < 2e-16 ***
#dist_forest         -3.511e-04  3.369e-05 -10.421  < 2e-16 ***
#dist_pond           -5.058e-05  1.012e-05  -4.998 5.95e-07 ***
#dist_river          -5.085e-05  8.267e-06  -6.150 8.24e-10 ***
#length_forest_edges  5.027e-04  5.586e-05   9.000  < 2e-16 ***
#length_hedges        5.942e-04  4.473e-05  13.285  < 2e-16 ***
#length_roads         2.193e-04  9.883e-05   2.219   0.0265 *  
#per_housing          1.398e-02  7.301e-04  19.152  < 2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1210117)
#
#    Null deviance: 826.84  on 5774  degrees of freedom
#Residual deviance: 697.87  on 5767  degrees of freedom
#AIC: 4202.7
#
#Number of Fisher Scoring iterations: 2

aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         10 4202.499
#reduced.model.present  9 4202.694

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "model.present"


#########################
#Lacerta bilineata####
#########################
lb=data.sp[data.sp@data$species=="Lacerta_bilineata",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(lb)[,7:8]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.2132  -0.3808   0.1619   0.3016   0.9101  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          6.445e-01  2.003e-02  32.174  < 2e-16 ***
#dist_forest         -7.018e-04  4.827e-05 -14.539  < 2e-16 ***
#dist_pond           -8.710e-05  1.726e-05  -5.047 4.75e-07 ***
#dist_river          -2.217e-05  1.247e-05  -1.778   0.0755 .  
#length_forest_edges  9.051e-04  9.414e-05   9.614  < 2e-16 ***
#length_hedges        1.124e-03  7.542e-05  14.905  < 2e-16 ***
#length_roads         2.006e-04  1.756e-04   1.143   0.2533    
#per_housing          9.425e-03  1.840e-03   5.123 3.18e-07 ***
#per_pastures         5.456e-05  2.201e-04   0.248   0.8042    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1715865)
#
#    Null deviance: 682.24  on 3146  degrees of freedom
#Residual deviance: 538.44  on 3138  degrees of freedom
#AIC: 3394.7
#
#Number of Fisher Scoring iterations: 2

sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.2167  -0.3843   0.1620   0.2998   0.8943  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          6.388e-01  1.831e-02  34.890  < 2e-16 ***
#dist_forest         -7.102e-04  4.799e-05 -14.799  < 2e-16 ***
#dist_pond           -9.843e-05  1.585e-05  -6.209 6.05e-10 ***
#length_forest_edges  9.284e-04  9.317e-05   9.965  < 2e-16 ***
#length_hedges        1.154e-03  7.193e-05  16.049  < 2e-16 ***
#per_housing          1.020e-02  1.735e-03   5.876 4.64e-09 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1716606)
#
#    Null deviance: 682.24  on 3146  degrees of freedom
#Residual deviance: 539.19  on 3141  degrees of freedom
#AIC: 3393
#
#Number of Fisher Scoring iterations: 2

aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         10 3394.671
#reduced.model.present  7 3393.037

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


#########################
#Natrix helvetica####
#########################
nh=data.sp[data.sp@data$species=="Natrix_helvetica",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(nh)[,7:8]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.1271  -0.4446   0.1828   0.3365   1.0996  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          6.320e-01  2.321e-02  27.226  < 2e-16 ***
#dist_forest         -6.588e-04  5.272e-05 -12.496  < 2e-16 ***
#dist_pond           -6.529e-05  1.865e-05  -3.501 0.000471 ***
#dist_river          -5.923e-05  1.360e-05  -4.354 1.39e-05 ***
#length_forest_edges  5.954e-04  1.191e-04   4.999 6.15e-07 ***
#length_hedges        9.078e-04  9.667e-05   9.391  < 2e-16 ***
#length_roads         6.140e-04  1.961e-04   3.131 0.001760 ** 
#per_housing          1.342e-02  1.850e-03   7.255 5.29e-13 ***
#per_pastures        -1.675e-04  2.704e-04  -0.619 0.535685    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1920134)
#
#    Null deviance: 618.61  on 2621  degrees of freedom
#Residual deviance: 501.73  on 2613  degrees of freedom
#AIC: 3125.1
#
#Number of Fisher Scoring iterations: 2

sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.1261  -0.4441   0.1777   0.3347   1.1020  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          6.270e-01  2.174e-02  28.843  < 2e-16 ***
#dist_forest         -6.589e-04  5.272e-05 -12.498  < 2e-16 ***
#dist_pond           -6.352e-05  1.843e-05  -3.447 0.000575 ***
#dist_river          -5.836e-05  1.353e-05  -4.314 1.66e-05 ***
#length_forest_edges  6.004e-04  1.188e-04   5.054 4.64e-07 ***
#length_hedges        8.921e-04  9.329e-05   9.563  < 2e-16 ***
#length_roads         6.258e-04  1.951e-04   3.207 0.001359 ** 
#per_housing          1.360e-02  1.827e-03   7.444 1.32e-13 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1919681)
#
#   Null deviance: 618.61  on 2621  degrees of freedom
#Residual deviance: 501.80  on 2614  degrees of freedom
#AIC: 3123.5
#
#Number of Fisher Scoring iterations: 2

aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         10 3125.100
#reduced.model.present  9 3123.485

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


##################
#Natrix maura####
##################
nm=data.sp[data.sp@data$species=="Natrix_maura",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(nm)[,7:8]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-0.96149  -0.22166  -0.07954   0.22026   1.13076  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          2.389e-01  2.478e-02   9.641  < 2e-16 ***
#dist_forest         -2.347e-04  5.320e-05  -4.412 1.11e-05 ***
#dist_pond           -7.504e-06  1.926e-05  -0.390    0.697    
#dist_river          -1.120e-04  1.378e-05  -8.130 9.53e-16 ***
#length_forest_edges  1.805e-03  1.283e-04  14.064  < 2e-16 ***
#length_hedges        6.837e-04  1.103e-04   6.201 7.40e-10 ***
#length_roads         3.274e-04  2.269e-04   1.443    0.149    
#per_housing          1.884e-02  2.263e-03   8.326  < 2e-16 ***
#per_pastures        -2.059e-03  3.122e-04  -6.595 6.06e-11 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1257047)
#
#    Null deviance: 273.78  on 1376  degrees of freedom
#Residual deviance: 171.96  on 1368  degrees of freedom
#AIC: 1063.1
#
#Number of Fisher Scoring iterations: 2

sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-0.95260  -0.22041  -0.08126   0.22363   1.13237  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          2.380e-01  2.315e-02  10.284  < 2e-16 ***
#dist_forest         -2.347e-04  5.307e-05  -4.422 1.06e-05 ***
#dist_river          -1.132e-04  1.262e-05  -8.968  < 2e-16 ***
#length_forest_edges  1.812e-03  1.281e-04  14.141  < 2e-16 ***
#length_hedges        7.157e-04  1.080e-04   6.626 4.96e-11 ***
#per_housing          1.971e-02  2.182e-03   9.030  < 2e-16 ***
#per_pastures        -2.080e-03  3.083e-04  -6.749 2.19e-11 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1257244)
#
#    Null deviance: 273.78  on 1376  degrees of freedom
#Residual deviance: 172.24  on 1370  degrees of freedom
#AIC: 1061.3
#
#Number of Fisher Scoring iterations: 2

aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         10 1063.078
#reduced.model.present  8 1061.305

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


###########################
#Hierophis viridiflavus####
###########################
hv=data.sp[data.sp@data$species=="Hierophis_viridiflavus",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(hv)[,7:8]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.1776  -0.1389   0.1465   0.2808   0.7131  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          6.768e-01  1.726e-02  39.200  < 2e-16 ***
#dist_forest         -3.492e-04  3.723e-05  -9.380  < 2e-16 ***
#dist_pond            2.279e-05  1.318e-05   1.728 0.083982 .  
#dist_river          -3.397e-05  9.325e-06  -3.643 0.000273 ***
#length_forest_edges  6.407e-04  8.261e-05   7.756 1.10e-14 ***
#length_hedges        7.036e-04  7.028e-05  10.012  < 2e-16 ***
#length_roads         1.027e-03  1.232e-04   8.337  < 2e-16 ***
#per_housing          1.079e-02  1.110e-03   9.722  < 2e-16 ***
#per_pastures        -1.283e-03  2.243e-04  -5.721 1.13e-08 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1596917)
#
#    Null deviance: 756.39  on 4104  degrees of freedom
#Residual deviance: 654.10  on 4096  degrees of freedom
#AIC: 4129.8
#
#Number of Fisher Scoring iterations: 2

sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.1810  -0.1405   0.1466   0.2764   0.7179  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          6.895e-01  1.561e-02  44.168  < 2e-16 ***
#dist_forest         -3.476e-04  3.723e-05  -9.339  < 2e-16 ***
#dist_river          -2.816e-05  8.699e-06  -3.237  0.00122 ** 
#length_forest_edges  6.437e-04  8.261e-05   7.792 8.31e-15 ***
#length_hedges        6.884e-04  6.974e-05   9.871  < 2e-16 ***
#length_roads         1.019e-03  1.231e-04   8.276  < 2e-16 ***
#per_housing          1.069e-02  1.109e-03   9.644  < 2e-16 ***
#per_pastures        -1.345e-03  2.215e-04  -6.073 1.37e-09 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1597692)
#
#    Null deviance: 756.39  on 4104  degrees of freedom
#Residual deviance: 654.57  on 4097  degrees of freedom
#AIC: 4130.8
#
#Number of Fisher Scoring iterations: 2

aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         10 4129.811
#reduced.model.present  9 4130.804

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "model.present"


#########################
#Zamenis longissimus####
#########################
zl=data.sp[data.sp@data$species=="Zamenis_longissimus",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(zl)[,7:8]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.1148  -0.4301   0.1510   0.3772   1.1881  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          5.009e-01  2.506e-02  19.989  < 2e-16 ***
#dist_forest         -6.474e-04  5.306e-05 -12.201  < 2e-16 ***
#dist_pond            8.048e-05  1.800e-05   4.472 8.11e-06 ***
#dist_river          -1.614e-05  1.305e-05  -1.237  0.21636    
#length_forest_edges  8.641e-04  1.269e-04   6.811 1.23e-11 ***
#length_hedges        6.220e-04  1.150e-04   5.408 7.03e-08 ***
#length_roads         9.493e-04  2.107e-04   4.505 6.97e-06 ***
#per_housing          1.500e-02  1.786e-03   8.400  < 2e-16 ***
#per_pastures        -1.056e-03  3.213e-04  -3.288  0.00102 ** 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.2007628)
#
#    Null deviance: 570.45  on 2327  degrees of freedom
#Residual deviance: 465.57  on 2319  degrees of freedom
#AIC: 2879.7
#
#Number of Fisher Scoring iterations: 2

sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.1136  -0.4310   0.1515   0.3876   1.1923  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          4.930e-01  2.424e-02  20.339  < 2e-16 ***
#dist_forest         -6.562e-04  5.259e-05 -12.478  < 2e-16 ***
#dist_pond            7.405e-05  1.723e-05   4.298 1.80e-05 ***
#length_forest_edges  8.571e-04  1.268e-04   6.762 1.71e-11 ***
#length_hedges        6.409e-04  1.140e-04   5.621 2.12e-08 ***
#length_roads         9.321e-04  2.103e-04   4.433 9.74e-06 ***
#per_housing          1.518e-02  1.780e-03   8.530  < 2e-16 ***
#per_pastures        -1.017e-03  3.197e-04  -3.180  0.00149 ** 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.2008086)
#
#    Null deviance: 570.45  on 2327  degrees of freedom
#Residual deviance: 465.88  on 2320  degrees of freedom
#AIC: 2879.2
#
#Number of Fisher Scoring iterations: 2

aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         10 2879.651
#reduced.model.present  9 2879.186

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"


########################
#Vipera aspis####
########################
va=data.sp[data.sp@data$species=="Vipera_aspis",] #select the data for the focal species

presvals=raster::extract(rasters.selected,as.data.frame(va)[,7:8]) #extract the values from the rasters. Columns 7 and 8 correspond to lat and lon values.
set.seed(1963) #setting random seed
backgr=randomPoints(rasters.selected,1000) #background values
absvals=raster::extract(rasters.selected,backgr) #absence values
#créer une dataframe avec les valeurs présence/absence
pb=c(rep(1,nrow(presvals)),rep(0,nrow(absvals))) #create a vector of 0 and 1 that correspond to the number of rows of presvals and absvals
sdmdata.present=data.frame(cbind(pb,rbind(presvals,absvals))) #create the dataframe
sdmdata.present=sdmdata.present[complete.cases(sdmdata.present),]

#Generalized linear model (GLM)
model.present=glm(pb~.,data=sdmdata.present) #build the model, where the response is the occurrence of the focal species and all environmental variables are included as explanatory variables
summary(model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-0.8416  -0.2991  -0.1528   0.4189   1.1059  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          2.284e-01  2.800e-02   8.157 7.61e-16 ***
#dist_forest         -2.709e-04  6.004e-05  -4.512 6.97e-06 ***
#dist_pond           -7.381e-05  2.171e-05  -3.401 0.000691 ***
#dist_river          -1.477e-05  1.573e-05  -0.939 0.347923    
#length_forest_edges  9.702e-04  1.648e-04   5.886 4.96e-09 ***
#length_hedges        7.107e-04  1.258e-04   5.649 1.95e-08 ***
#length_roads         1.448e-03  2.637e-04   5.491 4.74e-08 ***
#per_housing          1.295e-02  2.840e-03   4.558 5.62e-06 ***
#per_pastures         4.210e-04  3.471e-04   1.213 0.225398    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1718083)
#
#    Null deviance: 285.71  on 1399  degrees of freedom
#Residual deviance: 238.99  on 1391  degrees of freedom
#AIC: 1518.1
#
#Number of Fisher Scoring iterations: 2

sign.var=summary(model.present)$coeff[-1,4]<0.05
sign.var=names(sign.var)[sign.var==T]
sign.var
reduced.sdmdata.present=subset(sdmdata.present, select=c("pb", sign.var))

#build another GLM only with significant variables
reduced.model.present=glm(pb~.,data=reduced.sdmdata.present)
summary(reduced.model.present) #print summary
#Call:
#glm(formula = pb ~ ., data = reduced.sdmdata.present)
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-0.8451  -0.2959  -0.1539   0.4280   1.1287  
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          2.332e-01  2.600e-02   8.969  < 2e-16 ***
#dist_forest         -2.745e-04  5.987e-05  -4.584 4.96e-06 ***
#dist_pond           -8.738e-05  1.920e-05  -4.552 5.79e-06 ***
#length_forest_edges  9.783e-04  1.648e-04   5.936 3.68e-09 ***
#length_hedges        7.767e-04  1.184e-04   6.562 7.45e-11 ***
#length_roads         1.389e-03  2.613e-04   5.317 1.23e-07 ***
#per_housing          1.260e-02  2.818e-03   4.473 8.35e-06 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1718906)
#
#    Null deviance: 285.71  on 1399  degrees of freedom
#Residual deviance: 239.44  on 1393  degrees of freedom
#AIC: 1516.8
#
#Number of Fisher Scoring iterations: 2

aic=AIC(model.present,reduced.model.present)
aic
#                      df      AIC
#model.present         10 1518.073
#reduced.model.present  8 1516.755

selected.model=row.names(aic)[which.min(aic$AIC)] #create an object that contains the name of the model with the smallest AIC value
selected.model #print the name of the model
#[1] "reduced.model.present"