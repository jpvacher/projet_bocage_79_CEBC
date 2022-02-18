######################################################
#Adjust SDM models in Biomod2
#script written by Jean-Pierre Vacher 27 January 2022
#updated 18 February 2022
######################################################
#call libraries####
library(here)
library(rJava)
library(biomod2)
library(ROCR)

### Read the data####
herpmod=read.csv("DATA_SPECIES_79/data_species_79_L93_for_analysis.csv") #read the csv file where all the geographic data for all species are stored
#### A column with 1: presence and 0: absence should be added to the file
#### here, no absence, so only 1 in this column

#TRITURUS MARMORATUS####
tm.mod=herpmod[herpmod$species=="Triturus_marmoratus",] #select the species of interest
tm.mod$Triturus_marmoratus=1 #add a column with occurrence coded as 1

myrespname="Triturus_marmoratus" ### name of the species selected for the analysis
myresp=as.numeric(tm.mod[,myrespname]) ### select presence-absence data of the selected species
myrespxy= tm.mod[,c("lon","lat")] ### select lon and lat 

### import rasters of variables####
list.raster=(list.files("env_variables_log_reduced", full.names=T,pattern=".asc"))
myexpl=raster::stack(list.raster)


setwd("biomod2") #set the working directory where the files will be stored. As we used the package “here”, the current working directory should be just above the biomod2 directory
#### create a data object for biomod
mybiomoddata=BIOMOD_FormatingData(resp.var=myresp,
              expl.var=myexpl,
              resp.xy=myrespxy,
              resp.name=myrespname,
              PA.nb.rep=3, ## number of PA data sets
              PA.nb.absences = 1000, ## number of pseudoabsence points
              PA.strategy = "random")  ## mode of selection of PA

mybiomoddata
plot(mybiomoddata)

### MODEL
#store the file maxent.jar within the working directory, or provide the path to the maxent.jar file in the options
#Maxent citation: Steven J. Phillips, Miroslav Dudík, Robert E. Schapire. [Internet] Maxent software for modeling species niches and distributions (Version 3.4.1). Available from url: http://biodiversityinformatics.amnh.org/open_source/maxent/. Accessed on 2022-1-27.
mybiomodoption <- BIOMOD_ModelingOptions(GAM = list("s_smother",k=4),GLM = list(type="simple", test="AIC"),RF = list(do_classif=T, ntree=500, mtry="default", nodesize=5, maxnodes=NULL), GBM = list(n.trees=1000), MARS = list(type="simple"), FDA = list(method="mars")) ## options pour chaque modèle, si on veut les modifier il faut le faire dans cet objet
#check options here https://rdrr.io/rforge/biomod2/man/BIOMOD_ModelingOptions.html
#options for maxent: MAXENT.Phillips = list(betamultiplier=10, )
mybiomodelout <- BIOMOD_Modeling(mybiomoddata,models = c("GAM","GLM","RF","GBM","MARS","FDA"), #names of the model we want to test
                                    models.options = mybiomodoption,
                                    NbRunEval=5,   #number of runs per model/PA data set
                                    DataSplit=70,  #% of data used to calibrate the model
                                    Prevalence=0.5,
                                    VarImport=3,    # number of permutations to estimate the importance of variables
                                    models.eval.meth = c('TSS','ROC'), #statistics used to valuate models
                                    SaveObj = TRUE,
                                    rescal.all.models = TRUE, 
                                    do.full.models = FALSE,
                                    modeling.id = paste(myrespname,"FirstModeling",sep="")) ###modeling id = name of the model


##### Model evaluation 
mybiomodevals=get_evaluations(mybiomodelout)
dimnames(mybiomodevals)
### Get roc and tss values
mybiomodevals["ROC","Testing.data",,,]
mybiomodevals["TSS","Testing.data",,,]

#graphical abstract
models_scores_graph(mybiomodelout, by="models", metrics=c("ROC","TSS"), main="T. marmoratus")

#store the scores in a dataframe
eval.df<-t(data.frame(RUN_1=mybiomodevals[,1,,1,],
                      RUN_2=mybiomodevals[,1,,2,],
                      RUN_3=mybiomodevals[,1,,3,],
                      RUN_4=mybiomodevals[,1,,4,],
                      RUN_5=mybiomodevals[,1,,5,]))


######ROC AND AUC CURVES
######  http://rpubs.com/dgeorges/421347 (open with edge...)
form.dat <- get_formal_data(mybiomodelout, 'resp.var')
length(form.dat)
form.dat[is.na(form.dat)]=0

pred.val <- get_predictions(mybiomodelout) 
dim(pred.val)

calib.lines <- get_calib_lines(mybiomodelout)
dim(calib.lines)

model.comb <- 
  expand.grid(
    mod = dimnames(pred.val)[[2]],
    cv = dimnames(pred.val)[[3]],
    pa = dimnames(pred.val)[[4]],
    stringsAsFactors = FALSE
  )
model.comb ## check the numbers associated with the models 

mod.roc <-
  lapply(
    1:nrow(model.comb),
    function(i){
      mod <- model.comb$mod[i]
      cv <- model.comb$cv[i]
      pa <- model.comb$pa[i]
      
      eval.lines <- !calib.lines[, paste0('_', cv), paste0('_', pa)]
      
      resp <- form.dat[eval.lines]
      pred <- pred.val[eval.lines, mod, cv, pa] / 1000
      
      pROC::roc(resp, pred)
      
    }
  )


par(mfrow = c(4,4)) 
lapply(mod.roc, plot)
## 1 2 order of the plots ##
## 3 4 
## 5 6 

### variables contributions
get_variables_importance(mybiomodelout)

### get the response curves of the variables
# cherage models for which we want the response curves
myresponscurves <- BIOMOD_LoadModels(mybiomodelout, models='GAM')


response.plot2(models=myresponscurves,
               Data=get_formal_data(mybiomodelout,"expl.var"),
               show.variables =get_formal_data(mybiomodelout,"expl.var.names")
               ,do.bivariate = FALSE,
               fixed.var.metric = 'mean',
               save.file="no",
               name="response_curve",
               ImageSize=480,
               plot=TRUE)


### PROJECTION of the models
mybiomodproj <- BIOMOD_Projection(modeling.output = mybiomodelout,
                                  new.env = myexpl,
                                  proj.name ='current',
                                  selected.models ='all',
                                  binary.meth ='ROC',
                                  compress ='xz',
                                  clamping.mask = T,
                                  output.format ='.grd')
### summary 
mybiomodproj

### list of objects to create 
list.files("Triturus.marmoratus/proj_current")

#draw graphics by selection
#plot(mybiomodproj,str.grep='GAM') # to get plot associated to the method

### get the maps 
mycurrentproj=get_predictions(mybiomodproj)
mycurrentproj
plot(mycurrentproj)

#### build a consensus model
mybiomodem <- BIOMOD_EnsembleModeling(modeling.output = mybiomodelout,
                                      chosen.models ='all',
                                      em.by='all',
                                      eval.metric = c('TSS'),
                                      eval.metric.quality.threshold = c(0.5),
                                      prob.mean = T, ### mean of the models
                                      prob.cv = T, ### coefficient of variation
                                      prob.ci = T,  ### get max and min
                                      prob.ci.alpha = 0.05,
                                      prob.median = F, ## median model
                                      committee.averaging = F, 
                                      prob.mean.weight = F, ### mean of the modelswith weight attributed according to evaluation score
                                      prob.mean.weight.decay ='proportional',
                                      VarImport=3)

### check the models included in the consensus model
mybiomodem

### evaluation scores
get_evaluations(mybiomodem)
get_variables_importance(mybiomodem)
###create a table for evaluations
xx <-get_evaluations(mybiomodelout, as.data.frame=T)
yy <- get_evaluations(mybiomodem, as.data.frame=T)
rbind(xx,yy)


### Projection of the consensus model
mybiomodef <- BIOMOD_EnsembleForecasting(EM.output = mybiomodem,
                                         projection.output = mybiomodproj)


plot(mybiomodef) 

par(mfrow=c(1,1))
plot(mybiomodef@proj@val@layers[[1]]) 


#HYLA ARBOREA####
ha.mod=herpmod[herpmod$species=="Hyla_arborea",] #select the species of interest
ha.mod$Hyla_arborea=1 #add a column with occurrence coded as 1

myrespname="Hyla_arborea" ### name of the species selected for the analysis
myresp=as.numeric(ha.mod[,myrespname]) ### select presence-absence data of the selected species
myrespxy= ha.mod[,c("lon","lat")] ### select lon and lat 

### import rasters of variables####
list.raster=(list.files("env_variables_log_reduced", full.names=T,pattern=".asc"))
myexpl=raster::stack(list.raster)


setwd("biomod2") #set the working directory where the files will be stored. As we used the package “here”, the current working directory should be just above the biomod2 directory
#### create a data object for biomod
mybiomoddata=BIOMOD_FormatingData(resp.var=myresp,
                                  expl.var=myexpl,
                                  resp.xy=myrespxy,
                                  resp.name=myrespname,
                                  PA.nb.rep=3, ## number of PA data sets
                                  PA.nb.absences = 1000, ## number of pseudoabsence points
                                  PA.strategy = "random")  ## mode of selection of PA

mybiomoddata
plot(mybiomoddata)

### MODEL
#store the file maxent.jar within the working directory, or provide the path to the maxent.jar file in the options
#Maxent citation: Steven J. Phillips, Miroslav Dudík, Robert E. Schapire. [Internet] Maxent software for modeling species niches and distributions (Version 3.4.1). Available from url: http://biodiversityinformatics.amnh.org/open_source/maxent/. Accessed on 2022-1-27.
mybiomodoption <- BIOMOD_ModelingOptions(GAM = list("s_smother",k=4),GLM = list(type="simple", test="AIC"),RF = list(do_classif=T, ntree=500, mtry="default", nodesize=5, maxnodes=NULL), GBM = list(n.trees=1000), MARS = list(type="simple"), FDA = list(method="mars")) ## options pour chaque modèle, si on veut les modifier il faut le faire dans cet objet
#check options here https://rdrr.io/rforge/biomod2/man/BIOMOD_ModelingOptions.html
#options for maxent: MAXENT.Phillips = list(betamultiplier=10, )
mybiomodelout <- BIOMOD_Modeling(mybiomoddata,models = c("GAM","GLM","RF","GBM","MARS","FDA"), #names of the model we want to test
                                 models.options = mybiomodoption,
                                 NbRunEval=5,   #number of runs per model/PA data set
                                 DataSplit=70,  #% of data used to calibrate the model
                                 Prevalence=0.5,
                                 VarImport=3,    # number of permutations to estimate the importance of variables
                                 models.eval.meth = c('TSS','ROC'), #statistics used to valuate models
                                 SaveObj = TRUE,
                                 rescal.all.models = TRUE, 
                                 do.full.models = FALSE,
                                 modeling.id = paste(myrespname,"FirstModeling",sep="")) ###modeling id = name of the model


##### Model evaluation 
mybiomodevals=get_evaluations(mybiomodelout)
dimnames(mybiomodevals)
### Get roc and tss values
mybiomodevals["ROC","Testing.data",,,]
mybiomodevals["TSS","Testing.data",,,]

#graphical abstract
models_scores_graph(mybiomodelout, by="models", metrics=c("ROC","TSS"), main="H. arborea")

#store the scores in a dataframe
eval.df<-t(data.frame(RUN_1=mybiomodevals[,1,,1,],
                      RUN_2=mybiomodevals[,1,,2,],
                      RUN_3=mybiomodevals[,1,,3,],
                      RUN_4=mybiomodevals[,1,,4,],
                      RUN_5=mybiomodevals[,1,,5,]))


######ROC AND AUC CURVES
######  http://rpubs.com/dgeorges/421347 (open with edge...)
form.dat <- get_formal_data(mybiomodelout, 'resp.var')
length(form.dat)
form.dat[is.na(form.dat)]=0

pred.val <- get_predictions(mybiomodelout) 
dim(pred.val)

calib.lines <- get_calib_lines(mybiomodelout)
dim(calib.lines)

model.comb <- 
  expand.grid(
    mod = dimnames(pred.val)[[2]],
    cv = dimnames(pred.val)[[3]],
    pa = dimnames(pred.val)[[4]],
    stringsAsFactors = FALSE
  )
model.comb ## check the numbers associated with the models 

mod.roc <-
  lapply(
    1:nrow(model.comb),
    function(i){
      mod <- model.comb$mod[i]
      cv <- model.comb$cv[i]
      pa <- model.comb$pa[i]
      
      eval.lines <- !calib.lines[, paste0('_', cv), paste0('_', pa)]
      
      resp <- form.dat[eval.lines]
      pred <- pred.val[eval.lines, mod, cv, pa] / 1000
      
      pROC::roc(resp, pred)
      
    }
  )


par(mfrow = c(4,4)) 
lapply(mod.roc, plot)
## 1 2 order of the plots ##
## 3 4 
## 5 6 

### variables contributions
get_variables_importance(mybiomodelout)

### get the response curves of the variables
# cherage models for which we want the response curves
myresponscurves <- BIOMOD_LoadModels(mybiomodelout, models='GAM')


response.plot2(models=myresponscurves,
               Data=get_formal_data(mybiomodelout,"expl.var"),
               show.variables =get_formal_data(mybiomodelout,"expl.var.names")
               ,do.bivariate = FALSE,
               fixed.var.metric = 'mean',
               save.file="no",
               name="response_curve",
               ImageSize=480,
               plot=TRUE)


### PROJECTION of the models
mybiomodproj <- BIOMOD_Projection(modeling.output = mybiomodelout,
                                  new.env = myexpl,
                                  proj.name ='current',
                                  selected.models ='all',
                                  binary.meth ='ROC',
                                  compress ='xz',
                                  clamping.mask = T,
                                  output.format ='.grd')
### summary 
mybiomodproj

### list of objects to create 
list.files("Hyla.arborea/proj_current")

#draw graphics by selection
#plot(mybiomodproj,str.grep='GAM') # to get plot associated to the method

### get the maps 
mycurrentproj=get_predictions(mybiomodproj)
mycurrentproj
plot(mycurrentproj)

#### build a consensus model
mybiomodem <- BIOMOD_EnsembleModeling(modeling.output = mybiomodelout,
                                      chosen.models ='all',
                                      em.by='all',
                                      eval.metric = c('TSS'),
                                      eval.metric.quality.threshold = c(0.5),
                                      prob.mean = T, ### mean of the models
                                      prob.cv = T, ### coefficient of variation
                                      prob.ci = T,  ### get max and min
                                      prob.ci.alpha = 0.05,
                                      prob.median = F, ## median model
                                      committee.averaging = F, 
                                      prob.mean.weight = F, ### mean of the modelswith weight attributed according to evaluation score
                                      prob.mean.weight.decay ='proportional',
                                      VarImport=3)

### check the models included in the consensus model
mybiomodem

### evaluation scores
get_evaluations(mybiomodem)
get_variables_importance(mybiomodem)
###create a table for evaluations
xx <-get_evaluations(mybiomodelout, as.data.frame=T)
yy <- get_evaluations(mybiomodem, as.data.frame=T)
rbind(xx,yy)


### Projection of the consensus model
mybiomodef <- BIOMOD_EnsembleForecasting(EM.output = mybiomodem,
                                         projection.output = mybiomodproj)


plot(mybiomodef) 

par(mfrow=c(1,1))
plot(mybiomodef@proj@val@layers[[1]]) 


#LACERTA BILINEATA####
lb.mod=herpmod[herpmod$species=="Lacerta_bilineata",] #select the species of interest
lb.mod$Lacerta_bilineata=1 #add a column with occurrence coded as 1

myrespname="Lacerta_bilineata" ### name of the species selected for the analysis
myresp=as.numeric(lb.mod[,myrespname]) ### select presence-absence data of the selected species
myrespxy= lb.mod[,c("lon","lat")] ### select lon and lat 

### import rasters of variables####
list.raster=(list.files("env_variables_log_reduced", full.names=T,pattern=".asc"))
myexpl=raster::stack(list.raster)


setwd("biomod2") #set the working directory where the files will be stored. As we used the package “here”, the current working directory should be just above the biomod2 directory
#### create a data object for biomod
mybiomoddata=BIOMOD_FormatingData(resp.var=myresp,
                                  expl.var=myexpl,
                                  resp.xy=myrespxy,
                                  resp.name=myrespname,
                                  PA.nb.rep=3, ## number of PA data sets
                                  PA.nb.absences = 1000, ## number of pseudoabsence points
                                  PA.strategy = "random")  ## mode of selection of PA

mybiomoddata
plot(mybiomoddata)

### MODEL
#store the file maxent.jar within the working directory, or provide the path to the maxent.jar file in the options
#Maxent citation: Steven J. Phillips, Miroslav Dudík, Robert E. Schapire. [Internet] Maxent software for modeling species niches and distributions (Version 3.4.1). Available from url: http://biodiversityinformatics.amnh.org/open_source/maxent/. Accessed on 2022-1-27.
mybiomodoption <- BIOMOD_ModelingOptions(GAM = list("s_smother",k=4),GLM = list(type="simple", test="AIC"),RF = list(do_classif=T, ntree=500, mtry="default", nodesize=5, maxnodes=NULL), GBM = list(n.trees=1000), MARS = list(type="simple"), FDA = list(method="mars")) ## options pour chaque modèle, si on veut les modifier il faut le faire dans cet objet
#check options here https://rdrr.io/rforge/biomod2/man/BIOMOD_ModelingOptions.html
#options for maxent: MAXENT.Phillips = list(betamultiplier=10, )
mybiomodelout <- BIOMOD_Modeling(mybiomoddata,models = c("GAM","GLM","RF","GBM","MARS","FDA"), #names of the model we want to test
                                 models.options = mybiomodoption,
                                 NbRunEval=5,   #number of runs per model/PA data set
                                 DataSplit=70,  #% of data used to calibrate the model
                                 Prevalence=0.5,
                                 VarImport=3,    # number of permutations to estimate the importance of variables
                                 models.eval.meth = c('TSS','ROC'), #statistics used to valuate models
                                 SaveObj = TRUE,
                                 rescal.all.models = TRUE, 
                                 do.full.models = FALSE,
                                 modeling.id = paste(myrespname,"FirstModeling",sep="")) ###modeling id = name of the model


##### Model evaluation 
mybiomodevals=get_evaluations(mybiomodelout)
dimnames(mybiomodevals)
### Get roc and tss values
mybiomodevals["ROC","Testing.data",,,]
mybiomodevals["TSS","Testing.data",,,]

#graphical abstract
models_scores_graph(mybiomodelout, by="models", metrics=c("ROC","TSS"), main="L. bilineata")

#store the scores in a dataframe
eval.df<-t(data.frame(RUN_1=mybiomodevals[,1,,1,],
                      RUN_2=mybiomodevals[,1,,2,],
                      RUN_3=mybiomodevals[,1,,3,],
                      RUN_4=mybiomodevals[,1,,4,],
                      RUN_5=mybiomodevals[,1,,5,]))

######ROC AND AUC CURVES
######  http://rpubs.com/dgeorges/421347 (open with edge...)
form.dat <- get_formal_data(mybiomodelout, 'resp.var')
length(form.dat)
form.dat[is.na(form.dat)]=0

pred.val <- get_predictions(mybiomodelout) 
dim(pred.val)

calib.lines <- get_calib_lines(mybiomodelout)
dim(calib.lines)

model.comb <- 
  expand.grid(
    mod = dimnames(pred.val)[[2]],
    cv = dimnames(pred.val)[[3]],
    pa = dimnames(pred.val)[[4]],
    stringsAsFactors = FALSE
  )
model.comb ## check the numbers associated with the models 

mod.roc <-
  lapply(
    1:nrow(model.comb),
    function(i){
      mod <- model.comb$mod[i]
      cv <- model.comb$cv[i]
      pa <- model.comb$pa[i]
      
      eval.lines <- !calib.lines[, paste0('_', cv), paste0('_', pa)]
      
      resp <- form.dat[eval.lines]
      pred <- pred.val[eval.lines, mod, cv, pa] / 1000
      
      pROC::roc(resp, pred)
      
    }
  )


par(mfrow = c(4,4)) 
lapply(mod.roc, plot)
## 1 2 order of the plots ##
## 3 4 
## 5 6 

### variables contributions
get_variables_importance(mybiomodelout)

### get the response curves of the variables
# cherage models for which we want the response curves
myresponscurves <- BIOMOD_LoadModels(mybiomodelout, models='RF')


response.plot2(models=myresponscurves,
               Data=get_formal_data(mybiomodelout,"expl.var"),
               show.variables =get_formal_data(mybiomodelout,"expl.var.names")
               ,do.bivariate = FALSE,
               fixed.var.metric = 'mean',
               save.file="no",
               name="response_curve",
               ImageSize=480,
               plot=TRUE)


### PROJECTION of the models
mybiomodproj <- BIOMOD_Projection(modeling.output = mybiomodelout,
                                  new.env = myexpl,
                                  proj.name ='current',
                                  selected.models ='all',
                                  binary.meth ='ROC',
                                  compress ='xz',
                                  clamping.mask = T,
                                  output.format ='.grd')
### summary 
mybiomodproj

### list of objects to create 
list.files("Lacerta.bilineata/proj_current")

#draw graphics by selection
#plot(mybiomodproj,str.grep='GAM') # to get plot associated to the method

### get the maps 
mycurrentproj=get_predictions(mybiomodproj)
mycurrentproj
plot(mycurrentproj)

#### build a consensus model
mybiomodem <- BIOMOD_EnsembleModeling(modeling.output = mybiomodelout,
                                      chosen.models ='all',
                                      em.by='all',
                                      eval.metric = c('TSS'),
                                      eval.metric.quality.threshold = c(0.5),
                                      prob.mean = T, ### mean of the models
                                      prob.cv = T, ### coefficient of variation
                                      prob.ci = T,  ### get max and min
                                      prob.ci.alpha = 0.05,
                                      prob.median = F, ## median model
                                      committee.averaging = F, 
                                      prob.mean.weight = F, ### mean of the modelswith weight attributed according to evaluation score
                                      prob.mean.weight.decay ='proportional',
                                      VarImport=3)

### check the models included in the consensus model
mybiomodem

### evaluation scores
get_evaluations(mybiomodem)
get_variables_importance(mybiomodem)
###create a table for evaluations
xx <-get_evaluations(mybiomodelout, as.data.frame=T)
yy <- get_evaluations(mybiomodem, as.data.frame=T)
rbind(xx,yy)


### Projection of the consensus model
mybiomodef <- BIOMOD_EnsembleForecasting(EM.output = mybiomodem,
                                         projection.output = mybiomodproj)


plot(mybiomodef) 

par(mfrow=c(1,1))
plot(mybiomodef@proj@val@layers[[1]]) 


#VIPERA ASPIS####
va.mod=herpmod[herpmod$species=="Vipera_aspis",] #select the species of interest
va.mod$Vipera_aspis=1 #add a column with occurrence coded as 1

myrespname="Vipera_aspis" ### name of the species selected for the analysis
myresp=as.numeric(va.mod[,myrespname]) ### select presence-absence data of the selected species
myrespxy= va.mod[,c("lon","lat")] ### select lon and lat 

### import rasters of variables####
list.raster=(list.files("env_variables_log_reduced", full.names=T,pattern=".asc"))
myexpl=raster::stack(list.raster)


setwd("biomod2") #set the working directory where the files will be stored. As we used the package “here”, the current working directory should be just above the biomod2 directory
#### create a data object for biomod
mybiomoddata=BIOMOD_FormatingData(resp.var=myresp,
                                  expl.var=myexpl,
                                  resp.xy=myrespxy,
                                  resp.name=myrespname,
                                  PA.nb.rep=3, ## number of PA data sets
                                  PA.nb.absences = 1000, ## number of pseudoabsence points
                                  PA.strategy = "random")  ## mode of selection of PA

mybiomoddata
plot(mybiomoddata)

### MODEL
#store the file maxent.jar within the working directory, or provide the path to the maxent.jar file in the options
#Maxent citation: Steven J. Phillips, Miroslav Dudík, Robert E. Schapire. [Internet] Maxent software for modeling species niches and distributions (Version 3.4.1). Available from url: http://biodiversityinformatics.amnh.org/open_source/maxent/. Accessed on 2022-1-27.
mybiomodoption <- BIOMOD_ModelingOptions(GAM = list("s_smother",k=4),GLM = list(type="simple", test="AIC"),RF = list(do_classif=T, ntree=500, mtry="default", nodesize=5, maxnodes=NULL), GBM = list(n.trees=1000), MARS = list(type="simple"), FDA = list(method="mars")) ## options pour chaque modèle, si on veut les modifier il faut le faire dans cet objet
#check options here https://rdrr.io/rforge/biomod2/man/BIOMOD_ModelingOptions.html
#options for maxent: MAXENT.Phillips = list(betamultiplier=10, )
mybiomodelout <- BIOMOD_Modeling(mybiomoddata,models = c("GAM","GLM","RF","GBM","MARS","FDA"), #names of the model we want to test
                                 models.options = mybiomodoption,
                                 NbRunEval=5,   #number of runs per model/PA data set
                                 DataSplit=70,  #% of data used to calibrate the model
                                 Prevalence=0.5,
                                 VarImport=3,    # number of permutations to estimate the importance of variables
                                 models.eval.meth = c('TSS','ROC'), #statistics used to valuate models
                                 SaveObj = TRUE,
                                 rescal.all.models = TRUE, 
                                 do.full.models = FALSE,
                                 modeling.id = paste(myrespname,"FirstModeling",sep="")) ###modeling id = name of the model


##### Model evaluation 
mybiomodevals=get_evaluations(mybiomodelout)
dimnames(mybiomodevals)
### Get roc and tss values
mybiomodevals["ROC","Testing.data",,,]
mybiomodevals["TSS","Testing.data",,,]

#graphical abstract
models_scores_graph(mybiomodelout, by="models", metrics=c("ROC","TSS"), main="V. aspis")

#store the scores in a dataframe
eval.df<-t(data.frame(RUN_1=mybiomodevals[,1,,1,],
                      RUN_2=mybiomodevals[,1,,2,],
                      RUN_3=mybiomodevals[,1,,3,],
                      RUN_4=mybiomodevals[,1,,4,],
                      RUN_5=mybiomodevals[,1,,5,]))


######ROC AND AUC CURVES
######  http://rpubs.com/dgeorges/421347 (open with edge...)
form.dat <- get_formal_data(mybiomodelout, 'resp.var')
length(form.dat)
form.dat[is.na(form.dat)]=0

pred.val <- get_predictions(mybiomodelout) 
dim(pred.val)

calib.lines <- get_calib_lines(mybiomodelout)
dim(calib.lines)

model.comb <- 
  expand.grid(
    mod = dimnames(pred.val)[[2]],
    cv = dimnames(pred.val)[[3]],
    pa = dimnames(pred.val)[[4]],
    stringsAsFactors = FALSE
  )
model.comb ## check the numbers associated with the models 

mod.roc <-
  lapply(
    1:nrow(model.comb),
    function(i){
      mod <- model.comb$mod[i]
      cv <- model.comb$cv[i]
      pa <- model.comb$pa[i]
      
      eval.lines <- !calib.lines[, paste0('_', cv), paste0('_', pa)]
      
      resp <- form.dat[eval.lines]
      pred <- pred.val[eval.lines, mod, cv, pa] / 1000
      
      pROC::roc(resp, pred)
      
    }
  )


par(mfrow = c(4,4)) 
lapply(mod.roc, plot)
## 1 2 order of the plots ##
## 3 4 
## 5 6 

### variables contributions
get_variables_importance(mybiomodelout)

### get the response curves of the variables
# cherage models for which we want the response curves
myresponscurves <- BIOMOD_LoadModels(mybiomodelout, models='GAM')


response.plot2(models=myresponscurves,
               Data=get_formal_data(mybiomodelout,"expl.var"),
               show.variables =get_formal_data(mybiomodelout,"expl.var.names")
               ,do.bivariate = FALSE,
               fixed.var.metric = 'mean',
               save.file="no",
               name="response_curve",
               ImageSize=480,
               plot=TRUE)


### PROJECTION of the models
mybiomodproj <- BIOMOD_Projection(modeling.output = mybiomodelout,
                                  new.env = myexpl,
                                  proj.name ='current',
                                  selected.models ='all',
                                  binary.meth ='ROC',
                                  compress ='xz',
                                  clamping.mask = T,
                                  output.format ='.grd')
### summary 
mybiomodproj

### list of objects to create 
list.files("Vipera.aspis/proj_current")

#draw graphics by selection
#plot(mybiomodproj,str.grep='GAM') # to get plot associated to the method

### get the maps 
mycurrentproj=get_predictions(mybiomodproj)
mycurrentproj
plot(mycurrentproj)

#### build a consensus model
mybiomodem <- BIOMOD_EnsembleModeling(modeling.output = mybiomodelout,
                                      chosen.models ='all',
                                      em.by='all',
                                      eval.metric = c('TSS'),
                                      eval.metric.quality.threshold = c(0.5),
                                      prob.mean = T, ### mean of the models
                                      prob.cv = T, ### coefficient of variation
                                      prob.ci = T,  ### get max and min
                                      prob.ci.alpha = 0.05,
                                      prob.median = F, ## median model
                                      committee.averaging = F, 
                                      prob.mean.weight = F, ### mean of the modelswith weight attributed according to evaluation score
                                      prob.mean.weight.decay ='proportional',
                                      VarImport=3)

### check the models included in the consensus model
mybiomodem

### evaluation scores
get_evaluations(mybiomodem)
get_variables_importance(mybiomodem)
###create a table for evaluations
xx <-get_evaluations(mybiomodelout, as.data.frame=T)
yy <- get_evaluations(mybiomodem, as.data.frame=T)
rbind(xx,yy)


### Projection of the consensus model
mybiomodef <- BIOMOD_EnsembleForecasting(EM.output = mybiomodem,
                                         projection.output = mybiomodproj)


plot(mybiomodef) 

par(mfrow=c(1,1))
plot(mybiomodef@proj@val@layers[[1]]) 
