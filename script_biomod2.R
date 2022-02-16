######################################################
#Adjust SDM models in Biomod2
#script written by Jean-Pierre Vacher 27 January 2022
#updated 27 January 2022
######################################################
#call libraries
library(here)
library(rJava)
library(biomod2)
library(ROCR)

### Read the data 
herpmod=read.csv("DATA_SPECIES_79/data_species_79_L93_for_analysis.csv") #read the csv file where all the geographic data for all species are stored
#### A column with 1: presence and 0: absence should be added to the file
#### here, no absence, so only 1 in this column
#Rana temporaria
rt.mod=herpmod[herpmod$species=="Rana_temporaria",]
rt.mod$Rana_temporaria=1

myrespname="Rana_temporaria" ### name of the species selected for the analysis
myresp=as.numeric(rt.mod[,myrespname]) ### select presence-absence data of the selected species
myrespxy= rt.mod[,c("lon","lat")] ### select lon and lat 

### import rasters of variables
list.raster=(list.files("env_variables", full.names=T,pattern=".asc"))
myexpl=raster::stack(list.raster)


setwd("biomod2") #set the working directory where the files will be stored. As we used the package “here”, the current working directory should be just above the biomod2 directory
#### create a data object for biomod
mybiomoddata=BIOMOD_FormatingData(resp.var=myresp,
              expl.var=myexpl,
              resp.xy=myrespxy,
              resp.name=myrespname,
              PA.nb.rep=1, ## number of PA data sets
              PA.nb.absences = 1000, ## number of pseudoabsence points
              PA.strategy = "random")  ## mode of selection of PA

mybiomoddata
plot(mybiomoddata)

### MODEL
#penser à mettre le maxent.jar dans le dossier de travail, ou mettre le chemin d'accès dans les options ci dessous
#Maxent citation: Steven J. Phillips, Miroslav Dudík, Robert E. Schapire. [Internet] Maxent software for modeling species niches and distributions (Version 3.4.1). Available from url: http://biodiversityinformatics.amnh.org/open_source/maxent/. Accessed on 2022-1-27.
#mybiomodoption <- BIOMOD_ModelingOptions(GAM=list("s_smoother",k=4),MAXENT.Phillips = list(betamultiplier=10, )) ## options pour chaque modèle, si on veut les modifier il faut le faire dans cet objet
mybiomodoption <- BIOMOD_ModelingOptions(GAM=list("s_smoother",k=4),RF = list(do_classif=T, ntree=500, mtry="default", nodesize=5, maxnodes=NULL),MAXENT.Phillips = list(betamultiplier=10)) ## options pour chaque modèle, si on veut les modifier il faut le faire dans cet objet
#check options here https://rdrr.io/rforge/biomod2/man/BIOMOD_ModelingOptions.html

mybiomodelout <- BIOMOD_Modeling(mybiomoddata,models = c("GAM","RF","MAXENT.Phillips"), #nom des modèles que l'on veut utiliser
                                    models.options = mybiomodoption,
                                    NbRunEval=5,   #nombre de runs par modèle/jeu de donné de PA
                                    DataSplit=70,  #% de données utilisés pour calibrer le modèle
                                    Prevalence=0.5,
                                    VarImport=3,    # nombre de permutations pour estimer l'importance des variables
                                    models.eval.meth = c('TSS','ROC'), #statistiques utilisées pour l'évaluation du modèle
                                    SaveObj = TRUE,
                                    rescal.all.models = TRUE, 
                                    do.full.models = FALSE,
                                    modeling.id = paste(myrespname,"FirstModeling",sep="")) ###modeling id = nom du modèle


##### Evaluation des modèles 
mybiomodevals=get_evaluations(mybiomodelout)
dimnames(mybiomodevals)
### pour avoir les roc
mybiomodevals["ROC","Testing.data",,,]
mybiomodevals["TSS","Testing.data",,,]


###### COURBES ROC ET AUC
######  http://rpubs.com/dgeorges/421347 (ouvrir avec edge...)
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
model.comb ## pour voir les numeros associés aux modèles 

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


par(mfrow = c(4,5)) 
lapply(mod.roc, plot)
## 1 2 ordre des plots ##
## 3 4 
## 5 6 

### contributions des variables
get_variables_importance(mybiomodelout)

### pour avoir les courbes de réponse des variables
# on charge les modèles desquels on veut les courbres de réponse
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


### PROJECTION des modèles
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

### liste des objets créer 
list.files("Rana.temporaria/proj_current")

#faire des graphiques par sélection
#plot(mybiomodproj,str.grep='GAM') # pour avoir les plots associés à la méthode

### pour récupérer la map projeter 
mycurrentproj=get_predictions(mybiomodproj)
mycurrentproj
plot(mycurrentproj)

#### build a consensus model
mybiomodem <- BIOMOD_EnsembleModeling(modeling.output = mybiomodelout,
                                      chosen.models ='all',
                                      em.by='all',
                                      eval.metric = c('TSS'),
                                      eval.metric.quality.threshold = c(0.5),
                                      prob.mean = T, ### moyenne des modèles
                                      prob.cv = T, ### coeficient de variation
                                      prob.ci = T,  ### pour avoir max et min
                                      prob.ci.alpha = 0.05,
                                      prob.median = F, ## modèle médian
                                      committee.averaging = F, 
                                      prob.mean.weight = F, ### moyenne des modèles avec poids attribué en fonction du score d'évaluation
                                      prob.mean.weight.decay ='proportional',
                                      VarImport=3)

### voir les modèles inclus dans le modèle consensus
mybiomodem

### scores d'évaluation
get_evaluations(mybiomodem)
get_variables_importance(mybiomodem)
###creer un tableau des évaluations
xx <-get_evaluations(mybiomodelout, as.data.frame=T)
yy <- get_evaluations(mybiomodem, as.data.frame=T)
rbind(xx,yy)


### Projection du modèle consensus
mybiomodef <- BIOMOD_EnsembleForecasting(EM.output = mybiomodem,
                                         projection.output = mybiomodproj)


plot(mybiomodef) 

par(mfrow=c(1,1))
plot(mybiomodef@proj@val@layers[[1]]) 


