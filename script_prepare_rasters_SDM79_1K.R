######################################################
#PREPARE ENVIRONMENT VARIABLES FOR SDM BOCAGE
#script written by Jean-Pierre Vacher 17 February 2022
#updated 17 February 2022
######################################################

#charge the packages####
x=c("here","sf","dismo","ggmap", "rgdal", "rgeos", "maptools", "dplyr", "tidyr", "tmap", "raster","mapdata","sp","spdep","colorRamps","ggplot2","gridExtra","virtualspecies", "usdm") #create a list with names of packages
lapply(x, library, character.only=TRUE) #loop that read all the packages from the list


##########################
#ENVIRONMENT VARIABLES####
##########################

#Deux-Sevres department####
#create a vector layer for the delineation of Deux-Sevres
#dpt=readOGR(dsn="departements-20180101.shp",layer="departements-20180101") #open the layer
#dpt=spTransform(dpt,CRS("+init=epsg:2154")) #assign Lambert 93 geographic projection
#DS=dpt[dpt$code_insee=="79",] #select Deux-Sèvres
#raster::shapefile(DS, file="contour_deux_sevres_L93.shp") #save layer in the working directory
#DS=spTransform(DS, CRS("+init=epsg:4326")) #assign WGS84 geographic projection
#raster::shapefile(DS, file="contour_deux_sevres_WGS84.shp") #save layer in the working directory
#DSWGS84=readOGR(dsn="contour_deux_sevres_WGS84.shp", layer="contour_deux_sevres_WGS84") #read the layer
DSL93=readOGR(dsn="GIS_layers/contour_deux_sevres_L93.shp", layer="contour_deux_sevres_L93") #read the layer

#For the measurment of distances, we will split the department into 16 areas to optimize the computing time. For this, we will create a polygon that covers the whole department, and we will split this polygon into 16 sub-polygons, into which we will perform the distance calculation later.
#First, we create a rectangle polygon that fits the extent of Deux Sevres department
poly=data.frame(cbind(lon=c(extent(DSL93)[1],extent(DSL93)[1],extent(DSL93)[2],extent(DSL93)[2]), lat=c(extent(DSL93)[3], extent(DSL93)[4], extent(DSL93)[4], extent(DSL93)[3]))) #create a data frame that contains the extent values of DSL93
coordinates(poly)=~lon+lat #create a vector layer by defining the coordinates
proj4string(poly)=CRS("+init=epsg:2154") #assign the geographic projection system
poly #it is an object of class SpatialPoints. We need to transform it to a polygon
poly=Polygon(poly) #create a polygon out of the points
poly=Polygons(list(poly),1) #add slots to the polygon
poly=SpatialPolygons(list(poly)) #transform to a spatial polygon
proj4string(poly)=CRS("+init=epsg:2154") #we assign a coordinate system that corresponds to WGS84 (epsg:4326)
poly #check what it looks like
#if we want to check graphically how it looks like
#plot(DSL93)
#plot(poly,add=T)
#now we will split the polygon into 16 sub-polygons
newext1=c(extent(poly)[1],extent(poly)[2],extent(poly)[3], extent(poly)[3]+((extent(poly)[4]-extent(poly)[3])/16)) #define a new extent that is the fifth of the actual poly extent.
newext2=c(extent(poly)[1],extent(poly)[2], extent(poly)[3]+((extent(poly)[4]-extent(poly)[3])/16), extent(poly)[3]+(2*(extent(poly)[4]-extent(poly)[3])/16)) #define a new extent that is the fourth of the actual poly extent.
newext3=c(extent(poly)[1],extent(poly)[2],extent(poly)[3]+(2*(extent(poly)[4]-extent(poly)[3])/16),(extent(poly)[3]+(3*(extent(poly)[4]-extent(poly)[3])/16))) #define a new extent that is the fourth of the actual poly extent
newext4=c(extent(poly)[1],extent(poly)[2],(extent(poly)[3]+(3*(extent(poly)[4]-extent(poly)[3])/16)),(extent(poly)[3]+(4*(extent(poly)[4]-extent(poly)[3])/16))) #define a new extent that is the fourth of the actual poly extent
newext5=c(extent(poly)[1],extent(poly)[2],(extent(poly)[3]+(4*(extent(poly)[4]-extent(poly)[3])/16)),(extent(poly)[3]+(5*(extent(poly)[4]-extent(poly)[3])/16))) #define a new extent that is the fourth of the actual poly extent
newext6=c(extent(poly)[1],extent(poly)[2],(extent(poly)[3]+(5*(extent(poly)[4]-extent(poly)[3])/16)),(extent(poly)[3]+(6*(extent(poly)[4]-extent(poly)[3])/16))) #define a new extent that is the fourth of the actual poly extent
newext7=c(extent(poly)[1],extent(poly)[2],(extent(poly)[3]+(6*(extent(poly)[4]-extent(poly)[3])/16)),(extent(poly)[3]+(7*(extent(poly)[4]-extent(poly)[3])/16))) #define a new extent that is the fourth of the actual poly extent
newext8=c(extent(poly)[1],extent(poly)[2],(extent(poly)[3]+(7*(extent(poly)[4]-extent(poly)[3])/16)),(extent(poly)[3]+(8*(extent(poly)[4]-extent(poly)[3])/16))) #define a new extent that is the fourth of the actual poly extent
newext9=c(extent(poly)[1],extent(poly)[2],(extent(poly)[3]+(8*(extent(poly)[4]-extent(poly)[3])/16)),(extent(poly)[3]+(9*(extent(poly)[4]-extent(poly)[3])/16))) #define a new extent that is the fourth of the actual poly extent
newext10=c(extent(poly)[1],extent(poly)[2],(extent(poly)[3]+(9*(extent(poly)[4]-extent(poly)[3])/16)),(extent(poly)[3]+(10*(extent(poly)[4]-extent(poly)[3])/16))) #define a new extent that is the fourth of the actual poly extent
newext11=c(extent(poly)[1],extent(poly)[2],(extent(poly)[3]+(10*(extent(poly)[4]-extent(poly)[3])/16)),(extent(poly)[3]+(11*(extent(poly)[4]-extent(poly)[3])/16))) #define a new extent that is the fourth of the actual poly extent
newext12=c(extent(poly)[1],extent(poly)[2],(extent(poly)[3]+(11*(extent(poly)[4]-extent(poly)[3])/16)),(extent(poly)[3]+(12*(extent(poly)[4]-extent(poly)[3])/16))) #define a new extent that is the fourth of the actual poly extent
newext13=c(extent(poly)[1],extent(poly)[2],(extent(poly)[3]+(12*(extent(poly)[4]-extent(poly)[3])/16)),(extent(poly)[3]+(13*(extent(poly)[4]-extent(poly)[3])/16))) #define a new extent that is the fourth of the actual poly extent
newext14=c(extent(poly)[1],extent(poly)[2],(extent(poly)[3]+(13*(extent(poly)[4]-extent(poly)[3])/16)),(extent(poly)[3]+(14*(extent(poly)[4]-extent(poly)[3])/16))) #define a new extent that is the fourth of the actual poly extent
newext15=c(extent(poly)[1],extent(poly)[2],(extent(poly)[3]+(14*(extent(poly)[4]-extent(poly)[3])/16)),(extent(poly)[3]+(15*(extent(poly)[4]-extent(poly)[3])/16))) #define a new extent that is the fourth of the actual poly extent
newext16=c(extent(poly)[1],extent(poly)[2],(extent(poly)[3]+(15*(extent(poly)[4]-extent(poly)[3])/16)),(extent(poly)[3]+(16*(extent(poly)[4]-extent(poly)[3])/16))) #define a new extent that is the fourth of the actual poly extent
poly1=crop(poly, newext1) #create a sub-polygon with the new extent 1
poly2=crop(poly, newext2) #create a sub-polygon with the new extent 2
poly3=crop(poly, newext3) #create a sub-polygon with the new extent 3
poly4=crop(poly, newext4) #create a sub-polygon with the new extent 4
poly5=crop(poly, newext5) #create a sub-polygon with the new extent 5
poly6=crop(poly, newext6) #create a sub-polygon with the new extent 6
poly7=crop(poly, newext7) #create a sub-polygon with the new extent 7
poly8=crop(poly, newext8) #create a sub-polygon with the new extent 8
poly9=crop(poly, newext9) #create a sub-polygon with the new extent 9
poly10=crop(poly, newext10) #create a sub-polygon with the new extent 10
poly11=crop(poly, newext11) #create a sub-polygon with the new extent 11
poly12=crop(poly, newext12) #create a sub-polygon with the new extent 12
poly13=crop(poly, newext13) #create a sub-polygon with the new extent 13
poly14=crop(poly, newext14) #create a sub-polygon with the new extent 14
poly15=crop(poly, newext15) #create a sub-polygon with the new extent 15
poly16=crop(poly, newext16) #create a sub-polygon with the new extent 16

#check graphically if it worked and how it looks like
plot(poly)
plot(poly1, add=T, col=alpha("yellow",0.7), border=F)
plot(poly2, add=T, col=alpha("red",0.7), border=F)
plot(poly3, add=T, col=alpha("green",0.7), border=F)
plot(poly4, add=T, col=alpha("blue",0.7), border=F)
plot(poly5, add=T, col=alpha("turquoise",0.7), border=F)
plot(poly6, add=T, col=alpha("pink",0.7), border=F)
plot(poly7, add=T, col=alpha("forestgreen",0.7), border=F)
plot(poly8, add=T, col=alpha("grey",0.7), border=F)
plot(poly9, add=T, col=alpha("violet",0.7), border=F)
plot(poly10, add=T, col=alpha("orange",0.7), border=F)
plot(poly11, add=T, col=alpha("red4",0.7), border=F)
plot(poly12, add=T, col=alpha("dodgerblue",0.7), border=F)
plot(poly13, add=T, col=alpha("yellow2",0.7), border=F)
plot(poly14, add=T, col=alpha("red1",0.7), border=F)
plot(poly15, add=T, col=alpha("orange1",0.7), border=F)
plot(poly16, add=T, col=alpha("grey2",0.7), border=F)

#Now, we intersect the sub-polygons with the department to get sub-regions in which we will calculate the distances
DS1=crop(DSL93, poly1)
DS2=crop(DSL93, poly2)
DS3=crop(DSL93, poly3)
DS4=crop(DSL93, poly4)
DS5=crop(DSL93, poly5)
DS6=crop(DSL93, poly6)
DS7=crop(DSL93, poly7)
DS8=crop(DSL93, poly8)
DS9=crop(DSL93, poly9)
DS10=crop(DSL93, poly10)
DS11=crop(DSL93, poly11)
DS12=crop(DSL93, poly12)
DS13=crop(DSL93, poly13)
DS14=crop(DSL93, poly14)
DS15=crop(DSL93, poly15)
DS16=crop(DSL93, poly16)

#check graphically if it worked and how it looks like
#plot(DSL93, border=F)
#plot(DS1,add=T, col=alpha("yellow",0.7))
#plot(DS2, add=T, col=alpha("red",0.7))
#plot(DS3, add=T, col=alpha("green",0.7))
#plot(DS4, add=T, col=alpha("blue",0.7))
#plot(DS5, add=T, col=alpha("turquoise",0.7))
#plot(DS6, add=T, col=alpha("violet",0.7))
#plot(DS7, add=T, col=alpha("dodgerblue",0.7))
#plot(DS8, add=T, col=alpha("orange",0.7))
#plot(DS9, add=T, col=alpha("forestgreen",0.7))
#plot(DS10, add=T, col=alpha("grey",0.7))
#plot(DS11, add=T, col=alpha("red3",0.7))
#plot(DS12, add=T, col=alpha("pink",0.7))
#plot(DS13, add=T, col=alpha("yellow",0.7))
#plot(DS14, add=T, col=alpha("red",0.7))
#plot(DS15, add=T, col=alpha("dodgerblue",0.7))

#OK, now we have sub-regions of the department that can be used for calculating distances


#################
#BUILD A GRID####
#################

grid=raster(DSL93) #create a raster out of the delineation of the Deux-Sèvres
res(grid)=1000 #1000 m resolution grid (cells 1000 m x 1000 m)

#transform the grid from raster to vector layer
#grid.poly=rasterToPolygons(grid)
#grid.poly@data$ID=c(1:length(grid.poly))
#grid.poly@data$layer=NULL
#grid.poly@data$area=raster::area(grid.poly)
#grid.poly
#raster::shapefile(grid.poly, file="GIS_layers/grid_DeuxSevres_1K.shp")
grid.poly=readOGR(dsn="GIS_layers/grid_DeuxSevres_1K.shp", layer="grid_DeuxSevres_1K") #read the grid as a vector

#calculate centroids of the cells of the grid
#centroids=coordinates(grid.poly) #extract coordinates of the cells of the grid
#centroids=as.data.frame(centroids) #transform the list into a dataframe
#colnames(centroids)=c("x","y") #rename columns
#coordinates(centroids)=~x+y #define the coordinates of the layer
#proj4string(centroids)=CRS("+init=epsg:2154") #assign a projection system that corresponds to Lambert 93
#centroids #check what it looks like
#raster::shapefile(centroids, file="centroids_grid1K_L93.shp") #save the layer
centroids=readOGR(dsn="GIS_layers/centroids_grid1K_L93.shp", layer="centroids_grid1K_L93") #read the centroids layer

#We split the centroid layer to match the splited department layer into 16 sub-regions
centroids1=crop(centroids, poly1)
centroids2=crop(centroids, poly2)
centroids3=crop(centroids, poly3)
centroids4=crop(centroids, poly4)
centroids5=crop(centroids, poly5)
centroids6=crop(centroids, poly6)
centroids7=crop(centroids, poly7)
centroids8=crop(centroids, poly8)
centroids9=crop(centroids, poly9)
centroids10=crop(centroids, poly10)
centroids11=crop(centroids, poly11)
centroids12=crop(centroids, poly12)
centroids13=crop(centroids, poly13)
centroids14=crop(centroids, poly14)
centroids15=crop(centroids, poly15)
centroids16=crop(centroids, poly16)


############################
#BUILD VARIABLES RASTERS####
############################

###########
#RIVERS####
###########
river=readOGR(dsn="GIS_layers/TOPAGE_2019_DeuxSevres.shp",layer="TOPAGE_2019_DeuxSevres") #build from the layer TOPAGE (version 2019)
river

#Length of river per hectare
intersect.river=raster::intersect(river, grid.poly)
intersect.river=st_as_sf(intersect.river)
#intersect.river=st_read("GIS_layers/intersect_river_100m.shp") #read the layer with sf
intersect.river$length=st_length(intersect.river) #calculate the length of each line
length=data.frame(cbind(ID=intersect.river$ID,length=as.numeric(as.character(intersect.river$length)))) #transform the field “length” as numeric
length$ID=as.factor(length$ID)
length=aggregate(length~ID, sum, data=length) #calculate the sum of dens of each line within each cell of the grid
df=data.frame(grid.poly@data$ID) #extract the attribute table of grid.poly as a data frame
colnames(df)="ID" #provide a column name
df2=merge(df, length, by="ID", all.x=T, sort=T) #merge both data frame by field ID
df2[is.na(df2)]=0 #replace na by 0
df2$ID=as.factor(df2$ID) #change the ID column to factor
df2$dens=(df2$length*10000)/1000000
grid.poly.river=grid.poly #create a new object that is the same as grid.poly
grid.poly.river@data$dens=df2$dens[match(df2$ID, grid.poly.river@data$ID)] #transfer the dens values of sf2 to grid.poly.rivers. The ones with no density values will remain at 0.
grid.poly.river@data$log_dens=log(grid.poly.river@data$dens+1) #transfer the dens values of sf2 to grid.poly.rivers. The ones with no density values will remain at 0.

r.dens=raster(grid.poly.river) #create a raster from the object grid.poly.rivers
res(r.dens)=1000 #resolution of the grid is cells of 1000 m x 1000 m
values(r.dens)= grid.poly.river@data$dens #assign the values to the grids of the raster
r.dens=mask(r.dens, DSL93) #crop the raster to the delineation of the Deux-Sèvres department
writeRaster(r.dens, "env_variables_1K/dens_river.asc", format="ascii") #save raster

#log transform
r.dens=raster(grid.poly.river) #create a raster from the object grid.poly.rivers
res(r.dens)=1000 #resolution of the grid is cells of 1000 m x 1000 m
values(r.dens)= grid.poly.river@data$log_dens #assign the values to the grids of the raster
r.dens=mask(r.dens, DSL93) #crop the raster to the delineation of the Deux-Sèvres department
writeRaster(r.dens, "env_variables_1K_log/dens_river.asc", format="ascii") #save raster

#Distance to rivers####
#dist.river=apply(gDistance(centroids, river,byid=T),2,min)
#grid.river2=grid
#values(grid.river2)=dist.river
#grid.river2[is.na(grid.river2[])]=0
#projection(grid.river2)=CRS("+init=epsg:2154")
#grid.river2=mask(grid.river2,DSL93)
#writeRaster(grid.river2,"env_variables_1K/dist_river.asc", format="ascii")


##########
#ROADS####
##########
road=readOGR(dsn="GIS_layers/ROUTE500_2020_DeuxSevres.shp",layer="ROUTE500_2020_DeuxSevres") #build from the layer “TRONCON_ROUTE”of the ROUTE500 layer (version 2020)
road

#Density of roads####
intersect.road=raster::intersect(road, grid.poly)
intersect.road=st_as_sf(intersect.road)
#intersect.road=st_read("GIS_layers/intersect_road_100m.shp") #read the layer with sf
intersect.road$length=st_length(intersect.road)
length=data.frame(cbind(ID=intersect.road$ID,length=as.numeric(as.character(intersect.road$length)))) #transform the field “length” as numeric
length$ID=as.factor(length$ID)
length=aggregate(length~ID, sum, data=length) #calculate the sum of dens of each line within each cell of the grid
df=data.frame(grid.poly@data$ID) #extract the attribute table of grid.poly as a data frame
colnames(df)="ID" #provide a column name
df2=merge(df, length, by="ID", all.x=T, sort=T) #merge both data frame by field ID
df2[is.na(df2)]=0 #replace na by 0
df2$ID=as.factor(df2$ID) #change the ID column to factor
df2$dens=(df2$length*10000)/1000000
grid.poly.road=grid.poly #create a new object that is the same as grid.poly
grid.poly.road@data$dens=df2$dens[match(df2$ID, grid.poly.road@data$ID)] #transfer the dens values of sf2 to grid.poly.rivers. The ones with no density values will remain at 0.
grid.poly.road@data$log_dens=log(grid.poly.road@data$dens+1)

r.road=raster(grid.poly.road)
res(r.road)=1000
values(r.road)=grid.poly.road@data$dens
r.road=mask(r.road, DSL93)
writeRaster(r.road, "env_variables_1K/dens_roads.asc", format="ascii") #save raster

#log-transformed
r.road=raster(grid.poly.road)
res(r.road)=1000
values(r.road)=grid.poly.road@data$log_dens
r.road=mask(r.road, DSL93)
writeRaster(r.road, "env_variables_1K_log/dens_roads.asc", format="ascii") #save raster


###########
#HEDGES####
###########

haie=readOGR(dsn="GIS_layers/DNSB-HAIES_1-0__SHP_L93_D079_2020-06-24/DNSB-HAIES/1_DONNEES_LIVRAISON_2020-06-24/HAIE-LINEAIRE.shp", layer="HAIE-LINEAIRE")
haie

#Density of hedges####
#intersect.haie=raster::intersect(haie, grid.poly)
#intersect.haie=st_as_sf(intersect.haie)
#Process the intersection between grid.poly and hedges in QGIS. Open both layers, then vector -> geoprocessing tools -> intersection. First layer is hedges, second layer is grid. Then save as “intersect_haie_1K.shp”
intersect.haie=st_read("GIS_layers/intersect_haie_1K.shp") #read the layer with sf
intersect.haie$length=st_length(intersect.haie)
length=data.frame(cbind(ID=intersect.haie$ID,length=as.numeric(as.character(intersect.haie$length)))) #transform the field “length” as numeric
length$ID=as.factor(length$ID)
length=aggregate(length~ID, sum, data=length) #calculate the sum of dens of each line within each cell of the grid
df=data.frame(grid.poly@data$ID) #extract the attribute table of grid.poly as a data frame
colnames(df)="ID" #provide a column name
df2=merge(df, length, by="ID", all.x=T, sort=T) #merge both data frame by field ID
df2[is.na(df2)]=0 #replace na by 0
df2$ID=as.factor(df2$ID) #change the ID column to factor
df2$dens=(df2$length*10000)/1000000 #calculte density per ha
grid.poly.hedge=grid.poly #create a new object that is the same as grid.poly
grid.poly.hedge@data$dens=df2$dens[match(df2$ID, grid.poly.hedge@data$ID)] #transfer the dens values of sf2 to grid.poly.rivers. The ones with no density values will remain at 0.
grid.poly.hedge@data$log_dens=log(grid.poly.hedge@data$dens+1)

r.hedge=raster(grid.poly.hedge)
res(r.hedge)=1000
values(r.hedge)=grid.poly.hedge@data$dens
r.hedge=mask(r.hedge, DSL93)
writeRaster(r.hedge, "env_variables_1K/dens_hedges.asc", format="ascii") #save raster

#log-transformed
r.hedge=raster(grid.poly.hedge)
res(r.hedge)=1000
values(r.hedge)=grid.poly.hedge@data$log_dens
r.hedge=mask(r.hedge, DSL93)
writeRaster(r.hedge, "env_variables_1K_log/dens_hedges.asc", format="ascii") #save raster


##############
#PASTURES####
##############

#rpg=readOGR(dsn="~/Documents/jpv/projet_CEBC/RPG_2-0_SHP_LAMB93_R75-2019/RPG/1_DONNEES_LIVRAISON_2019/RPG_2-0_SHP_LAMB93_R75-2019/PARCELLES_GRAPHIQUES_R75.shp", layer="PARCELLES_GRAPHIQUES_R75")
#rpg79=crop(rpg, DS) #crop to fit the study site (lighter file)
#raster::shapefile(rpg79, file="RPG_79.shp") #save the layer
rpg79=readOGR(dsn="GIS_layers/RPG_79.shp",layer="RPG_79")
#rpg79=spTransform(rpg79, CRS("+init=epsg:4326"))

pasture=rbind(rpg79[rpg79@data$CODE_GROUP=="18",], rpg79[rpg79@data$CODE_GROUP=="19",])
#raster::shapefile(pasture, file="GIS_layers/pasture79.shp")
#In QGIS, open the RPG79 layer, open query and enter the SQL code CODE_GROUP = ‘18’ OR CODE_GROUP = ‘19’, then open the grid_DeuxSevres_100m.shp, then menu vector -> Geoprocessing Tools -> intersect, and save the resulting layer as intersect_pastures_100m.shp

#intersect.pasture=raster::intersect(pasture, grid.poly)
#intersect.pasture=st_as_sf(intersect.haie)
intersect.pasture=st_read("GIS_layers/intersect_pasture_1K.shp") #read the layer with sf
intersect.pasture$area2=st_area(intersect.pasture) #
area=as.integer(as.character(intersect.pasture$area))
area2=as.integer(as.character(intersect.pasture$area2))
per.pasture=data.frame(cbind(ID=intersect.pasture$ID,per=round((area2*100)/area),digits=4)) #bind the ID column with another column where we calculate the percentage
per.pasture$ID=as.factor(per.pasture$ID)
per.pasture$per=as.numeric(as.character(per.pasture$per))
per.pasture=aggregate(per~ID, sum, data=per.pasture)
df=data.frame(ID=as.factor(grid.poly@data$ID))
colnames(df)="ID"
df2=merge(df, per.pasture, by="ID", all.x=T, sort=T)
df2[is.na(df2)]=0
df2$ID=as.factor(df2$ID)
grid.poly.pasture=grid.poly
grid.poly.pasture@data$per=df2$per[match(df2$ID, grid.poly.pasture@data$ID)]
grid.poly.pasture@data$log_per=log(grid.poly.pasture@data$per+1)

#transform the vector layer into a raster
grid.pasture=grid
res(grid.pasture)=1000
values(grid.pasture)=grid.poly.pasture@data$per
grid.pasture=mask(grid.pasture, DSL93)
writeRaster(grid.pasture,"env_variables_1K/per_pastures.asc", format="ascii")

#log-transformed
grid.pasture=grid
res(grid.pasture)=1000
values(grid.pasture)=grid.poly.pasture@data$log_per
grid.pasture=mask(grid.pasture, DSL93)
writeRaster(grid.pasture,"env_variables_1K_log/per_pastures.asc", format="ascii")


################
#CROP FIELDS####
################
#rpg79=readOGR(dsn="RPG_79.shp",layer="RPG_79")
#agri=rbind(rpg79[rpg79@data$CODE_GROUP=="1",], rpg79[rpg79@data$CODE_GROUP=="2",], rpg79[rpg79@data$CODE_GROUP=="3",], rpg79[rpg79@data$CODE_GROUP=="4",], rpg79[rpg79@data$CODE_GROUP=="5",], rpg79[rpg79@data$CODE_GROUP=="6",], rpg79[rpg79@data$CODE_GROUP=="7",])
#raster::shapefile(agri, file="GIS_layers/cultures79.shp")
#In QGIS, open the RPG79 layer, open query and enter the SQL code CODE_GROUP = ‘1’ OR CODE_GROUP = ‘2’ OR CODE_GROUP = ‘3’ OR CODE_GROUP = ‘4’ OR CODE_GROUP = ‘5’ OR CODE_GROUP = ‘6’ OR CODE_GROUP = ‘7’, then open the grid_DeuxSevres_100m.shp, then menu vector -> Geoprocessing Tools -> intersect, and save the resulting layer as intersect_pastures_100m.shp

intersect.agri=st_read("GIS_layers/intersect_agri_1K.shp")
intersect.agri$area2=st_area(intersect.agri)
area=as.integer(as.character(intersect.agri$area))
area2=as.integer(as.character(intersect.agri$area2))
per.agri=data.frame(cbind(ID=intersect.agri$ID,per=round(((area2*100)/area),digits=4)))
per.agri$ID=as.factor(per.agri$ID)
per.agri$per=as.numeric(as.character(per.agri$per))
per.agri=aggregate(per~ID, sum, data=per.agri)
df=data.frame(ID=as.factor(grid.poly@data$ID))
colnames(df)="ID"
df2=merge(df, per.agri, by="ID", all.x=T, sort=T)
df2[is.na(df2)]=0
df2$ID=as.factor(df2$ID)
grid.poly.agri=grid.poly
grid.poly.agri@data$per=df2$per[match(df2$ID, grid.poly.agri@data$ID)]
grid.poly.agri@data$log_per=log(grid.poly.agri@data$per+1)

#transform the vector layer into a raster
grid.agri=grid
res(grid.agri)=1000
values(grid.agri)=grid.poly.agri@data$per
grid.agri=mask(grid.agri, DSL93)
writeRaster(grid.agri,"env_variables_1K/per_agri.asc", format="ascii")

#log-transformed
grid.agri=grid
res(grid.agri)=1000
values(grid.agri)=grid.poly.agri@data$log_per
grid.agri=mask(grid.agri, DSL93)
writeRaster(grid.agri,"env_variables_1K_log/per_agri.asc", format="ascii")


##################
#FORESTS EDGES####
##################

#forest.edge=readOGR(dsn="lisieres_foret_L93.shp", layer="lisieres_foret_L93")
#shrub.edge=readOGR(dsn="lisieres_bosquet_L93.shp", layer="lisieres_bosquet_L93")
#proj4string(shrub.edge)=crs(forest.edge)
#forest.edge=raster::bind(forest.edge, shrub.edge)
#raster::shapefile(forest.edge, file="forets_bosquets_edges_79_L93.shp")

forest.edge=readOGR(dsn="GIS_layers/forets_bosquets_edges_79_L93.shp", layer="forets_bosquets_edges_79_L93")

#Density of forest edges
#intersect.forest.edge=raster::intersect(forest.edge, grid.poly)
#intersect.forest.edge=st_as_sf(intersect.forest.edge)
#or
#forest.edge=st_as_sf(forest.edge)
#grid.poly.red=st_as_sf(grid.poly)
#intersect.forest.edge=st_intersection(forest.edge, grid.poly)
#or
#forest.edge=st_read("forets_bosquets_edges_79_L93.shp")
#grid.poly=st_read("grid_DeuxSevres.shp")
#intersect.forest.edge=st_intersection(forest.edge, grid.poly.red)
#or process this operation in QGIS with the intersection function in the vector menu
intersect.forest.edge=st_read("GIS_layers/intersect_forest_edges_1K.shp")
intersect.forest.edge$length=st_length(intersect.forest.edge)
length=data.frame(cbind(ID=intersect.forest.edge$ID_2,length=intersect.forest.edge$length)) #transform the field “length” as numeric
str(length)
length$ID=as.factor(length$ID)
length$length=as.numeric(length$length)
length=aggregate(length~ID, sum, data=length) #calculate the sum of dens of each line within each cell of the grid
df=data.frame(grid.poly@data$ID) #extract the attribute table of grid.poly as a data frame
colnames(df)="ID" #provide a column name
df2=merge(df, length, by="ID", all.x=T, sort=T) #merge both data frame by field ID
df2[is.na(df2)]=0 #replace na by 0
df2$ID=as.factor(df2$ID) #change the ID column to factor
df2$dens=(df2$length*10000)/1000000 #calculte density per ha
grid.poly.forest.edge=grid.poly #create a new object that is the same as grid.poly
grid.poly.forest.edge@data$dens=df2$dens[match(df2$ID, grid.poly.forest.edge@data$ID)] #transfer the dens values of sf2 to grid.poly.rivers. The ones with no density values will remain at 0.
grid.poly.forest.edge@data$log_dens=log(grid.poly.forest.edge@data$dens+1)

r.forest.edge=raster(grid.poly.forest.edge)
res(r.forest.edge)=1000
values(r.forest.edge)=grid.poly.forest.edge@data$dens
r.forest.edge=mask(r.forest.edge, DSL93)
writeRaster(r.forest.edge, "env_variables_1K/dens_forest_edges.asc", format="ascii") #save raster

#log-transformed
r.forest.edge=raster(grid.poly.forest.edge)
res(r.forest.edge)=1000
values(r.forest.edge)=grid.poly.forest.edge@data$log_dens
r.forest.edge=mask(r.forest.edge, DSL93)
writeRaster(r.forest.edge, "env_variables_1K_log/dens_forest_edges.asc", format="ascii") #save raster



############
#FORESTS####
############

#forest=readOGR(dsn="GIS_layers/BDFORET_2-0__SHP_LAMB93_D079_2014-04-01/BDFORET/1_DONNEES_LIVRAISON/BDF_2-0_SHP_LAMB93_D079/FORMATION_VEGETALE.shp", layer="FORMATION_VEGETALE")
#shrub=readOGR(dsn="GIS_layers/bosquets_79_L93.shp", layer="bosquets_79_L93")
#proj4string(shrub)=crs(forest)
#forest=raster::bind(forest, shrub)
#raster::shapefile(forest, file="forets_bosquets_79_L93.shp")
#In QGIS, open the RPG79 layer, open query and enter the SQL code CODE_GROUP = ‘18’ OR CODE_GROUP = ‘19’, then open the grid_DeuxSevres_100m.shp, then menu vector -> Geoprocessing Tools -> intersect, and save the resulting layer as intersect_pastures_100m.shp

#percentage of forest####
intersect.forest=st_read("GIS_layers/intersect_forest_1K.shp") #read the layer with sf
intersect.forest$area2=st_area(intersect.forest) 
area=as.integer(as.character(intersect.forest$area))
area2=as.integer(as.character(intersect.forest$area2))
per.forest=data.frame(cbind(ID=intersect.forest$ID_2,per=round((area2*100)/area,digits=4))) #bind the ID column with another column where we calculate the percentage
per.forest$ID=as.factor(per.forest$ID)
per.forest$pour=as.numeric(as.character(per.forest$per))
per.forest=aggregate(per~ID, sum, data=per.forest)
df=data.frame(ID=as.factor(grid.poly@data$ID))
colnames(df)="ID"
df2=merge(df, per.forest, by="ID", all.x=T, sort=T)
df2[is.na(df2)]=0
df2$ID=as.factor(df2$ID)
grid.poly.forest=grid.poly
grid.poly.forest@data$per=df2$per[match(df2$ID, grid.poly.forest@data$ID)]
summary(grid.poly.forest@data$per) #check the data
grid.poly.forest@data$log_per=log(grid.poly.forest@data$per+1)

#transform the vector layer into a raster
grid.forest=grid
res(grid.forest)=1000
values(grid.forest)=grid.poly.forest@data$per
grid.forest=mask(grid.forest, DSL93)
writeRaster(grid.forest, "env_variables_1K/per_forest.asc", format="ascii") #save raster

#log-transformed
grid.forest=grid
res(grid.forest)=1000
values(grid.forest)=grid.poly.forest@data$log_per
grid.forest=mask(grid.forest, DSL93)
writeRaster(grid.forest, "env_variables_1K_log/per_forest.asc", format="ascii") #save raster

#grid.foret[is.na(grid.foret[])]=0
#grid.foret=mask(grid.foret,DSL93)
#writeRaster(grid.foret,"env_variables/per_forest.asc", format="ascii")

#distance to the nearest forest edge####
#Distance from centroids of the grid to the nearest forest edge
#first read the forest shape layer
forest=readOGR(dsn="GIS_layers/forets_bosquets_79_L93.shp", layer="forets_bosquets_79_L93")
forest=gBuffer(forest, width=0, byid=T) #we apply a 0m buffer aroung each geometry of forest to correct for geometry flaws
#forest=spTransform(forest, CRS("+init=epsg:4326"))

dist.forest=apply(gDistance(centroids, forest,byid=T),2,min)
dist.forest=data.frame(cbind(ID=rownames(data.frame(dist.forest)), dist=dist.forest))
dist.forest$ID=as.factor(dist.forest$ID)
dist.forest$dist=as.numeric(dist.forest$dist)
dist.forest$log_dist=log(dist.forest$dist+1)



grid.forest2=grid
values(grid.forest2)=dist.forest$dist
grid.forest2[is.na(grid.forest2[])]=0
#projection(grid.forest2)=CRS("+init=epsg:2154")
grid.forest2=mask(grid.forest2,DSL93)
writeRaster(grid.forest2,"env_variables_1K/dist_forest.asc", format="ascii")

#log-transformed
grid.forest2=grid
values(grid.forest2)=dist.forest$log_dist
grid.forest2[is.na(grid.forest2[])]=0
#projection(grid.forest2)=CRS("+init=epsg:2154")
grid.forest2=mask(grid.forest2,DSL93)
writeRaster(grid.forest2,"env_variables_1K_log/dist_forest.asc", format="ascii")


###########
#PONDS####
###########
pond=readOGR(dsn="GIS_layers/Donnees_paysages_79/MARE_ETANG_2002/MARE_79_2002.shp", layer="MARE_79_2002")
pond

dist.pond=apply(gDistance(centroids, pond,byid=T),2,min)
dist.pond=data.frame(cbind(ID=row.names(data.frame(dist.pond)), dist=dist.pond))
dist.pond$ID=as.factor(dist.pond$ID)
dist.pond$dist=as.numeric(dist.pond$dist)
dist.pond$log_dist=log(dist.pond$dist+1)

grid.pond=grid
values(grid.pond)=dist.pond$dist
grid.pond[is.na(grid.pond[])]=0
#projection(grid.pond)=CRS("+init=epsg:2154")
grid.pond=mask(grid.pond,DSL93)
writeRaster(grid.pond,"env_variables_1K/dist_pond.asc", format="ascii")

#log-transformed
grid.pond=grid
values(grid.pond)=dist.pond$log_dist
grid.pond[is.na(grid.pond[])]=0
#projection(grid.pond)=CRS("+init=epsg:2154")
grid.pond=mask(grid.pond,DSL93)
writeRaster(grid.pond,"env_variables_1K_log/dist_pond.asc", format="ascii")


#############
#HOUSINGS####
#############

#housing=readOGR(dsn="GIS_layers/BDTOPO_3-0_TOUSTHEMES_SHP_LAMB93_D079_2020-12-15/BDTOPO/1_DONNEES_LIVRAISON_2021-01-00019/BDT_3-0_SHP_LAMB93_D079-ED2020-12-15/BATI/BATIMENT.shp", layer="BATIMENT")
#housing=spTransform(housing, CRS("+init=epsg:4326"))

intersect.housing=st_read("GIS_layers/intersect_housing_1K.shp")
intersect.housing$area2=st_area(intersect.housing)
area=as.integer(as.character(intersect.housing$area))
area2=as.integer(as.character(intersect.housing$area2))
per.housing=data.frame(cbind(ID=intersect.housing$ID_2,per=round((area2*100)/area,digits=4))) #bind the ID column with another column where we calculate the percentage
per.housing$ID=as.factor(per.housing$ID)
per.housing$per=as.numeric(as.character(per.housing$per))
per.housing=aggregate(per~ID, sum, data=per.housing)
df=data.frame(ID=as.factor(grid.poly@data$ID))
colnames(df)="ID"
df2=merge(df, per.housing, by="ID", all.x=T, sort=T)
df2[is.na(df2)]=0
df2$ID=as.factor(df2$ID)
grid.poly.housing=grid.poly
grid.poly.housing@data$per=df2$per[match(df2$ID, grid.poly.housing@data$ID)]
summary(grid.poly.housing@data$per) #check the data
grid.poly.housing@data$log_per=log(grid.poly.housing@data$per+1)

#transform the vector layer into a raster
grid.housing=grid
res(grid.housing)=1000
values(grid.housing)=grid.poly.housing@data$per
grid.housing=mask(grid.housing,DSL93)
writeRaster(grid.housing,"env_variables_1K/per_housing.asc", format="ascii")

#log-transformed
grid.housing=grid
res(grid.housing)=1000
values(grid.housing)=grid.poly.housing@data$log_per
grid.housing=mask(grid.housing,DSL93)
writeRaster(grid.housing,"env_variables_1K_log/per_housing.asc", format="ascii")


######################
#VARIABLES SUMMARY####
######################

#prepare a table that summarizes the variables####
r.pasture=raster("env_variables_1K/per_pastures.asc") #open the raster file
r.forest=raster("env_variables_1K/per_forest.asc")
r.housing=raster("env_variables_1K/per_housing.asc")
r.road=raster("env_variables_1K/dens_roads.asc")
r.forest.edge=raster("env_variables_1K/dens_forest_edges.asc")
r.hedge=raster("env_variables_1K/dens_hedges.asc")
r.river=raster("env_variables_1K/dens_river.asc")
r.forest2=raster("env_variables_1K/dist_forest.asc")
r.pond=raster("env_variables_1K/dist_pond.asc")
#grid.poly=readOGR(dsn="GIS_layers/grid_DeuxSevres_1K.shp", layer="grid_DeuxSevres_1K") #read the grid as a vector
per.pasture=data.frame(na.omit(cbind(ID=grid.poly@data$ID,per.pasture=values(r.pasture)))) #bind the ID from grid.poly and the values of the raster layer, while discarding the NA values that correspond to the cells that are outside the limits of the department but within the extent
per.forest=data.frame(na.omit(cbind(ID=grid.poly@data$ID,per.forest=values(r.forest))))
per.housing=data.frame(na.omit(cbind(ID=grid.poly@data$ID,per.housing=values(r.housing))))
dens.road=data.frame(na.omit(cbind(ID=grid.poly@data$ID,dens.road=values(r.road))))
dens.forest.edge=data.frame(na.omit(cbind(ID=grid.poly@data$ID,dens.forest.edge=values(r.forest.edge))))
dens.hedge=data.frame(na.omit(cbind(ID=grid.poly@data$ID,dens.hedge=values(r.hedge))))
dens.river=data.frame(na.omit(cbind(ID=grid.poly@data$ID,dens.river=values(r.river))))
dist.forest=data.frame(na.omit(cbind(ID=grid.poly@data$ID,dist.forest=values(r.forest2))))
dist.pond=data.frame(na.omit(cbind(ID=grid.poly@data$ID,dist.pond=values(r.pond))))

env.var=data.frame(cbind(ID=per.pasture$ID,per.pasture=per.pasture$per.pasture,per.housing=per.housing$per.housing,per.forest=per.forest$per.forest, dens.road=dens.road$dens.road, dens.forest.edge=dens.forest.edge$dens.forest.edge, dens.hedge=dens.hedge$dens.hedge, dens.river=dens.river$dens.river, dist.forest=dist.forest$dist.forest, dist.pond=dist.pond$dist.pond)) #bind the data from raster into a single dataframe
str(env.var) #check how it looks like
summary(env.var)
write.table(env.var,file="env_variables_1K.txt", sep="\t") #save the table in the working directory
sum.env.var=data.frame(t(rbind(min=sapply(env.var[,-1], min), max=sapply(env.var[,-1], max)))) #summarize in a table the min and max values for all variables
write.table(sum.env.var,file="summary_env_variables_1K.txt", sep="\t") #save the table in the working directory

env.var=read.table("env_variables_1K.txt", h=T) #read the table
summary(env.var)

#Check graphically the distribution of the variables####
p1=ggplot(data=env.var, aes(x=per.pasture))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Pastures", x="Percentage", y="Frequency")+
  theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))
p2=ggplot(data=env.var, aes(x=per.forest))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Forests", x="Percentage", y="Frequency")+
  theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))
p3=ggplot(data=env.var, aes(x=per.housing))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Housings", x="Percentage", y="Frequency")+
  theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))
p4=ggplot(data=env.var, aes(x=dist.forest))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Forests", x="Distance (m)", y="Frequency")+
  theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))
p5=ggplot(data=env.var, aes(x=dist.pond))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Ponds", x="Distance (m)", y="Frequency")+
  theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))
p6=ggplot(data=env.var, aes(x=length.forest.edge))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Forest margins", x="Density (m/ha)", y="Frequency")+
  theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))
p7=ggplot(data=env.var, aes(x=length.hedge))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Hedges", x="Density (m/ha)", y="Frequency")+
  theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))
p8=ggplot(data=env.var, aes(x=length.river))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Rivers", x="Density (m/ha)", y="Frequency")+
  theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))
p9=ggplot(data=env.var, aes(x=length.road))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Roads", x="Density (m/ha)", y="Frequency")+theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))

jpeg(file="plot_hist_env_var_1K.jpg", width=17, height=17, res=300, units="cm")
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,top="Data not transformed (grid 100m)")
dev.off()

#Variable transformation####
#transform the variables into log
env.var.log=env.var
str(env.var.log)
cols=colnames(env.var.log[,c(2:10)])
env.var.log[cols]=lapply(env.var.log[cols]+1, log) #on transforme en log (en ajoutant +1 pour éviter log(0)) les variables de distance

#Check graphically the distribution of the variables####
p10=ggplot(data=env.var.log, aes(x=per.pasture))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Pastures", x="Log percentage", y="Frequency")+
  theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))
p11=ggplot(data=env.var.log, aes(x=per.forest))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Forests", x="Log percentage", y="Frequency")+
  theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))
p12=ggplot(data=env.var.log, aes(x=per.housing))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Housings", x="Log percentage", y="Frequency")+
  theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))
p13=ggplot(data=env.var.log, aes(x=dist.forest))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Forests", x="Log distance (m)", y="Frequency")+
  theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))
p14=ggplot(data=env.var.log, aes(x=dist.pond))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Ponds", x="Log distance (m)", y="Frequency")+
  theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))
p15=ggplot(data=env.var.log, aes(x=length.forest.edge))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Forest margins", x="Log density (m/ha)", y="Frequency")+
  theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))
p16=ggplot(data=env.var.log, aes(x=length.hedge))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Hedges", x="Log density (m/ha)", y="Frequency")+
  theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))
p17=ggplot(data=env.var.log, aes(x=length.river))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Rivers", x="Log density (m/ha)", y="Frequency")+
  theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))
p18=ggplot(data=env.var.log, aes(x=length.road))+
  geom_histogram(color="black", fill="grey", size=0.3)+
  labs(title="Roads", x="Log density (m/ha)", y="Frequency")+theme(plot.title=element_text(size=11),axis.text=element_text(size=7), axis.title=element_text(size=9))

jpeg(file="plot_hist_env_var_logtransf_1K.jpg", width=17, height=17, res=300, units="cm")
grid.arrange(p10,p11,p12,p13,p14,p15,p16,p17,p18,top="Data log-transformed (grid 100m)")
dev.off()

#Better distribution for distance variables

#########################
#FILTER SPECIES DATA####
#########################
data.sp=read.table("DATA_SPECIES_79/DONNEES_AMPHREP_79.txt", h=T) #read the table of data
#remove rare species as well as introduced species
data.sp=data.sp[data.sp $species!="Bombina_variegata",]
data.sp=data.sp[data.sp$species!="Coronella_austriaca",]
data.sp=data.sp[data.sp$species!="Emys_orbicularis",]
data.sp=data.sp[data.sp$species!="Hyla_meridionalis",]
data.sp=data.sp[data.sp$species!="Ichthyosaura_alpestris",]
data.sp=data.sp[data.sp$species!="Lissotriton_vulgaris",]
data.sp=data.sp[data.sp$species!="Pelophylax_lessonae",]
data.sp=data.sp[data.sp$species!="Trachemys_scripta",]
data.sp=data.sp[data.sp$species!="Xenopus_laevis",]
str(data.sp) #check what it looks like
data.sp$class=as.factor(data.sp$class)
data.sp$species=as.factor(data.sp$species)
data.sp$date=as.Date(data.sp$date)
length(unique(data.sp$species)) #check the number of species
#[1] 22
#transform the table into a SpatialPointsDataFrame
coordinates(data.sp)=~lon+lat
proj4string(data.sp)=CRS("+init=epsg:4326")
data.sp=spTransform(data.sp, CRS("+init=epsg:2154"))
raster::shapefile(data.sp, file="occurrence_22species_79.shp")
#List of species
#sort(unique(data.sp@data$species))
#[1] Alytes_obstetricans    Anguis_fragilis        Bufo_spinosus          Epidalea_calamita      Hierophis_viridiflavus
# [6] Hyla_arborea           Lacerta_bilineata      Lissotriton_helveticus Natrix_helvetica       Natrix_maura          
#[11] Pelodytes_punctatus    Pelophylax_esculentus  Pelophylax_ridibundus  Podarcis_muralis       Rana_dalmatina        
#[16] Rana_temporaria        Salamandra_salamandra  Triturus_blasii        Triturus_cristatus     Triturus_marmoratus   
#[21] Vipera_aspis           Zamenis_longissimus 

#Filter the occurrence data. We will filter the occurrence of each species so it is downgraded to a 100m x 100m grid of occurrences
#First, we create a raster out of the grid.poly object
grid=raster(grid.poly) #transfrom grid.poly into a raster
res(grid)=100 #define the resolution. 100 correspond to 100m x 100m grid

#then, for each species, we rasterize the occurrence points so it fits the grid. Then, we transform the raster into a SpatialPointsDataFrame. The points correspond to the centroids of the cells of grid.poly into which the species occurs.

ao=rasterize(data.sp[data.sp@data$species=="Alytes_obstetricans",], grid,field="species",function(x,...)length(unique(na.omit(x)))) #overall richness
ao2=rasterToPoints(ao, spatial=T)
ao2=spTransform(ao2, CRS("+init=epsg:2154"))
ao2@data$species="Alytes_obstetricans"
ao3=data.frame(cbind(species=ao2@data$species, lon=coordinates(ao2)[,1], lat=coordinates(ao2)[,2]))

af=rasterize(data.sp[data.sp@data$species=="Anguis_fragilis",], grid,field="species",function(x,...)length(unique(na.omit(x))))
af2=rasterToPoints(af, spatial=T)
af2=spTransform(af2, CRS("+init=epsg:2154"))
af2@data$species="Anguis_fragilis"
af3=data.frame(cbind(species=af2@data$species, lon=coordinates(af2)[,1], lat=coordinates(af2)[,2]))

bs=rasterize(data.sp[data.sp@data$species=="Bufo_spinosus",], grid,field="species",function(x,...)length(unique(na.omit(x))))
bs2=rasterToPoints(bs, spatial=T)
bs2=spTransform(bs2, CRS("+init=epsg:2154"))
bs2@data$species="Bufo_spinosus"
bs3=data.frame(cbind(species=bs2@data$species, lon=coordinates(bs2)[,1], lat=coordinates(bs2)[,2]))

ec=rasterize(data.sp[data.sp@data$species=="Epidalea_calamita",], grid,field="species",function(x,...)length(unique(na.omit(x))))
ec2=rasterToPoints(ec, spatial=T)
ec2=spTransform(ec2, CRS("+init=epsg:2154"))
ec2@data$species="Epidalea_calamita"
ec3=data.frame(cbind(species=ec2@data$species, lon=coordinates(ec2)[,1], lat=coordinates(ec2)[,2]))

hv=rasterize(data.sp[data.sp@data$species=="Hierophis_viridiflavus",], grid,field="species",function(x,...)length(unique(na.omit(x))))
hv2=rasterToPoints(hv, spatial=T)
hv2=spTransform(hv2, CRS("+init=epsg:2154"))
hv2@data$species="Hierophis_viridiflavus"
hv3=data.frame(cbind(species=hv2@data$species, lon=coordinates(hv2)[,1], lat=coordinates(hv2)[,2]))

ha=rasterize(data.sp[data.sp@data$species=="Hyla_arborea",], grid,field="species",function(x,...)length(unique(na.omit(x))))
ha2=rasterToPoints(ha, spatial=T)
ha2=spTransform(ha2, CRS("+init=epsg:2154"))
ha2@data$species="Hyla_arborea"
ha3=data.frame(cbind(species=ha2@data$species, lon=coordinates(ha2)[,1], lat=coordinates(ha2)[,2]))

lb=rasterize(data.sp[data.sp@data$species=="Lacerta_bilineata",], grid,field="species",function(x,...)length(unique(na.omit(x))))
lb2=rasterToPoints(lb, spatial=T)
lb2=spTransform(lb2, CRS("+init=epsg:2154"))
lb2@data$species="Lacerta_bilineata"
lb3=data.frame(cbind(species=lb2@data$species, lon=coordinates(lb2)[,1], lat=coordinates(lb2)[,2]))

lh=rasterize(data.sp[data.sp@data$species=="Lissotriton_helveticus",], grid,field="species",function(x,...)length(unique(na.omit(x))))
lh2=rasterToPoints(lh, spatial=T)
lh2=spTransform(lh2, CRS("+init=epsg:2154"))
lh2@data$species="Lissotriton_helveticus"
lh3=data.frame(cbind(species=lh2@data$species, lon=coordinates(lh2)[,1], lat=coordinates(lh2)[,2]))

nh=rasterize(data.sp[data.sp@data$species=="Natrix_helvetica",], grid,field="species",function(x,...)length(unique(na.omit(x))))
nh2=rasterToPoints(nh, spatial=T)
nh2=spTransform(nh2, CRS("+init=epsg:2154"))
nh2@data$species="Natrix_helvetica"
nh3=data.frame(cbind(species=nh2@data$species, lon=coordinates(nh2)[,1], lat=coordinates(nh2)[,2]))

nm=rasterize(data.sp[data.sp@data$species=="Natrix_maura",], grid,field="species",function(x,...)length(unique(na.omit(x))))
nm2=rasterToPoints(nm, spatial=T)
nm2=spTransform(nm2, CRS("+init=epsg:2154"))
nm2@data$species="Natrix_maura"
nm3=data.frame(cbind(species=nm2@data$species, lon=coordinates(nm2)[,1], lat=coordinates(nm2)[,2]))

pp=rasterize(data.sp[data.sp@data$species=="Pelodytes_punctatus",], grid,field="species",function(x,...)length(unique(na.omit(x))))
pp2=rasterToPoints(pp, spatial=T)
pp2=spTransform(pp2, CRS("+init=epsg:2154"))
pp2@data$species="Pelodytes_punctatus"
pp3=data.frame(cbind(species=pp2@data$species, lon=coordinates(pp2)[,1], lat=coordinates(pp2)[,2]))

pe=rasterize(data.sp[data.sp@data$species=="Pelophylax_esculentus",], grid,field="species",function(x,...)length(unique(na.omit(x))))
pe2=rasterToPoints(pe, spatial=T)
pe2=spTransform(pe2, CRS("+init=epsg:2154"))
pe2@data$species="Pelophylax_esculentus"
pe3=data.frame(cbind(species=pe2@data$species, lon=coordinates(pe2)[,1], lat=coordinates(pe2)[,2]))

pr=rasterize(data.sp[data.sp@data$species=="Pelophylax_ridibundus",], grid,field="species",function(x,...)length(unique(na.omit(x)))) 
pr2=rasterToPoints(pr, spatial=T)
pr2=spTransform(pr2, CRS("+init=epsg:2154"))
pr2@data$species="Pelophylax_ridibundus"
pr3=data.frame(cbind(species=pr2@data$species, lon=coordinates(pr2)[,1], lat=coordinates(pr2)[,2]))

pm=rasterize(data.sp[data.sp@data$species=="Podarcis_muralis",], grid,field="species",function(x,...)length(unique(na.omit(x)))) 
pm2=rasterToPoints(pm, spatial=T)
pm2=spTransform(pm2, CRS("+init=epsg:2154"))
pm2@data$species="Podarcis_muralis"
pm3=data.frame(cbind(species=pm2@data$species, lon=coordinates(pm2)[,1], lat=coordinates(pm2)[,2]))

rd=rasterize(data.sp[data.sp@data$species=="Rana_dalmatina",], grid,field="species",function(x,...)length(unique(na.omit(x)))) 
rd2=rasterToPoints(rd, spatial=T)
rd2=spTransform(rd2, CRS("+init=epsg:2154"))
rd2@data$species="Rana_dalmatina"
rd3=data.frame(cbind(species=rd2@data$species, lon=coordinates(rd2)[,1], lat=coordinates(rd2)[,2]))

rt=rasterize(data.sp[data.sp@data$species=="Rana_temporaria",], grid,field="species",function(x,...)length(unique(na.omit(x)))) 
rt2=rasterToPoints(rt, spatial=T)
rt2=spTransform(rt2, CRS("+init=epsg:2154"))
rt2@data$species="Rana_temporaria"
rt3=data.frame(cbind(species=rt2@data$species, lon=coordinates(rt2)[,1], lat=coordinates(rt2)[,2]))

ss=rasterize(data.sp[data.sp@data$species=="Salamandra_salamandra",], grid,field="species",function(x,...)length(unique(na.omit(x)))) 
ss2=rasterToPoints(ss, spatial=T)
ss2=spTransform(ss2, CRS("+init=epsg:2154"))
ss2@data$species="Salamandra_salamandra"
ss3=data.frame(cbind(species=ss2@data$species, lon=coordinates(ss2)[,1], lat=coordinates(ss2)[,2]))

tb=rasterize(data.sp[data.sp@data$species=="Triturus_blasii",], grid,field="species",function(x,...)length(unique(na.omit(x)))) 
tb2=rasterToPoints(tb, spatial=T)
tb2=spTransform(tb2, CRS("+init=epsg:2154"))
tb2@data$species="Triturus_blasii"
tb3=data.frame(cbind(species=tb2@data$species, lon=coordinates(tb2)[,1], lat=coordinates(tb2)[,2]))

tc=rasterize(data.sp[data.sp@data$species=="Triturus_cristatus",], grid,field="species",function(x,...)length(unique(na.omit(x)))) 
tc2=rasterToPoints(tc, spatial=T)
tc2=spTransform(tc2, CRS("+init=epsg:2154"))
tc2@data$species="Triturus_cristatus"
tc3=data.frame(cbind(species=tc2@data$species, lon=coordinates(tc2)[,1], lat=coordinates(tc2)[,2]))

tm=rasterize(data.sp[data.sp@data$species=="Triturus_marmoratus",], grid,field="species",function(x,...)length(unique(na.omit(x)))) 
tm2=rasterToPoints(tm, spatial=T)
tm2=spTransform(tm2, CRS("+init=epsg:2154"))
tm2@data$species="Triturus_marmoratus"
tm3=data.frame(cbind(species=tm2@data$species, lon=coordinates(tm2)[,1], lat=coordinates(tm2)[,2]))

va=rasterize(data.sp[data.sp@data$species=="Vipera_aspis",], grid,field="species",function(x,...)length(unique(na.omit(x))))
va2=rasterToPoints(va, spatial=T)
va2=spTransform(va2, CRS("+init=epsg:2154"))
va2@data$species="Vipera_aspis"
va3=data.frame(cbind(species=va2@data$species, lon=coordinates(va2)[,1], lat=coordinates(va2)[,2]))

zl=rasterize(data.sp[data.sp@data$species=="Zamenis_longissimus",], grid,field="species",function(x,...)length(unique(na.omit(x)))) 
zl2=rasterToPoints(zl, spatial=T)
zl2=spTransform(zl2, CRS("+init=epsg:2154"))
zl2@data$species="Zamenis_longissimus"
zl3=data.frame(cbind(species=zl2@data$species, lon=coordinates(zl2)[,1], lat=coordinates(zl2)[,2]))

grid.species=rbind(ao3, af3, bs3, ec3, hv3, ha3, lb3, lh3, nh3, nm3, pp3, pe3, pr3, pm3, rd3, rt3, ss3, tc3, tm3, tb3, va3, zl3) #build a data frame with all species and the coordinates of the centroids where they occur within the grid
grid.species$species=as.factor(grid.species$species)
grid.species$lon=as.numeric(grid.species$lon)
grid.species$lat=as.numeric(grid.species$lat)
str(grid.species)
#write.csv(grid.species, file="DATA_SPECIES_79/data_species_grid100m_79.csv", row.names=F) #save the table as a csv file
coordinates(grid.species)=~lon+lat
proj4string(grid.species)=CRS("+init=epsg:2154")
grid.species
#raster::shapefile(grid.species, file="DATA_SPECIES_79/grid_occurrence_species_100m_79.shp")

