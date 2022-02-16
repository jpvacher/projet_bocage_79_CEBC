######################################################
#PREPARE ENVIRONMENT VARIABLES FOR SDM BOCAGE
#script written by Jean-Pierre Vacher 4 January 2022
#updated 16 February 2022
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

#########################################
#Work with the INPN grid - not necessary
#########################################

#Ajuster la grille INPN 1x1 km aux Deux-Sèvres sans tronquer les cellules périphériques 
#france=st_read("~/Documents/jpv/projet_CEBC/grid_france_L93_1x1.shp")
#ds=st_read("~/Documents/jpv/projet_CEBC/contour_deux_sevres_L93.shp")
#int=st_intersection(france, ds)
#grid=france[france$ID %in% int$ID,]
#grid=as_Spatial(grid)
#raster::shapefile(grid, file="grille_deux_sevres_L93_1x1_not_crop.shp")

#grid 100m x 100m
#grid.poly=readOGR(dsn="grid_DeuxSevres_100m_INPN.shp", layer="grid_DeuxSevres_100m_INPN")
#grid.poly@data$area=raster::area(grid.poly) #on ajoute un champ surface qui donne la surface de chaque polygone
#grid.poly@data$ID2=c(1:length(grid.poly)) #on ajoute un champ ID qui attribue un ID unique à chaque cellule
#grid.poly.red=grid.poly[grid.poly@data$area>10000,] #on exlut toutes les mailles tronquées qui se trouvent en marge
#raster::shapefile(grid.poly.red, file="grid_DeuxSevres_100m_trim.shp") #on sauvegarde la couche
#grid.poly.red=readOGR(dsn="grid_DeuxSevres_trim.shp", layer="grid_DeuxSevres_trim")
#raster::shapefile(grid.poly, "grid_DeuxSevres_100m_INPN.shp", overwrite=T) #on sauvegarde la couche avec l'information de surface des polygones et le nouvel ID

#Grille 250m x 250m
#grid.poly=readOGR(dsn="grid_DeuxSevres_250m_INPN.shp", layer="grid_DeuxSevres_250m_INPN")
#grid.poly@data$area=raster::area(grid.poly) #on ajoute un champ surface qui donne la surface de chaque polygone
#grid.poly@data$ID2=c(1:length(grid.poly)) #on ajoute un champ ID qui attribue un ID unique à chaque cellule
#raster::shapefile(grid.poly, "grid_DeuxSevres_250m_INPN.shp", overwrite=T) #on sauvegarde la couche avec l'information de surface des polygones et le nouvel ID
#grid.poly=readOGR(dsn="grid_DeuxSevres_250m_INPN.shp", layer="grid_DeuxSevres_250m_INPN")

#Grille 500m x 500m
#grid.poly=readOGR(dsn="grid_DeuxSevres_500m_INPN.shp", layer="grid_DeuxSevres_500m_INPN")
#grid.poly@data$area=raster::area(grid.poly) #on ajoute un champ surface qui donne la surface de chaque polygone


#################
#BUILD A GRID####
#################

grid=raster(DSL93) #create a raster out of the delineation of the Deux-Sèvres
res(grid)=100 #100 m resolution grid (cells 100 m x 100 m)

#transform the grid from raster to vector layer
#grid.poly=rasterToPolygons(grid)
#grid.poly@data$ID=c(1:length(grid.poly))
#grid.poly@data$layer=NULL
#grid.poly@data$area=raster::area(grid.poly)
#grid.poly
#raster::shapefile(grid.poly, file="grid_DeuxSevres_100m.shp")
grid.poly=readOGR(dsn="GIS_layers/grid_DeuxSevres_100m.shp", layer="grid_DeuxSevres_100m") #read the grid as a vector

#calculate centroids of the cells of the grid
#centroids=coordinates(grid.poly) #extract coordinates of the cells of the grid
#centroids=as.data.frame(centroids) #transform the list into a dataframe
#colnames(centroids)=c("x","y") #rename columns
#coordinates(centroids)=~x+y #define the coordinates of the layer
#proj4string(centroids)=CRS("+init=epsg:2154") #assign a projection system that corresponds to Lambert 93
#centroids #check what it looks like
#raster::shapefile(centroids, file="centroids_grid100m_L93.shp") #save the layer
centroids=readOGR(dsn="GIS_layers/centroids_grid100m_L93.shp", layer="centroids_grid100m_L93") #read the centroids layer

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
#intersect.river=raster::intersect(river, grid.poly)
#intersect.river=st_as_sf(intersect.river)
#or
#river=st_as_sf(river)
#grid.poly=st_as_sf(grid.poly)
#intersect.river=st_intersection(river, grid.poly)
#or
#river=st_read("TOPAGE_2020_DeuxSevres.shp")
#grid.poly=st_read("grid_DeuxSevres.shp")
#intersect.river=st_intersection(river, grid.poly)
#or process this operation in QGIS with the intersection function in the vector menu, it's more efficient than in R

intersect.river=st_read("GIS_layers/intersect_river_100m.shp") #read the layer with sf
intersect.river$length=st_length(intersect.river) #calculate the length of each line
length=data.frame(cbind(ID=intersect.river$ID,length=as.numeric(as.character(intersect.river$length)))) #transform the field “length” as numeric
length$ID=as.factor(length$ID)
length=aggregate(length~ID, sum, data=length) #calculate the sum of dens of each line within each cell of the grid
df=data.frame(grid.poly@data$ID) #extract the attribute table of grid.poly as a data frame
colnames(df)="ID" #provide a column name
df2=merge(df, length, by="ID", all.x=T, sort=T) #merge both data frame by field ID
df2[is.na(df2)]=0 #replace na by 0
df2$ID=as.factor(df2$ID) #change the ID column to factor
grid.poly.river=grid.poly #create a new object that is the same as grid.poly
grid.poly.river@data$length=df2$length[match(df2$ID, grid.poly.river@data$ID)] #transfer the dens values of sf2 to grid.poly.rivers. The ones with no density values will remain at 0.
grid.poly.river@data$log_length=log(grid.poly.river@data$length+1) #transfer the dens values of sf2 to grid.poly.rivers. The ones with no density values will remain at 0.

r.length=raster(grid.poly.river) #create a raster from the object grid.poly.rivers
res(r.length)=100 #resolution of the grid is cells of 100m x 100 m
values(r.length)= grid.poly.river@data$length #assign the values to the grids of the raster
r.length=mask(r.length, DSL93) #crop the raster to the delineation of the Deux-Sèvres department
writeRaster(r.length, "env_variables/length_river.asc", format="ascii") #save raster

#log transform
r.length=raster(grid.poly.river) #create a raster from the object grid.poly.rivers
res(r.length)=100 #resolution of the grid is cells of 100m x 100 m
values(r.length)= grid.poly.river@data$log_length #assign the values to the grids of the raster
r.length=mask(r.length, DSL93) #crop the raster to the delineation of the Deux-Sèvres department
writeRaster(r.length, "env_variables_log/length_river.asc", format="ascii") #save raster

length.river=grid.poly.river@data$length

#Distance to rivers####
#We fine 16 sub-layers of forest that encompass each zone of the previous sub-layers of centroids
river1=raster::bind(crop(river, DS1), crop(river, DS2))
river2=raster::bind(crop(river, DS1),crop(river, DS2), crop(river, DS3))
river3=raster::bind(crop(river, DS2),crop(river, DS3), crop(river, DS4))
river4=raster::bind(crop(river, DS3), crop(river, DS4), crop(river, DS5))
river5=raster::bind(crop(river, DS4), crop(river, DS5), crop(river, DS6))
river6=raster::bind(crop(river, DS5), crop(river, DS6), crop(river, DS7))
river7=raster::bind(crop(river, DS6), crop(river, DS7), crop(river, DS8))
river8=raster::bind(crop(river, DS7), crop(river, DS8), crop(river, DS9))
river9=raster::bind(crop(river, DS8), crop(river, DS9), crop(river, DS10))
river10=raster::bind(crop(river, DS9), crop(river, DS10), crop(river, DS11))
river11=raster::bind(crop(river, DS10), crop(river, DS11), crop(river, DS12))
river12=raster::bind(crop(river, DS11), crop(river, DS12), crop(river, DS13))
river13=raster::bind(crop(river, DS12), crop(river, DS13), crop(river, DS14))
river14=raster::bind(crop(river, DS13), crop(river, DS14), crop(river, DS15))
river15=raster::bind(crop(river, DS14), crop(river, DS15), crop(river, DS16))
river16=raster::bind(crop(river, DS15), crop(river, DS16))

#centroids=spTransform(centroids, CRS("+init=epsg:2154"))
#foret=spTransform(foret, CRS("+init=epsg:2154"))
dist.river1=apply(gDistance(centroids1, river1,byid=T),2,min)
dist.river2=apply(gDistance(centroids2, river2,byid=T),2,min)
dist.river3=apply(gDistance(centroids3, river3,byid=T),2,min)
dist.river4=apply(gDistance(centroids4, river4,byid=T),2,min)
dist.river5=apply(gDistance(centroids5, river5,byid=T),2,min)
dist.river6=apply(gDistance(centroids6, river6,byid=T),2,min)
dist.river7=apply(gDistance(centroids7, river7,byid=T),2,min)
dist.river8=apply(gDistance(centroids8, river8,byid=T),2,min)
dist.river9=apply(gDistance(centroids9, river9,byid=T),2,min)
dist.river10=apply(gDistance(centroids10, river10,byid=T),2,min)
dist.river11=apply(gDistance(centroids11, river11,byid=T),2,min)
dist.river12=apply(gDistance(centroids12, river12,byid=T),2,min)
dist.river13=apply(gDistance(centroids13, river13,byid=T),2,min)
dist.river14=apply(gDistance(centroids14, river14,byid=T),2,min)
dist.river15=apply(gDistance(centroids15, river15,byid=T),2,min)
dist.river16=apply(gDistance(centroids16, river16,byid=T),2,min)
dist.river=c(dist.river16,dist.river15,dist.river14,dist.river13,dist.river12,dist.river11,dist.river10,dist.river9,dist.river8,dist.river7,dist.river6,dist.river5,dist.river4,dist.river3,dist.river2,dist.river1)
dist.riverb=data.frame(ID=names(dist.river), dist=dist.river)
dist.river=dist.riverb[unique(dist.riverb$ID),]
dist.river=as.vector(dist.river$dist)

grid.river2=grid
values(grid.river2)=dist.river
grid.river2[is.na(grid.river2[])]=0
#projection(grid.river2)=CRS("+init=epsg:2154")
grid.river2=mask(grid.river2,DSL93)
writeRaster(grid.river2,"env_variables/dist_river.asc", format="ascii")


##########
#ROADS####
##########
road=readOGR(dsn="GIS_layers/ROUTE500_2020_DeuxSevres.shp",layer="ROUTE500_2020_DeuxSevres") #build from the layer “TRONCON_ROUTE”of the ROUTE500 layer (version 2020)
road

#Density of roads####
#intersect.road=raster::intersect(road, grid.poly)
#intersect.road=st_as_sf(intersect.road)
#or
#road=st_as_sf(road)
#grid.poly=st_as_sf(grid.poly)
#intersect.road=st_intersection(road, grid.poly)
#or
#road=st_read("ROUTE500_2020_DeuxSevres.shp")
#grid.poly=st_read("grid_DeuxSevres.shp")
#intersect.road=st_intersection(road, grid.poly)
#or process this operation in QGIS with the intersection function in the vector menu

intersect.road=st_read("GIS_layers/intersect_road_100m.shp") #read the layer with sf
intersect.road$length=st_length(intersect.road)
length=data.frame(cbind(ID=intersect.road$ID,length=as.numeric(as.character(intersect.road$length)))) #transform the field “length” as numeric
length$ID=as.factor(length$ID)
length=aggregate(length~ID, sum, data=length) #calculate the sum of dens of each line within each cell of the grid
df=data.frame(grid.poly@data$ID) #extract the attribute table of grid.poly as a data frame
colnames(df)="ID" #provide a column name
df2=merge(df, length, by="ID", all.x=T, sort=T) #merge both data frame by field ID
df2[is.na(df2)]=0 #replace na by 0
df2$ID=as.factor(df2$ID) #change the ID column to factor
grid.poly.road=grid.poly #create a new object that is the same as grid.poly
grid.poly.road@data$length=df2$length[match(df2$ID, grid.poly.road@data$ID)] #transfer the dens values of sf2 to grid.poly.rivers. The ones with no density values will remain at 0.
grid.poly.road@data$log_length=log(grid.poly.road@data$length+1)

r.road=raster(grid.poly.road)
res(r.road)=100
values(r.road)=grid.poly.road@data$length
r.road=mask(r.road, DSL93)
writeRaster(r.road, "env_variables/length_roads.asc", format="ascii") #save raster

#log-transformed
r.road=raster(grid.poly.road)
res(r.road)=100
values(r.road)=grid.poly.road@data$log_length
r.road=mask(r.road, DSL93)
writeRaster(r.road, "env_variables_log/length_roads.asc", format="ascii") #save raster


###########
#HEDGES####
###########

haie=readOGR(dsn="GIS_layers/DNSB-HAIES_1-0__SHP_L93_D079_2020-06-24/DNSB-HAIES/1_DONNEES_LIVRAISON_2020-06-24/HAIE-LINEAIRE.shp", layer="HAIE-LINEAIRE")
haie
#haie=spTransform(haie, CRS("+init=epsg:4326"))

#Density of hedges####
#intersect.haie=raster::intersect(haie, grid.poly)
#intersect.haie=st_as_sf(intersect.haie)
#or
#haie=st_as_sf(haie)
#grid.poly.red=st_as_sf(grid.poly)
#intersect.haie=st_intersection(haie, grid.poly)
#or
#haie=st_read("haie79_L93.shp")
#grid.poly=st_read("grid_DeuxSevres.shp")
#intersect.haie=st_intersection(haie, grid.poly.red)
#or process this operation in QGIS with the intersection function in the vector menu
intersect.haie=st_read("GIS_layers/intersect_haie_100m.shp") #read the layer with sf
intersect.haie$length=st_length(intersect.haie)
length=data.frame(cbind(ID=intersect.haie$ID,length=as.numeric(as.character(intersect.haie$length)))) #transform the field “length” as numeric
length$ID=as.factor(length$ID)
length=aggregate(length~ID, sum, data=length) #calculate the sum of dens of each line within each cell of the grid
df=data.frame(grid.poly@data$ID) #extract the attribute table of grid.poly as a data frame
colnames(df)="ID" #provide a column name
df2=merge(df, length, by="ID", all.x=T, sort=T) #merge both data frame by field ID
df2[is.na(df2)]=0 #replace na by 0
df2$ID=as.factor(df2$ID) #change the ID column to factor
grid.poly.hedge=grid.poly #create a new object that is the same as grid.poly
grid.poly.hedge@data$length=df2$length[match(df2$ID, grid.poly.hedge@data$ID)] #transfer the dens values of sf2 to grid.poly.rivers. The ones with no density values will remain at 0.
grid.poly.hedge@data$log_length=log(grid.poly.hedge@data$length+1)

r.hedge=raster(grid.poly.hedge)
res(r.hedge)=100
values(r.hedge)=grid.poly.hedge@data$length
r.hedge=mask(r.hedge, DSL93)
writeRaster(r.hedge, "env_variables/length_hedges.asc", format="ascii") #save raster

#log-transformed
r.hedge=raster(grid.poly.hedge)
res(r.hedge)=100
values(r.hedge)=grid.poly.hedge@data$log_length
r.hedge=mask(r.hedge, DSL93)
writeRaster(r.hedge, "env_variables_log/length_hedges.asc", format="ascii") #save raster


##############
#PASTURES####
##############

#rpg=readOGR(dsn="~/Documents/jpv/projet_CEBC/RPG_2-0_SHP_LAMB93_R75-2019/RPG/1_DONNEES_LIVRAISON_2019/RPG_2-0_SHP_LAMB93_R75-2019/PARCELLES_GRAPHIQUES_R75.shp", layer="PARCELLES_GRAPHIQUES_R75")
#rpg79=crop(rpg, DS) #crop to fit the study site (lighter file)
#raster::shapefile(rpg79, file="RPG_79.shp") #save the layer
rpg79=readOGR(dsn="GIS_layers/RPG_79.shp",layer="RPG_79")
#rpg79=spTransform(rpg79, CRS("+init=epsg:4326"))

#pasture=rbind(rpg79[rpg79@data$CODE_GROUP=="18",], rpg79[rpg79@data$CODE_GROUP=="19",])
#raster::shapefile(pasture, file="pasture79.shp")
#In QGIS, open the RPG79 layer, open query and enter the SQL code CODE_GROUP = ‘18’ OR CODE_GROUP = ‘19’, then open the grid_DeuxSevres_100m.shp, then menu vector -> Geoprocessing Tools -> intersect, and save the resulting layer as intersect_pastures_100m.shp

intersect.pasture=st_read("GIS_layers/intersect_pastures_100m.shp") #read the layer with sf
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
res(grid.pasture)=100
values(grid.pasture)=grid.poly.pasture@data$per
grid.pasture=mask(grid.pasture, DSL93)
writeRaster(grid.pasture,"env_variables/per_pastures.asc", format="ascii")

#log-transformed
grid.pasture=grid
res(grid.pasture)=100
values(grid.pasture)=grid.poly.pasture@data$log_per
grid.pasture=mask(grid.pasture, DSL93)
writeRaster(grid.pasture,"env_variables_log/per_pastures.asc", format="ascii")


################
#CROP FIELDS####
################
#rpg79=readOGR(dsn="RPG_79.shp",layer="RPG_79")
#agri=rbind(rpg79[rpg79@data$CODE_GROUP=="1",], rpg79[rpg79@data$CODE_GROUP=="2",], rpg79[rpg79@data$CODE_GROUP=="3",], rpg79[rpg79@data$CODE_GROUP=="4",], rpg79[rpg79@data$CODE_GROUP=="5",], rpg79[rpg79@data$CODE_GROUP=="6",], rpg79[rpg79@data$CODE_GROUP=="7",])
#raster::shapefile(agri, file="cultures79.shp")
#In QGIS, open the RPG79 layer, open query and enter the SQL code CODE_GROUP = ‘1’ OR CODE_GROUP = ‘2’ OR CODE_GROUP = ‘3’ OR CODE_GROUP = ‘4’ OR CODE_GROUP = ‘5’ OR CODE_GROUP = ‘6’ OR CODE_GROUP = ‘7’, then open the grid_DeuxSevres_100m.shp, then menu vector -> Geoprocessing Tools -> intersect, and save the resulting layer as intersect_pastures_100m.shp

intersect.agri=st_read("GIS_layers/intersect_agri_100m.shp")
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
res(grid.agri)=100
values(grid.agri)=grid.poly.agri@data$per
grid.agri=mask(grid.agri, DSL93)
writeRaster(grid.agri,"env_variables/per_agri.asc", format="ascii")

#log-transformed
grid.agri=grid
res(grid.agri)=100
values(grid.agri)=grid.poly.agri@data$log_per
grid.agri=mask(grid.agri, DSL93)
writeRaster(grid.agri,"env_variables_log/per_agri.asc", format="ascii")


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
intersect.forest.edge=st_read("GIS_layers/intersect_forest_edges_100m.shp")
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
grid.poly.forest.edge=grid.poly #create a new object that is the same as grid.poly
grid.poly.forest.edge@data$length=df2$length[match(df2$ID, grid.poly.forest.edge@data$ID)] #transfer the dens values of sf2 to grid.poly.rivers. The ones with no density values will remain at 0.
grid.poly.forest.edge@data$log_length=log(grid.poly.forest.edge@data$length+1)

r.forest.edge=raster(grid.poly.forest.edge)
res(r.forest.edge)=100
values(r.forest.edge)=grid.poly.forest.edge@data$length
r.forest.edge=mask(r.forest.edge, DSL93)
writeRaster(r.forest.edge, "env_variables/length_forest_edges.asc", format="ascii") #save raster

#log-transformed
r.forest.edge=raster(grid.poly.forest.edge)
res(r.forest.edge)=100
values(r.forest.edge)=grid.poly.forest.edge@data$log_length
r.forest.edge=mask(r.forest.edge, DSL93)
writeRaster(r.forest.edge, "env_variables_log/length_forest_edges.asc", format="ascii") #save raster



############
#FORESTS####
############

#forest=readOGR(dsn="BDFORET_2-0__SHP_LAMB93_D079_2014-04-01/BDFORET/1_DONNEES_LIVRAISON/BDF_2-0_SHP_LAMB93_D079/FORMATION_VEGETALE.shp", layer="FORMATION_VEGETALE")
#shrub=readOGR(dsn="bosquets_79_L93.shp", layer="bosquets_79_L93")
#proj4string(shrub)=crs(forest)
#forest=raster::bind(forest, shrub)
#raster::shapefile(forest, file="forets_bosquets_79_L93.shp")

#In QGIS, open the RPG79 layer, open query and enter the SQL code CODE_GROUP = ‘18’ OR CODE_GROUP = ‘19’, then open the grid_DeuxSevres_100m.shp, then menu vector -> Geoprocessing Tools -> intersect, and save the resulting layer as intersect_pastures_100m.shp

#percentage of forest####
intersect.forest=st_read("GIS_layers/intersect_forest_100m.shp") #read the layer with sf
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
res(grid.forest)=100
values(grid.forest)=grid.poly.forest@data$per
grid.forest=mask(grid.forest, DSL93)
writeRaster(grid.forest, "env_variables/per_forest.asc", format="ascii") #save raster

#log-transformed
grid.forest=grid
res(grid.forest)=100
values(grid.forest)=grid.poly.forest@data$log_per
grid.forest=mask(grid.forest, DSL93)
writeRaster(grid.forest, "env_variables_log/per_forest.asc", format="ascii") #save raster

#grid.foret[is.na(grid.foret[])]=0
#grid.foret=mask(grid.foret,DSL93)
#writeRaster(grid.foret,"env_variables/per_forest.asc", format="ascii")

#distance to the nearest forest edge####
#Distance from centroids of the grid to the nearest forest edge
#first read the forest shape layer
forest=readOGR(dsn="GIS_layers/forets_bosquets_79_L93.shp", layer="forets_bosquets_79_L93")
forest=gBuffer(forest, width=0, byid=T) #we apply a 0m buffer aroung each geometry of forest to correct for geometry flaws
#forest=spTransform(forest, CRS("+init=epsg:4326"))

#We define 16 sub-layers of forest that encompass each zone of the previous sub-layers of centroids
forest1=raster::bind(crop(forest, DS1), crop(forest, DS2))
forest2=raster::bind(crop(forest, DS1),crop(forest, DS2), crop(forest, DS3))
forest3=raster::bind(crop(forest, DS2),crop(forest, DS3), crop(forest, DS4))
forest4=raster::bind(crop(forest, DS3), crop(forest, DS4), crop(forest, DS5))
forest5=raster::bind(crop(forest, DS4), crop(forest, DS5), crop(forest, DS6))
forest6=raster::bind(crop(forest, DS5), crop(forest, DS6), crop(forest, DS7))
forest7=raster::bind(crop(forest, DS6), crop(forest, DS7), crop(forest, DS8))
forest8=raster::bind(crop(forest, DS7), crop(forest, DS8), crop(forest, DS9))
forest9=raster::bind(crop(forest, DS8), crop(forest, DS9), crop(forest, DS10))
forest10=raster::bind(crop(forest, DS9), crop(forest, DS10), crop(forest, DS11))
forest11=raster::bind(crop(forest, DS10), crop(forest, DS11), crop(forest, DS12))
forest12=raster::bind(crop(forest, DS11), crop(forest, DS12), crop(forest, DS13))
forest13=raster::bind(crop(forest, DS12), crop(forest, DS13), crop(forest, DS14))
forest14=raster::bind(crop(forest, DS13), crop(forest, DS14), crop(forest, DS15))
forest15=raster::bind(crop(forest, DS14), crop(forest, DS15), crop(forest, DS16))
forest16=raster::bind(crop(forest, DS15), crop(forest, DS16))

#centroids=spTransform(centroids, CRS("+init=epsg:2154"))
#foret=spTransform(foret, CRS("+init=epsg:2154"))
dist.forest1=apply(gDistance(centroids1, forest1,byid=T),2,min)
dist.forest2=apply(gDistance(centroids2, forest2,byid=T),2,min)
dist.forest3=apply(gDistance(centroids3, forest3,byid=T),2,min)
dist.forest4=apply(gDistance(centroids4, forest4,byid=T),2,min)
dist.forest5=apply(gDistance(centroids5, forest5,byid=T),2,min)
dist.forest6=apply(gDistance(centroids6, forest6,byid=T),2,min)
dist.forest7=apply(gDistance(centroids7, forest7,byid=T),2,min)
dist.forest8=apply(gDistance(centroids8, forest8,byid=T),2,min)
dist.forest9=apply(gDistance(centroids9, forest9,byid=T),2,min)
dist.forest10=apply(gDistance(centroids10, forest10,byid=T),2,min)
dist.forest11=apply(gDistance(centroids11, forest11,byid=T),2,min)
dist.forest12=apply(gDistance(centroids12, forest12,byid=T),2,min)
dist.forest13=apply(gDistance(centroids13, forest13,byid=T),2,min)
dist.forest14=apply(gDistance(centroids14, forest14,byid=T),2,min)
dist.forest15=apply(gDistance(centroids15, forest15,byid=T),2,min)
dist.forest16=apply(gDistance(centroids16, forest16,byid=T),2,min)
dist.forest=c(dist.forest16,dist.forest15,dist.forest14,dist.forest13,dist.forest12,dist.forest11,dist.forest10,dist.forest9,dist.forest8,dist.forest7,dist.forest6,dist.forest5,dist.forest4,dist.forest3,dist.forest2,dist.forest1)
dist.forestb=data.frame(ID=names(dist.forest), dist=dist.forest)
dist.forest=dist.forestb[unique(dist.forestb$ID),]
dist.forest=as.vector(dist.forest$dist)
dist.forest.log=log(dist.forest+1)

grid.forest2=grid
values(grid.forest2)=dist.forest
grid.forest2[is.na(grid.forest2[])]=0
#projection(grid.forest2)=CRS("+init=epsg:2154")
grid.forest2=mask(grid.forest2,DSL93)
writeRaster(grid.forest2,"env_variables/dist_forest.asc", format="ascii")

#log-transformed
grid.forest2=grid
values(grid.forest2)=dist.forest.log
grid.forest2[is.na(grid.forest2[])]=0
#projection(grid.forest2)=CRS("+init=epsg:2154")
grid.forest2=mask(grid.forest2,DSL93)
writeRaster(grid.forest2,"env_variables_log/dist_forest.asc", format="ascii")


###########
#PONDS####
###########
pond=readOGR(dsn="GIS_layers/Donnees_paysages_79/MARE_ETANG_2002/MARE_79_2002.shp", layer="MARE_79_2002")
pond

#Define 4 sub-layers of forest that encompass each zone of the previous sub-layers of centroids
pond1=raster::bind(crop(pond, DS1), crop(pond, DS2))
pond2=raster::bind(crop(pond, DS1),crop(pond, DS2), crop(pond, DS3))
pond3=raster::bind(crop(pond, DS2),crop(pond, DS3), crop(pond, DS4))
pond4=raster::bind(crop(pond, DS3), crop(pond, DS4), crop(pond, DS5))
pond5=raster::bind(crop(pond, DS4), crop(pond, DS5), crop(pond, DS6))
pond6=raster::bind(crop(pond, DS5), crop(pond, DS6), crop(pond, DS7))
pond7=raster::bind(crop(pond, DS6), crop(pond, DS7), crop(pond, DS8))
pond8=raster::bind(crop(pond, DS7), crop(pond, DS8), crop(pond, DS9))
pond9=raster::bind(crop(pond, DS8), crop(pond, DS9), crop(pond, DS10))
pond10=raster::bind(crop(pond, DS9), crop(pond, DS10), crop(pond, DS11))
pond11=raster::bind(crop(pond, DS10), crop(pond, DS11), crop(pond, DS12))
pond12=raster::bind(crop(pond, DS11), crop(pond, DS12), crop(pond, DS13))
pond13=raster::bind(crop(pond, DS12), crop(pond, DS13), crop(pond, DS14))
pond14=raster::bind(crop(pond, DS13), crop(pond, DS14), crop(pond, DS15))
pond15=raster::bind(crop(pond, DS14), crop(pond, DS15), crop(pond, DS16))
pond16=raster::bind(crop(pond, DS15), crop(pond, DS16))

#distance to ponds####
dist.pond1=apply(gDistance(centroids1, pond1,byid=T),2,min)
dist.pond2=apply(gDistance(centroids2, pond2,byid=T),2,min)
dist.pond3=apply(gDistance(centroids3, pond3,byid=T),2,min)
dist.pond4=apply(gDistance(centroids4, pond4,byid=T),2,min)
dist.pond5=apply(gDistance(centroids5, pond5,byid=T),2,min)
dist.pond6=apply(gDistance(centroids6, pond6,byid=T),2,min)
dist.pond7=apply(gDistance(centroids7, pond7,byid=T),2,min)
dist.pond8=apply(gDistance(centroids8, pond8,byid=T),2,min)
dist.pond9=apply(gDistance(centroids9, pond9,byid=T),2,min)
dist.pond10=apply(gDistance(centroids10, pond10,byid=T),2,min)
dist.pond11=apply(gDistance(centroids11, pond11,byid=T),2,min)
dist.pond12=apply(gDistance(centroids12, pond12,byid=T),2,min)
dist.pond13=apply(gDistance(centroids13, pond13,byid=T),2,min)
dist.pond14=apply(gDistance(centroids14, pond14,byid=T),2,min)
dist.pond15=apply(gDistance(centroids15, pond15,byid=T),2,min)
dist.pond16=apply(gDistance(centroids16, pond16,byid=T),2,min)
dist.pond=c(dist.pond16,dist.pond15,dist.pond14,dist.pond13,dist.pond12,dist.pond11,dist.pond10,dist.pond9,dist.pond8,dist.pond7,dist.pond6,dist.pond5,dist.pond4,dist.pond3,dist.pond2,dist.pond1)
dist.pondb=data.frame(ID=names(dist.pond), dist=dist.pond)
dist.pond=dist.pondb[unique(dist.pondb$ID),]
dist.pond=as.vector(dist.pond$dist)
dist.pond.log=log(dist.pond+1)

grid.pond=grid
values(grid.pond)=dist.pond
grid.pond[is.na(grid.pond[])]=0
#projection(grid.pond)=CRS("+init=epsg:2154")
grid.pond=mask(grid.pond,DSL93)
writeRaster(grid.pond,"env_variables/dist_pond.asc", format="ascii")

#log-transformed
grid.pond=grid
values(grid.pond)=dist.pond.log
grid.pond[is.na(grid.pond[])]=0
#projection(grid.pond)=CRS("+init=epsg:2154")
grid.pond=mask(grid.pond,DSL93)
writeRaster(grid.pond,"env_variables_log/dist_pond.asc", format="ascii")


#############
#HOUSINGS####
#############

housing=readOGR(dsn="GIS_layers/BDTOPO_3-0_TOUSTHEMES_SHP_LAMB93_D079_2020-12-15/BDTOPO/1_DONNEES_LIVRAISON_2021-01-00019/BDT_3-0_SHP_LAMB93_D079-ED2020-12-15/BATI/BATIMENT.shp", layer="BATIMENT")
#housing=spTransform(housing, CRS("+init=epsg:4326"))

intersect.housing=st_read("GIS_layers/intersect_housing_100m.shp")
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
res(grid.housing)=100
values(grid.housing)=grid.poly.housing@data$per
grid.housing=mask(grid.housing,DSL93)
writeRaster(grid.housing,"env_variables/per_housing.asc", format="ascii")

#log-transformed
grid.housing=grid
res(grid.housing)=100
values(grid.housing)=grid.poly.housing@data$log_per
grid.housing=mask(grid.housing,DSL93)
writeRaster(grid.housing,"env_variables_log/per_housing.asc", format="ascii")


######################
#VARIABLES SUMMARY####
######################

#prepare a table that summarizes the variables####
r.pasture=raster("env_variables/per_pastures.asc") #open the raster file
#r.agri=raster("env_variables/per_agri.asc")
r.forest=raster("env_variables/per_forest.asc")
r.housing=raster("env_variables/per_housing.asc")
r.road=raster("env_variables/length_roads.asc")
r.forest.edge=raster("env_variables/length_forest_edges.asc")
r.hedge=raster("env_variables/length_hedges.asc")
r.river=raster("env_variables/length_river.asc")
r.forest2=raster("env_variables/dist_forest.asc")
r.pond=raster("env_variables/dist_pond.asc")
per.pasture=data.frame(na.omit(cbind(ID=grid.poly@data$ID,per.pasture=values(r.pasture)))) #bind the ID from grid.poly and the values of the raster layer, while discarding the NA values that correspond to the cells that are outside the limits of the department but within the extent
#per.agri=data.frame(na.omit(cbind(ID=grid.poly@data$ID,per.agri=values(r.agri))))
per.forest=data.frame(na.omit(cbind(ID=grid.poly@data$ID,per.forest=values(r.forest))))
per.housing=data.frame(na.omit(cbind(ID=grid.poly@data$ID,per.housing=values(r.housing))))
length.road=data.frame(na.omit(cbind(ID=grid.poly@data$ID,length.road=values(r.road))))
length.forest.edge=data.frame(na.omit(cbind(ID=grid.poly@data$ID,length.forest.edge=values(r.forest.edge))))
length.hedge=data.frame(na.omit(cbind(ID=grid.poly@data$ID,length.hedge=values(r.hedge))))
length.river=data.frame(na.omit(cbind(ID=grid.poly@data$ID,length.river=values(r.river))))
dist.forest=data.frame(na.omit(cbind(ID=grid.poly@data$ID,dist.forest=values(r.forest2))))
dist.pond=data.frame(na.omit(cbind(ID=grid.poly@data$ID,dist.pond=values(r.pond))))

env.var=data.frame(cbind(ID=per.pasture$ID,per.pasture=per.pasture$per.pasture,per.housing=per.housing$per.housing,per.forest=per.forest$per.forest, length.road=length.road$length.road, length.forest.edge=length.forest.edge$length.forest.edge, length.hedge=length.hedge$length.hedge, length.river=length.river$length.river, dist.forest=dist.forest$dist.forest, dist.pond=dist.pond$dist.pond)) #bind the data from raster into a single dataframe
str(env.var) #check how it looks like
summary(env.var)
write.table(env.var,file="env_variables.txt", sep="\t") #save the table in the working directory
sum.env.var=data.frame(t(rbind(min=sapply(env.var[,-1], min), max=sapply(env.var[,-1], max)))) #summarize in a table the min and max values for all variables
write.table(sum.env.var,file="summary_env_variables.txt", sep="\t") #save the table in the working directory

env.var=read.table("env_variables.txt", h=T) #read the table
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

jpeg(file="plot_hist_env_var_100m.jpg", width=17, height=17, res=300, units="cm")
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,top="Data not transformed (grid 100m)")

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

jpeg(file="plot_hist_env_var_logtransf_100m.jpg", width=17, height=17, res=300, units="cm")
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

