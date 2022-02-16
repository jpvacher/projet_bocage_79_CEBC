######################################################
#Setting up rasters for Chloe
#script written by Jean-Pierre Vacher 4 January 2022
#updated 27 January 2022
######################################################
#Charger les packages nécessaires
x=c("here","dismo","ggmap", "rgdal", "rgeos", "maptools", "dplyr", "tidyr", "tmap", "raster","mapdata","sp","spdep","colorRamps","ggplot2","gridExtra")
lapply(x, library, character.only=TRUE)

DSL93=readOGR(dsn="contour_deux_sevres_L93.shp", layer="contour_deux_sevres_L93") #read the layer
grid=raster(DSL93) #create a raster out of the delineation of the Deux-Sèvres
res(grid)=10 #10 m resolution grid (cells 10 m x 10 m)

#woodlands
forest=readOGR(dsn="BDFORET_2-0__SHP_LAMB93_D079_2014-04-01/BDFORET/1_DONNEES_LIVRAISON/BDF_2-0_SHP_LAMB93_D079/FORMATION_VEGETALE.shp", layer="FORMATION_VEGETALE")
bosquets=readOGR(dsn="bosquets_79_L93.shp", layer="bosquets_79_L93")
proj4string(bosquets)=crs(forest)
forest=raster::bind(forest, bosquets)
r.forest=rasterize(forest, grid)
r.forest[is.na(r.forest[])]=0
r.forest=mask(r.forest,DSL93)
writeRaster(r.forest,"woodland_raster_10m.asc", format="ascii")

#woodland edges
forest2=readOGR(dsn="lisieres_foret_L93.shp", layer="lisieres_foret_L93")
bosquets2=readOGR(dsn="lisieres_bosquet_L93.shp", layer="lisieres_bosquet_L93")
proj4string(bosquets2)=crs(forest2)
forest.lis=raster::bind(forest2, bosquets2)
r.forest.lis=rasterize(forest.lis, grid)
r.forest.lis[is.na(r.forest.lis[])]=0
r.forest.lis=mask(r.forest.lis,DSL93)
writeRaster(r.forest.lis,"woodland_edges_raster_10m.asc", format="ascii")

#pastures
rpg79=readOGR(dsn="RPG_79.shp",layer="RPG_79")
prairie=rbind(rpg79[rpg79@data$CODE_GROUP=="18",], rpg79[rpg79@data$CODE_GROUP=="19",])
r.prairie=rasterize(prairie, grid)
r.prairie[is.na(r.prairie[])]=0
r.prairie=mask(r.prairie,DSL93)
writeRaster(r.prairie,"prairie_raster_10m.asc", format="ascii")


#crop fields
#rpg79=readOGR(dsn="RPG_79.shp",layer="RPG_79")
agri=rbind(rpg79[rpg79@data$CODE_GROUP=="1",], rpg79[rpg79@data$CODE_GROUP=="2",], rpg79[rpg79@data$CODE_GROUP=="3",], rpg79[rpg79@data$CODE_GROUP=="4",], rpg79[rpg79@data$CODE_GROUP=="5",], rpg79[rpg79@data$CODE_GROUP=="6",], rpg79[rpg79@data$CODE_GROUP=="7",])
r.agri=rasterize(agri, grid)
r.agri[is.na(r.agri[])]=0
r.agri=mask(r.agri,DSL93)
writeRaster(r.agri,"agri_raster_10m.asc", format="ascii")

#roads
road=readOGR(dsn="ROUTE500_2020_DeuxSevres.shp",layer="ROUTE500_2020_DeuxSevres") #build from the layer “TRONCON_ROUTE”of the ROUTE500 layer (version 2020)
r.roads=rasterize(road, grid)
r.roads[is.na(r.roads[])]=0
r.roads=mask(r.roads,DSL93)
writeRaster(r.roads,"roads_raster_20m.asc", format="ascii")

#rivers
rivers=readOGR(dsn="TOPAGE_2019_DeuxSevres.shp",layer="TOPAGE_2019_DeuxSevres") #build from the layer TOPAGE (version 2019)
r.rivers=rasterize(rivers, grid)
r.rivers[is.na(r.rivers[])]=0
r.rivers=mask(r.rivers,DSL93)
writeRaster(r.rivers,"rivers_raster_10m.asc", format="ascii")

#hedges
haie=readOGR(dsn="DNSB-HAIES_1-0__SHP_L93_D079_2020-06-24/DNSB-HAIES/1_DONNEES_LIVRAISON_2020-06-24/HAIE-LINEAIRE.shp", layer="HAIE-LINEAIRE")
r.haie=rasterize(haie, grid)
r.haie[is.na(r.haie[])]=0
r.haie=mask(r.haie,DSL93)
writeRaster(r.haie,"hedges_raster_20m.asc", format="ascii")

#housing
bati=readOGR(dsn="~/Documents/jpv/projet_CEBC/BDTOPO_3-0_TOUSTHEMES_SHP_LAMB93_D079_2020-12-15/BDTOPO/1_DONNEES_LIVRAISON_2021-01-00019/BDT_3-0_SHP_LAMB93_D079-ED2020-12-15/BATI/BATIMENT.shp", "BATIMENT")
r.bati=rasterize(bati, grid)
r.bati[is.na(r.bati[])]=0
r.bati=mask(r.bati,DSL93)
writeRaster(r.bati,"bati_raster_10m.asc", format="ascii")
