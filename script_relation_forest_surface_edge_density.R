######################################################
#script written by Jean-Pierre Vacher 4 January 2022
#updated 8 February 2022
######################################################

#charge the packages
x=c("here","sf","ggplot2", "gridExtra") #create a list with names of packages
lapply(x, library, character.only=TRUE) #loop that read all the packages from the list

########################
#ENVIRONMENT VARIABLES
########################
#Prepare the hexagon grids
#Deux-Sevres department
#DSL93=st_read("GIS_layers/contour_deux_sevres_L93.shp")
#DSL93=readOGR(dsn="GIS_layers/contour_deux_sevres_L93.shp", layer="contour_deux_sevres_L93") #read the layer

#Build hexagons
#hex.grid100=st_make_grid(DSL93, cellsize=100, square=F) #create a hexagon polygon layer of “cellsize” m in size
#hex.grid100=as_Spatial(hex.grid100) #transform as spatial polygon
#id=sapply(slot(hex.grid100, "polygons"), function(x)slot(x,"ID")) #create a polygon slot that will contain the IDs of the hexagons
#hex.grid100.df=data.frame(ID=1:length(hex.grid100), row.names=id) #create a dataframe with the ID of the polygons
#hex.grid100=SpatialPolygonsDataFrame(hex.grid100, hex.grid100.df) #transform the hexagon layer into a SpatialPolygonsDataFrame, with the polygon IDs as an attribute table
#hex.grid100=crop(hex.grid100, DSL93) #crop to adjust to the department limits
#hex.grid100 #check what it looks like
#plot(hex.grid100) #check graphicaly
#raster::shapefile(hex.grid100, file="hexagon_grid_100m.shp")

#hex.grid500=st_make_grid(DSL93, cellsize=500, square=F) #create a hexagon polygon layer of “cellsize” m in size
#hex.grid500=as_Spatial(hex.grid500) #transform as spatial polygon
#id=sapply(slot(hex.grid500, "polygons"), function(x)slot(x,"ID")) #create a polygon slot that will contain the IDs of the hexagons
#hex.grid500.df=data.frame(ID=1:length(hex.grid500), row.names=id) #create a dataframe with the ID of the polygons
#hex.grid500=SpatialPolygonsDataFrame(hex.grid500, hex.grid500.df) #transform the hexagon layer into a SpatialPolygonsDataFrame, with the polygon IDs as an attribute table
#hex.grid500=crop(hex.grid500, DSL93) #crop to adjust to the department limits
#hex.grid500 #check what it looks like
#plot(hex.grid500) #check graphicaly
#raster::shapefile(hex.grid500, file="hexagon_grid_500m.shp")

#hex.grid1K=st_make_grid(DSL93, cellsize=1000, square=F) #create a hexagon polygon layer of “cellsize” m in size
#hex.grid1K=as_Spatial(hex.grid1K) #transform as spatial polygon
#id=sapply(slot(hex.grid1K, "polygons"), function(x)slot(x,"ID")) #create a polygon slot that will contain the IDs of the hexagons
#hex.grid1K.df=data.frame(ID=1:length(hex.grid1K), row.names=id) #create a dataframe with the ID of the polygons
#hex.grid1K=SpatialPolygonsDataFrame(hex.grid1K, hex.grid1K.df) #transform the hexagon layer into a SpatialPolygonsDataFrame, with the polygon IDs as an attribute table
#hex.grid1K=crop(hex.grid1K, DSL93) #crop to adjust to the department limits
#hex.grid1K #check what it looks like
#plot(hex.grid1K) #check graphicaly
#raster::shapefile(hex.grid1K, file="hexagon_grid_1K.shp")

#read the hexagon layers
hex.grid100=st_read("GIS_layers/hexagon_grid_100m.shp")
hex.grid500=st_read("GIS_layers/hexagon_grid_500m.shp")
hex.grid1K=st_read("GIS_layers/hexagon_grid_1K.shp")

#Intersect files are produced in QGIS with the intersect function in the vector menu.

#100 m
intersect.forest100=st_read("GIS_layers/intersect_hexagon_forest_100m.shp")
intersect.forest.edge100=st_read("GIS_layers/intersect_hexagon_forest_edges_100m.shp")
intersect.forest.edge100$length=st_length(intersect.forest.edge100)
df100=data.frame(ID_2=unique(intersect.forest.edge100$ID_2))
length.forest.edge=aggregate(length~ID_2, sum, data= intersect.forest.edge100)
df2=merge(df100, length.forest.edge, by="ID_2", all.x=T, sort=T)
df2[is.na(df2)]=0
df2$ID_2=as.factor(df2$ID_2)
df2$n.poly=as.integer(table(intersect.forest.edge100$ID_2))
df2$length=as.numeric(as.character(df2$length))

p1=ggplot(data=df2, aes(x=length, y=n.poly))+
	geom_point(size=0.5)+
	geom_smooth()+
	labs(x="Length of edges (m)", y="Number of polygons", title="Scale = 100 m")

#500 m
intersect.forest500=st_read("GIS_layers/intersect_hexagon_forest_500m.shp")
intersect.forest.edge500=st_read("GIS_layers/intersect_hexagon_forest_edges_500m.shp")
intersect.forest.edge500$length=st_length(intersect.forest.edge500)
df500=data.frame(ID_2=unique(intersect.forest.edge500$ID_2))
length.forest.edge=aggregate(length~ID_2, sum, data= intersect.forest.edge500)
df3=merge(df500, length.forest.edge, by="ID_2", all.x=T, sort=T)
df3[is.na(df3)]=0
df3$ID_2=as.factor(df3$ID_2)
df3$n.poly=as.integer(table(intersect.forest.edge500$ID_2))
df3$length=as.numeric(as.character(df3$length))

p2=ggplot(data=df3, aes(x=length, y=n.poly))+
	geom_point(size=0.5)+
	geom_smooth()+
	labs(x="Length of edges (m)", y="Number of polygons", title="Scale = 500 m")

#1000 m
intersect.forest1K=st_read("GIS_layers/intersect_hexagon_forest_1K.shp")
intersect.forest.edge1K=st_read("GIS_layers/intersect_hexagon_forest_edges_1K.shp")
intersect.forest.edge1K$length=st_length(intersect.forest.edge1K)
df1K=data.frame(ID_2=unique(intersect.forest.edge1K$ID_2))
length.forest.edge=aggregate(length~ID_2, sum, data= intersect.forest.edge1K)
df4=merge(df1K, length.forest.edge, by="ID_2", all.x=T, sort=T)
df4[is.na(df4)]=0
df4$ID_2=as.factor(df4$ID_2)
df4$n.poly=as.integer(table(intersect.forest.edge1K$ID_2))
df4$length=as.numeric(as.character(df4$length))

p3=ggplot(data=df4, aes(x=length, y=n.poly))+
	geom_point(size=0.5)+
	geom_smooth()+
	labs(x="Length of edges (m)", y="Number of polygons", title="Scale = 1000 m")

jpeg(file="plot_nPolygons_length_forest.jpg", width=17, height=17, res=300, units="cm")
grid.arrange(p1, p2, p3)
dev.off()