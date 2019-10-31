# extract NLCD 2011 raster values to points and MCP polygons

library(FedData)

mendota_poly <- polygon_from_extent(raster::extent(736176, 738675, 4062655, 4065235), proj4string='+proj=utm +datum=NAD83 +zone=10')



#+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs


#mendota FOWL 36
#738675.00 m E
#4065235.00 m N
#736176.00 m E
#4062655.00 m N

nw_poly <- polygon_from_extent(raster::extent(473169, 504015, 4582389, 5204089), proj4string='+proj=utm +datum=WGS84 +zone=10 +ellps=WGS84')

#OR/WA/CA
#473169.00 m E
#5204089.00 m N
#504015.00 m E
#4582389.00 m N

#make sure to change to zone 12
id_poly <- polygon_from_extent(raster::extent(353000, 354303, 4727172, 4728774), proj4string='+proj=utm +datum=WGS84 +zone=12 +ellps=WGS84')
plot(id_poly)

#ID
#354303.00 m E
#4728774.00 m N
#353000.00 m E
#4727172.00 m N


##download the NLCD raster and clip to the NW polygon (note: can only download 2011 with this function, not 2016)
#this takes a few mins because it's a large area

mend_nlcd <- get_nlcd(template = mendota_poly, label = 'mend_poly_NEW',  year = 2011, dataset = "landcover")
plot(mend_nlcd)

#this one takes some time
nw_nlcd <- get_nlcd(template = nw_poly, label = 'nw_poly2',  year = 2011, dataset = "landcover")
plot(nw_nlcd)

id_nlcd <- get_nlcd(template = id_poly, label = 'id_poly3',  year = 2011, dataset = "landcover")
plot(id_nlcd)

library(adehabitatHR)
library(sp)
#bring in points for FOWL 36

noha36.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/FOWL 36.csv"))

coordinates(noha36.df) <- c("x", "y")
proj4string(noha36.df) <- CRS("+proj=utm +zone=10 +datum=NAD83 +units=m")  ## for example
head(noha36.df)
plot(noha36.df)

newproj <- "+proj=utm +zone=10 +datum=NAD83 +units=m"
#nlcd_winter <- st_transform(Suisun_NLCD_winter, crs = the_crs)
nlcd_mend_new <- projectRaster(mend_nlcd, crs = newproj, method = "ngb" )

proj4string(noha36.df)
proj4string(nlcd_mend_new)

#bring in MCP for FOWL 36
library(adehabitatHR)
noha36.mcp <-mcp(noha36.df, percent = 100, unin = "m",unout = "km2")
proj4string(noha36.mcp)
plot(nlcd_mend_new)
plot(noha36.mcp, add = T)

library(sp)
clip36 <- mask(nlcd_mend_new, noha36.mcp)
plot(clip36)
str(clip36)
med.clip <- values(clip36)

final.clip <- med.clip[!is.na(med.clip)]

plot(noha36.mcp)
plot(clip36, add =T)

unique(values(mend_nlcd))

#then make sure they share the same projection
crs(mend_nlcd) <- crs(noha36.df)
mend_nlcd_proj <- spTransform(noha36.df, CRS(projection(mend_nlcd)))

proj4string(mend_nlcd_proj)
proj4string(noha36.df)

plot(noha36.mcp, add = TRUE)
plot(noha36.df, add = TRUE)

#extract raster values from points
#nlcd_extract_36 <- extract(mend_nlcd, noha36.df)
#head(mend_nlcd)
#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_36 <- SpatialPoints(noha36.df, proj4string = crs(nlcd_mend_new))
nlcd_extract_36 <- extract(mend_nlcd, nlcd_sp_36)
#nlcd_extract_36.2 <- extract(mend_nlcd, noha36.df)
head(nlcd_extract_36)
nlcd_extract_36[which(!is.na(nlcd_extract_36))]

plot(mend_nlcd, add = T)
plot(noha36.df, add = T)

# check to make sure they have the correct # of columns and rows
str(mend_nlcd)
str(noha36.df)

as.matrix(mend_nlcd)
# combine the raster cell probabilities with their coord pairs with landcover grid cells
final_36 <- cbind.data.frame(noha36.df, nlcd_extract_36)
head(final_36)

#above works, but returns all columns, including empty grid cells with NA and 0 values
# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_36 <- final_36[which(!is.na(final_36[,3])),]
head(final_36)

write.csv(final_36, file = "nlcd.36.csv")

#look at this thread for help
#https://gis.stackexchange.com/questions/60527/how-to-extract-values-from-rasters-at-location-of-points-in-r

#work through steps from example above with my data
library(FedData)
library(rgdal)
library(raster)
library(dismo)
library(sf)

#first read in the point table
#noha36.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/FOWL 36.csv"))

noha36 <- readOGR(dsn = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/FOWL_36_MCP/points.shp", layer = "points" )
plot(noha36)
proj4string(noha36)

#coordinates(noha36) <- c("long", "lat")
projection(noha36) <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")  

#this is the nlcd layer for the mendota area
mend_nlcd <- get_nlcd(template = mendota_poly, label = 'mend_poly',  year = 2011, dataset = "landcover")
plot(mend_nlcd)
proj4string(mend_nlcd)

#then make sure they share the same projection
mend_nlcd_noha36 <- spTransform(noha36, CRS(projection(mend_nlcd)))

proj4string(mend_nlcd_noha36)
proj4string(noha36)

#then set data frame and extract - NOT WORKING
data <- data.frame(coordinates(mend_nlcd_noha36),
                   extract(mend_nlcd, mend_nlcd_noha36))

data[which(!is.na(data[,3]))]

names(data) <- c("x", "y", "value")
data


## another example

PortionClassInPoly <- function(MySingleShape, MainRaster) {
  
  # Start with one shape from your list
  # Rasterize this shape according to the parameters of your background raster
  mini_rast <- rasterize(MySingleShape, MainRaster, fun='last')  
  # divide by itself - to make new set new raster values to 1 (NA values stays as NA). 
  mini_rast <- mini_rast/mini_rast
  # get the total number of cells for this shape (save for below)
  total_cells <- cellStats(mini_rast, 'sum')
  
  # multiply by original - to make a mask layer from your original raster
  my_cutout <- MainRaster*mini_rast
  # Count the number of cells for each discrete landcover class
  in_this_poly_unit <- freq(my_cutout)
  # Divide the cell count for each of these classes by the total 
  # number of cells in your current shape
  class_percents <- in_this_poly_unit[,2]/total_cells
  
  # make the function return a dataframe with
  # % landcover in each class for each shape
  output <- data.frame(myclass=in_this_poly_unit[,1], portion=class_percents)
  
  return(output)
}

PortionClassInPoly(MySingleShape=noha36.mcp, MainRaster=mend_nlcd)
  
  
