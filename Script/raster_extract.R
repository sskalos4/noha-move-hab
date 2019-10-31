# extract NLCD 2011 raster values to points and MCP polygons

library(FedData)

mendota_poly <- polygon_from_extent(raster::extent(736176, 738675, 4062655, 4065235), proj4string='+proj=utm +datum=NAD83 +zone=10')

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

######## make sure to change to zone 12 ########
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
nw_nlcd <- get_nlcd(template = nw_poly, label = 'nw_poly3',  year = 2011, dataset = "landcover")
plot(nw_nlcd)

id_nlcd <- get_nlcd(template = id_poly, label = 'id_poly4',  year = 2011, dataset = "landcover")
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

plot(nlcd_mend_new)
plot(noha36.df, add=T)

#bring in MCP for FOWL 36
library(adehabitatHR)
noha36.mcp <-mcp(noha36.df, percent = 100, unin = "m",unout = "km2")
proj4string(noha36.mcp)

plot(nlcd_mend_new)
plot(noha36.mcp, add = T)
points(noha36.df, add=T)

#extract raster values from points
#nlcd_extract_36 <- extract(mend_nlcd, noha36.df)
#head(mend_nlcd)
#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_36 <- SpatialPoints(noha36.df, proj4string = crs(nlcd_mend_new))
nlcd_extract_36 <- extract(nlcd_mend_new, nlcd_sp_36)
#nlcd_extract_36.2 <- extract(mend_nlcd, noha36.df)
head(nlcd_extract_36)
nlcd_extract_36[which(!is.na(nlcd_extract_36))]

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final_36 <- cbind.data.frame(noha36.df, nlcd_extract_36)
head(final_36)

#above works, but returns all columns, including empty grid cells with NA and 0 values

# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_36 <- final_36[which(!is.na(final_36[,3])),]
head(final_36)

write.csv(final_36, file = "nlcd.36.points.csv")

## raster values to mcp
clip36 <- mask(nlcd_mend_new, noha36.mcp)
plot(clip36)
str(clip36)
med.clip <- values(clip36)
final.clip <- med.clip[!is.na(med.clip)]
final.clip

write.csv(final.clip, file = "nlcd.36.mcp.csv")

plot(noha36.mcp)
plot(clip36, add =T)
plot(noha36.df, add=T, pch = 16)

