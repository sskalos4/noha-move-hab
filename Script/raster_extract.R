# extract NLCD 2011 raster values to points and MCP polygons

library(FedData)

#37 2018 AK
fowl37.18_poly <- polygon_from_extent(raster::extent(498251, 500034, 7664495, 7665673), proj4string='+proj=utm +datum=NAD83 +zone=4')

#500034.00 m E
#7665673.00 m N
#498251.00 m E
#7664495.00 m N

#37 2019 AK
fowl37.19_poly <- polygon_from_extent(raster::extent(551052,560122, 7662701, 7670275), proj4string='+proj=utm +datum=NAD83 +zone=4')

#560122.00 m E
#7670275.00 m N
#551052.00 m E
#7662701.00 m N

yakima_poly <- polygon_from_extent(raster::extent(687180, 689159, 5144476, 5146132), proj4string='+proj=utm +datum=NAD83 +zone=10')

#689159.00 m E
#5146132.00 m N
#687180.00 m E
#5144476.00 m N

mal_poly <- polygon_from_extent(raster::extent(339280, 353173, 4740418, 4775425), proj4string='+proj=utm +datum=NAD83 +zone=11')

#353173.00 m E
#4775425.00 m N
#339280.00 m E
#4740418.00 m N


hart_poly <- polygon_from_extent(raster::extent(260639, 283341, 4687467, 4699922), proj4string='+proj=utm +datum=NAD83 +zone=11')

#283341.00 m E
#4699922.00 m N
#260639.00 m E
#4687467.00 m N

mendota_poly <- polygon_from_extent(raster::extent(736176, 738675, 4062655, 4065235), proj4string='+proj=utm +datum=NAD83 +zone=10')

#mendota FOWL 36
#738675.00 m E
#4065235.00 m N
#736176.00 m E
#4062655.00 m N


#nw_poly <- polygon_from_extent(raster::extent(473169, 504015, 4582389, 5204089), #proj4string='+proj=utm +datum=NAD83 +zone=10')

#OR/WA/CA
#473169.00 m E
#5204089.00 m N
#504015.00 m E
#4582389.00 m N

######## make sure to change to zone 12 ########
id_poly <- polygon_from_extent(raster::extent(353000, 354303, 4727172, 4728774), proj4string='+proj=utm +datum=NAD83 +zone=12')
plot(id_poly)

#ID
#354303.00 m E
#4728774.00 m N
#353000.00 m E
#4727172.00 m N

tl_poly <- polygon_from_extent(raster::extent(594549, 602563, 4645522, 4651452), proj4string='+proj=utm +datum=NAD83 +zone=10')
plot(tl_poly)

#602563.42 m E
#4651452.47 m N
#594549.86 m E
#4645522.56 m N

##download the NLCD raster and clip to the NW polygon (note: can only download 2011 with this function, not 2016)
#this takes a few mins because it's a large area

mend_nlcd <- get_nlcd(template = mendota_poly, label = 'mend_poly_NEW',  year = 2011, dataset = "landcover")
plot(mend_nlcd)

#this takes too much time - break up into smaller polygons of interest below
#nw_nlcd <- get_nlcd(template = nw_poly, label = 'nw_poly_NEW',  year = 2011, dataset = "landcover")
#plot(nw_nlcd)

#tule lake
tl_nlcd<- get_nlcd(template = tl_poly, label = 'tl_poly',  year = 2011, dataset = "landcover")
plot(tl_nlcd)

#fowl 41 hart mountain nwr
hart_nlcd <- get_nlcd(template = hart_poly, label = 'hart_poly',  year = 2011, dataset = "landcover")
plot(hart_nlcd)

#malheur nwr
mal_nlcd <- get_nlcd(template = mal_poly, label = 'mal_poly2',  year = 2011, dataset = "landcover")
plot(mal_nlcd)

#yakima valley
yak_nlcd <- get_nlcd(template = yakima_poly, label = 'yak_poly',  year = 2011, dataset = "landcover")
plot(yak_nlcd)

#idaho
id_nlcd <- get_nlcd(template = id_poly, label = 'id_poly4',  year = 2011, dataset = "landcover")
plot(id_nlcd)

#Alaska 2018 FOWL 37
ak18_nlcd <- get_nlcd(template = fowl37.18_poly, label = 'ak18_poly',  year = 2011, dataset = "landcover")
plot(ak18_nlcd)

#Alaska 2019 FOWL 37
ak19_nlcd <- get_nlcd(template = fowl37.19_poly, label = 'ak19_poly',  year = 2011, dataset = "landcover")
plot(ak19_nlcd)

library(adehabitatHR)
library(sp)

#bring in points for FOWL 36 - mendota

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

#bring in points for FOWL 40 -  Idaho

noha40.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/FOWL 40.csv"))

coordinates(noha40.df) <- c("x", "y")
proj4string(noha40.df) <- CRS("+proj=utm +zone=12 +datum=NAD83 +units=m")  ## for example
head(noha40.df)
plot(noha40.df)

newproj <- "+proj=utm +zone=12 +datum=NAD83 +units=m"
#nlcd_winter <- st_transform(Suisun_NLCD_winter, crs = the_crs)
nlcd_id_new <- projectRaster(id_nlcd, crs = newproj, method = "ngb" )

proj4string(noha40.df)
proj4string(nlcd_id_new)

plot(nlcd_id_new)
plot(noha40.df, add=T)

#bring in MCP for FOWL 40
library(adehabitatHR)
noha40.mcp <-mcp(noha40.df, percent = 100, unin = "m",unout = "km2")
proj4string(noha40.mcp)

plot(nlcd_id_new)
plot(noha40.mcp, add = T)
points(noha40.df, add=T)

#extract raster values from points
#nlcd_extract_40 <- extract(id_nlcd, noha40.df)
#head(id_nlcd)
#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_40 <- SpatialPoints(noha40.df, proj4string = crs(nlcd_id_new))
nlcd_extract_40 <- extract(nlcd_id_new, nlcd_sp_40)
#nlcd_extract_40.2 <- extract(id_nlcd, noha40.df)
head(nlcd_extract_40)
nlcd_extract_40[which(!is.na(nlcd_extract_40))]

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final_40 <- cbind.data.frame(noha40.df, nlcd_extract_40)
head(final_40)

#above works, but returns all columns, including empty grid cells with NA and 0 values

# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_40 <- final_40[which(!is.na(final_40[,3])),]
head(final_40)

write.csv(final_40, file = "nlcd.40.points.csv")

## raster values to mcp
clip40 <- mask(nlcd_id_new, noha40.mcp)
plot(clip40)

med.clip <- values(clip40)
final.clip <- med.clip[!is.na(med.clip)]
final.clip

write.csv(final.clip, file = "nlcd.40.mcp.csv")

plot(noha40.mcp)
plot(clip40, add =T)
plot(noha40.df, add=T, pch = 16)

#bring in points for FOWL 23 -  Tule Lake, CA Zone 10

noha23.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 23.csv"))

coordinates(noha23.df) <- c("x", "y")
proj4string(noha23.df) <- CRS("+proj=utm +zone=10 +datum=NAD83 +units=m")  ## for example
head(noha23.df)
plot(noha23.df)

newproj <- "+proj=utm +zone=10 +datum=NAD83 +units=m"
#nlcd_winter <- st_transform(Suisun_NLCD_winter, crs = the_crs)
nlcd_tl_new <- projectRaster(tl_nlcd, crs = newproj, method = "ngb" )

proj4string(noha23.df)
proj4string(nlcd_tl_new)

plot(nlcd_tl_new)
plot(noha23.df, add=T)

#bring in MCP for FOWL 23
library(adehabitatHR)
noha23.mcp <-mcp(noha23.df, percent = 100, unin = "m",unout = "km2")
proj4string(noha23.mcp)

plot(nlcd_tl_new)
plot(noha23.mcp, add = T)
points(noha23.df, add=T)

#extract raster values from points
#nlcd_extract_23 <- extract(tl_nlcd, noha23.df)
#head(tl_nlcd)
#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_23 <- SpatialPoints(noha23.df, proj4string = crs(nlcd_tl_new))
nlcd_extract_23 <- extract(nlcd_tl_new, nlcd_sp_23)
#nlcd_extract_23.2 <- extract(tl_nlcd, noha23.df)
head(nlcd_extract_23)
nlcd_extract_23[which(!is.na(nlcd_extract_23))]

# combine the raster cell probabilities with their coord pairs with landcover grtl cells
final_23 <- cbind.data.frame(noha23.df, nlcd_extract_23)
head(final_23)

#above works, but returns all columns, including empty grtl cells with NA and 0 values

# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_23 <- final_23[which(!is.na(final_23[,3])),]
head(final_23)

write.csv(final_23, file = "nlcd.23.points.csv")

## raster values to mcp
clip23 <- mask(nlcd_tl_new, noha23.mcp)
plot(clip23)
str(clip23)
med.clip <- values(clip23)
final.clip <- med.clip[!is.na(med.clip)]
final.clip

write.csv(final.clip, file = "nlcd.23.mcp.csv")

plot(noha23.mcp)
plot(clip23, add =T)
plot(noha23.df, add=T, pch = 16)

#bring in points for DRIB 61 -  Yakima Valley, WA Zone 10

noha23.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/DRIB 61.csv"))

coordinates(noha23.df) <- c("x", "y")
proj4string(noha23.df) <- CRS("+proj=utm +zone=10 +datum=NAD83 +units=m")  ## for example
head(noha23.df)
plot(noha23.df)

proj4string(noha23.df)
proj4string(nlcd_yak_new)

plot(nlcd_yak_new)
plot(noha23.df, add=T)

#bring in points for DRIB 61 -  Yakima Valley, WA Zone 10

drib61.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/DRIB 61.csv"))

coordinates(drib61.df) <- c("x", "y")
proj4string(drib61.df) <- CRS("+proj=utm +zone=10 +datum=NAD83 +units=m")  ## for example
head(drib61.df)
plot(drib61.df)

newproj <- "+proj=utm +zone=10 +datum=NAD83 +units=m"
nlcd_yak_new <- projectRaster(yak_nlcd, crs = newproj, method = "ngb" )

proj4string(drib61.df)
proj4string(nlcd_yak_new)

plot(nlcd_yak_new)
plot(drib61.df, add=T)

#bring in MCP for DRIB 61

drib61.mcp <-mcp(drib61.df, percent = 100, unin = "m",unout = "km2")
proj4string(drib61.mcp)

plot(nlcd_yak_new)
plot(drib61.mcp, add = T)
points(drib61.df, add=T)

#extract raster values from points
#nlcd_extract_23 <- extract(yak_nlcd, drib61.df)
#head(yak_nlcd)
#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_61 <- SpatialPoints(drib61.df, proj4string = crs(nlcd_yak_new))
nlcd_extract_61 <- extract(nlcd_yak_new, nlcd_sp_61)
#nlcd_extract_23.2 <- extract(yak_nlcd, drib61.df)
head(nlcd_extract_61)
nlcd_extract_61[which(!is.na(nlcd_extract_61))]

# combine the raster cell probabilities with their coord pairs with landcover gryak cells
final_61 <- cbind.data.frame(drib61.df, nlcd_extract_61)
head(final_61)

#above works, but returns all columns, including empty gryak cells with NA and 0 values

# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_61 <- final_61[which(!is.na(final_61[,3])),]
head(final_61)

write.csv(final_61, file = "nlcd.61.points.csv")

## raster values to mcp
clip61 <- mask(nlcd_yak_new, drib61.mcp)
plot(clip61)
str(clip61)
med.clip <- values(clip61)
final.clip <- med.clip[!is.na(med.clip)]
final.clip

write.csv(final.clip, file = "nlcd.61.mcp.csv")

plot(drib61.mcp)
plot(clip61, add =T)
plot(drib61.df, add=T, pch = 16)

#bring in points for DRIB 65 -  Malheur NWR  Zone 11

drib65.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/DRIB 65.csv"))

coordinates(drib65.df) <- c("x", "y")
proj4string(drib65.df) <- CRS("+proj=utm +zone=11 +datum=NAD83 +units=m")  ## for example
head(drib65.df)
plot(drib65.df)

newproj <- "+proj=utm +zone=11 +datum=NAD83 +units=m"
nlcd_mal_new <- projectRaster(mal_nlcd, crs = newproj, method = "ngb" )

proj4string(drib65.df)
proj4string(nlcd_mal_new)

plot(nlcd_mal_new)
plot(drib65.df, add=T)

#bring in MCP for DRIB 65

drib65.mcp <-mcp(drib65.df, percent = 100, unin = "m",unout = "km2")
proj4string(drib65.mcp)

plot(nlcd_mal_new)
plot(drib65.mcp, add = T)
points(drib65.df, add=T)

#extract raster values from points
#nlcd_extract_23 <- extract(mal_nlcd, drib65.df)
#head(mal_nlcd)
#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_65 <- SpatialPoints(drib65.df, proj4string = crs(nlcd_mal_new))
nlcd_extract_65 <- extract(nlcd_mal_new, nlcd_sp_65)
#nlcd_extract_23.2 <- extract(mal_nlcd, drib65.df)
head(nlcd_extract_65)
nlcd_extract_65[which(!is.na(nlcd_extract_65))]

# combine the raster cell probabilities with their coord pairs with landcover grmal cells
final_65 <- cbind.data.frame(drib65.df, nlcd_extract_65)
head(final_65)

#above works, but returns all columns, including empty grmal cells with NA and 0 values

# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_65 <- final_65[which(!is.na(final_65[,3])),]
head(final_65)

write.csv(final_65, file = "nlcd.65.points.csv")

## raster values to mcp
clip65 <- mask(nlcd_mal_new, drib65.mcp)
plot(clip65)
str(clip65)
med.clip <- values(clip65)
final.clip <- med.clip[!is.na(med.clip)]
final.clip

write.csv(final.clip, file = "nlcd.65.mcp.csv")

plot(drib65.mcp)
plot(clip65, add =T)
plot(drib65.df, add=T, pch = 16)

#bring in points for FOWL 22 -  Malheur NWR  Zone 11
library(raster)
drib22.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/FOWL 22.csv"))

coordinates(drib22.df) <- c("x", "y")
proj4string(drib22.df) <- CRS("+proj=utm +zone=11 +datum=NAD83 +units=m")  ## for example
head(drib22.df)
plot(drib22.df)

newproj <- "+proj=utm +zone=11 +datum=NAD83 +units=m"
nlcd_mal_new <- projectRaster(mal_nlcd, crs = newproj, method = "ngb" )

proj4string(drib22.df)
proj4string(nlcd_mal_new)

plot(nlcd_mal_new)
plot(drib22.df, add=T)

#bring in MCP for DRIB 22
library(adehabitatHR)
drib22.mcp <-mcp(drib22.df, percent = 100, unin = "m",unout = "km2")
proj4string(drib22.mcp)

plot(nlcd_mal_new)
plot(drib22.mcp, add = T)
points(drib22.df, add=T)

#extract raster values from points
#nlcd_extract_23 <- extract(mal_nlcd, drib22.df)
#head(mal_nlcd)
#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_22 <- SpatialPoints(drib22.df, proj4string = crs(nlcd_mal_new))
nlcd_extract_22 <- extract(nlcd_mal_new, nlcd_sp_22)
#nlcd_extract_23.2 <- extract(mal_nlcd, drib22.df)
head(nlcd_extract_22)
nlcd_extract_22[which(!is.na(nlcd_extract_22))]

# combine the raster cell probabilities with their coord pairs with landcover grmal cells
final_22 <- cbind.data.frame(drib22.df, nlcd_extract_22)
head(final_22)

#above works, but returns all columns, including empty grmal cells with NA and 0 values

# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_22 <- final_22[which(!is.na(final_22[,3])),]
head(final_22)

write.csv(final_22, file = "nlcd.22.points.csv")

## raster values to mcp
clip22 <- mask(nlcd_mal_new, drib22.mcp)
plot(clip22)
#str(clip22)
med.clip <- values(clip22)
final.clip <- med.clip[!is.na(med.clip)]
final.clip

write.csv(final.clip, file = "nlcd.22.mcp.csv")

plot(drib22.mcp)
plot(clip22, add =T)
plot(drib22.df, add=T, pch = 16)

#bring in points for FOWL 41 -  Hart Mountain NWR  Zone 11

fowl41.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/FOWL 41.csv"))

coordinates(fowl41.df) <- c("x", "y")
proj4string(fowl41.df) <- CRS("+proj=utm +zone=11 +datum=NAD83 +units=m")  ## for example
head(fowl41.df)
plot(fowl41.df)

newproj <- "+proj=utm +zone=11 +datum=NAD83 +units=m"
nlcd_hart_new <- projectRaster(hart_nlcd, crs = newproj, method = "ngb" )

proj4string(fowl41.df)
proj4string(nlcd_hart_new)

plot(nlcd_hart_new)
plot(fowl41.df, add=T)

#bring in MCP for FOWL 41

fowl41.mcp <-mcp(fowl41.df, percent = 100, unin = "m",unout = "km2")
proj4string(fowl41.mcp)

plot(nlcd_hart_new)
plot(fowl41.mcp, add = T)
points(fowl41.df, add=T)

#extract raster values from points
#nlcd_extract_23 <- extract(hart_nlcd, fowl41.df)
#head(hart_nlcd)
#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_41 <- SpatialPoints(fowl41.df, proj4string = crs(nlcd_hart_new))
nlcd_extract_41 <- extract(nlcd_hart_new, nlcd_sp_41)
#nlcd_extract_23.2 <- extract(hart_nlcd, fowl41.df)
head(nlcd_extract_41)
nlcd_extract_41[which(!is.na(nlcd_extract_41))]

# combine the raster cell probabilities with their coord pairs with landcover grhart cells
final_41 <- cbind.data.frame(fowl41.df, nlcd_extract_41)
head(final_41)

#above works, but returns all columns, including empty grhart cells with NA and 0 values

# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_41 <- final_41[which(!is.na(final_41[,3])),]
head(final_41)

write.csv(final_41, file = "nlcd.41.points.csv")

## raster values to mcp
clip41 <- mask(nlcd_hart_new, fowl41.mcp)
plot(clip41)
str(clip41)
med.clip <- values(clip41)
final.clip <- med.clip[!is.na(med.clip)]
final.clip

write.csv(final.clip, file = "nlcd.41.mcp.csv")

plot(fowl41.mcp)
plot(clip41, add =T)
plot(fowl41.df, add=T, pch = 16)

#bring in points for FOWL 37 -  Alaska 2018 Zone 4

fowl37.18.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/FOWL 37.2018.csv"))

coordinates(fowl37.18.df) <- c("x", "y")
proj4string(fowl37.18.df) <- CRS("+proj=utm +zone=4 +datum=NAD83 +units=m")  ## for example
head(fowl37.18.df)
plot(fowl37.18.df)

newproj <- "+proj=utm +zone=4 +datum=NAD83 +units=m"
nlcd_ak18_new <- projectRaster(ak18_nlcd, crs = newproj, method = "ngb" )

proj4string(fowl37.18.df)
proj4string(nlcd_ak18_new)

plot(nlcd_ak18_new)
plot(fowl37.18.df, add=T)

#bring in MCP for FOWL 37

fowl37.18.mcp <-mcp(fowl37.18.df, percent = 100, unin = "m",unout = "km2")
proj4string(fowl37.18.mcp)

plot(nlcd_ak18_new)
plot(fowl37.18.mcp, add = T)
points(fowl37.18.df, add=T)

#extract raster values from points
#nlcd_extract_23 <- extract(ak18_nlcd, fowl37.18.df)
#head(ak18_nlcd)
#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_37 <- SpatialPoints(fowl37.18.df, proj4string = crs(nlcd_ak18_new))
nlcd_extract_37 <- extract(nlcd_ak18_new, nlcd_sp_37)
#nlcd_extract_23.2 <- extract(ak18_nlcd, fowl37.18.df)
head(nlcd_extract_37)
nlcd_extract_37[which(!is.na(nlcd_extract_37))]

# combine the raster cell probabilities with their coord pairs with landcover grak18 cells
final_37 <- cbind.data.frame(fowl37.18.df, nlcd_extract_37)
head(final_37)

#above works, but returns all columns, including empty grak18 cells with NA and 0 values

# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_37 <- final_37[which(!is.na(final_37[,3])),]
head(final_37)

write.csv(final_37, file = "nlcd.37.points.csv")

## raster values to mcp
clip37 <- mask(nlcd_ak18_new, fowl37.18.mcp)
plot(clip37)
str(clip37)
med.clip <- values(clip37)
final.clip <- med.clip[!is.na(med.clip)]
final.clip

write.csv(final.clip, file = "nlcd.37.mcp.csv")

plot(fowl37.18.mcp)
plot(clip37, add =T)
plot(fowl37.18.df, add=T, pch = 16)

#bring in points for FOWL 37 -  Alaska 2019 Zone 4

fowl37.19.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/FOWL 37.2019.csv"))

coordinates(fowl37.19.df) <- c("x", "y")
proj4string(fowl37.19.df) <- CRS("+proj=utm +zone=4 +datum=NAD83 +units=m")  ## for example
head(fowl37.19.df)
plot(fowl37.19.df)

newproj <- "+proj=utm +zone=4 +datum=NAD83 +units=m"
nlcd_ak19_new <- projectRaster(ak19_nlcd, crs = newproj, method = "ngb" )

proj4string(fowl37.19.df)
proj4string(nlcd_ak19_new)

plot(nlcd_ak19_new)
plot(fowl37.19.df, add=T)

#bring in MCP for FOWL 37

fowl37.19.mcp <-mcp(fowl37.19.df, percent = 100, unin = "m",unout = "km2")
proj4string(fowl37.19.mcp)

plot(nlcd_ak19_new)
plot(fowl37.19.mcp, add = T)
points(fowl37.19.df, add=T)

#extract raster values from points
#nlcd_extract_23 <- extract(ak19_nlcd, fowl37.19.df)
#head(ak19_nlcd)
#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_37.19 <- SpatialPoints(fowl37.19.df, proj4string = crs(nlcd_ak19_new))
nlcd_extract_37.19 <- extract(nlcd_ak19_new, nlcd_sp_37.19)
#nlcd_extract_23.2 <- extract(ak19_nlcd, fowl37.19.df)
head(nlcd_extract_37.19)
nlcd_extract_37.19[which(!is.na(nlcd_extract_37.19))]

# combine the raster cell probabilities with their coord pairs with landcover grak19 cells
final_37.19 <- cbind.data.frame(fowl37.19.df, nlcd_extract_37.19)
head(final_37.19)

#above works, but returns all columns, including empty grak19 cells with NA and 0 values

# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_37.19 <- final_37.19[which(!is.na(final_37.19[,3])),]
head(final_37.19)

write.csv(final_37.19, file = "nlcd.37.19.points.csv")

## raster values to mcp
clip37.19 <- mask(nlcd_ak19_new, fowl37.19.mcp)
plot(clip37.19)
str(clip37.19)
med.clip <- values(clip37.19)
final.clip <- med.clip[!is.na(med.clip)]
final.clip

write.csv(final.clip, file = "nlcd.37.19.mcp.csv")

plot(fowl37.19.mcp)
plot(clip37.19, add =T)
plot(fowl37.19.df, add=T, pch = 16)

