## Attempting to create dBBMM utilization distributions##

## start with importing NLCD raster layer using help from: https://cran.r-project.org/web/packages/FedData/FedData.pdf 

## first we need to create an extent for a raster to which the NLCD layer will be cropped to when downloaded (because we can't just import the entire NLCD layer for the entire US)

getwd()
#setwd("/Users/Shannon/Desktop/R_Forever/Dissertation")

install.packages('FedData')
library(FedData)

suisun_polygon <- polygon_from_extent(raster::extent(578701, 601868, 4215762, 4234962), proj4string='+proj=utm +datum=WGS84 +zone=10 +ellps=WGS84')

##download the NLCD raster and clip to the Suisun polygon (note: can only download 2011 with this function, not 2016)

Suisun_NLCD <- get_nlcd(template = suisun_polygon, label = 'suisun',  year = 2011, dataset = "landcover")

plot(Suisun_NLCD)
#View(Suisun_NLCD)
#print(Suisun_NLCD)

## ok let's try dBBMM with Mama breeding locations from 2018 

install.packages("move")
library(move)

#bring in file from Movebank
mama_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/SIMP 02.csv")
#mama_move <- move(x = "C:/Users/sskalos/Documents/noha-move-hab/Data/SIMP 02.csv")
show(mama_move)
n.locs(mama_move) # number of locations
head(timeLag(mama_move, units="mins")) # time difference between locations - this is misleading because her first few days were set to two hour locations, but were then changd to 1 hour locations; make sure to look at ALL time differences
head(timestamps(mama_move))

# burst the movestack object to exclude any loactions that are greater than 72 minutes apart (because a couple locations are 72 and not 60 mis) - this is to prevent calculations of bridges and motion variance overnight between the last location of the previous day and the first location of the next morning, which are typically 400+ mins

mama_bursted <- move::burst(mama_move, c('normal','long')[1+(timeLag(mama_move, units='mins')>72)])

#plot mama's locations
#par(mfcol=1:2)
plot(mama_bursted, type="o", col=3, lwd=2, pch=20, xlab="location_long",ylab="location_lat")

#plot(mama_move, add =TRUE)
#str(mama_move)

#plot mama's locations with ggmap over map layer just to see if it is geographically correct! - it is
install.packages("ggmap")
library(ggmap)
require(ggmap) #these packages are necessary to work with google maps
#require(mapproj)
mama_df <- as(mama_bursted, "data.frame")
m <- get_map(bbox(extent(mama_bursted)*1.1), source="stamen", zoom=12)
ggmap(m)+geom_path(data=mama_df, aes(x=location.long, y=location.lat))

# transform coordinates from lat lon, center = T is required for the dbbmm to operate properly according to Bart on the movebank help chat

#mama_trans <- spTransform(mama_move, center=T)

mama_bursted_trans <- spTransform(x = mama_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m', center = T)
proj4string(mama_bursted_trans)

#str(mama_move)
#plot(mama_trans2, add = TRUE)
#str(mama_trans2)
str(mama_bursted_trans)

# bring in raster UTM with NAD 83 projection from ArcMAP because above code isn't working

nlcd_utm <- raster("~/Desktop/R_Forever/RRF/Data/Raster_UTM/NLCD_UTM.tif")
#nlcd_utm <- raster("C:/Users/sskalos/Documents/noha-move-hab/Data/Raster_UTM/NLCD_UTM.tif")
plot(nlcd_utm)
str(nlcd_utm)
as.data.frame(nlcd_utm, xy = TRUE)

#current projection - this is not working - something happens in this proejction transformation and the landcover values are lost for some reason.

Suisun_nlcd_trans <- raster(nlcd_utm)
crs(Suisun_nlcd_trans) <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-121.910674 +lat_0=38.115666 "
str(nlcd_utm)
str(Suisun_nlcd_trans)

plot(Suisun_nlcd_trans)


#run the dbbmm function - NOTE only odd numbers for margin and window size will produce valid likelihood estimations for motion variance

# From Kranstauber et al: increasing the size of the window increases reliability of the motion variance estimation at the cost of missing short term changes in the variation parameter, and should be close to, but less than, 24. Increasing the size of the margin, in contrast, enhances the power to identify "weak" breakpoints at the cost of not detecting breakpoints within the margin. 

mama_dbbmm <- brownian.bridge.dyn(mama_bursted_trans, burstType = 'normal', raster = Suisun_nlcd_trans, location.error = 10, ext = .3, time.step = 72, margin = 3, window.size = 7) #location error is 10 m as per the transmitter specifications, extent is 30% of raster extent, time step is 60 mins becasue locations were approximately every hour, margin is 3 which is the minimum number of locations needed to calculate breakpoints a a leave-one-out approach, and window size is 7 because this is equivalent to 7 locations, which equals  7 hours and may be able to detect behavioral changes within this relatively short window.

## below are the UDs calculated from the dbbmm
mama_dbbmm_UD<-new(".UD",calc(mama_dbbmm, sum)) ## it works!!!
head(mama_dbbmm_UD)
str(mama_dbbmm_UD)
summary(mama_dbbmm_UD)

## get the UD raster layer??
mama_ud <- UDStack(mama_dbbmm)
head(mama_ud)
str(mama_ud)
summary(mama_ud)

#now plot the UD on the left and the actual movement path on the right
#I can't figure out how to change the map area such that the map area is zoomed in, but whatever
par(mfrow=c(1,2))
plot(mama_dbbmm_UD, xlab="longitude", ylab="latitude")
#zoom(mama_dbbmm_UD)
plot(mama_dbbmm_UD, xlab="longitude", ylab="latitude")
lines(mama_bursted_trans, center=TRUE, col=3, lwd=2)
#plot(mama_dbbmm, xlab="location_long", ylab="location_lat")
#points(spTransform(mama_bursted, center=TRUE), col=8)

# this plots the 50% and 95% UD contour lines around the UDs
plot(mama_dbbmm_UD, xlab="location_long", ylab="location_lat")
contour(mama_dbbmm_UD, levels=c(.5, .95), col=c(6,2), add=TRUE, lwd=2)

#get the area of the 95% UD - i think these areas are in meters
mama_cont95 <- getVolumeUD(mama_dbbmm_UD)
mama_cont95 <- mama_cont95<=.95
area95 <- sum(values(mama_cont95))
area95

#get the area of the 50% UD - i think these areas are in meters
mama_cont5 <- getVolumeUD(mama_dbbmm_UD)
mama_cont5 <- mama_cont5<=.5
area5 <- sum(values(mama_cont5))
area5

## Ok, now let's follow Brian's steps and convert the DBBMM object to a SpatialLineDataFrame

#dbbmm dataframe- keep this!
dbbmm.df <- as.data.frame(mama_dbbmm_UD, xy = TRUE)
dbbmm.df[which(dbbmm.df[,3] == 0), 3] <- NA # this removes the zeros from the raster file
mama_ud_raster <- rasterFromXYZ(dbbmm.df, crs = "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0", digits = 5 )

## write raster - this is returning an empty raster and we don't know why
#writeRaster(mama_ud_raster, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster.tif", overwrite = TRUE)

writeRaster(mama_ud_raster, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/mama_ud_raster.tif", overwrite = TRUE)


#writeRaster(mama_dbbmm, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster2.tif", overwrite = TRUE)

#writeRaster(mama_ud, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster3.tif", overwrite = TRUE)


#writeRaster(mama_dbbmm_UD, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster4.tif", overwrite = TRUE)

writeRaster(mama_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/mama_ud_raster4.tif", overwrite = TRUE)

#but it does create the raster appropriately and plots below (weird)
plot(mama_ud_raster)
str(mama_ud_raster) # and there are values within this new raster (extra weird)

#save contours to shapefiles - this works

#example from Bart on Movebank - it works
#require(move)
#example(brownian.bridge.dyn)
#cont<-raster2contour(dbbmm, level=c(.5,.95))
#writeOGR(cont, dsn = '.', layer = 'mycontLines', driver = "ESRI Shapefile")

## try with mama data - it works
require(move)
cont2 <-raster2contour(mama_ud, level=c(.5,.95))
writeOGR(cont2, dsn = '.', layer = 'mama_contour2', driver = "ESRI Shapefile")

#works with either the mama_ud raster (above) or the mama_dbbmm_UD raster below, but below seems to be correct as it contains one layer, not every dbbmm step
require(move)
cont3 <-raster2contour(mama_dbbmm_UD, level=c(.5,.95))
writeOGR(cont3, dsn = '.', layer = 'mama_contour3', driver = "ESRI Shapefile")

# saves the contours to a kml file
install.packages('plotKML')
require(plotKML)
kml(cont3)

#getwd()
#ploygon <- readOGR(dsn = "/Users/Shannon/Desktop/R_Forever/Dissertation/noha-move-hab/Output" ,layer = "mama_contour3")
#mama_WGS <- spTransform(cont3, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) 
#writeOGR(cont3, dsn = "mama_contours.kml", layer = "mama_WGS", driver = "KML", overwrite_layer = TRUE)

#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp <- SpatialPoints(dbbmm.df[,1:2], proj4string = crs(Suisun_nlcd_trans))
nlcd_extract <- extract(nlcd_utm, nlcd_sp)
head(nlcd_extract)
nlcd_extract[which(!is.na(nlcd_extract))]
#str(Suisun_NLCD)
str(Suisun_nlcd_trans)
str(mama_dbbmm_UD)
str(nlcd_utm)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(mama_dbbmm_UD)
library(scales)
plot(nlcd_utm, col = alpha("red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final <- cbind.data.frame(dbbmm.df, nlcd_extract)
head(final)

#above works, but returns all columns, including empty grid cells with NA and 0 values
# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final <- final[which(!is.na(final[,4])),]
head(final)

# for loop to calculate probabilities of use within each landcover types using the UDs
prob.vec <- rep(NA, length(unique(final[,4])))
unique.vec <- unique(final[,4])
tot.prob <- sum(final[,3])
for (i in 1:length(prob.vec)){
  prob.vec[i] <- sum(final[which(final[,4] == unique.vec[i]),3])/tot.prob
}
#check that the for loop worked and the probabilities sum to 1 - they do
sum(prob.vec)

#save the probability table for each landcover class - it works!
probs.cover.tables <- cbind(prob.vec, unique.vec)

#view the entire table
probs.cover.tables
write.csv(probs.cover.tables, file = "mama_landcover_probs_final.csv")

longi <- mama_bursted@data$utm.easting
lattitude <- mama_bursted@data$utm.northing

plot(nlcd_utm)
points(longi,lattitude, cex = .1, add = TRUE)
#zoom(nlcd_utm, drawExtent())
summary(longi)
summary(lattitude)

## Now let's do the above with the other females
## SIMP 01 Laureen

laureen_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/SIMP 01.csv")
show(laureen_move)
n.locs(laureen_move) # number of locations
head(timeLag(laureen_move, units="mins")) # time difference between locations - this is misleading because her first few days were set to two hour locations, but were then changd to 1 hour locations; make sure to look at ALL time differences
head(timestamps(laureen_move))
summary(timeLag(laureen_move, units="mins"))
laureen_lag <- timeLag(laureen_move, units="mins")
write.csv(laureen_lag, file = "laureen_lag.csv")

# burst the movestack object to exclude any loactions that are greater than 72 minutes apart (because a several locations are 62-88 mins and not 60 mis) - this is to prevent calculations of bridges and motion variance overnight between the last location of the previous day and the first location of the next morning, which are typically 400+ mins

laureen_bursted <- move::burst(laureen_move, c('normal','long')[1+(timeLag(laureen_move, units='mins')>88)])

#plot laureen's locations
#par(mfcol=1:2)
plot(laureen_bursted, type="o", col=3, lwd=2, pch=20, xlab="location_long",ylab="location_lat")

#plot laureen's locations with ggmap over map layer just to see if it is geographically correct! - it is
install.packages("ggmap")
library(ggmap)
require(ggmap) #these packages are necessary to work with google maps
#require(mapproj)
laureen_df <- as(laureen_bursted, "data.frame")
m <- get_map(bbox(extent(laureen_bursted)*1.1), source="stamen", zoom=12)
ggmap(m)+geom_path(data=laureen_df, aes(x=location.long, y=location.lat))

# transform coordinates from lat lon, center = T is required for the dbbmm to operate properly according to Bart on the movebank help chat


laureen_bursted_trans <- spTransform(x = laureen_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575', center = T)
proj4string(laureen_bursted_trans)

str(laureen_bursted)

# bring in raster UTM with NAD 83 projection from ArcMAP because above code isn't working

nlcd_utm_laureen <- raster("~/Desktop/R_Forever/RRF/Data/Raster_UTM/NLCD_UTM.tif")
#plot(nlcd_utm)
#str(nlcd_utm)
#as.data.frame(nlcd_utm, xy = TRUE)

#current projection - this is not working - something happens in this proejction transformation and the landcover values are lost for some reason.
str(nlcd_utm_laureen)
str(laureen_bursted_trans)

Suisun_nlcd_trans_laureen <- raster(nlcd_utm_laureen)
crs(Suisun_nlcd_trans_laureen) <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575 "

#now they match coordinate systems
str(laureen_bursted_trans)
str(Suisun_nlcd_trans_laureen)

laureen_dbbmm <- brownian.bridge.dyn(laureen_bursted_trans, burstType = 'normal', raster = Suisun_nlcd_trans_laureen, location.error = 10, ext = .3, time.step = 60, margin = 3, window.size = 7) #location error is 10 m as per the transmitter specifications, extent is 30% of raster extent, time step is 60 mins becasue locations were approximately every hour, margin is 3 which is the minimum number of locations needed to calculate breakpoints a a leave-one-out approach, and window size is 7 because this is equivalent to 7 locations, which equals  7 hours and may be able to detect behavioral changes within this relatively short window.

## below are the UDs calculated from the dbbmm
laureen_dbbmm_UD<-new(".UD",calc(laureen_dbbmm, sum)) ## it works!!!
head(laureen_dbbmm_UD)
str(laureen_dbbmm_UD)
summary(laureen_dbbmm_UD)

## get the UD raster layer?? - don't need this
laureen_ud <- UDStack(laureen_dbbmm)
summary(laureen_ud)

#now plot the UD on the left and the actual movement path on the right
#I can't figure out how to change the map area such that the map area is zoomed in, but whatever
par(mfrow=c(1,2))
plot(laureen_dbbmm_UD, xlab="longitude", ylab="latitude")
#zoom(laureen_dbbmm_UD)
plot(laureen_dbbmm_UD, xlab="longitude", ylab="latitude")
lines(laureen_bursted_trans, center=TRUE, col=3, lwd=2)
#plot(mama_dbbmm, xlab="location_long", ylab="location_lat")
#points(spTransform(mama_bursted, center=TRUE), col=8)

# this plots the 50% and 95% UD contour lines around the UDs
plot(laureen_dbbmm_UD, xlab="location_long", ylab="location_lat")
contour(laureen_dbbmm_UD, levels=c(.5, .95), col=c(6,2), add=TRUE, lwd=2)

#get the area of the 95% UD - i think these areas are in meters
laureen_cont95 <- getVolumeUD(laureen_dbbmm_UD)
laureen_cont95 <- laureen_cont95<=.95
area95 <- sum(values(laureen_cont95))
area95

#get the area of the 50% UD - i think these areas are in meters
laureen_cont5 <- getVolumeUD(laureen_dbbmm_UD)
laureen_cont5 <- laureen_cont5<=.5
area5 <- sum(values(laureen_cont5))
area5

## Ok, now let's follow Brian's steps and convert the DBBMM object to a SpatialLineDataFrame

#dbbmm dataframe- keep this!
laureen.dbbmm.df <- as.data.frame(laureen_dbbmm_UD, xy = TRUE)
laureen_ud_raster <- rasterFromXYZ(laureen.dbbmm.df, crs = "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0", digits = 5 )

## write raster - this is returning an empty raster and we don't know why
#writeRaster(laureen_ud_raster, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster.tif", overwrite = TRUE)

#writeRaster(laureen_dbbmm, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster2.tif", overwrite = TRUE)

#writeRaster(laureen_ud, "C:/Users/sskalos/Documents/noha-move-hab/Output/lauren_ud_raster3.tif", overwrite = TRUE)

writeRaster(laureen_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/laureen_ud_raster.tif", overwrite = TRUE)

#but it does create the raster appropriately and plots below (weird)
plot(laureen_ud_raster)
str(laureen_ud_raster) # and there are values within this new raster (extra weird)

#save contours to shapefiles - this works

#example from Bart on Movebank - it works
#require(move)
#example(brownian.bridge.dyn)
#cont<-raster2contour(dbbmm, level=c(.5,.95))
#writeOGR(cont, dsn = '.', layer = 'mycontLines', driver = "ESRI Shapefile")

## try with mama data - it works
require(move)
#cont2 <-raster2contour(laureen_ud, level=c(.5,.95))
#writeOGR(cont2, dsn = '.', layer = 'laureen_contour2', driver = "ESRI Shapefile")

#works with either the mama_ud raster (above) or the mama_dbbmm_UD raster below, but below seems to be correct as it contains one layer, not every dbbmm step
require(move)
cont_laureen <-raster2contour(laureen_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_laureen, dsn = '.', layer = 'laureen_contour', driver = "ESRI Shapefile")

#getwd()
#ploygon <- readOGR(dsn = "/Users/Shannon/Desktop/R_Forever/Dissertation/noha-move-hab/Output" ,layer = "mama_contour3")
#mama_WGS <- spTransform(cont3, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) 
#writeOGR(cont3, dsn = "mama_contours.kml", layer = "mama_WGS", driver = "KML", overwrite_layer = TRUE)

#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_laureen <- SpatialPoints(laureen.dbbmm.df[,1:2], proj4string = crs(Suisun_nlcd_trans_laureen))
nlcd_extract_laureen <- extract(nlcd_utm, nlcd_sp_laureen)
head(nlcd_extract_laureen)
nlcd_extract_laureen[which(!is.na(nlcd_extract_laureen))]

## check that both layers below contain the same # of columns and rows
str(Suisun_nlcd_trans_laureen)
str(laureen_dbbmm_UD)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(laureen_dbbmm_UD)
library(scales)
plot(nlcd_utm, col = alpha("red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final_laureen <- cbind.data.frame(laureen.dbbmm.df, nlcd_extract_laureen)
head(final_laureen)

#above works, but returns all columns, including empty grid cells with NA and 0 values
# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_laureen <- final_laureen[which(!is.na(final_laureen[,4])),]
head(final_laureen)

# for loop to calculate probabilities of use within each landcover types using the UDs
prob.vec <- rep(NA, length(unique(final_laureen[,4])))
unique.vec <- unique(final_laureen[,4])
tot.prob <- sum(final_laureen[,3])
for (i in 1:length(prob.vec)){
  prob.vec[i] <- sum(final_laureen[which(final_laureen[,4] == unique.vec[i]),3])/tot.prob
}
#check that the for loop worked and the probabilities sum to 1 - they do
sum(prob.vec)

#save the probability table for each landcover class - it works!
probs.cover.tables <- cbind(prob.vec, unique.vec)

#view the entire table
probs.cover.tables
write.csv(probs.cover.tables, file = "laureen_landcover_probs_final.csv")

longi <- laureen_bursted@data$utm.easting
lattitude <- laureen_bursted@data$utm.northing

plot(nlcd_utm)
points(longi,lattitude, cex = 1, add = TRUE)
#zoom(nlcd_utm, drawExtent())
summary(longi)
summary(lattitude)

## Now let's do the above with the other females
## SIMP 07 Joyce

joyce_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/SIMP 07.csv")
show(joyce_move)
n.locs(joyce_move) # number of locations
head(timeLag(joyce_move, units="mins")) # time difference between locations - this is misleading because her first few days were set to two hour locations, but were then changd to 1 hour locations; make sure to look at ALL time differences
head(timestamps(joyce_move))
summary(timeLag(joyce_move, units="mins"))
joyce_lag <- timeLag(joyce_move, units="mins")
write.csv(joyce_lag, file = "joyce_lag.csv")

# burst the movestack object to exclude any loactions that are greater than 72 minutes apart (because a several locations are 62 mins and not 60 mis) - this is to prevent calculations of bridges and motion variance overnight between the last location of the previous day and the first location of the next morning, which are typically 400+ mins

joyce_bursted <- move::burst(joyce_move, c('normal','long')[1+(timeLag(joyce_move, units='mins')>62)])

#plot joyce's locations
#par(mfcol=1:2)
plot(joyce_bursted, type="o", col=3, lwd=2, pch=20, xlab="location_long",ylab="location_lat")

#plot joyce's locations with ggmap over map layer just to see if it is geographically correct! - it is
install.packages("ggmap")
library(ggmap)
require(ggmap) #these packages are necessary to work with google maps
#require(mapproj)
joyce_df <- as(joyce_bursted, "data.frame")
m <- get_map(bbox(extent(joyce_bursted)*1.1), source="stamen", zoom=12)
ggmap(m)+geom_path(data=joyce_df, aes(x=location.long, y=location.lat))

# transform coordinates from lat lon, center = T is required for the dbbmm to operate properly according to Bart on the movebank help chat

joyce_bursted_trans <- spTransform(x = joyce_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0', center = T)
proj4string(joyce_bursted_trans)

str(joyce_bursted)

# bring in raster UTM with NAD 83 projection from ArcMAP because above code isn't working

nlcd_utm_joyce <- raster("~/Desktop/R_Forever/RRF/Data/Raster_UTM/NLCD_UTM.tif")
#plot(nlcd_utm)
#str(nlcd_utm)
#as.data.frame(nlcd_utm, xy = TRUE)

#current projection - this is not working - something happens in this proejction transformation and the landcover values are lost for some reason.
str(nlcd_utm_joyce)
str(joyce_bursted)

Suisun_nlcd_trans_joyce <- raster(nlcd_utm_joyce)
crs(Suisun_nlcd_trans_joyce) <- "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.003483 +lat_0=38.2062325"

#now they match
str(joyce_bursted_trans)
str(Suisun_nlcd_trans_joyce)

joyce_dbbmm <- brownian.bridge.dyn(joyce_bursted_trans, burstType = 'normal', raster = Suisun_nlcd_trans_joyce, location.error = 10, ext = .3, time.step = 60, margin = 3, window.size = 7) #location error is 10 m as per the transmitter specifications, extent is 30% of raster extent, time step is 60 mins becasue locations were approximately every hour, margin is 3 which is the minimum number of locations needed to calculate breakpoints a a leave-one-out approach, and window size is 7 because this is equivalent to 7 locations, which equals  7 hours and may be able to detect behavioral changes within this relatively short window.

## below are the UDs calculated from the dbbmm
joyce_dbbmm_UD<-new(".UD",calc(joyce_dbbmm, sum)) ## it works!!!
head(joyce_dbbmm_UD)
str(joyce_dbbmm_UD)
summary(joyce_dbbmm_UD)

## get the UD raster layer??
joyce_ud <- UDStack(joyce_dbbmm)
summary(joyce_ud)

#now plot the UD on the left and the actual movement path on the right
#I can't figure out how to change the map area such that the map area is zoomed in, but whatever
par(mfrow=c(1,2))
plot(joyce_dbbmm_UD, xlab="longitude", ylab="latitude")
#zoom(joyce_dbbmm_UD)
plot(joyce_dbbmm_UD, xlab="longitude", ylab="latitude")
lines(joyce_bursted_trans, center=TRUE, col=3, lwd=2)
#plot(mama_dbbmm, xlab="location_long", ylab="location_lat")
#points(spTransform(mama_bursted, center=TRUE), col=8)

# this plots the 50% and 95% UD contour lines around the UDs
plot(joyce_dbbmm_UD, xlab="location_long", ylab="location_lat")
contour(joyce_dbbmm_UD, levels=c(.5, .95), col=c(6,2), add=TRUE, lwd=2)

#get the area of the 95% UD - i think these areas are in meters
joyce_cont95 <- getVolumeUD(joyce_dbbmm_UD)
joyce_cont95 <- joyce_cont95<=.95
area95 <- sum(values(joyce_cont95))
area95

#get the area of the 50% UD - i think these areas are in meters
joyce_cont5 <- getVolumeUD(joyce_dbbmm_UD)
joyce_cont5 <- joyce_cont5<=.5
area5 <- sum(values(joyce_cont5))
area5

## Ok, now let's follow Brian's steps and convert the DBBMM object to a SpatialLineDataFrame

#dbbmm dataframe- keep this!
joyce.dbbmm.df <- as.data.frame(joyce_dbbmm_UD, xy = TRUE)
joyce_ud_raster <- rasterFromXYZ(joyce.dbbmm.df, crs = "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0", digits = 5 )

## write raster - this is returning an empty raster and we don't know why
#writeRaster(joyce_ud_raster, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster.tif", overwrite = TRUE)

#writeRaster(joyce_dbbmm, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster2.tif", overwrite = TRUE)

#writeRaster(joyce_ud, "C:/Users/sskalos/Documents/noha-move-hab/Output/lauren_ud_raster3.tif", overwrite = TRUE)

writeRaster(joyce_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/joyce_ud_raster.tif", overwrite = TRUE)

#but it does create the raster appropriately and plots below (weird)
plot(joyce_ud_raster)
str(joyce_ud_raster) # and there are values within this new raster (extra weird)

#save contours to shapefiles - this works

#example from Bart on Movebank - it works
#require(move)
#example(brownian.bridge.dyn)
#cont<-raster2contour(dbbmm, level=c(.5,.95))
#writeOGR(cont, dsn = '.', layer = 'mycontLines', driver = "ESRI Shapefile")

## try with mama data - it works
require(move)
#cont2 <-raster2contour(joyce_ud, level=c(.5,.95))
#writeOGR(cont2, dsn = '.', layer = 'joyce_contour2', driver = "ESRI Shapefile")

#works with either the mama_ud raster (above) or the mama_dbbmm_UD raster below, but below seems to be correct as it contains one layer, not every dbbmm step
require(move)
cont_joyce <-raster2contour(joyce_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_joyce, dsn = '.', layer = 'joyce_contour', driver = "ESRI Shapefile")

#getwd()
#ploygon <- readOGR(dsn = "/Users/Shannon/Desktop/R_Forever/Dissertation/noha-move-hab/Output" ,layer = "mama_contour3")
#mama_WGS <- spTransform(cont3, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) 
#writeOGR(cont3, dsn = "mama_contours.kml", layer = "mama_WGS", driver = "KML", overwrite_layer = TRUE)

#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_joyce <- SpatialPoints(joyce.dbbmm.df[,1:2], proj4string = crs(Suisun_nlcd_trans_joyce))
nlcd_extract_joyce <- extract(nlcd_utm, nlcd_sp_joyce)
head(nlcd_extract_joyce)
nlcd_extract_joyce[which(!is.na(nlcd_extract_joyce))]

# check to make sure they have the correct # of columns and rows
str(Suisun_nlcd_trans_joyce)
str(joyce_dbbmm_UD)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(joyce_dbbmm_UD)
library(scales)
plot(nlcd_utm, col = alpha("red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final_joyce <- cbind.data.frame(joyce.dbbmm.df, nlcd_extract_joyce)
head(final_joyce)

#above works, but returns all columns, including empty grid cells with NA and 0 values
# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_joyce <- final_joyce[which(!is.na(final_joyce[,4])),]
head(final_joyce)

# for loop to calculate probabilities of use within each landcover types using the UDs
prob.vec <- rep(NA, length(unique(final_joyce[,4])))
unique.vec <- unique(final_joyce[,4])
tot.prob <- sum(final_joyce[,3])
for (i in 1:length(prob.vec)){
  prob.vec[i] <- sum(final_joyce[which(final_joyce[,4] == unique.vec[i]),3])/tot.prob
}
#check that the for loop worked and the probabilities sum to 1 - they do
sum(prob.vec)

#save the probability table for each landcover class - it works!
probs.cover.tables <- cbind(prob.vec, unique.vec)

#view the entire table
probs.cover.tables
write.csv(probs.cover.tables, file = "joyce_landcover_probs_final.csv")

longi <- joyce_bursted@data$utm.easting
lattitude <- joyce_bursted@data$utm.northing

plot(nlcd_utm)
points(longi,lattitude, cex = 1, add = TRUE)
#zoom(nlcd_utm, drawExtent())
summary(longi)
summary(lattitude)

## Now let's do the above with the other females
## SIMP 10 Cinnamon

cinnamon_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/SIMP 10.csv")
show(cinnamon_move)
n.locs(cinnamon_move) # number of locations
head(timeLag(cinnamon_move, units="mins")) # time difference between locations - this is misleading because her first few days were set to two hour locations, but were then changd to 1 hour locations; make sure to look at ALL time differences
head(timestamps(cinnamon_move))
summary(timeLag(cinnamon_move, units="mins"))
cinnamon_lag <- timeLag(cinnamon_move, units="mins")
write.csv(cinnamon_lag, file = "cinnamon_lag.csv")

# burst the movestack object to exclude any loactions that are greater than 72 minutes apart (because a several locations are 74 mins and not 60 mis) - this is to prevent calculations of bridges and motion variance overnight between the last location of the previous day and the first location of the next morning, which are typically 400+ mins

cinnamon_bursted <- move::burst(cinnamon_move, c('normal','long')[1+(timeLag(cinnamon_move, units='mins')>74)])

#plot cinnamon's locations
#par(mfcol=1:2)
plot(cinnamon_bursted, type="o", col=3, lwd=2, pch=20, xlab="location_long",ylab="location_lat")

#plot cinnamon's locations with ggmap over map layer just to see if it is geographically correct! - it is
install.packages("ggmap")
library(ggmap)
require(ggmap) #these packages are necessary to work with google maps
#require(mapproj)
cinnamon_df <- as(cinnamon_bursted, "data.frame")
m <- get_map(bbox(extent(cinnamon_bursted)*1.1), source="stamen", zoom=12)
ggmap(m)+geom_path(data=cinnamon_df, aes(x=location.long, y=location.lat))

# transform coordinates from lat lon, center = T is required for the dbbmm to operate properly according to Bart on the movebank help chat


cinnamon_bursted_trans <- spTransform(x = cinnamon_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575', center = T)
proj4string(cinnamon_bursted_trans)

str(cinnamon_bursted)

# bring in raster UTM with NAD 83 projection from ArcMAP because above code isn't working

nlcd_utm_cinn <- raster("~/Desktop/R_Forever/RRF/Data/Raster_UTM/NLCD_UTM.tif")
#plot(nlcd_utm)
#str(nlcd_utm)
#as.data.frame(nlcd_utm, xy = TRUE)

#current projection - this is not working - something happens in this proejction transformation and the landcover values are lost for some reason.
str(nlcd_utm_cinn)
str(cinnamon_bursted_trans)

Suisun_nlcd_trans_cinnamon <- raster(nlcd_utm_cinn)
crs(Suisun_nlcd_trans_cinnamon) <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575 "

#now they match
str(cinnamon_bursted_trans)
str(Suisun_nlcd_trans_cinnamon)

cinnamon_dbbmm <- brownian.bridge.dyn(cinnamon_bursted_trans, burstType = 'normal', raster = Suisun_nlcd_trans_cinnamon, location.error = 10, ext = .3, time.step = 60, margin = 3, window.size = 7) #location error is 10 m as per the transmitter specifications, extent is 30% of raster extent, time step is 60 mins becasue locations were approximately every hour, margin is 3 which is the minimum number of locations needed to calculate breakpoints a a leave-one-out approach, and window size is 7 because this is equivalent to 7 locations, which equals  7 hours and may be able to detect behavioral changes within this relatively short window.

## below are the UDs calculated from the dbbmm
cinnamon_dbbmm_UD<-new(".UD",calc(cinnamon_dbbmm, sum)) ## it works!!!
head(cinnamon_dbbmm_UD)
str(cinnamon_dbbmm_UD)
summary(cinnamon_dbbmm_UD)

## get the UD raster layer??
cinnamon_ud <- UDStack(cinnamon_dbbmm)
summary(cinnamon_ud)

#now plot the UD on the left and the actual movement path on the right
#I can't figure out how to change the map area such that the map area is zoomed in, but whatever
par(mfrow=c(1,2))
plot(cinnamon_dbbmm_UD, xlab="longitude", ylab="latitude")
#zoom(cinnamon_dbbmm_UD)
plot(cinnamon_dbbmm_UD, xlab="longitude", ylab="latitude")
lines(cinnamon_bursted_trans, center=TRUE, col=3, lwd=2)
#plot(mama_dbbmm, xlab="location_long", ylab="location_lat")
#points(spTransform(mama_bursted, center=TRUE), col=8)

# this plots the 50% and 95% UD contour lines around the UDs
plot(cinnamon_dbbmm_UD, xlab="location_long", ylab="location_lat")
contour(cinnamon_dbbmm_UD, levels=c(.5, .95), col=c(6,2), add=TRUE, lwd=2)

#get the area of the 95% UD - i think these areas are in meters
cinnamon_cont95 <- getVolumeUD(cinnamon_dbbmm_UD)
cinnamon_cont95 <- cinnamon_cont95<=.95
area95 <- sum(values(cinnamon_cont95))
area95

#get the area of the 50% UD - i think these areas are in meters
cinnamon_cont5 <- getVolumeUD(cinnamon_dbbmm_UD)
cinnamon_cont5 <- cinnamon_cont5<=.5
area5 <- sum(values(cinnamon_cont5))
area5

## Ok, now let's follow Brian's steps and convert the DBBMM object to a SpatialLineDataFrame

#dbbmm dataframe- keep this!
cinnamon.dbbmm.df <- as.data.frame(cinnamon_dbbmm_UD, xy = TRUE)
cinnamon_ud_raster <- rasterFromXYZ(cinnamon.dbbmm.df, crs = "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0", digits = 5 )

## write raster - this is returning an empty raster and we don't know why
#writeRaster(cinnamon_ud_raster, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster.tif", overwrite = TRUE)

#writeRaster(cinnamon_dbbmm, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster2.tif", overwrite = TRUE)

#writeRaster(cinnamon_ud, "C:/Users/sskalos/Documents/noha-move-hab/Output/lauren_ud_raster3.tif", overwrite = TRUE)

writeRaster(cinnamon_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/cinnamon_ud_raster.tif", overwrite = TRUE)

#but it does create the raster appropriately and plots below (weird)
plot(cinnamon_ud_raster)
str(cinnamon_ud_raster) # and there are values within this new raster (extra weird)

#save contours to shapefiles - this works

#example from Bart on Movebank - it works
#require(move)
#example(brownian.bridge.dyn)
#cont<-raster2contour(dbbmm, level=c(.5,.95))
#writeOGR(cont, dsn = '.', layer = 'mycontLines', driver = "ESRI Shapefile")

## try with mama data - it works
require(move)
#cont2 <-raster2contour(cinnamon_ud, level=c(.5,.95))
#writeOGR(cont2, dsn = '.', layer = 'cinnamon_contour2', driver = "ESRI Shapefile")

#works with either the mama_ud raster (above) or the mama_dbbmm_UD raster below, but below seems to be correct as it contains one layer, not every dbbmm step
require(move)
cont_cinnamon <-raster2contour(cinnamon_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_cinnamon, dsn = '.', layer = 'cinnamon_contour', driver = "ESRI Shapefile")

#getwd()
#ploygon <- readOGR(dsn = "/Users/Shannon/Desktop/R_Forever/Dissertation/noha-move-hab/Output" ,layer = "mama_contour3")
#mama_WGS <- spTransform(cont3, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) 
#writeOGR(cont3, dsn = "mama_contours.kml", layer = "mama_WGS", driver = "KML", overwrite_layer = TRUE)

#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_cinnamon <- SpatialPoints(cinnamon.dbbmm.df[,1:2], proj4string = crs(Suisun_nlcd_trans_cinnamon))
nlcd_extract_cinnamon <- extract(nlcd_utm, nlcd_sp_cinnamon)
head(nlcd_extract_cinnamon)
nlcd_extract_cinnamon[which(!is.na(nlcd_extract_cinnamon))]

#check to make sure the have the same # of columns and rows
str(Suisun_nlcd_trans_cinnamon)
str(cinnamon_dbbmm_UD)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(cinnamon_dbbmm_UD)
library(scales)
plot(nlcd_utm, col = alpha("red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final_cinnamon <- cbind.data.frame(cinnamon.dbbmm.df, nlcd_extract_cinnamon)
head(final_cinnamon)

#above works, but returns all columns, including empty grid cells with NA and 0 values
# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_cinnamon <- final_cinnamon[which(!is.na(final_cinnamon[,4])),]
head(final_cinnamon)

# for loop to calculate probabilities of use within each landcover types using the UDs
prob.vec <- rep(NA, length(unique(final_cinnamon[,4])))
unique.vec <- unique(final_cinnamon[,4])
tot.prob <- sum(final_cinnamon[,3])
for (i in 1:length(prob.vec)){
  prob.vec[i] <- sum(final_cinnamon[which(final_cinnamon[,4] == unique.vec[i]),3])/tot.prob
}
#check that the for loop worked and the probabilities sum to 1 - they do
sum(prob.vec)

#save the probability table for each landcover class - it works!
probs.cover.tables <- cbind(prob.vec, unique.vec)

#view the entire table
probs.cover.tables
write.csv(probs.cover.tables, file = "cinnamon_landcover_probs_final.csv")

longi <- cinnamon_bursted@data$utm.easting
lattitude <- cinnamon_bursted@data$utm.northing

plot(nlcd_utm)
points(longi,lattitude, cex = 1, add = TRUE)
#zoom(nlcd_utm, drawExtent())
summary(longi)
summary(lattitude)


## Now let's do the above with the other females
## NOHA 30.1 Bre
library(move)
bre_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 30.1.csv")
show(bre_move)
n.locs(bre_move) # number of locations
head(timeLag(bre_move, units="mins")) # time difference between locations - this is misleading because her first few days were set to two hour locations, but were then changd to 1 hour locations; make sure to look at ALL time differences
head(timestamps(bre_move))
summary(timeLag(bre_move, units="mins"))
bre_lag <- timeLag(bre_move, units="mins")
write.csv(bre_lag, file = "bre_lag.csv")

# burst the movestack object to exclude any loactions that are greater than 31.5 minutes apart (because a several locations are 31.5 mins and not 30 mis) - this is to prevent calculations of bridges and motion variance overnight between the last location of the previous day and the first location of the next morning, which are typically 400+ mins

bre_bursted <- move::burst(bre_move, c('normal','long')[1+(timeLag(bre_move, units='mins')>31.5)])

#plot bre's locations
#par(mfcol=1:2)
plot(bre_bursted, type="o", col=3, lwd=2, pch=20, xlab="location_long",ylab="location_lat")

#plot bre's locations with ggmap over map layer just to see if it is geographically correct! - it is
install.packages("ggmap")
library(ggmap)
require(ggmap) #these packages are necessary to work with google maps
#require(mapproj)
bre_df <- as(bre_bursted, "data.frame")
m <- get_map(bbox(extent(bre_bursted)*1.1), source="stamen", zoom=12)
ggmap(m)+geom_path(data=bre_df, aes(x=location.long, y=location.lat))

# transform coordinates from lat lon, center = T is required for the dbbmm to operate properly according to Bart on the movebank help chat


bre_bursted_trans <- spTransform(x = bre_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575', center = T)
proj4string(bre_bursted_trans)

str(bre_bursted)

# bring in raster UTM with NAD 83 projection from ArcMAP because above code isn't working

nlcd_utm_cinn <- raster("~/Desktop/R_Forever/RRF/Data/Raster_UTM/NLCD_UTM.tif")
#plot(nlcd_utm)
#str(nlcd_utm)
#as.data.frame(nlcd_utm, xy = TRUE)

#current projection - this is not working - something happens in this proejction transformation and the landcover values are lost for some reason.
str(nlcd_utm_cinn)
str(bre_bursted_trans)

Suisun_nlcd_trans_bre <- raster(nlcd_utm_cinn)
crs(Suisun_nlcd_trans_bre) <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575 "

#now they match
str(bre_bursted_trans)
str(Suisun_nlcd_trans_bre)

bre_dbbmm <- brownian.bridge.dyn(bre_bursted_trans, burstType = 'normal', raster = Suisun_nlcd_trans_bre, location.error = 10, ext = .3, time.step = 60, margin = 3, window.size = 7) #location error is 10 m as per the transmitter specifications, extent is 30% of raster extent, time step is 60 mins becasue locations were approximately every hour, margin is 3 which is the minimum number of locations needed to calculate breakpoints a a leave-one-out approach, and window size is 7 because this is equivalent to 7 locations, which equals  7 hours and may be able to detect behavioral changes within this relatively short window.

## below are the UDs calculated from the dbbmm
bre_dbbmm_UD<-new(".UD",calc(bre_dbbmm, sum)) ## it works!!!
head(bre_dbbmm_UD)
str(bre_dbbmm_UD)
summary(bre_dbbmm_UD)

## get the UD raster layer??
bre_ud <- UDStack(bre_dbbmm)
summary(bre_ud)

#now plot the UD on the left and the actual movement path on the right
#I can't figure out how to change the map area such that the map area is zoomed in, but whatever
par(mfrow=c(1,2))
plot(bre_dbbmm_UD, xlab="longitude", ylab="latitude")
#zoom(bre_dbbmm_UD)
plot(bre_dbbmm_UD, xlab="longitude", ylab="latitude")
lines(bre_bursted_trans, center=TRUE, col=3, lwd=2)
#plot(mama_dbbmm, xlab="location_long", ylab="location_lat")
#points(spTransform(mama_bursted, center=TRUE), col=8)

# this plots the 50% and 95% UD contour lines around the UDs
plot(bre_dbbmm_UD, xlab="location_long", ylab="location_lat")
contour(bre_dbbmm_UD, levels=c(.5, .95), col=c(6,2), add=TRUE, lwd=2)

#get the area of the 95% UD - i think these areas are in meters
bre_cont95 <- getVolumeUD(bre_dbbmm_UD)
bre_cont95 <- bre_cont95<=.95
area95 <- sum(values(bre_cont95))
area95

#get the area of the 50% UD - i think these areas are in meters
bre_cont5 <- getVolumeUD(bre_dbbmm_UD)
bre_cont5 <- bre_cont5<=.5
area5 <- sum(values(bre_cont5))
area5

## Ok, now let's follow Brian's steps and convert the DBBMM object to a SpatialLineDataFrame

#dbbmm dataframe- keep this!
bre.dbbmm.df <- as.data.frame(bre_dbbmm_UD, xy = TRUE)
bre_ud_raster <- rasterFromXYZ(bre.dbbmm.df, crs = "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0", digits = 5 )

## write raster - this is returning an empty raster and we don't know why
#writeRaster(bre_ud_raster, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster.tif", overwrite = TRUE)

#writeRaster(bre_dbbmm, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster2.tif", overwrite = TRUE)

#writeRaster(bre_ud, "C:/Users/sskalos/Documents/noha-move-hab/Output/lauren_ud_raster3.tif", overwrite = TRUE)

writeRaster(bre_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/bre_ud_raster.tif", overwrite = TRUE)

#but it does create the raster appropriately and plots below (weird)
plot(bre_ud_raster)
str(bre_ud_raster) # and there are values within this new raster (extra weird)

#save contours to shapefiles - this works

#example from Bart on Movebank - it works
#require(move)
#example(brownian.bridge.dyn)
#cont<-raster2contour(dbbmm, level=c(.5,.95))
#writeOGR(cont, dsn = '.', layer = 'mycontLines', driver = "ESRI Shapefile")

## try with mama data - it works
require(move)
#cont2 <-raster2contour(bre_ud, level=c(.5,.95))
#writeOGR(cont2, dsn = '.', layer = 'bre_contour2', driver = "ESRI Shapefile")

#works with either the mama_ud raster (above) or the mama_dbbmm_UD raster below, but below seems to be correct as it contains one layer, not every dbbmm step
require(move)
cont_bre <-raster2contour(bre_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_bre, dsn = '.', layer = 'bre_contour', driver = "ESRI Shapefile")

#getwd()
#ploygon <- readOGR(dsn = "/Users/Shannon/Desktop/R_Forever/Dissertation/noha-move-hab/Output" ,layer = "mama_contour3")
#mama_WGS <- spTransform(cont3, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) 
#writeOGR(cont3, dsn = "mama_contours.kml", layer = "mama_WGS", driver = "KML", overwrite_layer = TRUE)

#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_bre <- SpatialPoints(bre.dbbmm.df[,1:2], proj4string = crs(Suisun_nlcd_trans_bre))
nlcd_extract_bre <- extract(nlcd_utm, nlcd_sp_bre)
head(nlcd_extract_bre)
nlcd_extract_bre[which(!is.na(nlcd_extract_bre))]

#check to make sure the have the same # of columns and rows
str(Suisun_nlcd_trans_bre)
str(bre_dbbmm_UD)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(bre_dbbmm_UD)
library(scales)
plot(nlcd_utm, col = alpha("red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final_bre <- cbind.data.frame(bre.dbbmm.df, nlcd_extract_bre)
head(final_bre)

#above works, but returns all columns, including empty grid cells with NA and 0 values
# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_bre <- final_bre[which(!is.na(final_bre[,4])),]
head(final_bre)

# for loop to calculate probabilities of use within each landcover types using the UDs
prob.vec <- rep(NA, length(unique(final_bre[,4])))
unique.vec <- unique(final_bre[,4])
tot.prob <- sum(final_bre[,3])
for (i in 1:length(prob.vec)){
  prob.vec[i] <- sum(final_bre[which(final_bre[,4] == unique.vec[i]),3])/tot.prob
}
#check that the for loop worked and the probabilities sum to 1 - they do
sum(prob.vec)

#save the probability table for each landcover class - it works!
probs.cover.tables <- cbind(prob.vec, unique.vec)

#view the entire table
probs.cover.tables
write.csv(probs.cover.tables, file = "bre_landcover_probs_final.csv")


