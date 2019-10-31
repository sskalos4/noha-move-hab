## Attempting to create dBBMM utilization distributions##

## start with importing NLCD raster layer using help from: https://cran.r-project.org/web/packages/FedData/FedData.pdf 

## first we need to create an extent for a raster to which the NLCD layer will be cropped to when downloaded (because we can't just import the entire NLCD layer for the entire US)

getwd()
#setwd("/Users/Shannon/Desktop/R_Forever/Dissertation")

install.packages('FedData')
library(FedData)

suisun_polygon <- polygon_from_extent(raster::extent(578701, 601868, 4215762, 4234962), proj4string='+proj=utm +datum=WGS84 +zone=10 +ellps=WGS84')

##download the NLCD raster and clip to the Suisun polygon (note: can only download 2011 with this function, not 2016)

Suisun_NLCD <- get_nlcd(template = suisun_polygon, label = 'suisun5',  year = 2011, dataset = "landcover")

plot(Suisun_NLCD)
summary(Suisun_NLCD)
#View(Suisun_NLCD)
#print(Suisun_NLCD)

## ok let's try dBBMM with Mama breeding locations from 2018 

install.packages("move")
library(move)

#bring in file from Movebank
mama_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/SIMP 02.csv") 
#mama_move <- move(x = "C:/Users/sskalos/Documents/noha-move-hab/Data/SIMP 02.csv")
show(mama_move) 
n.locs(mama_move)  # number of locations
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
summary(nlcd_utm)
#current projection - this is not working - something happens in this proejction transformation and the landcover values are lost for some reason.

Suisun_nlcd_trans <- raster(nlcd_utm)
crs(Suisun_nlcd_trans) <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-121.910674 +lat_0=38.115666 "
str(nlcd_utm)
str(Suisun_nlcd_trans)

plot(Suisun_nlcd_trans)


#run the dbbmm function - NOTE only odd numbers for margin and window size will produce valid likelihood estimations for motion variance

# From Kranstauber et al: increasing the size of the window increases reliability of the motion variance estimation at the cost of missing short term changes in the variation parameter, and should be close to, but less than, 24. Increasing the size of the margin, in contrast, enhances the power to identify "weak" breakpoints at the cost of not detecting breakpoints within the margin. 

mama_dbbmm <- brownian.bridge.dyn(mama_bursted_trans, burstType = 'normal', raster = Suisun_nlcd_trans, location.error = 10, ext = .3, time.step = 60, margin = 3, window.size = 7) #location error is 10 m as per the transmitter specifications, extent is 30% of raster extent, time step is 60 mins becasue locations were approximately every hour, margin is 3 which is the minimum number of locations needed to calculate breakpoints a a leave-one-out approach, and window size is 7 because this is equivalent to 7 locations, which equals  7 hours and may be able to detect behavioral changes within this relatively short window.

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
#dbbmm.df[which(dbbmm.df[,3] == 0), 3] <- NA # this removes the zeros from the raster file
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
cont5 <-raster2contour(mama_dbbmm_UD, level=c(.5,.95))
writeOGR(cont5, dsn = '.', layer = 'mama_contour5', driver = "ESRI Shapefile")

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

head(laureen.dbbmm.df)

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
## NOHA 626 Megan
library(move)
Meg_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 626.csv")
show(Meg_move)
n.locs(Meg_move) # number of locations
head(timeLag(Meg_move, units="mins")) # time difference between locations - this is misleading because her first few days were set to two hour locations, but were then changd to 1 hour locations; make sure to look at ALL time differences
head(timestamps(Meg_move))
summary(timeLag(Meg_move, units="mins"))
Meg_lag <- timeLag(Meg_move, units="mins")
write.csv(Meg_lag, file = "Meg_lag.csv")

# burst the movestack object to exclude any loactions that are greater than 32.3 minutes apart (because a several locations are 32.3 mins and not 30 mis) - this is to prevent calculations of bridges and motion variance overnight between the last location of the previous day and the first location of the next morning, which are typically 400+ mins

Meg_bursted <- move::burst(Meg_move, c('normal','long')[1+(timeLag(Meg_move, units='mins')>32.3)])

#plot Meg's locations
#par(mfcol=1:2)
plot(Meg_bursted, type="o", col=3, lwd=2, pch=20, xlab="location_long",ylab="location_lat")

#plot Meg's locations with ggmap over map layer just to see if it is geographically correct! - it is
install.packages("ggmap")
library(ggmap)
require(ggmap) #these packages are necessary to work with google maps
#require(mapproj)
Meg_df <- as(Meg_bursted, "data.frame")
m <- get_map(bbox(extent(Meg_bursted)*1.1), source="stamen", zoom=12)
ggmap(m)+geom_path(data=Meg_df, aes(x=location.long, y=location.lat))

# transform coordinates from lat lon, center = T is required for the dbbmm to operate properly according to Bart on the movebank help chat


Meg_bursted_trans <- spTransform(x = Meg_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575', center = T)
proj4string(Meg_bursted_trans)

str(Meg_bursted)

# bring in raster UTM with NAD 83 projection from ArcMAP because above code isn't working

nlcd_utm_Meg <- raster("~/Desktop/R_Forever/RRF/Data/Raster_UTM/NLCD_UTM.tif")
#plot(nlcd_utm)
#str(nlcd_utm)
#as.data.frame(nlcd_utm, xy = TRUE)

#current projection - this is not working - something happens in this proejction transformation and the landcover values are lost for some reason.
str(nlcd_utm_Meg)
str(Meg_bursted_trans)

Suisun_nlcd_trans_Meg <- raster(nlcd_utm_Meg)
crs(Suisun_nlcd_trans_Meg) <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575 "

#now they match
str(Meg_bursted_trans)
str(Suisun_nlcd_trans_Meg)

Meg_dbbmm <- brownian.bridge.dyn(Meg_bursted_trans, burstType = 'normal', raster = Suisun_nlcd_trans_Meg, location.error = 10, ext = .3, time.step = 30, margin = 3, window.size = 13) #location error is 10 m as per the transmitter specifications, extent is 30% of raster extent, time step is 30 mins becasue locations were approximately every hour, margin is 3 which is the minimum number of locations needed to calculate Megakpoints a a leave-one-out approach, and window size is 13 because this is equivalent to 13 locations, which equals  ~7 hours and may be able to detect behavioral changes within this relatively short window.

## below are the UDs calculated from the dbbmm
Meg_dbbmm_UD<-new(".UD",calc(Meg_dbbmm, sum)) ## it works!!!
head(Meg_dbbmm_UD)
str(Meg_dbbmm_UD)
summary(Meg_dbbmm_UD)

## get the UD raster layer?? 
#Meg_ud <- UDStack(Meg_dbbmm)
#summary(Meg_ud)

#now plot the UD on the left and the actual movement path on the right
#I can't figure out how to change the map area such that the map area is zoomed in, but whatever
par(mfrow=c(1,2))
plot(Meg_dbbmm_UD, xlab="longitude", ylab="latitude")
#zoom(Meg_dbbmm_UD)
plot(Meg_dbbmm_UD, xlab="longitude", ylab="latitude")
lines(Meg_bursted_trans, center=TRUE, col=3, lwd=2)
#plot(mama_dbbmm, xlab="location_long", ylab="location_lat")
#points(spTransform(mama_bursted, center=TRUE), col=8)

# this plots the 50% and 95% UD contour lines around the UDs
plot(Meg_dbbmm_UD, xlab="location_long", ylab="location_lat")
contour(Meg_dbbmm_UD, levels=c(.5, .95), col=c(6,2), add=TRUE, lwd=2)

#get the area of the 95% UD - i think these areas are in meters
Meg_cont95 <- getVolumeUD(Meg_dbbmm_UD)
Meg_cont95 <- Meg_cont95<=.95
area95 <- sum(values(Meg_cont95))
area95

#get the area of the 50% UD - i think these areas are in meters
Meg_cont5 <- getVolumeUD(Meg_dbbmm_UD)
Meg_cont5 <- Meg_cont5<=.5
area5 <- sum(values(Meg_cont5))
area5

## Ok, now let's follow Brian's steps and convert the DBBMM object to a SpatialLineDataFrame

#dbbmm dataframe- keep this!
Meg.dbbmm.df <- as.data.frame(Meg_dbbmm_UD, xy = TRUE)
Meg_ud_raster <- rasterFromXYZ(Meg.dbbmm.df, crs = "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0", digits = 5 )

## write raster - this is returning an empty raster and we don't know why
#writeRaster(Meg_ud_raster, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster.tif", overwrite = TRUE)

#writeRaster(Meg_dbbmm, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster2.tif", overwrite = TRUE)

#writeRaster(Meg_ud, "C:/Users/sskalos/Documents/noha-move-hab/Output/lauren_ud_raster3.tif", overwrite = TRUE)

writeRaster(Meg_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/Meg_ud_raster.tif", overwrite = TRUE)

#but it does create the raster appropriately and plots below (weird)
plot(Meg_ud_raster)
str(Meg_ud_raster) # and there are values within this new raster (extra weird)

#save contours to shapefiles - this works

#example from Bart on Movebank - it works
#require(move)
#example(brownian.bridge.dyn)
#cont<-raster2contour(dbbmm, level=c(.5,.95))
#writeOGR(cont, dsn = '.', layer = 'mycontLines', driver = "ESRI Shapefile")

## try with mama data - it works
require(move)
#cont2 <-raster2contour(Meg_ud, level=c(.5,.95))
#writeOGR(cont2, dsn = '.', layer = 'Meg_contour2', driver = "ESRI Shapefile")

#works with either the mama_ud raster (above) or the mama_dbbmm_UD raster below, but below seems to be correct as it contains one layer, not every dbbmm step
require(move)
cont_Meg <-raster2contour(Meg_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_Meg, dsn = '.', layer = 'Meg_contour', driver = "ESRI Shapefile")

#getwd()
#ploygon <- readOGR(dsn = "/Users/Shannon/Desktop/R_Forever/Dissertation/noha-move-hab/Output" ,layer = "mama_contour3")
#mama_WGS <- spTransform(cont3, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) 
#writeOGR(cont3, dsn = "mama_contours.kml", layer = "mama_WGS", driver = "KML", overwrite_layer = TRUE)

#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_Meg <- SpatialPoints(Meg.dbbmm.df[,1:2], proj4string = crs(Suisun_nlcd_trans_Meg))
nlcd_extract_Meg <- extract(nlcd_utm_Meg, nlcd_sp_Meg)
head(nlcd_extract_Meg)
nlcd_extract_Meg[which(!is.na(nlcd_extract_Meg))]

#check to make sure the have the same # of columns and rows
str(Suisun_nlcd_trans_Meg)
str(Meg_dbbmm_UD)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(Meg_dbbmm_UD)
library(scales)
plot(nlcd_utm_Meg, col = alpha("red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final_Meg <- cbind.data.frame(Meg.dbbmm.df, nlcd_extract_Meg)
head(final_Meg)

#above works, but returns all columns, including empty grid cells with NA and 0 values
# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_Meg <- final_Meg[which(!is.na(final_Meg[,4])),]
head(final_Meg)

# for loop to calculate probabilities of use within each landcover types using the UDs
prob.vec <- rep(NA, length(unique(final_Meg[,4])))
unique.vec <- unique(final_Meg[,4])
tot.prob <- sum(final_Meg[,3])
for (i in 1:length(prob.vec)){
  prob.vec[i] <- sum(final_Meg[which(final_Meg[,4] == unique.vec[i]),3])/tot.prob
}
#check that the for loop worked and the probabilities sum to 1 - they do
sum(prob.vec)

#save the probability table for each landcover class - it works!
probs.cover.tables <- cbind(prob.vec, unique.vec)

#view the entire table
probs.cover.tables
write.csv(probs.cover.tables, file = "Meg_landcover_probs_final.csv")

## Now let's do the above with the other females
## NOHA 627 Red

library(move)
#note: there were 10 records that were duplicate timestamps with no data because they were SMS messages. If you removed them manually from the .csv file and save, this changes the timestamp column permenantly such that the dates are lost and the move function cannot read the file in. You MUST go to movebank and manually mark each row as an outlier, then redownload the file without these outliers first

getDuplicatedTimestamps(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 627.csv")

Red_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 627.csv")

show(Red_move)
n.locs(Red_move) # number of locations
head(timeLag(Red_move, units="mins")) # time difference between locations - this is misleading because her first few days were set to two hour locations, but were then changd to 1 hour locations; make sure to look at ALL time differences
head(timestamps(Red_move))
summary(timeLag(Red_move, units="mins"))
Red_lag <- timeLag(Red_move, units="mins")
write.csv(Red_lag, file = "Red_lag.csv")

# burst the movestack object to exclude any loactions that are greater than 32.3 minutes apart (because a several locations are 31.4 mins and not 30 mis) - this is to prevent calculations of bridges and motion variance overnight between the last location of the previous day and the first location of the next morning, which are typically 400+ mins

Red_bursted <- move::burst(Red_move, c('normal','long')[1+(timeLag(Red_move, units='mins')>31.4)])

#plot Red's locations
#par(mfcol=1:2)
plot(Red_bursted, type="o", col=3, lwd=2, pch=20, xlab="location_long",ylab="location_lat")

#plot Red's locations with ggmap over map layer just to see if it is geographically correct! - it is
install.packages("ggmap")
library(ggmap)
require(ggmap) #these packages are necessary to work with google maps
#require(mapproj)
Red_df <- as(Red_bursted, "data.frame")
m <- get_map(bbox(extent(Red_bursted)*1.1), source="stamen", zoom=12)
ggmap(m)+geom_path(data=Red_df, aes(x=location.long, y=location.lat))

# transform coordinates from lat lon, center = T is required for the dbbmm to operate properly according to Bart on the movebank help chat


Red_bursted_trans <- spTransform(x = Red_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575', center = T)
proj4string(Red_bursted_trans)

str(Red_bursted)

# bring in raster UTM with NAD 83 projection from ArcMAP because above code isn't working

nlcd_utm_Red <- raster("~/Desktop/R_Forever/RRF/Data/Raster_UTM/NLCD_UTM.tif")
proj4string(nlcd_utm_Red)
#plot(nlcd_utm)
#str(nlcd_utm)
#as.data.frame(nlcd_utm, xy = TRUE)

#current projection - this is not working - something happens in this proejction transformation and the landcover values are lost for some reason.
str(nlcd_utm_Red)
str(Red_bursted_trans)

Suisun_nlcd_trans_Red <- raster(nlcd_utm_Red)
proj4string(Suisun_nlcd_trans_Red) <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575 "

#now they match
str(Red_bursted_trans)
str(Suisun_nlcd_trans_Red)

Red_dbbmm <- brownian.bridge.dyn(Red_bursted_trans, burstType = 'normal', raster = Suisun_nlcd_trans_Red, location.error = 10, ext = .3, time.step = 30, margin = 3, window.size = 13) #location error is 10 m as per the transmitter specifications, extent is 30% of raster extent, time step is 30 mins becasue locations were approximately every hour, margin is 3 which is the minimum number of locations needed to calculate breakpoints with a leave-one-out approach, and window size is 13 (must be odd) because this is equivalent to 14 locations, which equals ~ 7 hours (30 min locations) and may be able to detect behavioral changes within this relatively short window.

## below are the UDs calculated from the dbbmm
Red_dbbmm_UD<-new(".UD",calc(Red_dbbmm, sum)) ## it works!!!
head(Red_dbbmm_UD)
str(Red_dbbmm_UD)
summary(Red_dbbmm_UD)

## get the UD raster layer?? 
#Red_ud <- UDStack(Red_dbbmm)
#summary(Red_ud)

#now plot the UD on the left and the actual movement path on the right
#I can't figure out how to change the map area such that the map area is zoomed in, but whatever
par(mfrow=c(1,2))
plot(Red_dbbmm_UD, xlab="longitude", ylab="latitude")
#zoom(Red_dbbmm_UD)
plot(Red_dbbmm_UD, xlab="longitude", ylab="latitude")
lines(Red_bursted_trans, center=TRUE, col=3, lwd=2)
#plot(mama_dbbmm, xlab="location_long", ylab="location_lat")
#points(spTransform(mama_bursted, center=TRUE), col=8)

# this plots the 50% and 95% UD contour lines around the UDs
plot(Red_dbbmm_UD, xlab="location_long", ylab="location_lat")
contour(Red_dbbmm_UD, levels=c(.5, .95), col=c(6,2), add=TRUE, lwd=2)

#get the area of the 95% UD - i think these areas are in meters
Red_cont95 <- getVolumeUD(Red_dbbmm_UD)
Red_cont95 <- Red_cont95<=.95
area95 <- sum(values(Red_cont95))
area95

#get the area of the 50% UD - i think these areas are in meters
Red_cont5 <- getVolumeUD(Red_dbbmm_UD)
Red_cont5 <- Red_cont5<=.5
area5 <- sum(values(Red_cont5))
area5

## Ok, now let's follow Brian's steps and convert the DBBMM object to a SpatialLineDataFrame

#dbbmm dataframe- keep this!
Red.dbbmm.df <- as.data.frame(Red_dbbmm_UD, xy = TRUE)
Red_ud_raster <- rasterFromXYZ(Red.dbbmm.df, crs = "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0", digits = 5 )

## write raster - this is returning an empty raster and we don't know why
#writeRaster(Red_ud_raster, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster.tif", overwrite = TRUE)

#writeRaster(Red_dbbmm, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster2.tif", overwrite = TRUE)

#writeRaster(Red_ud, "C:/Users/sskalos/Documents/noha-move-hab/Output/lauren_ud_raster3.tif", overwrite = TRUE)

writeRaster(Red_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/Red_ud_raster.tif", overwrite = TRUE)

#but it does create the raster appropriately and plots below (weird)
plot(Red_ud_raster)
str(Red_ud_raster) # and there are values within this new raster (extra weird)

#save contours to shapefiles - this works

#example from Bart on Movebank - it works
#require(move)
#example(brownian.bridge.dyn)
#cont<-raster2contour(dbbmm, level=c(.5,.95))
#writeOGR(cont, dsn = '.', layer = 'mycontLines', driver = "ESRI Shapefile")

## try with mama data - it works
require(move)
#cont2 <-raster2contour(Red_ud, level=c(.5,.95))
#writeOGR(cont2, dsn = '.', layer = 'Red_contour2', driver = "ESRI Shapefile")

#works with either the mama_ud raster (above) or the mama_dbbmm_UD raster below, but below seems to be correct as it contains one layer, not every dbbmm step
require(move)
cont_Red <-raster2contour(Red_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_Red, dsn = '.', layer = 'Red_contour', driver = "ESRI Shapefile")

#getwd()
#ploygon <- readOGR(dsn = "/Users/Shannon/Desktop/R_Forever/Dissertation/noha-move-hab/Output" ,layer = "mama_contour3")
#mama_WGS <- spTransform(cont3, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) 
#writeOGR(cont3, dsn = "mama_contours.kml", layer = "mama_WGS", driver = "KML", overwrite_layer = TRUE)

#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_Red <- SpatialPoints(Red.dbbmm.df[,1:2], proj4string = crs(Suisun_nlcd_trans_Red))
nlcd_extract_Red <- extract(nlcd_utm_Red, nlcd_sp_Red)
head(nlcd_extract_Red)
nlcd_extract_Red[which(!is.na(nlcd_extract_Red))]

#check to make sure the have the same # of columns and rows
str(Suisun_nlcd_trans_Red)
str(Red_dbbmm_UD)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(Red_dbbmm_UD)
library(scales)
plot(nlcd_utm_Red, col = alpha("red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final_Red <- cbind.data.frame(Red.dbbmm.df, nlcd_extract_Red)
head(final_Red)

#above works, but returns all columns, including empty grid cells with NA and 0 values
# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_Red <- final_Red[which(!is.na(final_Red[,4])),]
head(final_Red)

# for loop to calculate probabilities of use within each landcover types using the UDs
prob.vec <- rep(NA, length(unique(final_Red[,4])))
unique.vec <- unique(final_Red[,4])
tot.prob <- sum(final_Red[,3])
for (i in 1:length(prob.vec)){
  prob.vec[i] <- sum(final_Red[which(final_Red[,4] == unique.vec[i]),3])/tot.prob
}
#check that the for loop worked and the probabilities sum to 1 - they do
sum(prob.vec)

#save the probability table for each landcover class - it works!
probs.cover.tables <- cbind(prob.vec, unique.vec)

#view the entire table
probs.cover.tables
write.csv(probs.cover.tables, file = "Red_landcover_probs_final.csv")

## Now let's do the above with the other females
## NOHA 628 Salty

library(move)

### Ok now let's calculate dBBMM for NOHA 628

#getDuplicatedTimestamps(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 628.csv")

Salty_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 628.csv")

show(Salty_move)
n.locs(Salty_move) # number of locations
head(timeLag(Salty_move, units="mins")) # time difference between locations - this is misleading because her first few days were set to two hour locations, but were then changd to 1 hour locations; make sure to look at ALL time differences
head(timestamps(Salty_move))
summary(timeLag(Salty_move, units="mins"))
Salty_lag <- timeLag(Salty_move, units="mins")
write.csv(Salty_lag, file = "Salty_lag.csv")

# burst the movestack object to exclude any loactions that are greater than 32.3 minutes apart (because a several locations are 31.4 mins and not 30 mis) - this is to prevent calculations of bridges and motion variance overnight between the last location of the previous day and the first location of the next morning, which are typically 400+ mins

Salty_bursted <- move::burst(Salty_move, c('normal','long')[1+(timeLag(Salty_move, units='mins')>31.4)])

#plot Salty's locations
#par(mfcol=1:2)
plot(Salty_bursted, type="o", col=3, lwd=2, pch=20, xlab="location_long",ylab="location_lat")

#plot Salty's locations with ggmap over map layer just to see if it is geographically correct! - it is
install.packages("ggmap")
library(ggmap)
require(ggmap) #these packages are necessary to work with google maps
#require(mapproj)
Salty_df <- as(Salty_bursted, "data.frame")
m <- get_map(bbox(extent(Salty_bursted)*1.1), source="stamen", zoom=12)
ggmap(m)+geom_path(data=Salty_df, aes(x=location.long, y=location.lat))

# transform coordinates from lat lon, center = T is requiSalty for the dbbmm to operate properly according to Bart on the movebank help chat


Salty_bursted_trans <- spTransform(x = Salty_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575', center = T)
proj4string(Salty_bursted_trans)

str(Salty_bursted)

# bring in raster UTM with NAD 83 projection from ArcMAP because above code isn't working

nlcd_utm_Salty <- raster("~/Desktop/R_Forever/RRF/Data/Raster_UTM/NLCD_UTM.tif")
#plot(nlcd_utm)
#str(nlcd_utm)
#as.data.frame(nlcd_utm, xy = TRUE)

#current projection - this is not working - something happens in this proejction transformation and the landcover values are lost for some reason.
str(nlcd_utm_Salty)
str(Salty_bursted_trans)

Suisun_nlcd_trans_Salty <- raster(nlcd_utm_Salty)
crs(Suisun_nlcd_trans_Salty) <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575 "

#now they match
str(Salty_bursted_trans)
str(Suisun_nlcd_trans_Salty)

Salty_dbbmm <- brownian.bridge.dyn(Salty_bursted_trans, burstType = 'normal', raster = Suisun_nlcd_trans_Salty, location.error = 10, ext = .3, time.step = 30, margin = 3, window.size = 13) #location error is 10 m as per the transmitter specifications, extent is 30% of raster extent, time step is 30 mins becasue locations were approximately every hour, margin is 3 which is the minimum number of locations needed to calculate breakpoints with a leave-one-out approach, and window size is 13 (must be odd) because this is equivalent to 14 locations, which equals ~ 7 hours (30 min locations) and may be able to detect behavioral changes within this relatively short window.

## below are the UDs calculated from the dbbmm
Salty_dbbmm_UD<-new(".UD",calc(Salty_dbbmm, sum)) ## it works!!!
head(Salty_dbbmm_UD)
str(Salty_dbbmm_UD)
summary(Salty_dbbmm_UD)

## get the UD raster layer?? 
#Salty_ud <- UDStack(Salty_dbbmm)
#summary(Salty_ud)

#now plot the UD on the left and the actual movement path on the right
#I can't figure out how to change the map area such that the map area is zoomed in, but whatever
par(mfrow=c(1,2))
plot(Salty_dbbmm_UD, xlab="longitude", ylab="latitude")
#zoom(Salty_dbbmm_UD)
plot(Salty_dbbmm_UD, xlab="longitude", ylab="latitude")
lines(Salty_bursted_trans, center=TRUE, col=3, lwd=2)
#plot(mama_dbbmm, xlab="location_long", ylab="location_lat")
#points(spTransform(mama_bursted, center=TRUE), col=8)

# this plots the 50% and 95% UD contour lines around the UDs
plot(Salty_dbbmm_UD, xlab="location_long", ylab="location_lat")
contour(Salty_dbbmm_UD, levels=c(.5, .95), col=c(6,2), add=TRUE, lwd=2)

#get the area of the 95% UD - i think these areas are in meters
Salty_cont95 <- getVolumeUD(Salty_dbbmm_UD)
Salty_cont95 <- Salty_cont95<=.95
area95 <- sum(values(Salty_cont95))
area95

#get the area of the 50% UD - i think these areas are in meters
Salty_cont5 <- getVolumeUD(Salty_dbbmm_UD)
Salty_cont5 <- Salty_cont5<=.5
area5 <- sum(values(Salty_cont5))
area5

## Ok, now let's follow Brian's steps and convert the DBBMM object to a SpatialLineDataFrame

#dbbmm dataframe- keep this!
Salty.dbbmm.df <- as.data.frame(Salty_dbbmm_UD, xy = TRUE)
Salty_ud_raster <- rasterFromXYZ(Salty.dbbmm.df, crs = "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0", digits = 5 )

## write raster - this is returning an empty raster and we don't know why
#writeRaster(Salty_ud_raster, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster.tif", overwrite = TRUE)

#writeRaster(Salty_dbbmm, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster2.tif", overwrite = TRUE)

#writeRaster(Salty_ud, "C:/Users/sskalos/Documents/noha-move-hab/Output/lauren_ud_raster3.tif", overwrite = TRUE)

writeRaster(Salty_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/Salty_ud_raster.tif", overwrite = TRUE)

#but it does create the raster appropriately and plots below (weird)
plot(Salty_ud_raster)
str(Salty_ud_raster) # and there are values within this new raster (extra weird)

#save contours to shapefiles - this works

#example from Bart on Movebank - it works
#require(move)
#example(brownian.bridge.dyn)
#cont<-raster2contour(dbbmm, level=c(.5,.95))
#writeOGR(cont, dsn = '.', layer = 'mycontLines', driver = "ESRI Shapefile")

## try with mama data - it works
#require(move)
#cont2 <-raster2contour(Salty_ud, level=c(.5,.95))
#writeOGR(cont2, dsn = '.', layer = 'Salty_contour2', driver = "ESRI Shapefile")

#works with either the mama_ud raster (above) or the mama_dbbmm_UD raster below, but below seems to be correct as it contains one layer, not every dbbmm step
#require(move)
cont_Salty <-raster2contour(Salty_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_Salty, dsn = '.', layer = 'Salty_contour', driver = "ESRI Shapefile", overwrite_layer = TRUE)

#getwd()
#ploygon <- readOGR(dsn = "/Users/Shannon/Desktop/R_Forever/Dissertation/noha-move-hab/Output" ,layer = "mama_contour3")
#mama_WGS <- spTransform(cont3, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) 
#writeOGR(cont3, dsn = "mama_contours.kml", layer = "mama_WGS", driver = "KML", overwrite_layer = TRUE)

#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_Salty <- SpatialPoints(Salty.dbbmm.df[,1:2], proj4string = crs(Suisun_nlcd_trans_Salty))
nlcd_extract_Salty <- extract(nlcd_utm_Salty, nlcd_sp_Salty)
head(nlcd_extract_Salty)
nlcd_extract_Salty[which(!is.na(nlcd_extract_Salty))]

#check to make sure the have the same # of columns and rows
str(Suisun_nlcd_trans_Salty)
str(Salty_dbbmm_UD)

# test to make sure it works - it does (Salty square represents the nlcd raster layer)
plot(Salty_dbbmm_UD)
library(scales)
plot(nlcd_utm_Salty, col = alpha("Red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final_Salty <- cbind.data.frame(Salty.dbbmm.df, nlcd_extract_Salty)
head(final_Salty)

#above works, but returns all columns, including empty grid cells with NA and 0 values
# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_Salty <- final_Salty[which(!is.na(final_Salty[,4])),]
head(final_Salty)

# for loop to calculate probabilities of use within each landcover types using the UDs
prob.vec <- rep(NA, length(unique(final_Salty[,4])))
unique.vec <- unique(final_Salty[,4])
tot.prob <- sum(final_Salty[,3])
for (i in 1:length(prob.vec)){
  prob.vec[i] <- sum(final_Salty[which(final_Salty[,4] == unique.vec[i]),3])/tot.prob
}
#check that the for loop worked and the probabilities sum to 1 - they do
sum(prob.vec)

#save the probability table for each landcover class - it works!
probs.cover.tables <- cbind(prob.vec, unique.vec)

#view the entire table
probs.cover.tables
write.csv(probs.cover.tables, file = "Salty_landcover_probs_final.csv")

## Now let's do the above with the other females
## NOHA 629 Marlana

library(move)

Marlana_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 629.csv")
show(Marlana_move)
n.locs(Marlana_move) # number of locations
head(timeLag(Marlana_move, units="mins")) # time difference between locations - this is misleading because her first few days were set to two hour locations, but were then changd to 1 hour locations; make sure to look at ALL time differences
head(timestamps(Marlana_move))
summary(timeLag(Marlana_move, units="mins"))
Marlana_lag <- timeLag(Marlana_move, units="mins")
write.csv(Marlana_lag, file = "Marlana_lag.csv")

# burst the movestack object to exclude any loactions that are greater than 32.3 minutes apart (because a several locations are 31.5 mins and not 30 mis) - this is to prevent calculations of bridges and motion variance overnight between the last location of the previous day and the first location of the next morning, which are typically 400+ mins

Marlana_bursted <- move::burst(Marlana_move, c('normal','long')[1+(timeLag(Marlana_move, units='mins')>31.5)])

#plot Marlana's locations
#par(mfcol=1:2)
plot(Marlana_bursted, type="o", col=3, lwd=2, pch=20, xlab="location_long",ylab="location_lat")

#plot Marlana's locations with ggmap over map layer just to see if it is geographically correct! - it is
install.packages("ggmap")
library(ggmap)
require(ggmap) #these packages are necessary to work with google maps
#require(mapproj)
Marlana_df <- as(Marlana_bursted, "data.frame")
m <- get_map(bbox(extent(Marlana_bursted)*1.1), source="stamen", zoom=12)
ggmap(m)+geom_path(data=Marlana_df, aes(x=location.long, y=location.lat))

# transform coordinates from lat lon, center = T is requiMarlana for the dbbmm to operate properly according to Bart on the movebank help chat


Marlana_bursted_trans <- spTransform(x = Marlana_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575', center = T)
proj4string(Marlana_bursted_trans)

str(Marlana_bursted)

# bring in raster UTM with NAD 83 projection from ArcMAP because above code isn't working

nlcd_utm_Marlana <- raster("~/Desktop/R_Forever/RRF/Data/Raster_UTM/NLCD_UTM.tif")
#plot(nlcd_utm)
#str(nlcd_utm)
#as.data.frame(nlcd_utm, xy = TRUE)

#current projection - this is not working - something happens in this proejction transformation and the landcover values are lost for some reason.
str(nlcd_utm_Marlana)
str(Marlana_bursted_trans)

Suisun_nlcd_trans_Marlana <- raster(nlcd_utm_Marlana)
crs(Suisun_nlcd_trans_Marlana) <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575 "

#now they match
str(Marlana_bursted_trans)
str(Suisun_nlcd_trans_Marlana)

Marlana_dbbmm <- brownian.bridge.dyn(Marlana_bursted_trans, burstType = 'normal', raster = Suisun_nlcd_trans_Marlana, location.error = 10, ext = .3, time.step = 30, margin = 3, window.size = 13) #location error is 10 m as per the transmitter specifications, extent is 30% of raster extent, time step is 30 mins becasue locations were approximately every hour, margin is 3 which is the minimum number of locations needed to calculate breakpoints with a leave-one-out approach, and window size is 13 (must be odd) because this is equivalent to 14 locations, which equals ~ 7 hours (30 min locations) and may be able to detect behavioral changes within this relatively short window.

## below are the UDs calculated from the dbbmm
Marlana_dbbmm_UD<-new(".UD",calc(Marlana_dbbmm, sum)) ## it works!!!
head(Marlana_dbbmm_UD)
str(Marlana_dbbmm_UD)
summary(Marlana_dbbmm_UD)

## get the UD raster layer?? 
#Marlana_ud <- UDStack(Marlana_dbbmm)
#summary(Marlana_ud)

#now plot the UD on the left and the actual movement path on the right
#I can't figure out how to change the map area such that the map area is zoomed in, but whatever
par(mfrow=c(1,2))
plot(Marlana_dbbmm_UD, xlab="longitude", ylab="latitude")
#zoom(Marlana_dbbmm_UD)
plot(Marlana_dbbmm_UD, xlab="longitude", ylab="latitude")
lines(Marlana_bursted_trans, center=TRUE, col=3, lwd=2)
#plot(mama_dbbmm, xlab="location_long", ylab="location_lat")
#points(spTransform(mama_bursted, center=TRUE), col=8)

# this plots the 50% and 95% UD contour lines around the UDs
plot(Marlana_dbbmm_UD, xlab="location_long", ylab="location_lat")
contour(Marlana_dbbmm_UD, levels=c(.5, .95), col=c(6,2), add=TRUE, lwd=2)

#get the area of the 95% UD - i think these areas are in meters
Marlana_cont95 <- getVolumeUD(Marlana_dbbmm_UD)
Marlana_cont95 <- Marlana_cont95<=.95
area95 <- sum(values(Marlana_cont95))
area95

#get the area of the 50% UD - i think these areas are in meters
Marlana_cont5 <- getVolumeUD(Marlana_dbbmm_UD)
Marlana_cont5 <- Marlana_cont5<=.5
area5 <- sum(values(Marlana_cont5))
area5

## Ok, now let's follow Brian's steps and convert the DBBMM object to a SpatialLineDataFrame

#dbbmm dataframe- keep this!
Marlana.dbbmm.df <- as.data.frame(Marlana_dbbmm_UD, xy = TRUE)
Marlana_ud_raster <- rasterFromXYZ(Marlana.dbbmm.df, crs = "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0", digits = 5 )

## write raster - this is returning an empty raster and we don't know why
#writeRaster(Marlana_ud_raster, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster.tif", overwrite = TRUE)

#writeRaster(Marlana_dbbmm, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster2.tif", overwrite = TRUE)

#writeRaster(Marlana_ud, "C:/Users/sskalos/Documents/noha-move-hab/Output/lauren_ud_raster3.tif", overwrite = TRUE)

writeRaster(Marlana_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/Marlana_ud_raster.tif", overwrite = TRUE)

#but it does create the raster appropriately and plots below (weird)
plot(Marlana_ud_raster)
str(Marlana_ud_raster) # and there are values within this new raster (extra weird)

#save contours to shapefiles - this works

#example from Bart on Movebank - it works
#require(move)
#example(brownian.bridge.dyn)
#cont<-raster2contour(dbbmm, level=c(.5,.95))
#writeOGR(cont, dsn = '.', layer = 'mycontLines', driver = "ESRI Shapefile")

## try with mama data - it works
#require(move)
#cont2 <-raster2contour(Marlana_ud, level=c(.5,.95))
#writeOGR(cont2, dsn = '.', layer = 'Marlana_contour2', driver = "ESRI Shapefile")

#works with either the mama_ud raster (above) or the mama_dbbmm_UD raster below, but below seems to be correct as it contains one layer, not every dbbmm step
#require(move)
cont_Marlana <-raster2contour(Marlana_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_Marlana, dsn = '.', layer = 'Marlana_contour', driver = "ESRI Shapefile", overwrite_layer = TRUE)

#getwd()
#ploygon <- readOGR(dsn = "/Users/Shannon/Desktop/R_Forever/Dissertation/noha-move-hab/Output" ,layer = "mama_contour3")
#mama_WGS <- spTransform(cont3, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) 
#writeOGR(cont3, dsn = "mama_contours.kml", layer = "mama_WGS", driver = "KML", overwrite_layer = TRUE)

#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_Marlana <- SpatialPoints(Marlana.dbbmm.df[,1:2], proj4string = crs(Suisun_nlcd_trans_Marlana))
nlcd_extract_Marlana <- extract(nlcd_utm_Marlana, nlcd_sp_Marlana)
head(nlcd_extract_Marlana)
nlcd_extract_Marlana[which(!is.na(nlcd_extract_Marlana))]

#check to make sure the have the same # of columns and rows
str(Suisun_nlcd_trans_Marlana)
str(Marlana_dbbmm_UD)

# test to make sure it works - it does (Marlana square represents the nlcd raster layer)
plot(Marlana_dbbmm_UD)
library(scales)
plot(nlcd_utm_Marlana, col = alpha("Red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final_Marlana <- cbind.data.frame(Marlana.dbbmm.df, nlcd_extract_Marlana)
head(final_Marlana)

#above works, but returns all columns, including empty grid cells with NA and 0 values
# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_Marlana <- final_Marlana[which(!is.na(final_Marlana[,4])),]
head(final_Marlana)

# for loop to calculate probabilities of use within each landcover types using the UDs
prob.vec <- rep(NA, length(unique(final_Marlana[,4])))
unique.vec <- unique(final_Marlana[,4])
tot.prob <- sum(final_Marlana[,3])
for (i in 1:length(prob.vec)){
  prob.vec[i] <- sum(final_Marlana[which(final_Marlana[,4] == unique.vec[i]),3])/tot.prob
}
#check that the for loop worked and the probabilities sum to 1 - they do
sum(prob.vec)

#save the probability table for each landcover class - it works!
probs.cover.tables <- cbind(prob.vec, unique.vec)

#view the entire table
probs.cover.tables
write.csv(probs.cover.tables, file = "Marlana_landcover_probs_final.csv")

## Now let's do the above with the other females
## NOHA 630 Jelly

library(move)
Jelly_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 630.csv")
show(Jelly_move)
n.locs(Jelly_move) # number of locations
head(timeLag(Jelly_move, units="mins")) # time difference between locations - this is misleading because her first few days were set to two hour locations, but were then changd to 1 hour locations; make sure to look at ALL time differences
head(timestamps(Jelly_move))
summary(timeLag(Jelly_move, units="mins"))
Jelly_lag <- timeLag(Jelly_move, units="mins")
write.csv(Jelly_lag, file = "Jelly_lag.csv")

# burst the movestack object to exclude any loactions that are greater than 32.3 minutes apart (because a several locations are 35 mins and not 30 mis) - this is to prevent calculations of bridges and motion variance overnight between the last location of the previous day and the first location of the next morning, which are typically 400+ mins

Jelly_bursted <- move::burst(Jelly_move, c('normal','long')[1+(timeLag(Jelly_move, units='mins')>35)])

#plot Jelly's locations
#par(mfcol=1:2)
plot(Jelly_bursted, type="o", col=3, lwd=2, pch=20, xlab="location_long",ylab="location_lat")

#plot Jelly's locations with ggmap over map layer just to see if it is geographically correct! - it is
#install.packages("ggmap")
library(ggmap)
require(ggmap) #these packages are necessary to work with google maps
#require(mapproj)
Jelly_df <- as(Jelly_bursted, "data.frame")
m <- get_map(bbox(extent(Jelly_bursted)*1.1), source="stamen", zoom=12)
ggmap(m)+geom_path(data=Jelly_df, aes(x=location.long, y=location.lat))

# transform coordinates from lat lon, center = T is requiJelly for the dbbmm to operate properly according to Bart on the movebank help chat


Jelly_bursted_trans <- spTransform(x = Jelly_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575', center = T)
proj4string(Jelly_bursted_trans)

str(Jelly_bursted)

# bring in raster UTM with NAD 83 projection from ArcMAP because above code isn't working

nlcd_utm_Jelly <- raster("~/Desktop/R_Forever/RRF/Data/Raster_UTM/NLCD_UTM.tif")
#plot(nlcd_utm)
#str(nlcd_utm)
#as.data.frame(nlcd_utm, xy = TRUE)

#current projection - this is not working - something happens in this proejction transformation and the landcover values are lost for some reason.
str(nlcd_utm_Jelly)
str(Jelly_bursted_trans)

Suisun_nlcd_trans_Jelly <- raster(nlcd_utm_Jelly)
crs(Suisun_nlcd_trans_Jelly) <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575 "

#now they match
str(Jelly_bursted_trans)
str(Suisun_nlcd_trans_Jelly)

Jelly_dbbmm <- brownian.bridge.dyn(Jelly_bursted_trans, burstType = 'normal', raster = Suisun_nlcd_trans_Jelly, location.error = 10, ext = .3, time.step = 30, margin = 3, window.size = 13) #location error is 10 m as per the transmitter specifications, extent is 30% of raster extent, time step is 30 mins becasue locations were approximately every hour, margin is 3 which is the minimum number of locations needed to calculate breakpoints with a leave-one-out approach, and window size is 13 (must be odd) because this is equivalent to 14 locations, which equals ~ 7 hours (30 min locations) and may be able to detect behavioral changes within this relatively short window.

## below are the UDs calculated from the dbbmm
Jelly_dbbmm_UD<-new(".UD",calc(Jelly_dbbmm, sum)) ## it works!!!
head(Jelly_dbbmm_UD)
str(Jelly_dbbmm_UD)
summary(Jelly_dbbmm_UD)

## get the UD raster layer?? 
#Jelly_ud <- UDStack(Jelly_dbbmm)
#summary(Jelly_ud)

#now plot the UD on the left and the actual movement path on the right
#I can't figure out how to change the map area such that the map area is zoomed in, but whatever
par(mfrow=c(1,2))
plot(Jelly_dbbmm_UD, xlab="longitude", ylab="latitude")
#zoom(Jelly_dbbmm_UD)
plot(Jelly_dbbmm_UD, xlab="longitude", ylab="latitude")
lines(Jelly_bursted_trans, center=TRUE, col=3, lwd=2)
#plot(mama_dbbmm, xlab="location_long", ylab="location_lat")
#points(spTransform(mama_bursted, center=TRUE), col=8)

# this plots the 50% and 95% UD contour lines around the UDs
plot(Jelly_dbbmm_UD, xlab="location_long", ylab="location_lat")
contour(Jelly_dbbmm_UD, levels=c(.5, .95), col=c(6,2), add=TRUE, lwd=2)

#get the area of the 95% UD - i think these areas are in meters
Jelly_cont95 <- getVolumeUD(Jelly_dbbmm_UD)
Jelly_cont95 <- Jelly_cont95<=.95
area95 <- sum(values(Jelly_cont95))
area95

#get the area of the 50% UD - i think these areas are in meters
Jelly_cont5 <- getVolumeUD(Jelly_dbbmm_UD)
Jelly_cont5 <- Jelly_cont5<=.5
area5 <- sum(values(Jelly_cont5))
area5

## Ok, now let's follow Brian's steps and convert the DBBMM object to a SpatialLineDataFrame

#dbbmm dataframe- keep this!
Jelly.dbbmm.df <- as.data.frame(Jelly_dbbmm_UD, xy = TRUE)
Jelly_ud_raster <- rasterFromXYZ(Jelly.dbbmm.df, crs = "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0", digits = 5 )

## write raster - this is returning an empty raster and we don't know why
#writeRaster(Jelly_ud_raster, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster.tif", overwrite = TRUE)

#writeRaster(Jelly_dbbmm, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster2.tif", overwrite = TRUE)

#writeRaster(Jelly_ud, "C:/Users/sskalos/Documents/noha-move-hab/Output/lauren_ud_raster3.tif", overwrite = TRUE)

writeRaster(Jelly_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/Jelly_ud_raster.tif", overwrite = TRUE)

#but it does create the raster appropriately and plots below (weird)
plot(Jelly_ud_raster)
str(Jelly_ud_raster) # and there are values within this new raster (extra weird)

#save contours to shapefiles - this works

#example from Bart on Movebank - it works
#require(move)
#example(brownian.bridge.dyn)
#cont<-raster2contour(dbbmm, level=c(.5,.95))
#writeOGR(cont, dsn = '.', layer = 'mycontLines', driver = "ESRI Shapefile")

## try with mama data - it works
#require(move)
#cont2 <-raster2contour(Jelly_ud, level=c(.5,.95))
#writeOGR(cont2, dsn = '.', layer = 'Jelly_contour2', driver = "ESRI Shapefile")

#works with either the mama_ud raster (above) or the mama_dbbmm_UD raster below, but below seems to be correct as it contains one layer, not every dbbmm step
#require(move)
cont_Jelly <-raster2contour(Jelly_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_Jelly, dsn = '.', layer = 'Jelly_contour', driver = "ESRI Shapefile", overwrite_layer = TRUE)

#getwd()
#ploygon <- readOGR(dsn = "/Users/Shannon/Desktop/R_Forever/Dissertation/noha-move-hab/Output" ,layer = "mama_contour3")
#mama_WGS <- spTransform(cont3, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) 
#writeOGR(cont3, dsn = "mama_contours.kml", layer = "mama_WGS", driver = "KML", overwrite_layer = TRUE)

#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_Jelly <- SpatialPoints(Jelly.dbbmm.df[,1:2], proj4string = crs(Suisun_nlcd_trans_Jelly))
nlcd_extract_Jelly <- extract(nlcd_utm_Jelly, nlcd_sp_Jelly)
head(nlcd_extract_Jelly)
nlcd_extract_Jelly[which(!is.na(nlcd_extract_Jelly))]

#check to make sure the have the same # of columns and rows
str(Suisun_nlcd_trans_Jelly)
str(Jelly_dbbmm_UD)

# test to make sure it works - it does (Jelly square represents the nlcd raster layer)
plot(Jelly_dbbmm_UD)
library(scales)
plot(nlcd_utm_Jelly, col = alpha("Red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final_Jelly <- cbind.data.frame(Jelly.dbbmm.df, nlcd_extract_Jelly)
head(final_Jelly)

#above works, but returns all columns, including empty grid cells with NA and 0 values
# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_Jelly <- final_Jelly[which(!is.na(final_Jelly[,4])),]
head(final_Jelly)

# for loop to calculate probabilities of use within each landcover types using the UDs
prob.vec <- rep(NA, length(unique(final_Jelly[,4])))
unique.vec <- unique(final_Jelly[,4])
tot.prob <- sum(final_Jelly[,3])
for (i in 1:length(prob.vec)){
  prob.vec[i] <- sum(final_Jelly[which(final_Jelly[,4] == unique.vec[i]),3])/tot.prob
}
#check that the for loop worked and the probabilities sum to 1 - they do
sum(prob.vec)

#save the probability table for each landcover class - it works!
probs.cover.tables <- cbind(prob.vec, unique.vec)

#view the entire table
probs.cover.tables
write.csv(probs.cover.tables, file = "Jelly_landcover_probs_final.csv")

## Now let's do the above with the other females
## NOHA 30.1 Bre

library(move)
Bre_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 30.1.csv")
show(Bre_move)
n.locs(Bre_move) # number of locations
head(timeLag(Bre_move, units="mins")) # time difference between locations - this is misleading because her first few days were set to two hour locations, but were then changd to 1 hour locations; make sure to look at ALL time differences
head(timestamps(Bre_move))
summary(timeLag(Bre_move, units="mins"))
Bre_lag <- timeLag(Bre_move, units="mins")
write.csv(Bre_lag, file = "Bre_lag.csv")

# burst the movestack object to exclude any loactions that are greater than 32.3 minutes apart (because a several locations are 31.5 mins and not 30 mis) - this is to prevent calculations of bridges and motion variance overnight between the last location of the previous day and the first location of the next morning, which are typically 400+ mins

Bre_bursted <- move::burst(Bre_move, c('normal','long')[1+(timeLag(Bre_move, units='mins')>31.5)])

#plot Bre's locations
#par(mfcol=1:2)
plot(Bre_bursted, type="o", col=3, lwd=2, pch=20, xlab="location_long",ylab="location_lat")

#plot Bre's locations with ggmap over map layer just to see if it is geographically correct! - it is
#install.packages("ggmap")
library(ggmap)
require(ggmap) #these packages are necessary to work with google maps
#require(mapproj)
Bre_df <- as(Bre_bursted, "data.frame")
m <- get_map(bbox(extent(Bre_bursted)*1.1), source="stamen", zoom=12)
ggmap(m)+geom_path(data=Bre_df, aes(x=location.long, y=location.lat))

# transform coordinates from lat lon, center = T is requiBre for the dbbmm to operate properly according to Bart on the movebank help chat

str(Bre_bursted)
Bre_bursted_trans <- spTransform(x = Bre_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575', center = T)
proj4string(Bre_bursted_trans)



# bring in raster UTM with NAD 83 projection from ArcMAP because above code isn't working

nlcd_utm_Bre <- raster("~/Desktop/R_Forever/RRF/Data/Raster_UTM/NLCD_UTM.tif")
#plot(nlcd_utm)
#str(nlcd_utm)
#as.data.frame(nlcd_utm, xy = TRUE)

#current projection - this is not working - something happens in this proejction transformation and the landcover values are lost for some reason.
crs(nlcd_utm_Bre)
crs(nlcd_new)
crs(Bre_bursted_trans)

Suisun_nlcd_trans_Bre <- raster(nlcd_utm_Bre)
crs(Suisun_nlcd_trans_Bre) <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575 "

#now they match
str(Bre_bursted_trans)
str(Suisun_nlcd_trans_Bre)
crs(nlcd_utm_Bre)

plot(Bre_bursted_trans, add = TRUE)
plot(nlcd_utm_Bre)

Bre_dbbmm <- brownian.bridge.dyn(Bre_bursted_trans, burstType = 'normal', raster = Suisun_nlcd_trans_Bre, location.error = 10, ext = .3, time.step = 30, margin = 3, window.size = 13) #location error is 10 m as per the transmitter specifications, extent is 30% of raster extent, time step is 30 mins becasue locations were approximately every hour, margin is 3 which is the minimum number of locations needed to calculate breakpoints with a leave-one-out approach, and window size is 13 (must be odd) because this is equivalent to 14 locations, which equals ~ 7 hours (30 min locations) and may be able to detect behavioral changes within this relatively short window.

## below are the UDs calculated from the dbbmm
Bre_dbbmm_UD<-new(".UD",calc(Bre_dbbmm, sum)) ## it works!!!
head(Bre_dbbmm_UD)
str(Bre_dbbmm_UD)
summary(Bre_dbbmm_UD)

## get the UD raster layer?? 
#Bre_ud <- UDStack(Bre_dbbmm)
#summary(Bre_ud)

#now plot the UD on the left and the actual movement path on the right
#I can't figure out how to change the map area such that the map area is zoomed in, but whatever
par(mfrow=c(1,2))
plot(Bre_dbbmm_UD, xlab="longitude", ylab="latitude")
#zoom(Bre_dbbmm_UD)
plot(Bre_dbbmm_UD, xlab="longitude", ylab="latitude")
lines(Bre_bursted_trans, center=TRUE, col=3, lwd=2)
#plot(mama_dbbmm, xlab="location_long", ylab="location_lat")
#points(spTransform(mama_bursted, center=TRUE), col=8)

# this plots the 50% and 95% UD contour lines around the UDs
plot(Bre_dbbmm_UD, xlab="location_long", ylab="location_lat")
contour(Bre_dbbmm_UD, levels=c(.5, .95), col=c(6,2), add=TRUE, lwd=2)

#get the area of the 95% UD - i think these areas are in meters
Bre_cont95 <- getVolumeUD(Bre_dbbmm_UD)
Bre_cont95 <- Bre_cont95<=.95
area95 <- sum(values(Bre_cont95))
area95

#get the area of the 50% UD - i think these areas are in meters
Bre_cont5 <- getVolumeUD(Bre_dbbmm_UD)
Bre_cont5 <- Bre_cont5<=.5
area5 <- sum(values(Bre_cont5))
area5

## Ok, now let's follow Brian's steps and convert the DBBMM object to a SpatialLineDataFrame

#dbbmm dataframe- keep this!
Bre.dbbmm.df <- as.data.frame(Bre_dbbmm_UD, xy = TRUE)
Bre_ud_raster <- rasterFromXYZ(Bre.dbbmm.df, crs = "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0", digits = 5 )

## write raster - this is returning an empty raster and we don't know why
#writeRaster(Bre_ud_raster, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster.tif", overwrite = TRUE)

#writeRaster(Bre_dbbmm, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster2.tif", overwrite = TRUE)

#writeRaster(Bre_ud, "C:/Users/sskalos/Documents/noha-move-hab/Output/lauren_ud_raster3.tif", overwrite = TRUE)

writeRaster(Bre_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/Bre_ud_raster.tif", overwrite = TRUE)

#but it does create the raster appropriately and plots below (weird)
plot(Bre_ud_raster)
str(Bre_ud_raster) # and there are values within this new raster (extra weird)

#save contours to shapefiles - this works

#example from Bart on Movebank - it works
#require(move)
#example(brownian.bridge.dyn)
#cont<-raster2contour(dbbmm, level=c(.5,.95))
#writeOGR(cont, dsn = '.', layer = 'mycontLines', driver = "ESRI Shapefile")

## try with mama data - it works
#require(move)
#cont2 <-raster2contour(Bre_ud, level=c(.5,.95))
#writeOGR(cont2, dsn = '.', layer = 'Bre_contour2', driver = "ESRI Shapefile")

#works with either the mama_ud raster (above) or the mama_dbbmm_UD raster below, but below seems to be correct as it contains one layer, not every dbbmm step
#require(move)
cont_Bre <-raster2contour(Bre_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_Bre, dsn = '.', layer = 'Bre_contour', driver = "ESRI Shapefile", overwrite_layer = TRUE)

#getwd()
#ploygon <- readOGR(dsn = "/Users/Shannon/Desktop/R_Forever/Dissertation/noha-move-hab/Output" ,layer = "mama_contour3")
#mama_WGS <- spTransform(cont3, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) 
#writeOGR(cont3, dsn = "mama_contours.kml", layer = "mama_WGS", driver = "KML", overwrite_layer = TRUE)

#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_Bre <- SpatialPoints(Bre.dbbmm.df[,1:2], proj4string = crs(Suisun_nlcd_trans_Bre))
nlcd_extract_Bre <- extract(nlcd_utm_Bre, nlcd_sp_Bre)
head(nlcd_extract_Bre)
nlcd_extract_Bre[which(!is.na(nlcd_extract_Bre))]

#check to make sure the have the same # of columns and rows
str(Suisun_nlcd_trans_Bre)
str(Bre_dbbmm_UD)

# test to make sure it works - it does (Bre square represents the nlcd raster layer)
plot(Bre_dbbmm_UD)
contour(Bre_dbbmm_UD, levels=c(.5, .95), col=c(6,2), add=TRUE, lwd=2)
library(scales)
plot(nlcd_utm_Bre, col = alpha("Red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final_Bre <- cbind.data.frame(Bre.dbbmm.df, nlcd_extract_Bre)
head(final_Bre)

#above works, but returns all columns, including empty grid cells with NA and 0 values
# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_Bre <- final_Bre[which(!is.na(final_Bre[,4])),]
head(final_Bre)

# for loop to calculate probabilities of use within each landcover types using the UDs
prob.vec <- rep(NA, length(unique(final_Bre[,4])))
unique.vec <- unique(final_Bre[,4])
tot.prob <- sum(final_Bre[,3])
for (i in 1:length(prob.vec)){
  prob.vec[i] <- sum(final_Bre[which(final_Bre[,4] == unique.vec[i]),3])/tot.prob
}
#check that the for loop worked and the probabilities sum to 1 - they do
sum(prob.vec)

#save the probability table for each landcover class - it works!
probs.cover.tables <- cbind(prob.vec, unique.vec)

#view the entire table
probs.cover.tables
write.csv(probs.cover.tables, file = "Bre_landcover_probs_final.csv")

## need to make a new polygon because dBBMM for Bre is outside of the current polygon nlcd_utm

library(FedData)

suisun_polygon_new <- polygon_from_extent(raster::extent(573476, 609771, 4212972, 4233434), proj4string='+proj=utm +datum=WGS84 +zone=10 +ellps=WGS84')


#suisun_polygon_new <- polygon_from_extent(raster::extent(573476, 609771, 4212972, 4233434), proj4string="+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-121.910674 +lat_0=38.115666 ")

#test
#sac_polygon_new <- polygon_from_extent(raster::extent(604115, 633352, 4268919, 4283557), proj4string="+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-121.910674 +lat_0=38.115666 ")

plot(suisun_polygon_new)

##download the NLCD raster and clip to the Suisun polygon (note: can only download 2011 with this function, not 2016)

Suisun_NLCD_new <- get_nlcd(template = suisun_polygon_new, label = 'suisun4',  year = 2011, dataset = "landcover")
summary(Suisun_NLCD_new)
str(Suisun_NLCD_new)

proj4string(Suisun_NLCD_new)
proj4string(suisun_polygon_new)
proj4string(nlcd_utm)
proj4string(Bre_bursted_trans)

Bre_bursted <- move::burst(Bre_move, c('normal','long')[1+(timeLag(Bre_move, units='mins')>31.5)])

#plot Bre's locations
#par(mfcol=1:2)
plot(Bre_bursted, type="o", col=3, lwd=2, pch=20, xlab="location_long",ylab="location_lat")

#plot Bre's locations with ggmap over map layer just to see if it is geographically correct! - it is
#install.packages("ggmap")
library(ggmap)
require(ggmap) #these packages are necessary to work with google maps
#require(mapproj)
Bre_df <- as(Bre_bursted, "data.frame")
m <- get_map(bbox(extent(Bre_bursted)*1.1), source="stamen", zoom=12)
ggmap(m)+geom_path(data=Bre_df, aes(x=location.long, y=location.lat))

# transform coordinates from lat lon, center = T is requiBre for the dbbmm to operate properly according to Bart on the movebank help chat

str(Bre_bursted)
Bre_bursted_trans <- spTransform(x = Bre_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575', center = T)
proj4string(Bre_bursted_trans)

#Sac_NLCD_new <- get_nlcd(template = sac_polygon_new, label = 'sac',  year = 2011, dataset = "landcover")

library(sf)
library(raster)
#r <- raster(suisun_polygon_new)
#r <- setValues(r, 1:ncell(r))
newproj <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575"
#nlcd_new <- st_transform(Suisun_NLCD_new, crs = the_crs)
nlcd_new <- projectRaster(Suisun_NLCD_new, crs = newproj)
proj4string(nlcd_new)
proj4string(Bre_bursted_trans)
str(nlcd_new)
str(Bre_bursted_trans)

#Suisun_NLCD_new <- spTransform(Suisun_NLCD,CRSobj = suisun_polygon)

# transform NLCD raster to same projection as the move object above

#current projection - this is not working
Suisun_nlcd_trans_Bre2 <- raster(Suisun_NLCD_new)
crs(Suisun_nlcd_trans) <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-121.910674 +lat_0=38.115666 "
str(Suisun_nlcd_trans_Bre2)
str(Bre_bursted_trans)
plot(nlcd_new)

Bre_dbbmm <- brownian.bridge.dyn(Bre_bursted_trans, burstType = 'normal', raster = nlcd_new, location.error = 10, ext = .3, time.step = 30, margin = 3, window.size = 13) #location error is 10 m as per the transmitter specifications, extent is 30% of raster extent, time step is 30 mins becasue locations were approximately every hour, margin is 3 which is the minimum number of locations needed to calculate breakpoints with a leave-one-out approach, and window size is 13 (must be odd) because this is equivalent to 14 locations, which equals ~ 7 hours (30 min locations) and may be able to detect behavioral changes within this relatively short window.

## below are the UDs calculated from the dbbmm
Bre_dbbmm_UD<-new(".UD",calc(Bre_dbbmm, sum)) ## it works!!!
head(Bre_dbbmm_UD)
str(Bre_dbbmm_UD)
summary(Bre_dbbmm_UD)

#dbbmm dataframe- keep this!
Bre.dbbmm.df <- as.data.frame(Bre_dbbmm_UD, xy = TRUE)

#using the nlcd_utm raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp_Bre <- SpatialPoints(Bre.dbbmm.df[,1:2], proj4string = crs(nlcd_new))
nlcd_extract_Bre <- extract(Suisun_NLCD_new, nlcd_sp_Bre)
head(nlcd_extract_Bre)
nlcd_extract_Bre[which(!is.na(nlcd_extract_Bre))]

#check to make sure the have the same # of columns and rows
str(nlcd_new)
str(Bre_dbbmm_UD)

# test to make sure it works - it does (Bre square represents the nlcd raster layer)
plot(Bre_dbbmm_UD)
contour(Bre_dbbmm_UD, levels=c(.5, .95), col=c(6,2), add=TRUE, lwd=2)
library(scales)
plot(nlcd_new, col = alpha("Red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final_Bre2 <- cbind.data.frame(Bre.dbbmm.df, nlcd_extract_Bre)
head(final_Bre2)

#above works, but returns all columns, including empty grid cells with NA and 0 values
# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
final_Bre2 <- final_Bre2[which(!is.na(final_Bre2[,4])),]
head(final_Bre2)

# for loop to calculate probabilities of use within each landcover types using the UDs
prob.vec <- rep(NA, length(unique(final_Bre2[,4])))
unique.vec <- unique(final_Bre2[,4])
tot.prob <- sum(final_Bre2[,3])
for (i in 1:length(prob.vec)){
  prob.vec[i] <- sum(final_Bre2[which(final_Bre2[,4] == unique.vec[i]),3])/tot.prob
}
#check that the for loop worked and the probabilities sum to 1 - they do
sum(prob.vec)

install.packages("sdmvspecies")
library(sdmvspecies)
prob.scale <- rescale(prob.vec, to = c(0,1), from = range(prob.vec, na.rm = TRUE))
sum(prob.scale)

#save the probability table for each landcover class - it works!
probs.cover.tables <- cbind(prob.vec, unique.vec)

#view the entire table
probs.cover.tables
write.csv(probs.cover.tables, file = "Bre_landcover_probs_final2.csv")

#get the area of the 95% UD - i think these areas are in meters
Bre_cont95 <- getVolumeUD(Bre_dbbmm_UD)
Bre_cont95 <- Bre_cont95<=.95
area95 <- sum(values(Bre_cont95))
area95

#get the area of the 50% UD - i think these areas are in meters
Bre_cont5 <- getVolumeUD(Bre_dbbmm_UD)
Bre_cont5 <- Bre_cont5<=.5
area5 <- sum(values(Bre_cont5))
area5

