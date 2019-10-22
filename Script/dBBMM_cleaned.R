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
#mama_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/SIMP 02.csv")
mama_move <- move(x = "C:/Users/sskalos/Documents/noha-move-hab/Data/SIMP 02.csv")
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
str(mama_bursted)

# bring in raster UTM with NAD 83 projection from ArcMAP because above code isn't working

nlcd_utm <- raster("C:/Users/sskalos/Documents/noha-move-hab/Data/Raster_UTM/NLCD_UTM.tif")
plot(nlcd_utm)
str(nlcd_utm)
as.data.frame(nlcd_utm, xy = TRUE)

#current projection - this is not working - something happens in this proejction transformation and the landcover values are lost for some reason.

Suisun_nlcd_trans <- raster(nlcd_utm)
crs(Suisun_nlcd_trans) <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-121.910674 +lat_0=38.115666 "
str(nlcd_utm)
str(Suisun_nlcd_trans)

plot(Suisun_nlcd_trans)
plot(Suisun_NLCD)

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
mama_ud_raster <- rasterFromXYZ(dbbmm.df, crs = "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0", digits = 5 )

## write raster - this is returning an empty raster and we don't know why
writeRaster(mama_ud_raster, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster.tif", overwrite = TRUE)

writeRaster(mama_dbbmm, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster2.tif", overwrite = TRUE)

writeRaster(mama_ud, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster3.tif", overwrite = TRUE)


writeRaster(mama_dbbmm_UD, "C:/Users/sskalos/Documents/noha-move-hab/Output/mama_ud_raster4.tif", overwrite = TRUE)

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

#orks with either the mama_ud raster (above) or the mama_dbbmm_UD raster below
require(move)
cont3 <-raster2contour(mama_dbbmm_UD, level=c(.5,.95))
writeOGR(cont3, dsn = '.', layer = 'mama_contour3', driver = "ESRI Shapefile")

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
#str(Suisun_nlcd_trans)
#str(mama_dbbmm_UD)

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

