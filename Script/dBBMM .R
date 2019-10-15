## Attempting to create dBBMM utilization distributions##

## start with importing NLCD raster layer using help from: https://cran.r-project.org/web/packages/FedData/FedData.pdf 

## first we need to create an extent for a raster to which the NLCD layer will be cropped to when downloaded (because we can't just import the entire NLCD layer for the entire US)

getwd()
#setwd("/Users/Shannon/Desktop/R_Forever/Dissertation")

install.packages('FedData')
library(FedData)

suisun_polygon <- polygon_from_extent(raster::extent(578701, 601868, 4215762, 4234962), proj4string='+proj=utm +datum=WGS84 +zone=10')

#get_nlcd(template = suisun_polygon, label, year = 2016, dataset = "landcover",
#         raw.dir = "./RAW/NLCD", extraction.dir = paste0("./EXTRACTIONS/",
#         label, "/NLCD"), raster.options = c("COMPRESS=DEFLATE", "ZLEVEL=9",
#         "INTERLEAVE=BAND"), force.redo = F)

Suisun_NLCD <- get_nlcd(template = suisun_polygon, label = 'suisun',  year = 2011, dataset = "landcover")

plot(Suisun_NLCD)
#View(Suisun_NLCD)
print(Suisun_NLCD)

## DON'T DO BELOW CODE INITIALLY - SKIP TO JUST TESTING MAMA DATA BELOW
## import GPS locations and convert lat/long to utm and create a spatial points file

library(sp)
library(rgdal)

breeding_2018 <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/RRF/Data/All_Breeding_2018_Nesting_Period.csv"))
coordinates(breeding_2018) <- c("location.long", "location.lat")
proj4string(breeding_2018) <- CRS("+proj=longlat +datum=WGS84 + zone10")  ## for example

breeding_2018_WGS <- spTransform(breeding_2018, CRS("+proj=utm +zone=10 ellps=WGS84 + zone10" ))
breeding_2018_WGS

## For a SpatialPoints object rather than a SpatialPointsDataFrame, just do: 
breeding_2018_sp <- as(breeding_2018_WGS, "SpatialPoints")
breeding_2018_sp


plot(Suisun_NLCD)
plot(breeding_2018_sp, pch = 19, color = "black", add = TRUE)


## trying to bring in 2016 NLCD image below and it doesn't work
library(raster)
library(sp)
library(rgdal)

#nlcd2016 <- brick("~/Desktop/R_Forever/RRF/Data/NLCD_2016_Land_Cover_L48_20190424.tiff")

nlcd2016 <- stack("~/Desktop/R_Forever/RRF/Data/NLCD_2016_Land_Cover_L48_20190424.tiff")

crs(nlcd2016) <- "+proj=utm +datum=WGS84 +zone=10 +ellps=WGS84 +towgs84=0,0,0 " 
suisun_nlcd_extent <- spTransform(suisun_polygon, CRS(proj4string(nlcd2016)))
suisun_nlcd_crop <- crop(nlcd2016, suisun_polygon) # this doesn't work and I don't know why

nlcd2016
suisun_polygon

plot(suisun_nlcd_crop)

## ok let's try dBBMM with the Mama from 2018 

install.packages("move")
library(move)

#mama <- as.data.frame(read.csv("~/Desktop/R_Forever/RRF/Data/SIMP02NestingPeriod.csv"))

#coordinates(mama) <- c("lon", "lat")  
#proj4string(mama) <- CRS("+proj=longlat +datum=WGS84")
#mama
#plot(mama, pch = 19)

#bring in file from Movebank
mama_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/SIMP 02.csv")
show(mama_move)
n.locs(mama_move) # number of locations
head(timeLag(mama_move, units="mins")) # time difference between locations
head(timestamps(mama_move))

# burst the movestack object to exclude any loactions that are greater than 72 minutes apart (because a couple locations are 72 and not 60 mis) - this is to prevent calculations of bridges and motion variance overnight between the last location of the previous day and the first location of the next morning, which are typically 400+ mins

mama_bursted <- move::burst(mama_move, c('normal','long')[1+(timeLag(mama_move, units='mins')>72)])

#plot mama's locations
#par(mfcol=1:2)
plot(mama_bursted, type="o", col=3, lwd=2, pch=20, xlab="location_long",ylab="location_lat")

#plot mama's locations with ggmap over map layer just to see if it is geographically correct! - it is
require(ggmap) #these packages are necessary to work with google maps
#require(mapproj)
mama_df <- as(mama_bursted, "data.frame")
m <- get_map(bbox(extent(mama_bursted)*1.1), source="stamen", zoom=12)
ggmap(m)+geom_path(data=mama_df, aes(x=location.long, y=location.lat))


# transform coordinates from lat lon, center = T is required for the dbbmm to operate properly according to Bart on the movebank help chat

#mama_trans <- spTransform(mama_move, center=T)
mama_trans2 <- spTransform(x = mama_bursted, CRSobj = "+proj=aeqd + ellps=WGS84", center = T)
proj4string(mama_trans2)

# transform NLCD raster to same projection as the move object above

#current projection
r <- raster(Suisun_NLCD)
crs(r) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

new_crs <- "+proj=aeqd +ellps=WGS84 +lon_0=-121.910674 +lat_0=38.115666"

#transformed projection

Suisun_NLCD_trans <- projectRaster(from = r, crs = new_crs )
proj4string(Suisun_NLCD_trans)

#run the dbbmm function - NOTE only odd numbers for margin and window size will produce valid likelihood estimations for motion variance

# From Kranstauber et al: increasing the size of the window increases reliability of the motion variance estimation at the cost of missing short term changes in the variation parameter, and should be close to, but less than, 24. Increasing the size of the margin, in contrast, enhances the power to identify "weak" breakpoints at the cost of not detecting breakpoints within the margin. 

mama_dbbmm <- brownian.bridge.dyn(mama_trans2, burstType = 'normal', raster = Suisun_NLCD_trans, location.error = 10, ext = .3, time.step = 60, margin = 3, window.size = 7) #location error is 10 m as per the transmitter specifications, extent is 30% of raster extent, time step is 60 mins becasue locations were approximately every hour, margin is 3 which is the minimum number of locations needed to calculate breakpoints a a leave-one-out approach, and window size is 7 because this is equivalent to 7 locations, which equals  7 hours and may be able to detect behavioral changes within this relatively short window.

#try without the raster - doesn't work, don't worry about it because you have a raster set
#mama_dbbmm <- brownian.bridge.dyn(mama_trans2, burstType = 'normal', dimSize = 30, location.error = 10, ext = .3, time.step = 60, margin = 13)

## below are the UDs calculated from the dbbmm
mama_dbbmm_UD<-new(".UD",calc(mama_dbbmm, sum)) ## it works!!!

## get the UD raster layer??
mama_ud <- UDStack(mama_dbbmm)
head(mama_ud)

#now plot the UD on the left and the actual movement path on the right
#I can't figure out how to change the map area such that the map area is zoomed in, but whatever
par(mfrow=c(1,2))
plot(mama_dbbmm_UD, xlab="longitude", ylab="latitude")
plot(mama_dbbmm_UD, xlab="longitude", ylab="latitude")
lines(spTransform(mama_bursted, center=TRUE), col=3, lwd=2)
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

## save the objects above
#save(x=mama_dbbmm, file="~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/mama.RData")
#rm(mama_dbbmm)
#load(file="~/Desktop/R_Forever/Dissertaion/noha-move-hab/Output/mama.RData")

#save(x=mama_dbbmm_UD, file="~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/mama.RData")
#rm(mama_dbbmm_UD)

# saves the contours to a kml file
install.packages('plotKML')
require(plotKML)
kml(mama_move)

## extract the motion variance values from the dbbmm - this defines how diffusive or how irregular the path of an animal is
mama_var <- getMotionVariance(mama_dbbmm) ## with a DBBMM object

# convert a move object to adehabitat compatible object
mama_ade <- move2ade(mama_move)

#segment lengths and speed of track

head(seglength(mama_ade))
head(speed(mama_move))

#convert raster UD to contour to import into ArcMap or ArcGIS

#example from Bart on Movebank - it works
require(move)
example(brownian.bridge.dyn)
cont<-raster2contour(dbbmm, level=c(.5,.95))
writeOGR(cont, dsn = '.', layer = 'mycontLines', driver = "ESRI Shapefile")

## try with mama data - it works
require(move)
cont<-raster2contour(mama_dbbmm, level=c(.5,.95))
writeOGR(cont, dsn = '.', layer = 'mama_contour2', driver = "ESRI Shapefile")

require(move)
cont2 <-raster2contour(mama_ud, level=c(.5,.95))
writeOGR(cont2, dsn = '.', layer = 'mama_contour', driver = "ESRI Shapefile")

install.packages("GISTools")
library(GISTools)
library(rgdal)
mama_points_new <- writeOGR(mama_bursted, dsn = "tempdir", layer = "mama_points_new", driver = "ESRI Shapefile")

# try to save contour as a kml file - does not work
require(move)
cont<-raster2contour(mama_dbbmm_UD, level=c(.5,.95))
writeOGR(cont, dsn = '.', layer = 'mama_contour_kml', driver = "KML")

## troubleshooting the lower y grid not large enough error from above

#I had a similar error message, and changing the extent or dimension did not help. I only succeeded in changing the error message to "Upper x grid not large enough" or "Lower y grid not large enough" or "Upper y grid not large enough." Bart Kranstauber was very helpful, immediately catching that my data had a long delay between two locations (over 40 days); here's his response to my e-mail, including the solution:


#The breaking down of the code most likely relates to changes in the move package, calculations can now be done much quicker for much larger rasters due to some optimizations. I took some time to investigate the problem relates to the fact that you have a large data gap (80 days, more than half the tracking period) for which a large movement variance was estimated because of fast movements before and after. This means that for a long period the location is very uncertain. This segment should probably be omited from the calculation this can be done like this: 

#bursted<-move::burst(porc.move.ae, c('normal','long')[1+(timeLag(porc.move.ae, units='days')>40)])

#porc_dbbmm <- brownian.bridge.dyn(bursted, burstType='normal', dimSize=1500, location.error=20, ext=1.4, time.step=10, margin=13)
#porc_dbbmm<-new(".UD",calc(porc_dbbmm, sum))
#plot(porc_dbbmm)

## now let's us moveud to extract uds around each movement step
## Note don't need this now that I was able to calculate motion variance for each step - not true?

install.packages("remotes")
remotes::install_github("bacollier/moveud")
library(moveud)

## another way to save the UDs

rangeud <- 1104
ts <- 72

mama_ud <- move.forud(mama_dbbmm, range.subset = 10:12, ts = .2, ras = Suisun_NLCD_trans, le = 10, lev = c(50, 95), crs = "+proj=aeqd + ellps=WGS84", path = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output", name = "mama_ud")

