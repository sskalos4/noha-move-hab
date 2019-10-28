#Test Run: create new dbbmm on mama data with new Suisun layer that fits ALL breeding female locations (had to make a new layer for Bre because some locations were slightly outside of the polygon border). If the probs are the same for mama as they were for the original polygon, then I should be able to calculate cumulative probability from all of the UDs (which will all need to be re-run with this new layer).

library(FedData)

suisun_polygon_new <- polygon_from_extent(raster::extent(573476, 609771, 4212972, 4233434), proj4string='+proj=utm +datum=WGS84 +zone=10 +ellps=WGS84')

plot(suisun_polygon_new)

##download the NLCD raster and clip to the Suisun polygon (note: can only download 2011 with this function, not 2016)

Suisun_NLCD_new <- get_nlcd(template = suisun_polygon_new, label = 'suisun_new',  year = 2011, dataset = "landcover")
summary(Suisun_NLCD_new)
str(Suisun_NLCD_new)

library(move)

#bring in file from Movebank
mama_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/SIMP 02.csv") 
mama_bursted <- move::burst(mama_move, c('normal','long')[1+(timeLag(mama_move, units='mins')>72)])

mama_bursted_trans <- spTransform(x = mama_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m', center = T)

proj4string(Suisun_NLCD_new) #this raster is incorrect
proj4string(suisun_polygon_new) # this polygon is incorrect
proj4string(mama_bursted_trans) # need all layers to match this projection

#matching projections below
library(sf)
library(raster)
r <- raster(suisun_polygon_new)
#r <- setValues(r, 1:ncell(r))
newproj <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-121.910674 +lat_0=38.115666"
nlcd_new <- projectRaster(Suisun_NLCD_new, crs = newproj)

#now they match, but still need the Suisun_nlcd_trans_mama layer to match, too
proj4string(nlcd_new)
proj4string(mama_bursted_trans)

#current projection - this is not working and we don't need this layer below anyway because the nlcd_new and bursted_trans files above match already

#Suisun_nlcd_trans_mama <- raster(Suisun_NLCD_new)
#crs(Suisun_nlcd_trans_mama) <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-121.910674 +lat_0=38.115666 "
#plot(Suisun_nlcd_trans_mama)

#now these match, too
#proj4string(Suisun_nlcd_trans_mama)
#proj4string(mama_bursted_trans)

mama_dbbmm <- brownian.bridge.dyn(mama_bursted_trans, burstType = 'normal', raster = nlcd_new, location.error = 10, ext = .3, time.step = 60, margin = 3, window.size = 7) #location error is 10 m as per the transmitter specifications, extent is 30% of raster extent, time step is 60 mins becasue locations were approximately every hour, margin is 3 which is the minimum number of locations needed to calculate breakpoints a a leave-one-out approach, and window size is 7 because this is equivalent to 7 locations, which equals  7 hours and may be able to detect behavioral changes within this relatively short window.

## below are the UDs calculated from the dbbmm
mama_dbbmm_UD<-new(".UD",calc(mama_dbbmm, sum)) ## it works!!!

#get the area of the 95% UD - i think these areas are in meters?? NO these are the number of cells
mama_cont95 <- getVolumeUD(mama_dbbmm_UD)
mama_cont95 <- mama_cont95<=.95
area95 <- sum(values(mama_cont95))
area95

#get the area of the 50% UD - i think these areas are in meters
mama_cont5 <- getVolumeUD(mama_dbbmm_UD)
mama_cont5 <- mama_cont5<=.5
area5 <- sum(values(mama_cont5))
area5

#dbbmm dataframe- keep this!
dbbmm.df <- as.data.frame(mama_dbbmm_UD, xy = TRUE)

#save UD raster
writeRaster(mama_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/mama_ud_raster_new.tif", overwrite = TRUE)

#save contours
cont_new <-raster2contour(mama_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_new, dsn = '.', layer = 'mama_contour_new', driver = "ESRI Shapefile", overwrite_layer =  TRUE)

#using the nlcd_new raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
library(sp)
nlcd_sp <- SpatialPoints(dbbmm.df[,1:2], proj4string = crs(nlcd_new))
nlcd_extract <- extract(Suisun_NLCD_new, nlcd_sp)
head(nlcd_extract)
nlcd_extract[which(!is.na(nlcd_extract))]

#check that the columns and rows match - they do
str(mama_dbbmm_UD)
str(nlcd_new)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(mama_dbbmm_UD)
library(scales)
plot(nlcd_new, col = alpha("red", .5), add = TRUE)

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
write.csv(probs.cover.tables, file = "mama_landcover_probs_final_new.csv")


## Laureen SIMP 1 with new suisun raster
library(move)

#bring in file from Movebank
laureen_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/SIMP 01.csv") 
laureen_bursted <- move::burst(laureen_move, c('normal','long')[1+(timeLag(laureen_move, units='mins')>88)])

laureen_bursted_trans <- spTransform(x = laureen_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m', center = T)

proj4string(Suisun_NLCD_new) #this raster is incorrect
proj4string(suisun_polygon_new) # this polygon is incorrect
proj4string(laureen_bursted_trans) # need all layers to match this projection

#matching projections below
library(sf)
library(raster)
r <- raster(suisun_polygon_new)
#r <- setValues(r, 1:ncell(r))
newproj <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.0374075 +lat_0=38.2021575"
nlcd_new_laureen <- projectRaster(Suisun_NLCD_new, crs = newproj)

#now they match, but still need the Suisun_nlcd_trans_laureen layer to match, too
proj4string(nlcd_new_laureen)
proj4string(laureen_bursted_trans)


laureen_dbbmm <- brownian.bridge.dyn(laureen_bursted_trans, burstType = 'normal', raster = nlcd_new_laureen, location.error = 10, ext = .3, time.step = 60, margin = 3, window.size = 7) 

## below are the UDs calculated from the dbbmm
laureen_dbbmm_UD<-new(".UD",calc(laureen_dbbmm, sum)) ## it works!!!

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

#dbbmm dataframe- keep this!
laureen.dbbmm.df <- as.data.frame(laureen_dbbmm_UD, xy = TRUE)

#save UD raster
writeRaster(laureen_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/laureen_ud_raster_new.tif", overwrite = TRUE)

#save contours
cont_new <-raster2contour(laureen_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_new, dsn = '.', layer = 'laureen_contour_new', driver = "ESRI Shapefile", overwrite_layer =  TRUE)

#using the nlcd_new raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp <- SpatialPoints(laureen.dbbmm.df[,1:2], proj4string = crs(nlcd_new_laureen))
nlcd_extract <- extract(Suisun_NLCD_new, nlcd_sp)
head(nlcd_extract)
nlcd_extract[which(!is.na(nlcd_extract))]

#check that the columns and rows match - they do
str(laureen_dbbmm_UD)
str(nlcd_new_laureen)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(laureen_dbbmm_UD)
library(scales)
plot(nlcd_new, col = alpha("red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final <- cbind.data.frame(laureen.dbbmm.df, nlcd_extract)
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
write.csv(probs.cover.tables, file = "laureen_landcover_probs_final_new.csv")


## Joyce SIMP 7 with new suisun raster
library(move)

#bring in file from Movebank
joyce_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/SIMP 07.csv") 
joyce_bursted <- move::burst(joyce_move, c('normal','long')[1+(timeLag(joyce_move, units='mins')>62)])

joyce_bursted_trans <- spTransform(x = joyce_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m', center = T)

proj4string(Suisun_NLCD_new) #this raster is incorrect
proj4string(suisun_polygon_new) # this polygon is incorrect
proj4string(joyce_bursted_trans) # need all layers to match this projection

#matching projections below
library(sf)
library(raster)
r <- raster(suisun_polygon_new)
#r <- setValues(r, 1:ncell(r))
newproj <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.003483 +lat_0=38.2062325"
nlcd_new_joyce <- projectRaster(Suisun_NLCD_new, crs = newproj)

#now they match, but still need the Suisun_nlcd_trans_joyce layer to match, too
proj4string(nlcd_new_joyce)
proj4string(joyce_bursted_trans)


joyce_dbbmm <- brownian.bridge.dyn(joyce_bursted_trans, burstType = 'normal', raster = nlcd_new_joyce, location.error = 10, ext = .3, time.step = 60, margin = 3, window.size = 7) 

## below are the UDs calculated from the dbbmm
joyce_dbbmm_UD<-new(".UD",calc(joyce_dbbmm, sum)) ## it works!!!

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

#dbbmm dataframe- keep this!
joyce.dbbmm.df <- as.data.frame(joyce_dbbmm_UD, xy = TRUE)

#save UD raster
writeRaster(joyce_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/joyce_ud_raster_new.tif", overwrite = TRUE)

#save contours
cont_new <-raster2contour(joyce_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_new, dsn = '.', layer = 'joyce_contour_new', driver = "ESRI Shapefile", overwrite_layer =  TRUE)

#using the nlcd_new raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp <- SpatialPoints(joyce.dbbmm.df[,1:2], proj4string = crs(nlcd_new_joyce))
nlcd_extract <- extract(Suisun_NLCD_new, nlcd_sp)
head(nlcd_extract)
nlcd_extract[which(!is.na(nlcd_extract))]

#check that the columns and rows match - they do
str(joyce_dbbmm_UD)
str(nlcd_new_joyce)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(joyce_dbbmm_UD)
library(scales)
plot(nlcd_new, col = alpha("red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final <- cbind.data.frame(joyce.dbbmm.df, nlcd_extract)
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
write.csv(probs.cover.tables, file = "joyce_landcover_probs_final_new.csv")

## Cinnamon SIMP 10 with new suisun raster
library(move)

#bring in file from Movebank
cinnamon_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/SIMP 10.csv") 
cinnamon_bursted <- move::burst(cinnamon_move, c('normal','long')[1+(timeLag(cinnamon_move, units='mins')>74)])

cinnamon_bursted_trans <- spTransform(x = cinnamon_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m', center = T)

proj4string(Suisun_NLCD_new) #this raster is incorrect
proj4string(suisun_polygon_new) # this polygon is incorrect
proj4string(cinnamon_bursted_trans) # need all layers to match this projection

#matching projections below
library(sf)
library(raster)
r <- raster(suisun_polygon_new)
#r <- setValues(r, 1:ncell(r))
newproj <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-121.872816 +lat_0=38.1182245"
nlcd_new_cinnamon <- projectRaster(Suisun_NLCD_new, crs = newproj)

#now they match, but still need the Suisun_nlcd_trans_cinnamon layer to match, too
proj4string(nlcd_new_cinnamon)
proj4string(cinnamon_bursted_trans)


cinnamon_dbbmm <- brownian.bridge.dyn(cinnamon_bursted_trans, burstType = 'normal', raster = nlcd_new_cinnamon, location.error = 10, ext = .3, time.step = 60, margin = 3, window.size = 7) 

## below are the UDs calculated from the dbbmm
cinnamon_dbbmm_UD<-new(".UD",calc(cinnamon_dbbmm, sum)) ## it works!!!

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

#dbbmm dataframe- keep this!
cinnamon.dbbmm.df <- as.data.frame(cinnamon_dbbmm_UD, xy = TRUE)

#save UD raster
writeRaster(cinnamon_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/cinnamon_ud_raster_new.tif", overwrite = TRUE)

#save contours
cont_new <-raster2contour(cinnamon_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_new, dsn = '.', layer = 'cinnamon_contour_new', driver = "ESRI Shapefile", overwrite_layer =  TRUE)

#using the nlcd_new raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp <- SpatialPoints(cinnamon.dbbmm.df[,1:2], proj4string = crs(nlcd_new_cinnamon))
nlcd_extract <- extract(Suisun_NLCD_new, nlcd_sp)
head(nlcd_extract)
nlcd_extract[which(!is.na(nlcd_extract))]

#check that the columns and rows match - they do
str(cinnamon_dbbmm_UD)
str(nlcd_new_cinnamon)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(cinnamon_dbbmm_UD)
library(scales)
plot(nlcd_new, col = alpha("red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final <- cbind.data.frame(cinnamon.dbbmm.df, nlcd_extract)
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
write.csv(probs.cover.tables, file = "cinnamon_landcover_probs_final_new.csv")

## Megan NOHA 626 with new suisun raster
library(move)

#bring in file from Movebank
megan_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 626.csv") 
megan_bursted <- move::burst(megan_move, c('normal','long')[1+(timeLag(megan_move, units='mins')>32.3)])

megan_bursted_trans <- spTransform(x = megan_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m', center = T)

proj4string(Suisun_NLCD_new) #this raster is incorrect
proj4string(suisun_polygon_new) # this polygon is incorrect
proj4string(megan_bursted_trans) # need all layers to match this projection

#matching projections below
library(sf)
library(raster)
r <- raster(suisun_polygon_new)
#r <- setValues(r, 1:ncell(r))
newproj <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-121.970436096191 +lat_0=38.1442394256592"
nlcd_new_megan <- projectRaster(Suisun_NLCD_new, crs = newproj)

#now they match, but still need the Suisun_nlcd_trans_megan layer to match, too
proj4string(nlcd_new_megan)
proj4string(megan_bursted_trans)


megan_dbbmm <- brownian.bridge.dyn(megan_bursted_trans, burstType = 'normal', raster = nlcd_new_megan, location.error = 10, ext = .3, time.step = 60, margin = 3, window.size = 7) 

## below are the UDs calculated from the dbbmm
megan_dbbmm_UD<-new(".UD",calc(megan_dbbmm, sum)) ## it works!!!

#get the area of the 95% UD - i think these areas are in meters
megan_cont95 <- getVolumeUD(megan_dbbmm_UD)
megan_cont95 <- megan_cont95<=.95
area95 <- sum(values(megan_cont95))
area95

#get the area of the 50% UD - i think these areas are in meters
megan_cont5 <- getVolumeUD(megan_dbbmm_UD)
megan_cont5 <- megan_cont5<=.5
area5 <- sum(values(megan_cont5))
area5

#dbbmm dataframe- keep this!
megan.dbbmm.df <- as.data.frame(megan_dbbmm_UD, xy = TRUE)

#save UD raster
writeRaster(megan_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/megan_ud_raster_new.tif", overwrite = TRUE)

#save contours
cont_new <-raster2contour(megan_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_new, dsn = '.', layer = 'megan_contour_new', driver = "ESRI Shapefile", overwrite_layer =  TRUE)

#using the nlcd_new raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp <- SpatialPoints(megan.dbbmm.df[,1:2], proj4string = crs(nlcd_new_megan))
nlcd_extract <- extract(Suisun_NLCD_new, nlcd_sp)
head(nlcd_extract)
nlcd_extract[which(!is.na(nlcd_extract))]

#check that the columns and rows match - they do
str(megan_dbbmm_UD)
str(nlcd_new_megan)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(megan_dbbmm_UD)
library(scales)
plot(nlcd_new, col = alpha("red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final <- cbind.data.frame(megan.dbbmm.df, nlcd_extract)
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
write.csv(probs.cover.tables, file = "megan_landcover_probs_final_new.csv")

## Redhouse NOHA 627 with new suisun raster
library(move)

#bring in file from Movebank
Redhouse_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 627.csv") 
Redhouse_bursted <- move::burst(Redhouse_move, c('normal','long')[1+(timeLag(Redhouse_move, units='mins')>31.4)])

Redhouse_bursted_trans <- spTransform(x = Redhouse_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m', center = T)

proj4string(Suisun_NLCD_new) #this raster is incorrect
proj4string(suisun_polygon_new) # this polygon is incorrect
proj4string(Redhouse_bursted_trans) # need all layers to match this projection

#matching projections below
library(sf)
library(raster)
r <- raster(suisun_polygon_new)
#r <- setValues(r, 1:ncell(r))
newproj <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-121.944164276123 +lat_0=38.1145782470703"
nlcd_new_Redhouse <- projectRaster(Suisun_NLCD_new, crs = newproj)

#now they match, but still need the Suisun_nlcd_trans_Redhouse layer to match, too
proj4string(nlcd_new_Redhouse)
proj4string(Redhouse_bursted_trans)


Redhouse_dbbmm <- brownian.bridge.dyn(Redhouse_bursted_trans, burstType = 'normal', raster = nlcd_new_Redhouse, location.error = 10, ext = .3, time.step = 60, margin = 3, window.size = 7) 

## below are the UDs calculated from the dbbmm
Redhouse_dbbmm_UD<-new(".UD",calc(Redhouse_dbbmm, sum)) ## it works!!!

#get the area of the 95% UD - i think these areas are in meters
Redhouse_cont95 <- getVolumeUD(Redhouse_dbbmm_UD)
Redhouse_cont95 <- Redhouse_cont95<=.95
area95 <- sum(values(Redhouse_cont95))
area95

#get the area of the 50% UD - i think these areas are in meters
Redhouse_cont5 <- getVolumeUD(Redhouse_dbbmm_UD)
Redhouse_cont5 <- Redhouse_cont5<=.5
area5 <- sum(values(Redhouse_cont5))
area5

#dbbmm dataframe- keep this!
Redhouse.dbbmm.df <- as.data.frame(Redhouse_dbbmm_UD, xy = TRUE)

#save UD raster
writeRaster(Redhouse_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/Redhouse_ud_raster_new.tif", overwrite = TRUE)

#save contours
cont_new <-raster2contour(Redhouse_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_new, dsn = '.', layer = 'Redhouse_contour_new', driver = "ESRI Shapefile", overwrite_layer =  TRUE)

#using the nlcd_new raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp <- SpatialPoints(Redhouse.dbbmm.df[,1:2], proj4string = crs(nlcd_new_Redhouse))
nlcd_extract <- extract(Suisun_NLCD_new, nlcd_sp)
head(nlcd_extract)
nlcd_extract[which(!is.na(nlcd_extract))]

#check that the columns and rows match - they do
str(Redhouse_dbbmm_UD)
str(nlcd_new_Redhouse)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(Redhouse_dbbmm_UD)
library(scales)
plot(nlcd_new, col = alpha("red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final <- cbind.data.frame(Redhouse.dbbmm.df, nlcd_extract)
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
write.csv(probs.cover.tables, file = "Redhouse_landcover_probs_final_new.csv")

## Salty NOHA 628 with new suisun raster
library(move)

#bring in file from Movebank
Salty_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 628.csv") 
Salty_bursted <- move::burst(Salty_move, c('normal','long')[1+(timeLag(Salty_move, units='mins')>31.4)])

Salty_bursted_trans <- spTransform(x = Salty_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m', center = T)

proj4string(Suisun_NLCD_new) #this raster is incorrect
proj4string(suisun_polygon_new) # this polygon is incorrect
proj4string(Salty_bursted_trans) # need all layers to match this projection

#matching projections below
library(sf)
library(raster)
r <- raster(suisun_polygon_new)
#r <- setValues(r, 1:ncell(r))
newproj <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-121.990695953369 +lat_0=38.1896324157715"
nlcd_new_Salty <- projectRaster(Suisun_NLCD_new, crs = newproj)

#now they match, but still need the Suisun_nlcd_trans_Salty layer to match, too
proj4string(nlcd_new_Salty)
proj4string(Salty_bursted_trans)


Salty_dbbmm <- brownian.bridge.dyn(Salty_bursted_trans, burstType = 'normal', raster = nlcd_new_Salty, location.error = 10, ext = .3, time.step = 60, margin = 3, window.size = 7) 

## below are the UDs calculated from the dbbmm
Salty_dbbmm_UD<-new(".UD",calc(Salty_dbbmm, sum)) ## it works!!!

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

#dbbmm dataframe- keep this!
Salty.dbbmm.df <- as.data.frame(Salty_dbbmm_UD, xy = TRUE)

#save UD raster
writeRaster(Salty_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/Salty_ud_raster_new.tif", overwrite = TRUE)

#save contours
cont_new <-raster2contour(Salty_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_new, dsn = '.', layer = 'Salty_contour_new', driver = "ESRI Shapefile", overwrite_layer =  TRUE)

#using the nlcd_new raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp <- SpatialPoints(Salty.dbbmm.df[,1:2], proj4string = crs(nlcd_new_Salty))
nlcd_extract <- extract(Suisun_NLCD_new, nlcd_sp)
head(nlcd_extract)
nlcd_extract[which(!is.na(nlcd_extract))]

#check that the columns and rows match - they do
str(Salty_dbbmm_UD)
str(nlcd_new_Salty)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(Salty_dbbmm_UD)
library(scales)
plot(nlcd_new, col = alpha("red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final <- cbind.data.frame(Salty.dbbmm.df, nlcd_extract)
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
write.csv(probs.cover.tables, file = "Salty_landcover_probs_final_new.csv")

## Marlana NOHA 629 with new suisun raster
library(move)

#bring in file from Movebank
Marlana_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 629.csv") 
Marlana_bursted <- move::burst(Marlana_move, c('normal','long')[1+(timeLag(Marlana_move, units='mins')>31.5)])

Marlana_bursted_trans <- spTransform(x = Marlana_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m', center = T)

proj4string(Suisun_NLCD_new) #this raster is incorrect
proj4string(suisun_polygon_new) # this polygon is incorrect
proj4string(Marlana_bursted_trans) # need all layers to match this projection

#matching projections below
library(sf)
library(raster)
r <- raster(suisun_polygon_new)
#r <- setValues(r, 1:ncell(r))
newproj <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-121.940395355225 +lat_0=38.1387825012207"
nlcd_new_Marlana <- projectRaster(Suisun_NLCD_new, crs = newproj)

#now they match, but still need the Suisun_nlcd_trans_Marlana layer to match, too
proj4string(nlcd_new_Marlana)
proj4string(Marlana_bursted_trans)


Marlana_dbbmm <- brownian.bridge.dyn(Marlana_bursted_trans, burstType = 'normal', raster = nlcd_new_Marlana, location.error = 10, ext = .3, time.step = 60, margin = 3, window.size = 7) 

## below are the UDs calculated from the dbbmm
Marlana_dbbmm_UD<-new(".UD",calc(Marlana_dbbmm, sum)) ## it works!!!

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

#dbbmm dataframe- keep this!
Marlana.dbbmm.df <- as.data.frame(Marlana_dbbmm_UD, xy = TRUE)

#save UD raster
writeRaster(Marlana_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/Marlana_ud_raster_new.tif", overwrite = TRUE)

#save contours
cont_new <-raster2contour(Marlana_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_new, dsn = '.', layer = 'Marlana_contour_new', driver = "ESRI Shapefile", overwrite_layer =  TRUE)

#using the nlcd_new raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp <- SpatialPoints(Marlana.dbbmm.df[,1:2], proj4string = crs(nlcd_new_Marlana))
nlcd_extract <- extract(Suisun_NLCD_new, nlcd_sp)
head(nlcd_extract)
nlcd_extract[which(!is.na(nlcd_extract))]

#check that the columns and rows match - they do
str(Marlana_dbbmm_UD)
str(nlcd_new_Marlana)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(Marlana_dbbmm_UD)
library(scales)
plot(nlcd_new, col = alpha("red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final <- cbind.data.frame(Marlana.dbbmm.df, nlcd_extract)
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
write.csv(probs.cover.tables, file = "Marlana_landcover_probs_final_new.csv")

## Jelly NOHA 630 with new suisun raster
library(move)

#bring in file from Movebank
Jelly_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 630.csv") 
Jelly_bursted <- move::burst(Jelly_move, c('normal','long')[1+(timeLag(Jelly_move, units='mins')>35)])

Jelly_bursted_trans <- spTransform(x = Jelly_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m', center = T)

proj4string(Suisun_NLCD_new) #this raster is incorrect
proj4string(suisun_polygon_new) # this polygon is incorrect
proj4string(Jelly_bursted_trans) # need all layers to match this projection

#matching projections below
library(sf)
library(raster)
r <- raster(suisun_polygon_new)
#r <- setValues(r, 1:ncell(r))
newproj <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-122.000186920166 +lat_0=38.1836967468262"
nlcd_new_Jelly <- projectRaster(Suisun_NLCD_new, crs = newproj)

#now they match, but still need the Suisun_nlcd_trans_Jelly layer to match, too
proj4string(nlcd_new_Jelly)
proj4string(Jelly_bursted_trans)


Jelly_dbbmm <- brownian.bridge.dyn(Jelly_bursted_trans, burstType = 'normal', raster = nlcd_new_Jelly, location.error = 10, ext = .3, time.step = 60, margin = 3, window.size = 7) 

## below are the UDs calculated from the dbbmm
Jelly_dbbmm_UD<-new(".UD",calc(Jelly_dbbmm, sum)) ## it works!!!

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

#dbbmm dataframe- keep this!
Jelly.dbbmm.df <- as.data.frame(Jelly_dbbmm_UD, xy = TRUE)

#save UD raster
writeRaster(Jelly_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/Jelly_ud_raster_new.tif", overwrite = TRUE)

#save contours
cont_new <-raster2contour(Jelly_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_new, dsn = '.', layer = 'Jelly_contour_new', driver = "ESRI Shapefile", overwrite_layer =  TRUE)

#using the nlcd_new raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp <- SpatialPoints(Jelly.dbbmm.df[,1:2], proj4string = crs(nlcd_new_Jelly))
nlcd_extract <- extract(Suisun_NLCD_new, nlcd_sp)
head(nlcd_extract)
nlcd_extract[which(!is.na(nlcd_extract))]

#check that the columns and rows match - they do
str(Jelly_dbbmm_UD)
str(nlcd_new_Jelly)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(Jelly_dbbmm_UD)
library(scales)
plot(nlcd_new, col = alpha("red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final <- cbind.data.frame(Jelly.dbbmm.df, nlcd_extract)
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
write.csv(probs.cover.tables, file = "Jelly_landcover_probs_final_new.csv")

## Bre NOHA 30.1 with new suisun raster
library(move)

#bring in file from Movebank
Bre_move <- move(x = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 30.1.csv") 
Bre_bursted <- move::burst(Bre_move, c('normal','long')[1+(timeLag(Bre_move, units='mins')>31.5)])

Bre_bursted_trans <- spTransform(x = Bre_bursted, CRSobj = '+proj=utm +zone=10 +datum=NAD83 +units=m', center = T)

proj4string(Suisun_NLCD_new) #this raster is incorrect
proj4string(suisun_polygon_new) # this polygon is incorrect
proj4string(Bre_bursted_trans) # need all layers to match this projection

#matching projections below
library(sf)
library(raster)
r <- raster(suisun_polygon_new)
#r <- setValues(r, 1:ncell(r))
newproj <- "+proj=utm +zone=10 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0 +lon_0=-121.889591217041 +lat_0=38.088960647583"
nlcd_new_Bre <- projectRaster(Suisun_NLCD_new, crs = newproj)

#now they match, but still need the Suisun_nlcd_trans_Bre layer to match, too
proj4string(nlcd_new_Bre)
proj4string(Bre_bursted_trans)


Bre_dbbmm <- brownian.bridge.dyn(Bre_bursted_trans, burstType = 'normal', raster = nlcd_new_Bre, location.error = 10, ext = .3, time.step = 60, margin = 3, window.size = 7) 

## below are the UDs calculated from the dbbmm
Bre_dbbmm_UD<-new(".UD",calc(Bre_dbbmm, sum)) ## it works!!!

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

#dbbmm dataframe- keep this!
Bre.dbbmm.df <- as.data.frame(Bre_dbbmm_UD, xy = TRUE)

#save UD raster
writeRaster(Bre_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/Bre_ud_raster_new.tif", overwrite = TRUE)

#save contours
cont_new <-raster2contour(Bre_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_new, dsn = '.', layer = 'Bre_contour_new', driver = "ESRI Shapefile", overwrite_layer =  TRUE)

#using the nlcd_new raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
nlcd_sp <- SpatialPoints(Bre.dbbmm.df[,1:2], proj4string = crs(nlcd_new_Bre))
nlcd_extract <- extract(Suisun_NLCD_new, nlcd_sp)
head(nlcd_extract)
nlcd_extract[which(!is.na(nlcd_extract))]

#check that the columns and rows match - they do
str(Bre_dbbmm_UD)
str(nlcd_new_Bre)

# test to make sure it works - it does (red square represents the nlcd raster layer)
plot(Bre_dbbmm_UD)
library(scales)
plot(nlcd_new, col = alpha("red", .5), add = TRUE)

# combine the raster cell probabilities with their coord pairs with landcover grid cells
final <- cbind.data.frame(Bre.dbbmm.df, nlcd_extract)
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
write.csv(probs.cover.tables, file = "Bre_landcover_probs_final_new.csv")


## example data for calculating cumulative raster values and rescaling to 1 - use this to sum all breedig female' UD values to calculate cumulative probability for each landcover class

r <- raster(ncol=10, nrow=10)
r1 <- init(r, fun=runif)
r2 <- init(r, fun=runif)
r3 <- overlay(r1, r2, fun=function(x,y){return(x+y)})
head(r1)
head(r2)
head(r3)

0.43208508+0.08343650
0.1293262+0.5808329
add <- sum(r1[1,10], r2[1,10])
add
summary(r3)
rnew <- rescale(r3)
summary(rnew)

# combine the raster cell probabilities with their coord pairs with landcover grid cells for TWO FEMALES

mama_final <- cbind.data.frame(dbbmm.df, nlcd_extract)
laureen_final <- cbind.data.frame(laureen.dbbmm.df, nlcd_extract)
head(mama_final)
head(laureen_final)

#above works, but returns all columns, including empty grid cells with NA and 0 values
# below code removes NA in the 4th column (the landcover column) and returns only columns with landcover values 
mama_final <- mama_final[which(!is.na(mama_final[,4])),]
laureen_final <- laureen_final[which(!is.na(laureen_final[,4])),]
head(mama_final)
head(laureen_final)
combined <- rbind(mama_final, laureen_final)
head(combined)


#add the probability columns (layer) together for duplicate rows (x,y, and nlcd_extract)
#[,-2] means all columns except the first 2
#add.combine <- aggregate(combined[,-2], list(layer=combined[,3]), FUN = sum)

#let's try sing tidyverse
library(tidyverse)
library(dplyr)
#add.combine <- inner_join(mama_final, laureen_final)
str(mama_final)
str(laureen_final)
str(combined)

combined %>%
  group_by(x,y) %>%
  summarise_all(sum) %>%
  data.frame() -> new.combined # so that new.combined can further be used, if needed

#alternative approach
combined2 <- inner_join(mama_final, laureen_final, by = "x")

newdf
new.combined
# for loop to calculate probabilities of use within each landcover types using the UDs
prob.vec.1 <- rep(NA, length(unique(mama_final[,4])))
unique.vec <- unique(mama_final[,4])
tot.prob <- sum(mama_final[,3])
for (i in 1:length(prob.vec.1)){
  prob.vec.1[i] <- sum(mama_final[which(mama_final[,4] == unique.vec[i]),3])/tot.prob
}

mama_sum <- sum(prob.vec.1)

# for loop to calculate probabilities of use within each landcover types using the UDs
prob.vec.2 <- rep(NA, length(unique(laureen_final[,4])))
unique.vec <- unique(laureen_final[,4])
tot.prob <- sum(laureen_final[,3])
for (i in 1:length(prob.vec.2)){
  prob.vec.2[i] <- sum(laureen_final[which(laureen_final[,4] == unique.vec[i]),3])/tot.prob
}

#check that the for loop worked and the probabilities sum to 1 - they do

laureen_sum <- sum(prob.vec.2)

laureen_sum
mama_sum
prob.vec.1
prob.vec.2

vec.add <- prob.vec.1+prob.vec.2
sum(vec.add)

#rescale probabilities - these aren't working correctly and are normalizing the entire total to 1 such that when similar landcover class probabilities are combined, the sum of the probabilities is greater than 1. We want the cumulative probability across all 15 landcover classes to be 1, not each one scaled such that one equals 1 and all others are less than one - this results in a sum prob >1.

#vec.add.scale <- (1/(prob.vec.1+prob.vec.2))
#vec.cum <- (prob.vec.1*vec.add.scale)+(prob.vec.2*vec.add.scale)
#sum(vec.cum)

#scaled.cum.prob <- function(x){((x-min(x))/(max(x)-min(x)))}
#scaled.cum.prob <- function(x){(for(i in:length(x))/sum(x)}
#cum.prob <- scaled.cum.prob(vec.add)
#sum(cum.prob)

#save the probability table for each landcover class 
probs.cover.tables.1 <- cbind(prob.vec.1, unique.vec)
probs.cover.tables.2 <- cbind(prob.vec.2, unique.vec)
probs.cover.tables.cum <- cbind(cum.prob, unique.vec)

#view the entire table
probs.cover.tables.cum
write.csv(probs.cover.tables.cum, file = "cumulative_landcover_probs_final.csv")
