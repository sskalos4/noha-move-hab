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

#dbbmm dataframe- keep this!
dbbmm.df <- as.data.frame(mama_dbbmm_UD, xy = TRUE)

#save UD raster
writeRaster(mama_dbbmm_UD, "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/mama_ud_raster_new.tif", overwrite = TRUE)

#save contours
cont_new <-raster2contour(mama_dbbmm_UD, level=c(.5,.95))
writeOGR(cont_new, dsn = '.', layer = 'mama_contour_new', driver = "ESRI Shapefile", overwrite_layer =  TRUE)

#using the nlcd_new raster lines up with the correct number of rows and columns from our dbbmm dataframe, and we can extract the landcover values (finally!)
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
add.combine <-

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

#rescale probabilities

#vec.add.scale <- (1/(prob.vec.1+prob.vec.2))
#vec.cum <- (prob.vec.1*vec.add.scale)+(prob.vec.2*vec.add.scale)
#sum(vec.cum)

#scaled.cum.prob <- function(x){((x-min(x))/(max(x)-min(x)))}
#scaled.cum.prob <- function(x){(for(i in:length(x))/sum(x)}
#cum.prob <- scaled.cum.prob(vec.add)
#sum(cum.prob)

#save the probability table for each landcover class - it works!
probs.cover.tables.1 <- cbind(prob.vec.1, unique.vec)
probs.cover.tables.2 <- cbind(prob.vec.2, unique.vec)
probs.cover.tables.cum <- cbind(cum.prob, unique.vec)

#view the entire table
probs.cover.tables.cum
write.csv(probs.cover.tables.cum, file = "cumulative_landcover_probs_final.csv")
