# re-running KDE for wintering females in Suisun, as well as trying to see if dBBMM will work for these individuals, instead

### KDE example with mama data ####
getwd()
library(adehabitatHR)
library(rgdal)
library(maptools)
library(move)

coords <- read.csv("/Users/Shannon/Desktop/R_Forever/Dissertation/noha-move-hab/Data/mama_utms.csv", header = TRUE)
coords.xy <- coords[c("x", "y")]
coords.spatialpoints <- SpatialPoints(coords.xy)
proj4string(coords.spatialpoints) = CRS("+init=epsg:32610")

ud <- kernelUD(coords.spatialpoints, h = "href", grid = 30)
image(ud)

ud <- kernelUD(coords.spatialpoints, h = "href", grid = 30)
image(ud)

ver99 <- getverticeshr(ud, 99, unout = "km2")
ver95 <- getverticeshr(ud, 95, unout = "km2")
ver50 <- getverticeshr(ud, 50, unout = "km2")
plot(ver99, col = "yellow")
plot(ver95, col = "orange", add = TRUE)
plot(ver50, col = "red", add = TRUE)

getvolumeUD(ud)
kernel.area(ud, percent = seq(20, 95, by = 5),
            unin = c("m", "km"),
            unout = c("ha", "km2", "m2"), standardize = FALSE)

ud.df <- as.data.frame.estUD(ud)

ver99kml <- spTransform(ver99, CRS("+proj=longlat +datum=WGS84"))
writeOGR(obj=ver99kml, "UD99.kml" , layer="ver99kml", driver="KML")

ver95kml <- spTransform(ver95, CRS("+proj=longlat +datum=WGS84"))
writeOGR(obj=ver95kml, "UD95.kml" , layer="ver95kml", driver="KML")

ver50kml <- spTransform(ver50, CRS("+proj=longlat +datum=WGS84"))
writeOGR(obj=ver50kml, "UD50.kml" , layer="ver50kml", driver="KML")



### SIMP 2 Post Fledge (July 31st to January 31st) ###
coords2post <- read.csv("/Users/Shannon/Desktop/R_Forever/Dissertation/noha-move-hab/Data/SIMP02_postfledgejanuary.csv", header = TRUE)
coords2post.xy <- coords2post[c("x", "y")]
coords2post.spatialpoints <- SpatialPoints(coords2post.xy)
proj4string(coords2post.spatialpoints) = CRS("+proj=longlat +datum=WGS84")
coords2post.utm <- spTransform(coords2post.spatialpoints, CRS("+init=epsg:32610"))
ud2post <- kernelUD(coords2post.spatialpoints,  h = "href", grid = 30)
image(ud2post)

ver99_2post <- getverticeshr(ud2post, 99, unout = "km2")
ver95_2post <- getverticeshr(ud2post, 95, unout = "km2")
ver50_2post <- getverticeshr(ud2post, 50, unout = "km2")
plot(ver99_2post, col = "yellow")
plot(ver95_2post, col = "orange", add = TRUE)
plot(ver50_2post, col = "red", add = TRUE)

ver99_2post <- spTransform(ver99_2post, CRS("+proj=longlat +datum=WGS84"))
#writeOGR(obj=ver99_2post, "SIMP02post_UD99.kml" , layer="ver99_2post", driver="KML")
writeOGR(obj=ver99_2post, "SIMP02post_UD99.shp" , layer="ver99_2post", driver="ESRI Shapefile")

ver95_2post <- spTransform(ver95_2post, CRS("+proj=longlat +datum=WGS84"))
#writeOGR(obj=ver95kml_2post, "SIMP02post_UD95.kml" , layer="ver95kml_2post", driver="KML")
writeOGR(obj=ver95_2post, "SIMP02post_UD95.shp" , layer="ver95_2post", driver="ESRI Shapefile")

ver50_2post <- spTransform(ver50_2post, CRS("+proj=longlat +datum=WGS84"))
#writeOGR(obj=ver50kml_2post, "SIMP02post_UD50.kml" , layer="ver50kml_2post", driver="KML")
writeOGR(obj=ver50_2post, "SIMP02post_UD50.shp" , layer="ver50_2post", driver="ESRI Shapefile")

