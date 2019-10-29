# migration distances and MCPs

getwd()
library(trajr)

#my data - female NOHA 27

path27 <- read.csv(file = "Data/NOHA_27.csv")
trj27 <- TrajFromCoords(path27)
path27
plot(trj27)

TrajLength(trj27)
TrajKm27 <-(TrajLength(trj27))/1000
TrajKm27

#my data - female NOHA 29


path29 <- read.csv(file = "Data/NOHA_29.csv")
trj29 <- TrajFromCoords(path29)
path29
plot(trj29)

TrajLength(trj29)
TrajKm29 <-(TrajLength(trj29))/1000
TrajKm29

## MCPs for 2018 migrant breeding females
library(adehabitatHR)
library(sp)

noha61.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/DRIB 61.csv"))

coordinates(noha61.df) <- c("x", "y")
proj4string(noha61.df) <- CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")  ## for example
head(noha61.df)

noha61.mcp <-mcp(noha61.df, percent = 100, unin = "m",unout = "km2")
plot(noha61.mcp)
plot(noha61.df, add = TRUE)

mcp.area(noha61.df, percent = 100, unin = "m", unout = "km2")

## FOWL 40
noha40.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/FOWL 40.csv"))

coordinates(noha40.df) <- c("x", "y")
proj4string(noha40.df) <- CRS("+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs")  ## for example
head(noha40.df)

noha40.mcp <-mcp(noha40.df, percent = 100, unin = "m",unout = "km2")
plot(noha40.mcp)
plot(noha40.df, add = TRUE)

mcp.area(noha40.df, percent = 100, unin = "m", unout = "km2")

## FOWL 22

noha22.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/FOWL 22.csv"))

coordinates(noha22.df) <- c("x", "y")
proj4string(noha22.df) <- CRS("+proj=utm +zone=11 +datum=WGS84 +units=m +no_defs")  ## for example
head(noha22.df)

noha22.mcp <-mcp(noha22.df, percent = 100, unin = "m",unout = "km2")
plot(noha22.mcp)
plot(noha22.df, add = TRUE)

mcp.area(noha22.df, percent = 100, unin = "m", unout = "km2")


## FOWL 41

noha41.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/FOWL 41.csv"))

coordinates(noha41.df) <- c("x", "y")
proj4string(noha41.df) <- CRS("+proj=utm +zone=11 +datum=WGS84 +units=m +no_defs")  ## for example
head(noha41.df)

noha41.mcp <-mcp(noha41.df, percent = 100, unin = "m",unout = "km2")
plot(noha41.mcp)
plot(noha41.df, add = TRUE)

mcp.area(noha41.df, percent = 100, unin = "m", unout = "km2")

## FOWL 36

noha36.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/FOWL 36.csv"))

coordinates(noha36.df) <- c("x", "y")
proj4string(noha36.df) <- CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")  ## for example
head(noha36.df)

noha36.mcp <-mcp(noha36.df, percent = 100, unin = "m",unout = "km2")
plot(noha36.mcp)
plot(noha36.df, add = TRUE)

mcp.area(noha36.df, percent = 100, unin = "m", unout = "km2")

## DRIB 65

noha65.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/DRIB 65.csv"))

coordinates(noha65.df) <- c("x", "y")
proj4string(noha65.df) <- CRS("+proj=utm +zone=11 +datum=WGS84 +units=m +no_defs")  ## for example
head(noha65.df)

noha65.mcp <-mcp(noha65.df, percent = 100, unin = "m",unout = "km2")
plot(noha65.mcp)
plot(noha65.df, add = TRUE)

mcp.area(noha65.df, percent = 100, unin = "m", unout = "km2")

## NOHA 23

noha23.df <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/NOHA 23.csv"))

coordinates(noha23.df) <- c("x", "y")
proj4string(noha23.df) <- CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")  ## for example
head(noha23.df)

noha23.mcp <-mcp(noha23.df, percent = 100, unin = "m",unout = "km2")
plot(noha23.mcp)
plot(noha23.df, add = TRUE)

mcp.area(noha23.df, percent = 100, unin = "m", unout = "km2")


library(rgdal)
writeOGR(noha22.mcp, dsn = ".", layer = "NOHA 22 MCP", driver = "ESRI Shapefile", overwrite_layer = TRUE )

writeOGR(noha23.mcp, dsn = ".", layer = "NOHA 23 MCP", driver = "ESRI Shapefile", overwrite_layer = TRUE )

writeOGR(noha40.mcp, dsn = ".", layer = "NOHA 40 MCP", driver = "ESRI Shapefile", overwrite_layer = TRUE )

writeOGR(noha41.mcp, dsn = ".", layer = "NOHA 41 MCP", driver = "ESRI Shapefile", overwrite_layer = TRUE )

writeOGR(noha36.mcp, dsn = ".", layer = "NOHA 36 MCP", driver = "ESRI Shapefile", overwrite_layer = TRUE )

writeOGR(noha61.mcp, dsn = ".", layer = "NOHA 61 MCP", driver = "ESRI Shapefile", overwrite_layer = TRUE )

writeOGR(noha65.mcp, dsn = ".", layer = "NOHA 65 MCP", driver = "ESRI Shapefile", overwrite_layer = TRUE )
