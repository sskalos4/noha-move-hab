#third order hab selection comparing hab at points (used) to hab within mcp home range (available)

library(adehabitatHS)

mcp.hab <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/mcp.csv"))

point.hab <- as.data.frame(read.csv(file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Data/points.csv"))

hab.list <- list(mcp = mcp.hab, point = point.hab)

hab.selc <- compana(hab.list$point[,-c(1:2)], hab.list$mcp[,-c(1:2)], test = "randomisation", rnv = 0.01, nrep = 500, alpha = 0.1)

hab.selc2 <- compana(hab.list$point[,-c(1:2)], hab.list$mcp[,-c(1:2)])

hab.selc2
hab.selc

hab.selc$rank
hab.selc$rm

hab.selc2$rank
hab.selc2$rm

hab.selc$test
hab.selc$type.test
hab.selc$rmnb
## The ranking matrix:
print(hab.selc$rm, quote = FALSE)

##example data from compana function

data("pheasant")
pheasant$studyarea
pheasant$mcp
pheasant$locs
print(pheasant)
str(pheasant)
compana

write.csv(hab.selc$rm)

## Third order habitat selection
## (relocation within home range)
## We remove the second column
## (as in the paper of Aebischer et al.)

data("squirrel")
squirrel
squiana3 <- compana(squirrel$locs[,-2], squirrel$mcp[,-2])
squiana3
## The ranking matrix:
print(squiana3$rm, quote = FALSE)
squiana3$test
squiana3$type.test

