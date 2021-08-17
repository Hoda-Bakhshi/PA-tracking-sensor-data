# The linear interpolation code was developed by:
# @author Jed Long
# @references
# Long, JA (2015) Kinematic interpolation of movement data. \emph{International Journal of Geographical Information Science}. DOI: 10.1080/13658816.2015.1081909. 
# and adapted by Author: Lindsey Conrow 

###############################################
###functions to support linear interpolation###
###############################################

#read preprocessed GPS file without NA values (with whatever function, but you need
#Lat, Long, DateTime for subsequent functions)
#
gps<- readGPSPT2("<whateverfile.csv") 


##updateDiffTime - output is vector of time difference between fixes,
##recalculates elapsed time between points (in seconds)
##based on the first valid fix each day,
##good for finding breaks in reception once NAs are removed
##note: since it splits based on date, it doesn't account for
##overnight activity.
updateDiffTime <- function(gps, timeCol = "DateTime"){
  gps <- split(gps, as.Date(gps[,timeCol]))
  timediff <- lapply(gps,function(gps)  c(0, diff(as.numeric(gps[,timeCol]))))
  timediff <- unlist(timediff)
  return(timediff)
}

gps$newDiffTime <- updateDiffTime(gps, timeCol = "DateTime")
which(gps$newDiffTime>1)
str(gps)

#run this function on gps dataframe with NAs removed and updated diff times
#returns a dataframe of whatever size gaps with preceding
#and relevant rows
#"gapsize" adjusts the time in seconds between fixes that you want to call a gap
#timeCol = recalculated difftime column

findGaps <- function(gps,timeCol = "Timediff", gapSize = 2){
  row.names(gps) <- NULL
  rows1 <- which(gps[,timeCol] >= gapSize) - 1
  rows2 <- which(gps[,timeCol] >= gapSize) 
  point1 <- gps[row.names(gps) %in% rows1,]
  point2 <- gps[row.names(gps)%in% rows2,]
  gapslist <- rbind(point1,point2)
  
  #gapslist <- gapslist[ order(as.numeric(row.names(gapslist))),]
  
  gapslist <- gapslist[ order(as.numeric(gapslist$list)),]
  
  return(list(point2,gapslist))
}

ms <- findGaps(gps, timeCol = "newDiffTime", gapSize = 2) #gaps >5 mins
gapsDF <- as.data.frame(ms[2])  

#generates a list of vectors of time stamps between gaps
#to be used in the linear2 function to generate interpolated fixes
#(future: make time interval customizable with by ?)
getTimeSlice <- function(gapsDF, dateCol = "DateTime"){
  require("dplyr", character.only = T, quietly = TRUE, warn.conflicts = FALSE)
  s <- gapsDF[seq(1,nrow(gapsDF),2),]
  sTime <- s[,dateCol]
  e <- gapsDF[seq(2,nrow(gapsDF),2),]
  eTime <- e[,dateCol]
  
  bTimes <- data.frame(sTime,eTime)
  
  timeSlice <- bTimes %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
       timeList = list(seq.POSIXt(sTime,eTime, by = "sec"))
    )
  
  return(timeSlice[,3])
}

which(bTimes[,2]<bTimes[,1])
t.slice <- getTimeSlice(gapsDF)

#function to interpolate - uses the DF of gaps and time slice list
linear2 <- function(gapsDF, t.slice, dateCol = "DateTime", 
                    lat = 'Latitude',lon = 'Longitude'){
  tzone <- attr(gapsDF[,dateCol],'tzone')
  s <- gapsDF[seq(1,nrow(gapsDF),2),]
  e <- gapsDF[seq(2,nrow(gapsDF),2),]
  x1 <- cbind(s[,lat],s[,lon]) #get first point
  x2 <- cbind(e[,lat],e[,lon]) #get second point
  t1 <- s[,dateCol] #get start time
  t2 <- e[,dateCol] #get end time
  t <- as.numeric(t2 - t1, units = 'secs') #subtract times (seconds)
  ff <- data.frame(t1,t.slice)
  t.s <- mapply('-',ff[,2],ff[,1],SIMPLIFY=FALSE)
  rat <- mapply('/',t.s,t,SIMPLIFY=FALSE)#t.s/t get rate of change
  x <- mapply('*',rat,(x2[,2] - x1[,2]),SIMPLIFY=FALSE) #mutiply rate by difference in longs
  y <- mapply('*',rat,(x2[,1] - x1[,1]),SIMPLIFY=FALSE) #mutiply rate by difference in lats
  bx <- mapply('+',x,x1[,2],SIMPLIFY=FALSE) #add point1 long
  by <- mapply('+',y,x1[,1],SIMPLIFY=FALSE) #add point1 lat
  options(digits = 12)
  bx <- unlist(bx)
  by <- unlist(by)
  zz <- unname(do.call('c', unlist(t.slice,recursive = F)))
  zz <- format(zz,tz = tzone, usetz = T)
  xy <- data.frame(DateTime = zz,Latitude = by,Longitude = bx)
  return(xy)
}

newpoints <- linear2(gapsDF,t.slice)

#### do sometiong here
###combining original and new fixes
library(plyr)
ff <- rbind.fill(gps,newpoints) #merge original gps & new interpolated points
ff <- ff[order(ff$DateTime),] #sort by datetime
ff= ff[!duplicated(ff$DateTime),]
ff$uniqueID <-c(seq.int(nrow(ff))) #make new uniqueID for easy sorting
rownames(ff) <- ff$uniqueID #change row names to match new order
ff<- as.data.frame(ff)

##
##functions to calculate new values between fixes, after merging w/ original GPS
##
#get distance between sucessive observations
getDist <- function(gps,lat = "Latitude",lon = "Longitude"){
  toRadians = pi/180
  latOrig <- gps[1:(nrow(gps) - 1), lat]
  latDest <- gps[2:nrow(gps), lat]
  latDistance = (latOrig - latDest) * toRadians
  lonOrig <- gps[1:(nrow(gps) - 1), lon]
  lonDest <- gps[2:nrow(gps), lon]
  lonDistance = (lonOrig - lonDest) * toRadians
  sindLat <- sin(latDistance/2)
  sindLon <- sin(lonDistance/2)
  a <- sindLat^2 + sindLon^2 * cos(latOrig * toRadians) * cos(latDest * 
                                                                toRadians)
  
  distcalc <- c(0,(6371000 * (2 * atan2(sqrt(a), sqrt(1 - a)))))
  return(distcalc)
}

##updateDist - recalculates distance between points
##based on the first observation each day
##possibly doesn't change anything 
##uses getDist() function, should be computed
##after smoothing method
updateDist <- function(gps, timeCol = "DateTime",
                       lat = "Latitude", lon = "Longitude"){
  x <- split(gps, as.Date(gps[,timeCol]))
  distdiff<- lapply(x,function(x)  getDist(x))
  distdiff <- unlist(distdiff)
  return(distdiff)
}

##calculate updateDiffTime again

ff$newDist<- updateDist(gps=ff, timeCol = "DateTime",
                         lat = "Latitude", lon = "Longitude")

ff$newTimediff <- updateDiffTime(ff)

##update speed
##
##
ff$newSpeedms2 <- ifelse(is.nan(ff$newDist/ff$newTimediff),0,ff$newDist/ff$newTimediff)




