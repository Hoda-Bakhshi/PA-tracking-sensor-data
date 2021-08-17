
#@Author: Germán Carrillo
# https://geotux.tuxfamily.org/en/2012/02/20/snapping-points-to-lines-in-r-en/
# The functions snapPointsToLines, nearestPointOnSegment and nearestPointOnLine
# have been included into the R package called maptools (v. 0.8-16).
######################################################
#@Adapted by: Hoda Allahbakhshi


# Data= labeled preprocessed GPS data
devices<- c("chest","Knee","Lhip-front","Lpocket","Rhip-side")

#######################################################
p=Data

require("raster")
require("sp")
require("rgdal")

for (i in 1:length(Data))
{
  
  for (j in seq_along(devices))
  {
    test=Data[[i]][[4]][[j]]
    rr<- SpatialPointsDataFrame(test[,c(3,4)],test)
    crs<- "+init=epsg:2056"
    proj4string(rr) = crs
    
    test2=as.data.frame(c(1:nrow(test)))
    test2$Lat=test$Latitude
    test2$Long=test$Longitude
    test2$Time=test$Time
    test2$Date=test$Date
    test2$pm10=NULL
    test2$TimeTotal=test$list
    test2$t=test$DateTime
    test3=test2[,2:7]
    x=test3
    m = matrix(c(x[,2],x[,1]),ncol=2)
    df = data.frame(x$t, x[,5],x[,6])
    
    p[[i]][[j]] = SpatialPointsDataFrame(coords=m,data=df,proj4string=CRS("+init=epsg:2056"))
    p[[i]][[j]] = spTransform( p[[i]][[j]],crs)
    
  }
  
}

# read OSM road data 
l = readOGR(dsn="osm data", layer="osmroad")

proj4string(l) = crs
l = spTransform(l,crs)


nearestPointOnSegment = function(s, p){
  ap = c(p[1] - s[1,1], p[2] - s[1,2])
  ab = c(s[2,1] - s[1,1], s[2,2] - s[1,2])
  t = sum(ap*ab) / sum(ab*ab)
  t = ifelse(t<0,0,ifelse(t>1,1,t))
  x = s[1,1] + ab[1] * t 
  y = s[1,2] + ab[2] * t
  c(x, y, (x-p[1])^2 + (y-p[2])^2)  # Return nearest point and distance
}
nearestPointOnLine = function(coordsLine, coordsPoints){
  nearest_points = vapply(2:nrow(coordsLine), 
                          function(x) 
                            nearestPointOnSegment(coordsLine[(x-1):x,], coordsPoints),
                          FUN.VALUE=c(0,0,0))
  
  # Return coordinates of the nearest point in this line  
  nearest_points[1:2, which.min(nearest_points[3,])]  
}

snapPointsToLines <- function(points, lines, maxDist=NA, withAttrs=TRUE) {
  
  require("rgeos")
  
  if (!is.na(maxDist)) {
    w = gWithinDistance(points, lines, dist=maxDist, byid=TRUE)
    validPoints = apply(w,2,any)
    validLines = apply(w,1,any)
    points = points[validPoints,]
    lines =  lines[validLines,]
  }

  d = gDistance(points, lines, byid=TRUE) 
  nearest_line_index = apply(d, 2, which.min) # Position of each nearest line in lines object 
  
  coordsLines = coordinates(lines)  
  coordsPoints = coordinates(points)  
  
  # Get coordinates of nearest points lying on nearest lines
  mNewCoords = vapply(1:length(points), 
                      function(x) 
                        nearestPointOnLine(coordsLines[[nearest_line_index[x]]][[1]], 
                                           coordsPoints[x,]), FUN.VALUE=c(0,0))
  
  # Recover lines' Ids (Ids and index differ if maxDist is given)
  if (!is.na(maxDist)) nearest_line_id = as.numeric(rownames(d)[nearest_line_index])+1 
  else nearest_line_id = nearest_line_index 
  
  # Create data frame and sp points
  if (withAttrs) df = cbind(points@data, nearest_line_id) 
  else df = data.frame(nearest_line_id, row.names=names(nearest_line_index))
  
  SpatialPointsDataFrame(coords=t(mNewCoords), data=df, 
                         proj4string=CRS(proj4string(points)))
}

##########################################

library(tictoc)
tic()
  for (i in 1:length(Data))

  {
    tic()
    for (j in seq_along(devices))
   { snappedPoints[[i]][[j]] = snapPointsToLines(p[[i]][[j]], s=l, 5, FALSE) 
   snappedPointsFinal<- snappedPoints
    }
    message("This is step:", i)
    toc()
  }
toc()

Map_matched_GPS_Points<-snappedPointsFinal
