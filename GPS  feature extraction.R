#@author: Hoda Allahbakhshi
#@references: Allahbakhshi, H., Conrow, L., Naimi, B., & Weibel, R. (2020). Using accelerometer and GPS data for real-life physical activity type detection. Sensors, 20(3), 588.


#The function extracts multiple features relevant from GPS data for classification of PA types. 
# The extracted features are:
#Average speed and elevation difference 

library('dplyr')
library("rgdal")



#Data: load labeled GPS data
devices<- c("chest","Knee","Lhip-front","Lpocket","Rhip-side")

# Interpolate GPS missing values
# first run the relevant function in interpolation_Linear.R

for (i in seq_along(Data)){ 
  
  ff <- as.list(seq_along(devices))
  names(ff ) <- devices
  
  stop<- as.POSIXct(Data[[i]][[1]][1,7])
  attr(stop, "tzone") <- "GMT"
  start<-as.POSIXct(Data[[i]][[1]][nrow(Data[[i]][[1]]),6])
  attr(start, "tzone") <- "GMT"
  gps_nrow<- as.numeric(abs(difftime(start,stop,units = "secs"))*1+1)
  
  SE_GPS<-Start_End(x=Data[[i]][[2]])
  

  for (j in seq_along(devices)) {
    
    gps <- SE_GPS[[j]]
    #segment's start time
    SgpsT<-Data[[i]][[2]][[j]]$DateTime[1]
    #segment's end time
    EgpsT<-Data[[i]][[2]][[j]]$DateTime[nrow(gps)]
    
    gpsTimeLength<-as.numeric(abs(difftime(SgpsT,EgpsT,units = "secs"))*1 +1)
    
    
    if (nrow(gps)==as.numeric(gps_nrow))
    {
      ff[[j]] <- gps
      ff[[j]]$newDiffTime <- updateDiffTime(ff[[j]], timeCol = "DateTime")
      library(plyr)
      ff[[j]] <- ff[[j]][order(ff[[j]]$DateTime),] #sort by datetime
      ff[[j]]= ff[[j]][!duplicated(ff[[j]]$DateTime),]
      ff[[j]]$uniqueID <-c(seq.int(nrow(ff[[j]]))) #make new uniqueID for easy sorting
      rownames(ff[[j]]) <- ff[[j]]$uniqueID #change row names to match new order
      ff[[j]]<- as.data.frame(ff[[j]])
      ##calculate updateDiffTime again
      ff[[j]]$newDist<- updateDist(ff[[j]], timeCol = "DateTime",lat = "Latitude", lon = "Longitude")
      ff[[j]]$newTimediff <- updateDiffTime(ff[[j]])
      ff[[j]]$newSpeedms2 <- ifelse(is.nan(ff[[j]]$newDist/ff[[j]]$newTimediff),0,ff[[j]]$newDist/ff[[j]]$newTimediff)
      
    } else if (nrow(gps)!=as.numeric(gps_nrow)){
      gps <- rmNAGPS(gps)
      gps$newDiffTime <- updateDiffTime(gps, timeCol = "DateTime")
      ms <- findGaps(gps, timeCol = "newDiffTime", gapSize = 2) #gaps >5 mins
      gapsDF <- as.data.frame(ms[2])  
      
      if (nrow(gapsDF)==0){
        ff[[j]] <- gps
        ff[[j]]$newDiffTime <- updateDiffTime(ff[[j]], timeCol = "DateTime")
        library(plyr)
        ff[[j]] <- ff[[j]][order(ff[[j]]$DateTime),] #sort by datetime
        ff[[j]]= ff[[j]][!duplicated(ff[[j]]$DateTime),]
        ff[[j]]$uniqueID <-c(seq.int(nrow(ff[[j]]))) #make new uniqueID for easy sorting
        rownames(ff[[j]]) <- ff[[j]]$uniqueID #change row names to match new order
        ff[[j]]<- as.data.frame(ff[[j]])
        ##calculate updateDiffTime again
        ff[[j]]$newDist<- updateDist(ff[[j]], timeCol = "DateTime",lat = "Latitude", lon = "Longitude")
        ff[[j]]$newTimediff <- updateDiffTime(ff[[j]])
        ff[[j]]$newSpeedms2 <- ifelse(is.nan(ff[[j]]$newDist/ff[[j]]$newTimediff),0,ff[[j]]$newDist/ff[[j]]$newTimediff)
        
      }
      else if (nrow(gapsDF)> 0){
        
        t.slice <- getTimeSlice(gapsDF)
        newpoints <- linear2(gapsDF,t.slice)
        
        library(plyr)
        ff[[j]] <- rbind.fill(gps,newpoints) #merge original gps & new interpolated points
        ff[[j]] <- ff[[j]][order(ff[[j]]$DateTime),] #sort by datetime
        ff[[j]]= ff[[j]][!duplicated(ff[[j]]$DateTime),]
        ff[[j]]$uniqueID <-c(seq.int(nrow(ff[[j]]))) #make new uniqueID for easy sorting
        rownames(ff[[j]]) <- ff[[j]]$uniqueID #change row names to match new order
        ff[[j]]<- as.data.frame(ff[[j]])
        ##calculate updateDiffTime again
        ff[[j]]$newDist<- updateDist(gps=ff[[j]], timeCol = "DateTime",lat = "Latitude", lon = "Longitude")
        ff[[j]]$newTimediff <- updateDiffTime(ff[[j]])
        ##update speed
        
        ff[[j]]$newSpeedms2 <- ifelse(is.nan(ff[[j]]$newDist/ff[[j]]$newTimediff),0,ff[[j]]$newDist/ff[[j]]$newTimediff)
        
      }
    }
    
  }
  
  Data[[i]]<- rlist::list.append(Data[[i]],ff)
  
  gps<- NULL
  ms <- NULL
  gapsDF <-NULL
  t.slice <- NULL
  newpoints <- NULL
  
}


# To derive elevation difference feature, first map-match the interpolated GPS data
#Afterward, use the map-matched GPS coordinates to derive an elevation value from DEM in ArcGIS software v.10.6.1
# using "Extract value to points" tool in order to assign an elevation value to each GPS point

#  read elevation data
#Nuser:number of participants

elev_data=as.list(c(1:Nuser))
for (i in 1:Nuser)
{
  elev_data[[i]] <- append(elev_data[[i]],as.list(c(1:4)))}


for (i in 1:Nuser)
{
  ind=seq(1,150,by=5)
  
  if (i!=Nuser)
    
    k=(ind[i]):(ind[i+1]-1)
  
  if (i==Nuser)
    k=(ind[i]):(length(files))
  
  for(j in seq_along(devices))
  
    elev_data[[i]][[j]]=readOGR(dsn=files[[k[j]]],layer=gsub('.{4}$', '',sub(".*/", "", files[k[j]])), stringsAsFactors=default.stringsAsFactors())
}


elev_data2=as.list(c(1:Nuser))
for(i in 1:Nuser)
  elev_data2[[i]]=append(elev_data_Semi2[[i]],as.list(c(1:4))) 

for (i in 1:Nuser) {
  for (j in seq_along(devices))
    elev_data[[i]][[j]]<-elev_data[[i]][[j]]@data
}


#Data: load labeled GPS data
#GPS_Data_with_elv <-Data

for (i in 1:Nuser) {
  
  for (j in seq_along(devices))
  {
    
    GPS_Data_with_elv[[i]][[4]][[j]]$elevation<-elev_data2[[i]][[j]]$RASTERVALU
    
    
  }
}

####################################################################################
