
# Extract accelerometer features

Final_features <- as.list(seq_along(x))

names(Final_features ) <- seq_along(x)

for (i in seq_along(x)){ 
  tic()
  Final_features[[i]]<- Get_Final_features(x[[i]][[3]])
  toc()
  i
}
gc()



# To combine the GPS and accelerometer data, they had to be time-matched together. 
#I matched the GPS timestamps with the first and last timestamps of accelerometer segments 
#rounded down to the nearest second to combine the GPS data with the accelerometer data



for (i in 1:Nuser){ 
  
  for (j in seq_along(devices)){
    op <- options(digits.secs = 3)
    Final_features[[i]][[j]]$minT<- as.POSIXct(Final_features[[i]][[j]]$minT, origin="1970-01-01",tz='GMT')
    Final_features[[i]][[j]]$maxT<- as.POSIXct(Final_features[[i]][[j]]$maxT,origin="1970-01-01",tz='GMT')
    
    matchtest=as.data.frame (match_inds(GPS_Data_with_elv[[i]][[4]][[j]]$DateTime,Final_features[[i]][[j]]$minT))
    matchtest2=as.data.frame (match_inds(GPS_Data_with_elv[[i]][[4]][[j]]$DateTime,Final_features[[i]][[j]]$maxT))
    mtch=cbind(matchtest,matchtest2)
    
    matchtestraw=as.data.frame (match_inds( GPS_Data_with_elv[[i]][[4]][[j]]$DateTime,Final_features[[i]][[j]]$minT))
    matchtestraw2=as.data.frame (match_inds( GPS_Data_with_elv[[i]][[4]][[j]]$DateTime,Final_features[[i]][[j]]$maxT))
    mtch2=cbind(matchtestraw,matchtestraw2)
    
    a= GPS_Data_with_elv[[i]][[4]][[j]]$newSpeedms2[mtch2[,1]]
    b= GPS_Data_with_elv[[i]][[4]][[j]]$newSpeedms2[mtch2[,2]]
    c=(a+b)/2
    
    e= GPS_Data_with_elv[[i]][[4]][[j]]$elevation[mtch[,2]]-All_Senior_Users_Final_Semi_Data_with_elv_map_matched[[i]][[4]][[j]]$elevation[mtch[,1]]
    
    Final_features[[i]][[j]]$Speed <- c
    Final_features[[i]][[j]]$Elevation_diff <- e
  }
}




