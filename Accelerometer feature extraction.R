#@author: Hoda Allahbakhshi
#@references: Allahbakhshi, H., Conrow, L., Naimi, B., & Weibel, R. (2020). Using accelerometer and GPS data for real-life physical activity type detection. Sensors, 20(3), 588.

#The code extracts multiple features relevant from accelerometer data for classification of PA types. 
# The extracted features are:

# Time domain features: 
# mean, standard deviation and range of three axes and total acceleration,
# correlation among three axes,
# kurtosis, skewness and average absolute difference of three axes,
# number of observations falling within each of 10 bins of the three axes, 
# time interval between local peaks and number of peaks of three axes.

# Frequency domain features using FFT:
# power spectral density, energy of the signal,
# mean of the first three dominant frequencies,
# amplitude of the first three dominant frequencies of three axes and total acceleration

library('seewave')
library('signal')
library('plyr')
library('reshape2')
library('dplyr')
library('stats')
library('e1071')
library("data.table")
library("zoo")
library("pracma")
library("DescTools")


# segmentation

#Method: overlapping fixed size windowing technique
# window size = 2 seconds= 2*50(sampling rate)=100
windowsize=100
# set 50% overalp
ov=0.5
# x=labeled accelerometer data

# Time domain features:

binFreqTable <- function(x,bins,windowsize) {
  freq = hist(x, breaks=bins, include.lowest=TRUE, plot=FALSE)
  ranges = paste(head(freq$breaks,-1), freq$breaks[-1], sep=" - ")
  return(list(range = ranges, frequency = (freq$counts)/windowsize))
}

feature <- function(x,windowsize,ov){
  # start and end time of the segment
  minT<- rollapply((x$DateTime),windowsize,min,by=ov*windowsize)
  maxT<- rollapply((x$DateTime),windowsize,max,by=ov*windowsize)
  # mean of three axes and total acceleration - 4x
  ax <- rollapply((x$Acc_X),windowsize,mean,by=ov*windowsize)
  ay <- rollapply((x$Acc_Y),windowsize,mean,by=ov*windowsize)
  az <- rollapply((x$Acc_Z),windowsize,mean,by=ov*windowsize)
  aSVM<-rollapply((x$SVM),windowsize,mean,by=ov*windowsize)
  # Standard deviation of three axes and total acceleration - 4x
  sx <- rollapply((x$Acc_X),windowsize,sd,by=ov*windowsize)
  sy <- rollapply((x$Acc_Y),windowsize,sd,by=ov*windowsize)
  sz <- rollapply((x$Acc_Z),windowsize,sd,by=ov*windowsize)
  sSVM <- rollapply((x$SVM),windowsize,sd,by=ov*windowsize)
  # Range of three axes and total acceleration, correlation - 4x
  range_x <- rollapply(x$Acc_X,windowsize,range,by=ov*windowsize)
  range_y <- rollapply(x$Acc_Y,windowsize,range,by=ov*windowsize)
  range_z <- rollapply(x$Acc_Z,windowsize,range,by=ov*windowsize)
  range_SVM <- rollapply(x$SVM,windowsize,range,by=ov*windowsize)
  rx <- abs(range_x[,2] - range_x[,1])
  ry <- abs(range_y[,2] - range_y[,1])
  rz <- abs(range_z[,2] - range_z[,1])
  rSVM <- abs(range_SVM[,2] - range_SVM[,1])
  # correlation among three axes- 3x
  r <- cbind(x$Acc_X,x$Acc_Y)
  s <- cbind(x$Acc_X,x$Acc_Z)
  t <- cbind(x$Acc_Z,x$Acc_Y)
  corxy <-  rollapply(r,windowsize,function(r) cor(r[,1],r[,2],use = "everything",method = c("pearson")),by=ov*windowsize,by.column=FALSE)
  corxz <-  rollapply(s,windowsize,function(s) cor(s[,1],s[,2],use = "everything",method = c("pearson")),by=ov*windowsize,by.column=FALSE)
  coryz <-  rollapply(t,windowsize,function(t) cor(t[,1],t[,2],use = "everything",method = c("pearson")),by=ov*windowsize,by.column=FALSE)
  # kurtosis, skewness of three axes- 6x
  skewness_x <- rollapply((x$Acc_X),windowsize,skewness,na.rm = FALSE,by=ov*windowsize)
  kurtosis_x <- rollapply((x$Acc_X),windowsize,kurtosis,na.rm = FALSE,by=ov*windowsize)
  skewness_y <- rollapply((x$Acc_Y),windowsize,skewness,na.rm = FALSE,by=ov*windowsize)
  kurtosis_y <- rollapply((x$Acc_Y),windowsize,kurtosis,na.rm = FALSE,by=ov*windowsize)
  skewness_z <- rollapply((x$Acc_Z),windowsize,skewness,na.rm = FALSE,by=ov*windowsize)
  kurtosis_z <- rollapply((x$Acc_Z),windowsize,kurtosis,na.rm = FALSE,by=ov*windowsize)
  # average absolute difference of three axes and total acceleration - 4x
  aadx  <- rollapply((x$Acc_X),windowsize,MeanAD,by=ov*windowsize)
  aady  <- rollapply((x$Acc_Y),windowsize,MeanAD,by=ov*windowsize)
  aadz <- rollapply((x$Acc_Y),windowsize,MeanAD,by=ov*windowsize)
  aadSVM <- rollapply((x$SVM),windowsize,MeanAD,by=ov*windowsize)
  # number of observations falling within each of 10 bins of the three axes - 30 x
  bins_x=seq(min(x$Acc_X),max(x$Acc_X),by=((max(x$Acc_X) - min(x$Acc_X))/10))
  bins_y=seq(min(x$Acc_Y),max(x$Acc_Y),by=((max(x$Acc_Y) - min(x$Acc_Y))/10))
  bins_z=seq(min(x$Acc_Z),max(x$Acc_Z),by=((max(x$Acc_Z) - min(x$Acc_Z))/10))
  bin_distributionx <- rollapply(x$Acc_X,windowsize,function(x) binFreqTable(x,bins_x,windowsize),by=ov*windowsize,by.column=FALSE)
  bin_distributionx <- data.frame(matrix(unlist(bin_distributionx[,2]), nrow=length(bin_distributionx[,2]), byrow=T))
  bin_distributiony <- rollapply(x$Acc_Y,windowsize,function(x) binFreqTable(x,bins_y,windowsize),by=ov*windowsize,by.column=FALSE)
  bin_distributiony <- data.frame(matrix(unlist(bin_distributiony[,2]), nrow=length(bin_distributiony[,2]), byrow=T))
  bin_distributionz <- rollapply(x$Acc_Z,windowsize,function(x) binFreqTable(x,bins_z,windowsize),by=ov*windowsize,by.column=FALSE)
  bin_distributionz <- data.frame(matrix(unlist(bin_distributionz[,2]), nrow=length(bin_distributionz[,2]), byrow=T))
 # time interval between local peaks and number of peaks of three axes - 6 x
  PeaksX <-  rollapply((x$Acc_X),windowsize,function (x) find_peaks(x),by=ov*windowsize,by.column=FALSE)
  PeaksX<-PeaksX[,1]
  tt=c()
  diffTime=c()
  T_btw_Peaks=c()
  count=c()
  for( i in 1:length(PeaksX)){
    tt<-x$DateTime[PeaksX[[i]]]
    diffTime<- abs(tt - lag(tt, default =tt[1]))
    T_btw_Peaks[[i]]=mean( diffTime)
    count[[i]]<-length(PeaksX[[i]])}
  T_btw_Peaksx=T_btw_Peaks
  PeakNOX=count
  
  PeaksY <-  rollapply((x$Acc_Y),windowsize,function (x) find_peaks(x),by=ov*windowsize,by.column=FALSE)
  PeaksY<-PeaksY[,1]
  tt=c()
  diffTime=c()
  T_btw_Peaks=c()
  count=c()
  for( i in 1:length(PeaksY)){
    tt<- x$DateTime[PeaksY[[i]]]
    diffTime <- abs(tt - lag(tt, default =tt[1]))
    T_btw_Peaks[[i]]=mean( diffTime)
    count[[i]]<-length(PeaksY[[i]])}
  T_btw_Peaksy=T_btw_Peaks
  PeakNOY=count
  
  PeaksZ <-  rollapply((x$Acc_Z),windowsize,function (x) find_peaks2(x),by=ov*windowsize,by.column=FALSE)
  PeaksZ<-PeaksZ[,1]
  tt=c()
  diffTime=c()
  T_btw_Peaks=c()
  count=c()
  for( i in 1:length(PeaksZ)){
    tt<- x$DateTime[PeaksZ[[i]]]
    diffTime <- abs(tt - lag(tt, default =tt[1]))
    T_btw_Peaks[[i]]=mean( diffTime)
    count[[i]]<-length(PeaksZ[[i]])}
  T_btw_Peaksz=T_btw_Peaks
  PeakNOZ=count

  Re <- data.frame(minT,maxT,
                   ax,ay,az,aSVM,
                   sx,sy,sz,sSVM,
                   skewness_x,kurtosis_x,skewness_y,kurtosis_y,skewness_z,kurtosis_z,
                   corxy,corxz,coryz,
                   rx,ry,rz,rSVM,
                   aadx,aady,aadz,aadSVM,
                   bin_distributionx,bin_distributiony,bin_distributionz,
                   T_btw_Peaksx,T_btw_Peaksy,T_btw_Peaksz,
                   PeakNOX,PeakNOY,PeakNOZ)  

  arr <- c('minT','maxT',
           'Avgaccx','Avgaccy','Avgaccz','AvgSVM',
           'Stdx','Stdy','Stdz','StdSVM','skewnessx','kurtosisx','skewnessy','kurtosisy','skewnessz','kurtosisz',
           'corxy','corxz','coryz',
           'Rangex','Rangey','Rangez','RangeSVM',
           'Avgabsdiffx','Avgabsdiffy','Avgabsdiffz','AvgabsdiffSVM',
           "bin1x","bin2x","bin3x","bin4x","bin5x","bin6x","bin7x","bin8x","bin9x","bin10x",
           "bin1y","bin2y","bin3y","bin4y","bin5y","bin6y","bin7y","bin8y","bin9y","bin10y",
           "bin1z","bin2z","bin3z","bin4z","bin5z","bin6z","bin7z","bin8z","bin9z","bin10z",
           "T_btw_Peaksx","T_btw_Peaksy","T_btw_Peaksz",
           "PeakNOX","PeakNOY","PeakNOZ")
  colnames(Re) <- arr 
  return(Re)}
#Frequency domain features
windowsize=100
ov=0.5
# x=labeled accelerometer data
dfeatures <- function (x,windowsize,f) { 
  #energy of the signal
  E_f <- sum(abs(x)*abs(x))/windowsize 
  X <- x[((windowsize/2) + 1):windowsize]
  X <- data.frame(X)
  colnames(X) <- "Index"
  indexes <- sort(abs(X$Index), decreasing = TRUE)
  #  amplitude of the first three dominant frequencies of signal
  Adf1=10*log10(indexes[1] ) 
  Adf2=10*log10(indexes[2])
  Adf3=10*log10(indexes[3])
  f <- data.frame(f)
  colnames(f) <- "frequency"
  f1 <- f[((windowsize/2)+1):nrow(f),]
  df1 <- f1[which(abs(X) == indexes[1])] 
  df2 <- f1[which(abs(X) == indexes[2])]
  df3 <- f1[which(abs(X) == indexes[3])] 
  # power spectral density
  PSD <- X
  PSD <- abs(PSD)*abs(PSD)
  PSD <- PSD/windowsize
  PSD <- 2*PSD 
  PSD <- data.frame(PSD)
  colnames(PSD ) <- "PSD "
  APSD <- 10*log10(PSD$PSD[1])
  #mean of the first three dominant frequencies
  meandf=mean(c(df1,df2,df3))
  a=cbind(E_f,Adf1,Adf2,Adf3,meandf,APSD)
  arr <- c('ef ','Adf1','Adf2','Adf3','meandf','APSD')
  colnames(a)<- arr 
  return(a)
} 

Fourier1 <- function(x) {  
  Fourie <- fftshift((fft(x))) 
  return( Fourie)}

extractfeatures <- function(x,windowsize,ov) { 
  final <- as.list(seq(floor(nrow(x)/windowsize)))
  finalf <- as.list(seq(floor(nrow(x)/windowsize)))
  xf <- as.list(seq(floor(nrow(x)/windowsize)))
  finalfx <- as.list(seq(floor(nrow(x)/windowsize)))
  xfx <- as.list(seq(floor(nrow(x)/windowsize)))
  finalfy <- as.list(seq(floor(nrow(x)/windowsize)))
  xfy <- as.list(seq(floor(nrow(x)/windowsize)))
  finalfz <- as.list(seq(floor(nrow(x)/windowsize)))
  xfz <- as.list(seq(floor(nrow(x)/windowsize)))
  Ts <- as.numeric(x$DateTime[2]-x$DateTime[1])
  ##AXIS GENERATION.
  t  = c((0:(windowsize - 1))*Ts)  ##Time axis, [s].
  Fs  = 1/Ts 
  f  <- seq((-0.5)*Fs , (0.5 - 1/windowsize )*Fs  ,by= (1/windowsize )*Fs )##Frequency axis, [Hz].
  # Removing the DC component and considering the DC as the mean
  ax <- rollapply((x$Acc_X),windowsize, function(x){abs(mean(x) - x)},by=ov*windowsize)
  ay <- rollapply((x$Acc_Y),windowsize, function(x){abs(mean(x) - x)},by=ov*windowsize)
  az <- rollapply((x$Acc_Z),windowsize, function(x){abs(mean(x) - x)},by=ov*windowsize)
  aSVM<-rollapply((x$SVM),windowsize, function(x){abs(mean(x) - x)},by=ov*windowsize)
  xf<-  apply(aSVM,1,Fourier1)
  xf<-t(xf)
  xfx <-apply(ax,1,Fourier1)
  xfx<-t(xfx)
  xfy <- apply(ay,1,Fourier1)
  xfy<-t(xfy)
  xfz <- apply(az,1,Fourier1)
  xfz<-t(xfz)
  
  for (i in 1:nrow(xf)){
    finalf[[i]]  <- dfeatures(xf[i,],windowsize,f=f)}
  for (i in 1:nrow(xfx)){
    finalfx[[i]]  <- dfeatures(xfx[i,],windowsize,f=f)}
  for (i in 1:nrow(xfy)){
    finalfy[[i]]  <- dfeatures(xfy[i,],windowsize,f=f)}
  for (i in 1:nrow(xfz)){
    finalfz[[i]]  <- dfeatures(xfz[i,],windowsize,f=f)}
  for (i in 1:nrow(xf)){
    final[[i]]=cbind(finalf[[i]],finalfx[[i]],finalfy[[i]],finalfz[[i]])
  }
  arr <- c('ef ','Adf1','Adf2','Adf3','meandf','APSD',
           'efx ','Adf1x','Adf2x','Adf3x','meandfx','APSDx',
           'efy','Adf1y','Adf2y','Adf3y','meandfy','APSDy',
           'efz','Adf1z','Adf2z','Adf3z','meandfz','APSDz')
  for (i in 1:length(final)){
    colnames( final[[i]])<- arr 
  }
  return(final)
}

#######################################################
# Calculate all time and frequency domain features from accelerometer data
Window=100
Get_Final_features<- function(x)
  
{Labelled_data <- as.list(seq_along(x))
names(Labelled_data ) <- devices
Time_features <- as.list(seq_along(x))
names(Time_features ) <- devices
frequency_feature<- as.list(seq_along(x))
names(frequency_feature ) <- devices
freqfeature <- as.list(seq_along(x))
names(freqfeature ) <- devices
Feature_label <- as.list(seq_along(x))
names( Feature_label ) <- devices
Final_features <- as.list(seq_along(x))
names(Final_features ) <- devices

for (i in seq_along(x))
  
{Labelled_data[[i]]<-x[[i]]
Labelled_data[[i]]$Segment_id <- (seq(nrow(x[[i]])) - 1) %/% Window
Time_features[[i]] <- feature(x[[i]],Window,ov)

Feature_label[[i]]=as.factor(x[[i]]$activity)
Feature_label[[i]]=rollapply(Feature_label[[i]],Window,table,by=ov*Window)
Feature_label[[i]]=as.data.frame(Feature_label[[i]])
Feature_label[[i]]=colnames(Feature_label[[i]])[(apply(Feature_label[[i]],MARGIN=1,which.max))]
Time_features[[i]]$activity<- Feature_label[[i]]
freqfeature[[i]] <- extractfeatures(x[[i]],Window,ov)
frequency_feature[[i]]  = do.call("rbind",freqfeature[[i]])
Final_features[[i]]<- cbind(Time_features[[i]],frequency_feature[[i]] )}
return(Final_features) 

}



