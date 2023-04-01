#setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/Sikorsky")
#setwd("/home/nvp19001/Final/Sikorsky")
setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")

load("Sikorsky-Bridgeport_merged.Rda")
#can call the data as Data
#head(Data)
#Wave Data from Buoy
load("44039_reg_merged.Rda")
#can call the data as merged
#head(merged)
library(mgcv)


#whichYear:2004-2013
#indSeason:Annual-Annual Data,Windy-Windy Season Data,Calm-Calm Season Data
data.in_CS<-function(whichYear,indSeason)
{
  library(lubridate)
  library(dplyr)
  #For Wave-height
  filled<-merged
  filled$Timestamp<-ymd_hms(filled$Timestamp)
  filled$Timestamp1<-floor_date(filled$Timestamp, "hour")
  Data1<-data.frame(seq(ymd_hms("2004-1-1 00:00:00"),ymd_hms(" 2013-12-31 23:00:00"),by="hours"))
  colnames(Data1)<-"TimeStamp_1"
  Data_New<-left_join(Data1,filled,by=c("TimeStamp_1"="Timestamp1"),multiple = "all")
  d1<-duplicated(Data_New$TimeStamp_1)
  Data_New_1<-Data_New[-which(d1==TRUE),]
  filled<-Data_New_1
  
  #For Wind Data
  Data1<-data.frame(seq(ymd_hms("2004-1-1 00:00:00"),ymd_hms(" 2013-12-31 23:00:00"),by="hours"))
  
  colnames(Data1)<-"TimeStamp_1"
  Data$Timestamp<-ymd_hms(paste(Data$year, Data$month, Data$day, Data$hour, Data$minute, Data$second, sep="-"))
  Data1<-data.frame(seq(ymd_hms("2004-1-1 00:00:00"),ymd_hms(" 2013-12-31 23:00:00"),by="hours"))
  colnames(Data1)<-"TimeStamp_1"
  Data1$event_ind<-rep(0,nrow(Data1))
  #Data1
  EventData <-read.csv("~/Dropbox/Namitha Pais/Marine Sciences/Storm_Data_Filtered_2004_2013.csv")
  EventData$Timestamp <- as.POSIXct(EventData$Timestamp , format = "%m/%d/%y %H:%M")
  EventData$Timestamp <-round(EventData$Timestamp , "hours")
  EventData$Timestamp  <- ymd_hms(format(EventData$Timestamp , "%Y-%m-%d %H:%M:%S"))
  EventData_sorted <- EventData[order(EventData$Timestamp), ]
  indices<-which(Data1$TimeStamp_1 %in% EventData_sorted$Timestamp)
  for(i in 1:length(indices))
  {
    Data1$event_ind[ (indices[i]-12):(indices[i]+12)]<- 1
  }
  Data1$event_ind<-as.factor(Data1$event_ind)
  Data$Timestamp1<-floor_date(Data$Timestamp, "hour")
  Data_New<-left_join(Data1,Data,by=c("TimeStamp_1"="Timestamp1"),multiple = "all")
  d1<-duplicated(Data_New$TimeStamp_1)
  Data_New_1<-Data_New[-which(d1==TRUE),]
  Data<-Data_New_1
  
  Data$YYYY<-year(Data$TimeStamp_1)
  Data$MM<-month(Data$TimeStamp_1)
  Data$DD<-day(Data$TimeStamp_1)
  Data$hh<-hour(Data$TimeStamp_1)
  Data$mm<-minute(Data$TimeStamp_1)
  Data$WVHT<-filled$WVHT
  Data$WDIR_B<-filled$WDIR
  Data$WSPD_B<-filled$WSPD
  
  #cor(Data$WSPD,Data$WSPD_B,use="pairwise.complete.obs")
  #cor(Data$WDIR,Data$WDIR_B,use="pairwise.complete.obs",method="kendall")
  
  #to be used on post-interpolation files
  #take 1 year
  library(imputeTS)
  data.in_CY = Data[Data$YYYY==whichYear,]
  data.in_CY_M<-data.in_CY[(data.in_CY$MM == 1 | data.in_CY$MM == 2 |data.in_CY$MM == 3|
                              data.in_CY$MM == 4|data.in_CY$MM == 5|data.in_CY$MM == 6|data.in_CY$MM == 7|
                              data.in_CY$MM == 8 |data.in_CY$MM == 9|data.in_CY$MM == 10),]
  data.in_PY=Data[Data$YYYY==(whichYear-1),]
  data.in_PY_M<-data.in_PY[(data.in_PY$MM == 11 | data.in_PY$MM == 12),]
  data.in<-rbind(data.in_PY_M,data.in_CY_M)
  data.in$WSPD[which(is.na(data.in$WSPD)==FALSE & is.na(data.in$WDIR)==TRUE)]<-NaN
  data.in$WDIR[which(is.na(data.in$WSPD)==TRUE & is.na(data.in$WDIR)==FALSE)]<-NaN
  data.in$mmReg = "00"
  data.in$regTimeString = paste(data.in$YYYY,"-",data.in$MM,"-",data.in$DD," ",data.in$hh,":",data.in$mmReg,sep="")
  data.in$regTimestamp = strptime(data.in$regTimeString, "%Y-%m-%d %H:%M", tz="GMT")
  summary(data.in$regTimestamp)
  #add Imperial feet, threshold
  
  data.in$SWHm = data.in$WVHT
  data.in$SWHft1 = 3.28084 * data.in$SWHm
  data.in$colSWHft<-ifelse(is.na(data.in$SWHft1)==TRUE,"yellow2","black")
  data.in$indImp<-is.na(data.in$SWHft1)
  #data.in$SWHft<-na_interpolation(data.in$SWHft1)
  CalmInd<-(data.in$MM == 4|data.in$MM == 5|data.in$MM == 6|data.in$MM == 7|
              data.in$MM == 8 |data.in$MM == 9|data.in$MM == 10)
  WindyInd<-(data.in$MM == 11 | data.in$MM == 12 |data.in$MM == 1 | data.in$MM == 2 |data.in$MM == 3)
  #indSWHft<-which(colnames(data.in)=="SWHft1")
  DataWindy<-data.in[WindyInd,]
  DataCalm<-data.in[CalmInd,]

  
  #Imputing (WSPD,WDIR) pairwise.
  set.seed(3)
  mWindy=rle(is.na(DataWindy$WSPD))
  mWindy1<-rep(ifelse(mWindy$values,mWindy$lengths,0),times = mWindy$lengths)
  ZWindy<-rep(0,nrow(DataWindy))
  ZWindy[which(mWindy1<=5  & mWindy1>0)]<-rnorm(length(which(mWindy1<=5 & mWindy1>0)),0,1)
  sdWSPDWindy<-sd(na.omit(DataWindy$WSPD))
  sdWDIRWindy<-sd(na.omit(DataWindy$WDIR))
  
  mCalm=rle(is.na(DataCalm$WSPD))
  mCalm1<-rep(ifelse(mCalm$values,mCalm$lengths,0),times = mCalm$lengths)
  ZCalm<-rep(0,nrow(DataCalm))
  ZCalm[which(mCalm1<=5  & mCalm1>0)]<-rnorm(length(which(mCalm1<=5 & mCalm1>0)),0,1)
  sdWSPDCalm<-sd(na.omit(DataCalm$WSPD))
  sdWDIRCalm<-sd(na.omit(DataCalm$WDIR))
  
  #Case1
  #if the NA length<=5 replace by Lin Int + c*sd*N(0,1)
  c<-1
  #Windy
  DataWindy$WSPD1<-DataWindy$WSPD
  DataWindy$WSPD<-na_interpolation(DataWindy$WSPD1,maxgap=5)+(c*sdWSPDWindy*ZWindy)
  DataWindy$WDIR1<-DataWindy$WDIR
  DataWindy$WDIR<-na_interpolation(DataWindy$WDIR1,maxgap=5)
  #Calm
  DataCalm$WSPD1<-DataCalm$WSPD
  DataCalm$WSPD<-na_interpolation(DataCalm$WSPD1,maxgap=5)+(c*sdWSPDCalm*ZCalm)
  DataCalm$WDIR1<-DataCalm$WDIR
  DataCalm$WDIR<-na_interpolation(DataCalm$WDIR1,maxgap=5)
  
  #Case2
  #if the NA length>5 replace by rs from Known pairs of(WSPD,WDIR)
  #Windy
  indW<-which(mWindy1>5)
  indW_0<-which(mWindy1==0)
  set.seed(5)
  sampleWindy<-sample(which(mWindy1==0),size=length(indW),replace=TRUE)
  DataWindy$WSPD[indW]<-DataWindy$WSPD[sampleWindy]
  DataWindy$WDIR[indW]<-DataWindy$WDIR[sampleWindy]
  DataWindy$WSPD<-ifelse(DataWindy$WSPD<=0.25,0.25,DataWindy$WSPD)
  DataWindy$WDIR<-ifelse(DataWindy$WDIR<=0,0,DataWindy$WDIR)
  #Calm
  indC<-which(mCalm1>5)
  indC_0<-which(mCalm1==0)
  set.seed(5)
  sampleCalm<-sample(which(mCalm1==0),size=length(indC),replace=TRUE)
  DataCalm$WSPD[indC]<-DataCalm$WSPD[sampleCalm]
  DataCalm$WDIR[indC]<-DataCalm$WDIR[sampleCalm]
  DataCalm$WSPD<-ifelse(DataCalm$WSPD<=0.25,0.25,DataCalm$WSPD)
  DataCalm$WDIR<-ifelse(DataCalm$WDIR<=0,0,DataCalm$WDIR)
  
  threshold=0
  #data.in$Flag = as.integer(data.in$SWHft > threshold)
  
  library(plyr)
  data.in$WSPD<-c(DataWindy$WSPD,DataCalm$WSPD)
  data.in$WDIR<-c(DataWindy$WDIR,DataCalm$WDIR)
  data.in$WDIR<-round_any(data.in$WDIR,10)
  
  
  #orthogonal decomposition of waves
  NoSoComponents = sin(data.in$WDIR/360*2*pi)
  EaWeComponents = cos(data.in$WDIR/360*2*pi)
  NoSoBaseline = mean(NoSoComponents)
  EaWeComponents = mean(EaWeComponents)
  TangentRatio = NoSoBaseline / EaWeComponents
  MeanDir = (atan(TangentRatio) %% (2*pi)) / (2*pi) * 360
  OrthoDir = (MeanDir+90) %% 360
  
  #http://www.movable-type.co.uk/scripts/latlong.html
  #40°53'22"N 73°45'51"W
  #41°15'25"N 71°49'52"W
  #Distance:	167.1 km (to 4 SF*)
  #Initial bearing:	075° 12′ 44″
  #Final bearing:	076° 28′ 57″
  #Midpoint:	41° 04′ 38″ N, 072° 48′ 01″ W
  
  #Using 76° bearing or 14° north of due east (90°)
  #Wind is measured in FROM direction...
  #From Southwest
  AlongBasisAngle = (76 + 180) %% 360
  #From Northwest
  AcrossBasisAngle = (AlongBasisAngle + 90) %% 360
  
  #compass degress to geometric degress
  # G = [(360-(90-(360-(C %% 360)) ) %% 360]
  data.in$geometricWDIR = (360-(90-(360-(data.in$WDIR %% 360)))) %% 360
  
  #adjust for LI Sound Orientation
  reordientedGeometric = (data.in$geometricWDIR - 14) %% 360
  acrossGeometric = sin(reordientedGeometric/360*2*pi)*data.in$WSPD
  alongGeometric = cos(reordientedGeometric/360*2*pi)*data.in$WSPD
  
  #the -sin() is to do with the unorthodox basis
  data.in$Minor = -(sin(reordientedGeometric/360*2*pi)*data.in$WSPD)
  data.in$Major = cos(reordientedGeometric/360*2*pi)*data.in$WSPD
  
  
  
  #Compass degrees and Fetch distances for Western Bouy & Central Bouy
  #Degrees at true not magnetic
  #distances are miles (eyeballed from FAA sectional with protractor)
  #WLIS 44040 --  40.956 N 73.580 W (40°57'21" N 73°34'48" W)
  #CLIS 44029 --  41.138 N 72.655 W (41°8'15" N 72°39'18" W)
  #
  # >>>  Note that the buoys seem to report both 0 and 360
  fetchVecMiles <- c(
    0 , 3 , 9 ,
    10 , 5 , 10 ,
    20 , 5 , 10 ,
    30 , 10 , 10 ,
    40 , 12 , 10 ,
    50 , 17 , 15 ,
    60 , 40 , 16 ,
    70 , 75 , 40 ,
    80 , 61 , 55 ,
    90 , 22 , 15 ,
    100 , 17 , 14 ,
    110 , 11 , 11 ,
    120 , 4 , 11 ,
    130 , 4 , 11 ,
    140 , 3 , 10 ,
    150 , 3 , 10 ,
    160 , 3 , 11 ,
    170 , 3 , 11 ,
    180 , 3 , 11 ,
    190 , 3 , 11 ,
    200 , 3 , 12 ,
    210 , 4 , 14 ,
    220 , 5 , 15 ,
    230 , 10 , 19 ,
    240 , 14 , 25 ,
    250 , 11 , 40 ,
    260 , 8 , 53 ,
    270 , 7 , 31 ,
    280 , 6 , 23 ,
    290 , 5 , 19 ,
    300 , 5 , 17 ,
    310 , 5 , 14 ,
    320 , 5 , 11 ,
    330 , 5 , 10 ,
    340 , 5 , 8 ,
    350 , 4 , 9 ,
    360 , 3 , 9
  )
  fetchTabMiles <- as.data.frame(matrix(data = fetchVecMiles, nrow = 37, byrow = TRUE))
  names(fetchTabMiles) <- c("Deg","wlisMiles","clisMiles")
  whichBuoy = "central"
  #need to assign "wlisMiles" or "clisMiles"
  #
  # also doing miles-to-meters conversions 1mi = 1609.34m
  if (whichBuoy == "central"){
    data.in$Fetch <- 1609.34*merge(data.in, fetchTabMiles, by.x = "WDIR", by.y = "Deg", all.x = TRUE)$clisMiles
  }
  if (whichBuoy == "western"){
    data.in$Fetch <- 1609.34*merge(data.in, fetchTabMiles, by.x = "WDIR", by.y = "Deg", all.x = TRUE)$wlisMiles
  }
  
  ##################
  
  #seasonal dummies
  
  data.in$JFM = as.integer(
    (data.in$MM == 1 | data.in$MM == 2 | data.in$MM == 3 ) )
  data.in$JAS = as.integer(
    (data.in$MM == 7 | data.in$MM == 8 | data.in$MM == 9 ) )
  data.in$JJA = as.integer(
    (data.in$MM == 6 | data.in$MM == 7 | data.in$MM == 8 ) )
  
  
  
  ###
  #add in AR(r) on response
  r = 1  #hours
  findLagLocation = function(testTime){
    lagTime = testTime-3600*r
    lagLoc = which(data.in$regTimestamp == lagTime)
    if (length(lagLoc)<=0) {lagLoc = NA}
    return(lagLoc)
  }
  
  #is the lag we want actually available?
  data.in$lagLoc = rep(NA,dim(data.in)[1])
  for (i in 1:length(data.in$lagLoc)){
    data.in$lagLoc[i] = findLagLocation(data.in$regTimestamp[i])
  }
  
  #data.in$rFlag = data.in$Flag[data.in$lagLoc]
  #1
  data.in$rWSPD = data.in$WSPD[data.in$lagLoc]
  data.in$rMinor = (data.in$Minor[data.in$lagLoc])
  data.in$rMajor = (data.in$Major[data.in$lagLoc])
  data.in$rFetch = data.in$Fetch[data.in$lagLoc]
  DataWindy<-data.in[WindyInd,]
  DataCalm<-data.in[CalmInd,]
  data.in$rWSPDind=as.factor(ifelse(data.in$rWSPD<=0.5,1,0))
  
  rWSPDsumWindy<-summary(DataWindy$rWSPD)[c(2,3,5)]
  rWSPDsumCalm<-summary(DataCalm$rWSPD)[c(2,3,5)]
  
  SWHftsumWindy<-c(median(na.omit(DataWindy$SWHft1[which(DataWindy$rWSPD<=rWSPDsumWindy[1])])),
                   median(na.omit(DataWindy$SWHft1[which(DataWindy$rWSPD<=rWSPDsumWindy[2] &   DataWindy$rWSPD>rWSPDsumWindy[1])])),
                   median(na.omit(DataWindy$SWHft1[which(DataWindy$rWSPD<=rWSPDsumWindy[3] &   DataWindy$rWSPD>rWSPDsumWindy[2])])),
                   median(na.omit(DataWindy$SWHft1[which(DataWindy$rWSPD>rWSPDsumWindy[3])])))
  SWHftsumCalm<-c(median(na.omit(DataCalm$SWHft1[which(DataCalm$rWSPD<=rWSPDsumCalm[1])])),
                  median(na.omit(DataCalm$SWHft1[which(DataCalm$rWSPD<=rWSPDsumCalm[2]& DataCalm$rWSPD>rWSPDsumCalm[1])])),
                  median(na.omit(DataCalm$SWHft1[which(DataCalm$rWSPD<=rWSPDsumCalm[3] & DataCalm$rWSPD>rWSPDsumCalm[2])])),
                  median(na.omit(DataCalm$SWHft1[which(DataCalm$rWSPD>rWSPDsumCalm[3])])))
  
  #Imputing Wave Height-SWHft
  #Windy
  data.in$SWHft<-data.in$SWHft1
  DataWindy$SWHft<-DataWindy$SWHft1
  DataCalm$SWHft<-DataCalm$SWHft1
  set.seed(353)
  ZWindySWhft<-rnorm(length(which(is.na(DataWindy$SWHft1)==TRUE)),0,1)
  ZCalmSWhft<-rnorm(length(which(is.na(DataCalm$SWHft1)==TRUE)),0,1)
  sdW1<-rep(0,length(DataWindy$SWHft1))
  sdW1[which(is.na(DataWindy$SWHft1)==TRUE)]<-ZWindySWhft
  indW1<-which(DataWindy$rWSPD<= rWSPDsumWindy[1])
  indW2<-which(DataWindy$rWSPD<= rWSPDsumWindy[2] &DataWindy$rWSPD>rWSPDsumWindy[1])
  indW3<-which(DataWindy$rWSPD<= rWSPDsumWindy[3] &DataWindy$rWSPD> rWSPDsumWindy[2])
  indW4<-which(DataWindy$rWSPD>rWSPDsumWindy[3])
  sdW<-0.03
  DataWindy$SWHft<-DataWindy$SWHft1
  DataWindy$SWHft[indW1]<-na_replace(DataWindy$SWHft[indW1],SWHftsumWindy[1])+(sdW*sdW1[indW1])
  DataWindy$SWHft[indW2]<-na_replace(DataWindy$SWHft[indW2],SWHftsumWindy[2])+(sdW*sdW1[indW2])      
  DataWindy$SWHft[indW3]<-na_replace(DataWindy$SWHft[indW3],SWHftsumWindy[3])+(sdW*sdW1[indW3])      
  DataWindy$SWHft[indW4]<-na_replace(DataWindy$SWHft[indW4],SWHftsumWindy[4])+(sdW*sdW1[indW4])    
  summary(DataWindy$SWHft)
  #sd(na.omit(DataWindy$SWHft1))
  #sd(DataWindy$SWHft)
  
  #Calm
  sdC1<-rep(0,length(DataCalm$SWHft1))
  sdC1[which(is.na(DataCalm$SWHft1)==TRUE)]<-ZCalmSWhft
  indC1<-which(DataCalm$rWSPD<= rWSPDsumCalm[1])
  indC2<-which(DataCalm$rWSPD<= rWSPDsumCalm[2] & DataCalm$rWSPD>rWSPDsumCalm[1])
  indC3<-which(DataCalm$rWSPD <= rWSPDsumCalm[3] & DataCalm$rWSPD > rWSPDsumCalm[2])
  indC4<-which(DataCalm$rWSPD > rWSPDsumCalm[3])
  sdC<-0.01
  DataCalm$SWHft<-DataCalm$SWHft1
  DataCalm$SWHft[indC1]<-na_replace(DataCalm$SWHft[indC1],SWHftsumCalm[1])+(sdC*sdC1[indC1])
  DataCalm$SWHft[indC2]<-na_replace(DataCalm$SWHft[indC2],SWHftsumCalm[2])+(sdC*sdC1[indC2])   
  DataCalm$SWHft[indC3]<-na_replace(DataCalm$SWHft[indC3],SWHftsumCalm[3])+(sdC*sdC1[indC3])       
  DataCalm$SWHft[indC4]<-na_replace(DataCalm$SWHft[indC4],SWHftsumCalm[4])+(sdC*sdC1[indC4])       
  summary(DataCalm$SWHft) 
  #sd(na.omit(DataCalm$SWHft1))
  #sd(DataCalm$SWHft)
  
  data.in$SWHft<-c(DataWindy$SWHft,DataCalm$SWHft)
  data.in$SWHft<-ifelse(data.in$SWHft<0,0,data.in$SWHft)
  #Months= format(ISOdatetime(2000,1:12,1,0,0,0),"%b")
  #Year1<-c(rep(2008,10),rep(2007,2))
  #par(mar = c(2, 2, 2, 2))
  #par(mfrow=c(3,4))
  #for(i in c(11,12,1,2,3,4,5,6,8,9,10))
  #{
  #plot(data.in$SWHft[which(data.in$MM==i)],main=paste(Months[i],Year1[i]),type="l")
  #lines(data.in$SWHft[which(data.in$MM==i)],col=data.in$colSWHft[which(data.in$MM==i)],pch=19,cex=0.5,type="p")
  #}
  
  
  data.in$rATMP = data.in$ATMP[data.in$lagLoc]
  data.in$rWTMP = data.in$WTMP[data.in$lagLoc]
  data.in$rTempDiff = data.in$TempDiff[data.in$lagLoc]
  data.in$idx1 = (1:dim(data.in)[1])
  data.in$idx2 = (1:dim(data.in)[1])
  data.in$idx3 = (1:dim(data.in)[1])
  data.in$idx4 = (1:dim(data.in)[1])
  data.in$idx5 = (1:dim(data.in)[1])
  data.in$idx6 = (1:dim(data.in)[1])
  data.in$idx7 = (1:dim(data.in)[1])
  data.in$idx8 = (1:dim(data.in)[1])
  data.in$idx9 = (1:dim(data.in)[1])            
  
  data.in$rWSPD2 = data.in$rWSPD^2
  #2
  #data.in$rFetch = data.in$Fetch[data.in$lagLoc]
  
  
  goda <- function(u,F){
    H <- 0.3/9.807*u^2*(1-(1+0.004*(sqrt(9.807*F/u^2)))^(-2))
    return (H)
  }
  
  #note that this is in meters not feet like data.in variables elsewhere
  data.in$rGoda <- rep(NA,length(data.in$rWSPD))
  for (t in 1:length(data.in$rWSPD)){
    data.in$rGoda[t] <- goda(data.in$rWSPD[t],data.in$rFetch[t])
  }
  
  
  #added manually
  data.in$rAbsMinor = abs(data.in$rMinor)
  data.in$rAbsMajor = abs(data.in$rMajor)
  data.in$rAbsInt = data.in$rAbsMinor * data.in$rAbsMajor
  
  data.in$rPosMinor = 1/2*(data.in$rMinor + data.in$rAbsMinor)
  data.in$rNegMinor = -1/2*(data.in$rMinor - data.in$rAbsMinor)
  
  data.in$rPosMajor = 1/2*(data.in$rMajor + data.in$rAbsMajor)
  data.in$rNegMajor = -1/2*(data.in$rMajor - data.in$rAbsMajor)
  
  data.keep <- data.in
  
  if ( indSeason=="Annual") {
    data.in <- data.in
  } else if ( indSeason=="Windy") {
    data.in <- data.in[(data.in$MM == 11 | data.in$MM == 12 |data.in$MM == 1|data.in$MM == 2|data.in$MM == 3 ),]
    
  } else if( indSeason=="Calm"){
    data.in <-data.in[(data.in$MM == 4 | data.in$MM == 5 |data.in$MM == 6|data.in$MM == 7|data.in$MM == 8|data.in$MM == 9|data.in$MM == 10),]
  } 
  
  data.ar.exercise <- data.in
  #define lagged wave height variables for lm() call
  data.ar.exercise$lag1SWHft = c(rep(NA,(0+1)),data.in$SWHft[(1:(length(data.in$SWHft)-1))])
  data.ar.exercise$lag2SWHft = c(rep(NA,(0+2)),data.in$SWHft[(1:(length(data.in$SWHft)-2))])
  data.ar.exercise$lag3SWHft = c(rep(NA,(0+3)),data.in$SWHft[(1:(length(data.in$SWHft)-3))])
  data.ar.exercise$lag4SWHft = c(rep(NA,(0+4)),data.in$SWHft[(1:(length(data.in$SWHft)-4))])
  data.ar.exercise$lag5SWHft = c(rep(NA,(0+5)),data.in$SWHft[(1:(length(data.in$SWHft)-5))])
  data.ar.exercise$lag6SWHft = c(rep(NA,(0+6)),data.in$SWHft[(1:(length(data.in$SWHft)-6))])
  data.ar.exercise$lag21SWHft = c(rep(NA,(0+21)),data.in$SWHft[(1:(length(data.in$SWHft)-21))])
  data.ar.exercise$lag22SWHft = c(rep(NA,(0+22)),data.in$SWHft[(1:(length(data.in$SWHft)-22))])
  data.ar.exercise$lag23SWHft = c(rep(NA,(0+23)),data.in$SWHft[(1:(length(data.in$SWHft)-23))])
  data.ar.exercise$lag24SWHft = c(rep(NA,(0+24)),data.in$SWHft[(1:(length(data.in$SWHft)-24))])
  
  data.ar.exercise$rWSPD2 = data.ar.exercise$rWSPD^2
  data.ar.exercise$rWSPD_ind9=as.numeric(data.ar.exercise$event_ind)-1
  data.ar.exercise$rWSPD_ind9[which(is.na(data.ar.exercise$rWSPD)==TRUE)]<-NA
  #data.ar.exercise$rWSPD_ind9<-rep(0,length(data.ar.exercise$rWSPD))
  #indy<-rle(is.na(ifelse(data.ar.exercise$rWSPD>=9,NA,0)))
  #indY1<-rep(ifelse(indy$values,indy$lengths,0),times = indy$lengths)
  #length(which(indY1>5))
  #data.ar.exercise$rWSPD_ind9<-ifelse(indY1>=2,1,0)
  #data.ar.exercise$rWSPD_ind9[which(is.na(data.ar.exercise$rWSPD)==TRUE)]<-NA
  if(indSeason=="Windy")
  {
    ind1<-which(data.ar.exercise$WDIR>=60 &data.ar.exercise$WDIR<=120)
    ind2<-which(data.ar.exercise$WDIR>=240 &data.ar.exercise$WDIR<=300)
    ind31<-which(data.ar.exercise$WDIR>=0 &data.ar.exercise$WDIR< 60)
    ind32<-which(data.ar.exercise$WDIR>120 &data.ar.exercise$WDIR<240)
    ind33<-which(data.ar.exercise$WDIR>300&data.ar.exercise$WDIR<=360)
    ind3<-c(ind31,ind32,ind33)
    data.ar.exercise$rWSPD_T<-rep(NA,nrow(data.ar.exercise))
    data1_rWSPD<-data.frame(data.ar.exercise$rWSPD[ind1])
    data2_rWSPD<-data.frame(data.ar.exercise$rWSPD[ind2])
    data3_rWSPD<-data.frame(data.ar.exercise$rWSPD[ind3])
    colnames(data1_rWSPD)=colnames(data2_rWSPD)=colnames(data3_rWSPD)="rWSPD"
    
    data.ar.exercise$rWSPD_T[ind1]<- 1.3503+(0.8429*data.ar.exercise$rWSPD[ind1])+(0.8338*data.ar.exercise$rWSPD_ind9[ind1])
    data.ar.exercise$rWSPD_T[ind2]<- 4.7836 +(0.6602 *data.ar.exercise$rWSPD[ind2])+(-0.6969*data.ar.exercise$rWSPD_ind9[ind2])
    data.ar.exercise$rWSPD_T[ind3]<- 3.68845 +(0.62530*data.ar.exercise$rWSPD[ind3])+(-0.08779*data.ar.exercise$rWSPD_ind9[ind3])
    
    

  }else{
    ind1<-which(data.ar.exercise$WDIR>=60 &data.ar.exercise$WDIR<=120)
    ind2<-which(data.ar.exercise$WDIR>=240 &data.ar.exercise$WDIR<=300)
    ind31<-which(data.ar.exercise$WDIR>=0 &data.ar.exercise$WDIR< 60)
    ind32<-which(data.ar.exercise$WDIR>120 &data.ar.exercise$WDIR<240)
    ind33<-which(data.ar.exercise$WDIR>300&data.ar.exercise$WDIR<=360)
    ind3<-c(ind31,ind32,ind33)
    data.ar.exercise$rWSPD_T<-rep(NA,nrow(data.ar.exercise))
    data1_rWSPD<-data.frame(data.ar.exercise$rWSPD[ind1])
    data2_rWSPD<-data.frame(data.ar.exercise$rWSPD[ind2])
    data3_rWSPD<-data.frame(data.ar.exercise$rWSPD[ind3])
    colnames(data1_rWSPD)=colnames(data2_rWSPD)=colnames(data3_rWSPD)="rWSPD"
    
    data.ar.exercise$rWSPD_T[ind1]<- 1.2087+(0.8209*data.ar.exercise$rWSPD[ind1])+(0.3530*data.ar.exercise$rWSPD_ind9[ind1])
    data.ar.exercise$rWSPD_T[ind2]<- 3.9107 +(0.4092 *data.ar.exercise$rWSPD[ind2])+(-0.5313*data.ar.exercise$rWSPD_ind9[ind2])
    data.ar.exercise$rWSPD_T[ind3]<- 1.8282 +(0.5289*data.ar.exercise$rWSPD[ind3])+(0.8880*data.ar.exercise$rWSPD_ind9[ind3])
    
    
  }
  data.ar.exercise$rWSPD_T<-ifelse(data.ar.exercise$rWSPD_T<=0.25,0.25,data.ar.exercise$rWSPD_T)
  data.ar.exercise$rWSPD2_T<-(data.ar.exercise$rWSPD_T)^2
  #MacLaurin terms
  data.ar.exercise$z2 <- data.ar.exercise$rWSPD2*(sqrt(data.ar.exercise$rFetch/data.ar.exercise$rWSPD2))^1
  data.ar.exercise$z3 <- data.ar.exercise$rWSPD2*(sqrt(data.ar.exercise$rFetch/data.ar.exercise$rWSPD2))^2
  data.ar.exercise$z4 <- data.ar.exercise$rWSPD2*(sqrt(data.ar.exercise$rFetch/data.ar.exercise$rWSPD2))^3
  
  data.ar.exercise$z2_T <- data.ar.exercise$rWSPD2_T*(sqrt(data.ar.exercise$rFetch/data.ar.exercise$rWSPD2_T))^1
  data.ar.exercise$z3_T <- data.ar.exercise$rWSPD2_T*(sqrt(data.ar.exercise$rFetch/data.ar.exercise$rWSPD2_T))^2
  data.ar.exercise$z4_T <- data.ar.exercise$rWSPD2_T*(sqrt(data.ar.exercise$rFetch/data.ar.exercise$rWSPD2_T))^3
  
  #define lagged wave height variables  for speed^2,z2,z3 and z4 for lm() call
  data.ar.exercise$lag1rWSPD2 = c(rep(NA,(0+1)),data.in$rWSPD2[(1:(length(data.in$rWSPD2)-1))])
  data.ar.exercise$lag2rWSPD2 = c(rep(NA,(0+2)),data.in$rWSPD2[(1:(length(data.in$rWSPD2)-2))])
  data.ar.exercise$lag3rWSPD2 = c(rep(NA,(0+3)),data.in$rWSPD2[(1:(length(data.in$rWSPD2)-3))])
  data.ar.exercise$lag4rWSPD2 = c(rep(NA,(0+4)),data.in$rWSPD2[(1:(length(data.in$rWSPD2)-4))])
  data.ar.exercise$lag5rWSPD2 = c(rep(NA,(0+5)),data.in$rWSPD2[(1:(length(data.in$rWSPD2)-5))])
  data.ar.exercise$lag6rWSPD2 = c(rep(NA,(0+6)),data.in$rWSPD2[(1:(length(data.in$rWSPD2)-6))])
  
  ##
  data.ar.exercise$lag1rMinor = c(rep(NA,(0+1)),data.in$rMinor[(1:(length(data.in$rMinor)-1))])
  data.ar.exercise$lag2rMinor = c(rep(NA,(0+2)),data.in$rMinor[(1:(length(data.in$rMinor)-2))])
  data.ar.exercise$lag3rMinor = c(rep(NA,(0+3)),data.in$rMinor[(1:(length(data.in$rMinor)-3))])
  data.ar.exercise$lag4rMinor = c(rep(NA,(0+4)),data.in$rMinor[(1:(length(data.in$rMinor)-4))])
  data.ar.exercise$lag5rMinor = c(rep(NA,(0+5)),data.in$rMinor[(1:(length(data.in$rMinor)-5))])
  data.ar.exercise$lag6rMinor = c(rep(NA,(0+6)),data.in$rMinor[(1:(length(data.in$rMinor)-6))])
  ##
  data.ar.exercise$lag1rMajor = c(rep(NA,(0+1)),data.in$rMajor[(1:(length(data.in$rMajor)-1))])
  data.ar.exercise$lag2rMajor = c(rep(NA,(0+2)),data.in$rMajor[(1:(length(data.in$rMajor)-2))])
  data.ar.exercise$lag3rMajor = c(rep(NA,(0+3)),data.in$rMajor[(1:(length(data.in$rMajor)-3))])
  data.ar.exercise$lag4rMajor = c(rep(NA,(0+4)),data.in$rMajor[(1:(length(data.in$rMajor)-4))])
  data.ar.exercise$lag5rMajor = c(rep(NA,(0+5)),data.in$rMajor[(1:(length(data.in$rMajor)-5))])
  data.ar.exercise$lag6rMajor = c(rep(NA,(0+6)),data.in$rMajor[(1:(length(data.in$rMajor)-6))])
  
  #
  data.ar.exercise$lag1z2 = c(rep(NA,(0+1)),data.ar.exercise$z2[(1:(length(data.ar.exercise$z2)-1))])
  data.ar.exercise$lag2z2 = c(rep(NA,(0+2)),data.ar.exercise$z2[(1:(length(data.ar.exercise$z2)-2))])
  data.ar.exercise$lag3z2 = c(rep(NA,(0+3)),data.ar.exercise$z2[(1:(length(data.ar.exercise$z2)-3))])
  data.ar.exercise$lag4z2 = c(rep(NA,(0+4)),data.ar.exercise$z2[(1:(length(data.ar.exercise$z2)-4))])
  data.ar.exercise$lag5z2 = c(rep(NA,(0+5)),data.ar.exercise$z2[(1:(length(data.ar.exercise$z2)-5))])
  data.ar.exercise$lag6z2 = c(rep(NA,(0+6)),data.ar.exercise$z2[(1:(length(data.ar.exercise$z2)-6))])
  ##
  data.ar.exercise$lag1z3 = c(rep(NA,(0+1)),data.ar.exercise$z3[(1:(length(data.ar.exercise$z3)-1))])
  data.ar.exercise$lag2z3 = c(rep(NA,(0+2)),data.ar.exercise$z3[(1:(length(data.ar.exercise$z3)-2))])
  data.ar.exercise$lag3z3 = c(rep(NA,(0+3)),data.ar.exercise$z3[(1:(length(data.ar.exercise$z3)-3))])
  data.ar.exercise$lag4z3 = c(rep(NA,(0+4)),data.ar.exercise$z3[(1:(length(data.ar.exercise$z3)-4))])
  data.ar.exercise$lag5z3 = c(rep(NA,(0+5)),data.ar.exercise$z3[(1:(length(data.ar.exercise$z3)-5))])
  data.ar.exercise$lag6z3 = c(rep(NA,(0+6)),data.ar.exercise$z3[(1:(length(data.ar.exercise$z3)-6))])
  ##
  data.ar.exercise$lag1z4 = c(rep(NA,(0+1)),data.ar.exercise$z4[(1:(length(data.ar.exercise$z4)-1))])
  data.ar.exercise$lag2z4 = c(rep(NA,(0+2)),data.ar.exercise$z4[(1:(length(data.ar.exercise$z4)-2))])
  data.ar.exercise$lag3z4 = c(rep(NA,(0+3)),data.ar.exercise$z4[(1:(length(data.ar.exercise$z4)-3))])
  data.ar.exercise$lag4z4 = c(rep(NA,(0+4)),data.ar.exercise$z4[(1:(length(data.ar.exercise$z4)-4))])
  data.ar.exercise$lag5z4 = c(rep(NA,(0+5)),data.ar.exercise$z4[(1:(length(data.ar.exercise$z4)-5))])
  data.ar.exercise$lag6z4 = c(rep(NA,(0+6)),data.ar.exercise$z4[(1:(length(data.ar.exercise$z4)-6))])
  return(data.ar.exercise)
}

save("data.in_CS", file="data_fn_CS.Rdata")



