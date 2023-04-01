library(lubridate)
library(dplyr)
#For Wave-height
#filled<-merged
#filled$Timestamp<-ymd_hms(filled$Timestamp)
#filled$Timestamp1<-floor_date(filled$Timestamp, "hour")
Data1<-data.frame(seq(ymd_hms("2004-1-1 00:00:00"),ymd_hms(" 2013-12-31 23:00:00"),by="hours"))
colnames(Data1)<-"TimeStamp"
head(Data1)
dim(Data1)
Data1$event_ind<-rep(0,nrow(Data1))
#Data1

EventData <-read.csv("~/Dropbox/Namitha Pais/Marine Sciences/Storm_Data_Filtered_2004_2013.csv")
EventData$Timestamp <- as.POSIXct(EventData$Timestamp , format = "%m/%d/%y %H:%M")
EventData$Timestamp <-round(EventData$Timestamp , "hours")
EventData$Timestamp  <- ymd_hms(format(EventData$Timestamp , "%Y-%m-%d %H:%M:%S"))
EventData_sorted <- EventData[order(EventData$Timestamp), ]
indices<-which(Data1$TimeStamp %in% EventData_sorted$Timestamp)
head(Data1$TimeStamp[(indices)])
head(EventData_sorted$Timestamp)

length(unique(EventData_sorted$Timestamp))
length(Data1$TimeStamp[(indices)])
for(i in 1:length(indices))
{
  Data1$event_ind[ (indices[i]-12):(indices[i]+12)]<- 1
}

Data1$event_ind<-as.factor(Data1$event_ind)
Data_event_ind<-Data1
Data_event_ind$event_ind<-as.factor(Data_event_ind$event_ind)
summary(Data_event_ind)



Y<-c()
X<-c()
X_WDIR<-c()
Y_WDIR<-c()
X_SWHft<-c()
Y_SWHft<-c()
event_ind<-c()
for(j in 2005:2013)
{
  setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
  load("44039_reg_merged.Rda")
  source("Data_Function_CS.R")
  Data_CS<-data.in_CS(j,"Windy")
  X<-c(X,Data_CS$rWSPD)
  X_WDIR<-c(X_WDIR,Data_CS$WDIR)
  X_SWHft<-c(X_SWHft,Data_CS$SWHft)
  event_ind<-c(event_ind,Data_CS$event_ind)
  setwd('/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs/Threshold_Model')
  load("44039_reg_merged.Rda")
  source("Data_Function.R")
  Data_Buoy<-data.in("central",j,"Windy")
  Y<-c(Y,Data_Buoy$rWSPD)
  Y_WDIR<-c(Y_WDIR,Data_Buoy$WDIR)
  Y_SWHft<-c(Y_SWHft,Data_Buoy$SWHft)
}

Y_ind<-as.factor(event_ind-1)
Y_ind[which(is.na(Y)==TRUE)]<-NA
summary(Y_ind)

ind1<-which(Y_WDIR>=60 &Y_WDIR<=120 )
ind11<-which(X_WDIR>=60 &X_WDIR<=120 )
a1<-intersect(ind1,ind11)
length(a1)
length(ind1)
summary(Y_SWHft[ind1])



ind2<-which(Y_WDIR>=240 &Y_WDIR<=300 )
ind21<-which(X_WDIR>=240 &X_WDIR<=300 )
a2<-intersect(ind2,ind21)
summary(Y_SWHft[ind2])

ind3_1<-which(Y_WDIR>=0 &Y_WDIR< 60)
ind3_2<-which(Y_WDIR>120 &Y_WDIR<240)
ind3_3<-which(Y_WDIR>300&Y_WDIR<=360)
ind3<-c(ind3_1,ind3_2,ind3_3)

ind31_1<-which(X_WDIR>=0 &X_WDIR< 60)
ind31_2<-which(X_WDIR>120 &X_WDIR<240)
ind31_3<-which(X_WDIR>300&X_WDIR<=360)
ind31<-c(ind31_1,ind31_2,ind31_3)
a3<-intersect(ind3,ind31)

#When WDIR in CS is from East to West(60-120)
plot(Y[ind1],X[ind1],pch=19,ylim=c(0,21),xlim=c(0,21),cex=0.5)

#When WDIR in CS  and Buoy is from East to West(60-120)
#plot(Y[a1],X[a1],pch=19,ylim=c(0,21),xlim=c(0,21),cex=0.5)


#When WDIR in CS is from West to East (240-300)
plot(Y[ind2],X[ind2],pch=19,ylim=c(0,21),xlim=c(0,21),cex=0.5)

#When WDIR in CS and Buoy is from West to East (240-300)
#plot(Y[a2],X[a2],pch=19,ylim=c(0,21),xlim=c(0,21),cex=0.5)


plot(X,Y)

#When WDIR in CS is from "Other"
plot(X[ind3],Y[ind3],pch=19,ylim=c(0,21),xlim=c(0,21),cex=0.5)

#When WDIR in CS and Buoy is from "Other" 
#plot(Y[a3],X[a3],pch=19,ylim=c(0,21),xlim=c(0,21),cex=0.5)

#CS_T= a+b*Cs+ c*
#Time 

library(mgcv)
library(ggplot2)
#summary(X[ind1])
#When WDIR in CS is from East to West(60-120)
Y1<-na.omit(Y[ind1])
X1<-na.omit(X[ind1])
Y_ind1<-na.omit(Y_ind[ind1])
Data1<-data.frame(cbind(Y1,X1,Y_ind1))

colnames(Data1)<-c("rWSPD_T","rWSPD","rWSPD_T_ind")
Model_Windy11<-lm(rWSPD_T ~ rWSPD+ rWSPD_T_ind,data=Data1)
summary(Model_Windy11)
ggplot(Data1, aes(rWSPD, rWSPD_T,rWSPD_T_ind) ) +
  geom_point() +
  #stat_smooth(method = gam, formula = y ~ s(x,k=5))+
  stat_smooth(method = lm, formula = y ~ x,color="darkgreen")

plot(Data1$rWSPD_T,Model_Windy11$fitted.values,pch=19,cex=0.3)

#getwd()
#saveRDS(Model_Windy1,"GAM_Windy1.RDA")

#predictions <- model %>% predict(test.data)
#When WDIR in CS is from West to East (240-300)
Y2<-na.omit(Y[ind2])
X2<-na.omit(X[ind2])
Y_ind2<-na.omit(Y_ind[ind2])
Data2<-data.frame(cbind(Y2,X2,Y_ind2))
colnames(Data2)<-c("rWSPD_T","rWSPD","rWSPD_T_ind")
Model_Windy21<-lm(rWSPD_T ~ rWSPD+ rWSPD_T_ind,data=Data2)
summary(Model_Windy21)
ggplot(Data2, aes(rWSPD, rWSPD_T) ) +
  geom_point() +
  stat_smooth(method = gam, formula = y ~ s(x,k=5))+
  stat_smooth(method = lm, formula = y ~ x,color="darkgreen")
plot(Data2$rWSPD_T,Model_Windy21$fitted.values,pch=19)

#unique(Y2[which(Y2<=5)])
#When WDIR in CS is from "Other"
Y3<-na.omit(Y[ind3])
X3<-na.omit(X[ind3])
Y_ind3<-na.omit(Y_ind[ind3])
Data3<-data.frame(cbind(Y3,X3,Y_ind3))
colnames(Data3)<-c("rWSPD_T","rWSPD","rWSPD_T_ind")
Model_Windy31<-lm(rWSPD_T ~ rWSPD+rWSPD_T_ind,data=Data3)
summary(Model_Windy31)
ggplot(Data3, aes(rWSPD, rWSPD_T)) +
  geom_point() +
  stat_smooth(method = gam, formula = y ~ s(x,k=5))+
  stat_smooth(method = lm, formula = y ~ x,color="darkgreen")
plot(Data3$rWSPD_T,Model_Windy31$fitted.values,pch=19)


summary(Model_Windy11)
Model_Windy11
summary(Model_Windy21)
Model_Windy21
#5.6491       0.6595      -1.5610  
summary(Model_Windy31)
#3.7458       0.6245      -0.1429  
Model_Windy31


Model_Windy11
Model_Windy21
Model_Windy31


Y<-c()
X<-c()
X_WDIR<-c()
Y_WDIR<-c()
X_SWHft<-c()
Y_SWHft<-c()
event_ind<-c()
for(j in 2005:2013)
{
  setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
  load("44039_reg_merged.Rda")
  source("Data_Function_CS.R")
  Data_CS<-data.in_CS(j,"Calm")
  X<-c(X,Data_CS$rWSPD)
  X_WDIR<-c(X_WDIR,Data_CS$WDIR)
  X_SWHft<-c(X_SWHft,Data_CS$SWHft)
  event_ind<-c(event_ind,Data_CS$event_ind)
  setwd('/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs/Threshold_Model')
  load("44039_reg_merged.Rda")
  source("Data_Function.R")
  Data_Buoy<-data.in("central",j,"Calm")
  Y<-c(Y,Data_Buoy$rWSPD)
  Y_WDIR<-c(Y_WDIR,Data_Buoy$WDIR)
  Y_SWHft<-c(Y_SWHft,Data_Buoy$SWHft)
}

Y_ind<-as.factor(event_ind-1)
Y_ind[which(is.na(Y)==TRUE)]<-NA
summary(Y_ind)

ind1<-which(Y_WDIR>=60 &Y_WDIR<=120 )
ind11<-which(X_WDIR>=60 &X_WDIR<=120 )
a1<-intersect(ind1,ind11)
length(a1)
length(ind1)
summary(Y_SWHft[ind1])



ind2<-which(Y_WDIR>=240 &Y_WDIR<=300 )
ind21<-which(X_WDIR>=240 &X_WDIR<=300 )
a2<-intersect(ind2,ind21)
summary(Y_SWHft[ind2])

ind3_1<-which(Y_WDIR>=0 &Y_WDIR< 60)
ind3_2<-which(Y_WDIR>120 &Y_WDIR<240)
ind3_3<-which(Y_WDIR>300&Y_WDIR<=360)
ind3<-c(ind3_1,ind3_2,ind3_3)

ind31_1<-which(X_WDIR>=0 &X_WDIR< 60)
ind31_2<-which(X_WDIR>120 &X_WDIR<240)
ind31_3<-which(X_WDIR>300&X_WDIR<=360)
ind31<-c(ind31_1,ind31_2,ind31_3)
a3<-intersect(ind3,ind31)

#When WDIR in CS is from East to West(60-120)
plot(Y[ind1],X[ind1],pch=19,ylim=c(0,21),xlim=c(0,21),cex=0.5)

#When WDIR in CS  and Buoy is from East to West(60-120)
#plot(Y[a1],X[a1],pch=19,ylim=c(0,21),xlim=c(0,21),cex=0.5)


#When WDIR in CS is from West to East (240-300)
plot(Y[ind2],X[ind2],pch=19,ylim=c(0,21),xlim=c(0,21),cex=0.5)

#When WDIR in CS and Buoy is from West to East (240-300)
#plot(Y[a2],X[a2],pch=19,ylim=c(0,21),xlim=c(0,21),cex=0.5)


plot(X,Y)

#When WDIR in CS is from "Other"
plot(X[ind3],Y[ind3],pch=19,ylim=c(0,21),xlim=c(0,21),cex=0.5)

#When WDIR in CS and Buoy is from "Other" 
#plot(Y[a3],X[a3],pch=19,ylim=c(0,21),xlim=c(0,21),cex=0.5)

#CS_T= a+b*Cs+ c*
#Time 

library(mgcv)
library(ggplot2)
#summary(X[ind1])
#When WDIR in CS is from East to West(60-120)
Y1<-na.omit(Y[ind1])
X1<-na.omit(X[ind1])
Y_ind1<-na.omit(Y_ind[ind1])
Data1<-data.frame(cbind(Y1,X1,Y_ind1))

colnames(Data1)<-c("rWSPD_T","rWSPD","rWSPD_T_ind")
Model_Calm11<-lm(rWSPD_T ~ rWSPD+ rWSPD_T_ind,data=Data1)
summary(Model_Calm11)
ggplot(Data1, aes(rWSPD, rWSPD_T,rWSPD_T_ind) ) +
  geom_point() +
  #stat_smooth(method = gam, formula = y ~ s(x,k=5))+
  stat_smooth(method = lm, formula = y ~ x,color="darkgreen")

plot(Data1$rWSPD_T,Model_Calm11$fitted.values,pch=19,cex=0.3)

#getwd()
#saveRDS(Model_Calm1,"GAM_Calm1.RDA")

#predictions <- model %>% predict(test.data)
#When WDIR in CS is from West to East (240-300)
Y2<-na.omit(Y[ind2])
X2<-na.omit(X[ind2])
Y_ind2<-na.omit(Y_ind[ind2])
Data2<-data.frame(cbind(Y2,X2,Y_ind2))
colnames(Data2)<-c("rWSPD_T","rWSPD","rWSPD_T_ind")
Model_Calm21<-lm(rWSPD_T ~ rWSPD+ rWSPD_T_ind,data=Data2)
summary(Model_Calm21)
ggplot(Data2, aes(rWSPD, rWSPD_T) ) +
  geom_point() +
  stat_smooth(method = gam, formula = y ~ s(x,k=5))+
  stat_smooth(method = lm, formula = y ~ x,color="darkgreen")
plot(Data2$rWSPD_T,Model_Calm21$fitted.values,pch=19)

#unique(Y2[which(Y2<=5)])
#When WDIR in CS is from "Other"
Y3<-na.omit(Y[ind3])
X3<-na.omit(X[ind3])
Y_ind3<-na.omit(Y_ind[ind3])
Data3<-data.frame(cbind(Y3,X3,Y_ind3))
colnames(Data3)<-c("rWSPD_T","rWSPD","rWSPD_T_ind")
Model_Calm31<-lm(rWSPD_T ~ rWSPD+rWSPD_T_ind,data=Data3)
summary(Model_Calm31)
ggplot(Data3, aes(rWSPD, rWSPD_T)) +
  geom_point() +
  stat_smooth(method = gam, formula = y ~ s(x,k=5))+
  stat_smooth(method = lm, formula = y ~ x,color="darkgreen")
plot(Data3$rWSPD_T,Model_Calm31$fitted.values,pch=19)


summary(Model_Calm11)
Model_Calm11
summary(Model_Calm21)
Model_Calm21
#5.6491       0.6595      -1.5610  
summary(Model_Calm31)
#3.7458       0.6245      -0.1429  
Model_Calm31


Model_Calm11
Model_Calm21
Model_Calm31




