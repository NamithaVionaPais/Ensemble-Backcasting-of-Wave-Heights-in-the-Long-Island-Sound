setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
load("44039_reg_merged.Rda")
source("Data_Function_CS.R")
library(fGarch)
library(chngpt)
library(Metrics)
#rmse(data$actual, data$predicted)
MonthWindy<-c("Nov","Dec","Jan","Feb","Mar")
MonthWindyNum<-c(11,12,1,2,3)
ED_TM_Windy<-matrix(nrow=5,ncol=9)
row.names(ED_TM_Windy)<-MonthWindy
colnames(ED_TM_Windy)<-2005:2013
for(j in 2008)
{
  setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
  load("44039_reg_merged.Rda")
  source("Data_Function_CS.R")
  Data<-data.in_CS(j,"Windy")
  #Model 1
  Data$lag1rWSPD2_T<-c(rep(NA,(0+1)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-1))])
  Data$lag2rWSPD2_T<-c(rep(NA,(0+2)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-2))])
  Data$lag3rWSPD2_T<-c(rep(NA,(0+3)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-3))])
  Data$lag4rWSPD2_T<-c(rep(NA,(0+4)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-4))])
  Data$rWSPDind_T<-as.factor(ifelse(Data$rWSPD<= 0.25,1,0))
  formula_1<-SWHft ~ lag1SWHft+lag2SWHft+lag3SWHft+lag4SWHft+
    lag1rWSPD2_T + lag2rWSPD2_T + lag3rWSPD2_T +lag4rWSPD2_T+
    z2_T+ z3_T +z4_T+
    I(lag1SWHft*rWSPD2_T)+
    I(lag3SWHft*lag2rWSPD2_T)+
    rWSPDind_T
  formula_2<-~ rWSPD2_T
  runchngpt<-chngptm(formula.1= formula_1, formula.2=formula_2 , family="gaussian", data=Data, 
                       type="segmented", var.type="bootstrap", weights=NULL,ci.bootstrap.size = 100)
  ResWindy<-residuals(runchngpt)
  ModelRes_Windy<-garchFit(~garch(1,2), data = ResWindy,trace=FALSE,include.mean = FALSE)
  Ind_Res_Windy<-ifelse(ResWindy<0,-1,1)
  sigW<-ModelRes_Windy@sigma.t
  RhatWindy<-(sigW*Ind_Res_Windy)
  #PredInt<-matrix(nrow=nrow(Data),ncol=3)
  #PredInt[1:5,]<-rep(0,15)
  #PredInt[6:nrow(Data),]<-predict(runchngpt,interval="predict")
  #PredInt[6:nrow(Data),1]<-predict(runchngpt)+RhatWindy
  PredInt<-c()
  PredInt[1:5]<-rep(0,5)
  PredInt[6:nrow(Data)]<-predict(runchngpt)+RhatWindy
  euclidean <- function(a, b)
  {
    return(sqrt(sum((a - b)^2))/length(a))
  }
  EDataTM=c()
  i=1
  EDataTM[i]<-rmse(Data$SWHft[which(Data$MM==MonthWindyNum[i])][-(1:5)],PredInt[which(Data$MM==MonthWindyNum[i])][-(1:5)])
  for(i in 2:length(MonthWindy))
  {
    EDataTM[i]<-rmse(Data$SWHft[which(Data$MM==MonthWindyNum[i])],PredInt[which(Data$MM==MonthWindyNum[i])])
  }
  ED_TM_Windy[,(j-2004)]<-EDataTM
}

MonthWindy<-c("Nov","Dec","Jan","Feb","Mar")
MonthWindyNum<-c(11,12,1,2,3)
j=2008
par(mar = c(2, 2, 2, 2))
par(mfrow=c(3,2))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:3, ylim=0:3)
legend("center",c("Observed-SWHft","Predicted-SWHft","NA-Indicator"),col=c("black","indianred3","yellow2"),pch=19)
i=1
plot(Data$TimeStamp_1[which(Data$MM==MonthWindyNum[i])][-(1:5)],Data$SWHft[which(Data$MM==MonthWindyNum[i])][-(1:5)],type="l",
     ylim=c(0,10),main =paste(MonthWindy[i],(j-1),sep="-"),xaxt='n')
lines(Data$TimeStamp_1[which(Data$MM==MonthWindyNum[i])][-(1:5)],
      Data$SWHft[which(Data$MM==MonthWindyNum[i])][-(1:5)],pch=19,cex=0.2,type="p",
      col=Data$colSWHft[which(Data$MM==MonthWindyNum[i])][-(1:5)])

lines(Data$TimeStamp_1[which(Data$MM==MonthWindyNum[i])][-(1:5)],PredInt[which(Data$MM==MonthWindyNum[i])][-(1:5)],
      col="indianred3")
lines(Data$TimeStamp_1[which(Data$MM==MonthWindyNum[i])][-(1:5)],PredInt[which(Data$MM==MonthWindyNum[i])][-(1:5)],
      col="indianred3",pch=20,cex=0.1,type="p")

i=2
plot(Data$TimeStamp_1[which(Data$MM==MonthWindyNum[i])],Data$SWHft[which(Data$MM==MonthWindyNum[i])],type="l",
     ylim=c(0,10),main =paste(MonthWindy[i],(j-1),sep="-"),xaxt='n')
lines(Data$TimeStamp_1[which(Data$MM==MonthWindyNum[i])],Data$SWHft[which(Data$MM==MonthWindyNum[i])],
      pch=19,cex=0.2,type="p",col=Data$colSWHft[which(Data$MM==MonthWindyNum[i])])
lines(Data$TimeStamp_1[which(Data$MM==MonthWindyNum[i])],PredInt[which(Data$MM==MonthWindyNum[i])],
      col="indianred3")
lines(Data$TimeStamp_1[which(Data$MM==MonthWindyNum[i])],PredInt[which(Data$MM==MonthWindyNum[i])],
      col="indianred3",pch=20,cex=0.1,type="p")

for(i in 3:length(MonthWindy))
{
  plot(Data$TimeStamp_1[which(Data$MM==MonthWindyNum[i])],Data$SWHft[which(Data$MM==MonthWindyNum[i])],type="l",
       ylim=c(0,10),main =paste(MonthWindy[i],j,sep="-"),xaxt='n')
  lines(Data$TimeStamp_1[which(Data$MM==MonthWindyNum[i])],Data$SWHft[which(Data$MM==MonthWindyNum[i])],
        pch=19,cex=0.2,type="p",col=Data$colSWHft[which(Data$MM==MonthWindyNum[i])])
  lines(Data$TimeStamp_1[which(Data$MM==MonthWindyNum[i])],PredInt[which(Data$MM==MonthWindyNum[i])],
        col="indianred3")
  lines(Data$TimeStamp_1[which(Data$MM==MonthWindyNum[i])],PredInt[which(Data$MM==MonthWindyNum[i])],
        col="indianred3",pch=20,cex=0.1,type="p")
}



MonthCalm<-c("Apr","May","Jun","Jul","Aug","Sep","Oct")
MonthCalmNum<-c(4,5,6,7,8,9,10)
ED_TM_Calm<-matrix(nrow=7,ncol=9)
row.names(ED_TM_Calm)<-MonthCalm
{
  setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
  load("44039_reg_merged.Rda")
  source("Data_Function_CS.R")
  Data<-data.in_CS(j,"Calm")
  #Model 1
  Data$lag1rWSPD2_T<-c(rep(NA,(0+1)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-1))])
  Data$lag2rWSPD2_T<-c(rep(NA,(0+2)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-2))])
  Data$lag3rWSPD2_T<-c(rep(NA,(0+3)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-3))])
  Data$lag4rWSPD2_T<-c(rep(NA,(0+4)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-4))])
  Data$rWSPDind_T<-as.factor(ifelse(Data$rWSPD<= 0.25,1,0))
  formula_1<-SWHft ~ lag1SWHft+lag2SWHft+lag3SWHft+lag4SWHft+
    lag1rWSPD2_T + lag2rWSPD2_T + lag3rWSPD2_T +lag4rWSPD2_T+
    z2_T+ z3_T +z4_T+
    I(lag1SWHft*rWSPD2_T)+
    I(lag3SWHft*lag2rWSPD2_T)+
    rWSPDind_T
  formula_2<-~ rWSPD2_T
  runchngpt<-chngptm(formula.1= formula_1, formula.2=formula_2 , family="gaussian", data=Data, 
                     type="segmented", var.type="bootstrap", weights=NULL,ci.bootstrap.size = 100)
  ResCalm<-residuals(runchngpt)
  
  #ModelRes_Calm<-garchFit(~arma(6,0)+garch(1,2), data = ResCalm,trace=FALSE,include.mean = FALSE)
  ModelRes_Calm<-garchFit(~garch(1,2), data = ResCalm,trace=FALSE,include.mean = FALSE)
  Ind_Res_Calm<-ifelse(ResCalm<0,-1,1)
  sigC<-ModelRes_Calm@sigma.t
  RhatCalm<-(sigC*Ind_Res_Calm)
  #PredInt<-matrix(nrow=nrow(Data),ncol=3)
  #PredInt[1:4,]<-rep(0,12)
  #PredInt[5:nrow(Data),]<-predict(runchngpt,interval="predict")
  #PredInt[5:nrow(Data),1]<-predict(runchngpt)+RhatCalm
  PredInt<-c()
  PredInt[1:4]<-rep(0,4)
  PredInt[5:nrow(Data)]<-predict(runchngpt)+RhatCalm
  EDataTM=c()
  i=1
  EDataTM[i]<-rmse(Data$SWHft[which(Data$MM==MonthCalmNum[i])][-(1:5)],PredInt[which(Data$MM==MonthCalmNum[i])][-(1:5)])
  for(i in 2:length(MonthCalm))
  {
    EDataTM[i]<-rmse(Data$SWHft[which(Data$MM==MonthCalmNum[i])],PredInt[which(Data$MM==MonthCalmNum[i])])
  }
  ED_TM_Calm[,(j-2004)]<-EDataTM
}


MonthCalm<-c("Apr","May","Jun","Jul","Aug","Sep","Oct")
MonthCalmNum<-c(4,5,6,7,8,9,10)
par(mar = c(2, 2, 2, 2))
par(mfrow=c(4,2))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:3, ylim=0:3)
legend("center",c("Observed-SWHft","Predicted-SWHft","NA-Indicator"),col=c("black","indianred3","yellow2"),pch=19)
i=1
plot(Data$SWHft[which(Data$MM==MonthCalmNum[i])][-(1:4)],
     type="l",ylim=c(0,10),main =paste(MonthCalm[i],j,sep="-"),xaxt="n")
lines(Data$SWHft[which(Data$MM==MonthCalmNum[i])][-(1:4)],col=Data$colSWHft[which(Data$MM==MonthCalmNum[i])][-(1:5)],pch=19,cex=0.2,type="p")
lines(PredInt[which(Data$MM==MonthCalmNum[i])][-(1:4)],col="indianred3")
lines(PredInt[which(Data$MM==MonthCalmNum[i])][-(1:4)],col="indianred3",pch=20,cex=0.2)
for(i in 2:length(MonthCalm))
{
  plot(Data$SWHft[which(Data$MM==MonthCalmNum[i])],type="l",ylim=c(0,10),
       main =paste(MonthCalm[i],j,sep="-"),xaxt='n')
  lines(Data$SWHft[which(Data$MM==MonthCalmNum[i])],col=Data$colSWHft[which(Data$MM==MonthCalmNum[i])],pch=19,cex=0.2,type="p")
  lines(PredInt[which(Data$MM==MonthCalmNum[i])] ,col="indianred3")
  lines(PredInt[which(Data$MM==MonthCalmNum[i])],col="indianred3",pch=20,cex=0.2)
}


