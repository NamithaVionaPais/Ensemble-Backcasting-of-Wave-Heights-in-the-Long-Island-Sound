setwd('/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs/Threshold_Model')
load("44039_reg_merged.Rda")
source("Data_Function.R")
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
  Data<-data.in("central",j,"Windy")
  #Model 1
  runchngpt<-readRDS(paste('TM_',j,'_Windy_slurm','.RDS',sep=""))
  ResWindy<-residuals(runchngpt)
  ModelRes_Windy<-garchFit(~garch(1,2), data = ResWindy,trace=FALSE,include.mean = FALSE)
  Ind_Res_Windy<-ifelse(ResWindy<0,-1,1)
  sigW<-ModelRes_Windy@sigma.t
  RhatWindy<-(sigW*Ind_Res_Windy)
  #round(ModelRes_Windy@fit$matcoef[,1],4)
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
colnames(2008)
for(j in 2008)
{
  Data<-data.in("central",j,"Calm")
  runchngpt<-readRDS( paste('TM_',j,'_Calm_slurm','.RDS',sep=""))
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


