setwd('/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs/Threshold_Model')
load("44039_reg_merged.Rda")
source("Data_Function.R")
library(chngpt)
library(fGarch)
library(Metrics)
ED_Year_W=ED_Year_C=c()
for(j in 2005:2013)
{
  DataW<-data.in("central",j,"Windy")
  DataC<-data.in("central",j,"Calm")
  #Model 1
  runchngptW<-readRDS(paste('TM_',j,'_Windy_slurm','.RDS',sep=""))
  runchngptC<-readRDS(paste('TM_',j,'_Calm_slurm','.RDS',sep=""))
  ResWindy<-residuals(runchngptW)
  ResCalm<-residuals(runchngptC)
  ModelRes_Windy<-garchFit(~garch(1,2), data = ResWindy,trace=FALSE,include.mean = FALSE)
  ModelRes_Calm<-garchFit(~garch(1,2), data = ResCalm,trace=FALSE,include.mean = FALSE)
  Ind_Res_Windy<-ifelse(ResWindy<0,-1,1)
  Ind_Res_Calm<-ifelse(ResCalm<0,-1,1)
  sigW<-ModelRes_Windy@sigma.t
  sigC<-ModelRes_Calm@sigma.t
  RhatWindy<-(sigW*Ind_Res_Windy)
  RhatCalm<-(sigC*Ind_Res_Calm)
  PredIntW<-c()
  PredIntW[1:5]<-rep(0,5)
  PredIntW[6:nrow(DataW)]<-predict(runchngptW)+RhatWindy
  ED_Year_W[(j-2004)]<-round(rmse(PredIntW[6:nrow(DataW)],DataW$SWHft[-(1:5)]),4)
  PredIntC<-c()
  PredIntC[1:4]<-rep(0,4)
  PredIntC[5:nrow(DataC)]<-predict(runchngptC)+RhatCalm
  ED_Year_C[(j-2004)]<-round(rmse(PredIntC[5:nrow(DataC)],DataC$SWHft[-(1:4)]),4)
}
cbind(2005:2013,ED_Year_W,ED_Year_C)

