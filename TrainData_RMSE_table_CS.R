setwd('/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs/Threshold_Model')
load("44039_reg_merged.Rda")
source("Data_Function.R")
library(chngpt)
library(fGarch)
library(Metrics)
ED_Year_W=ED_Year_C=c()
for(j in 2005:2013)
{
  setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
  load("44039_reg_merged.Rda")
  source("Data_Function_CS.R")
  DataW<-data.in_CS(j,"Windy")
  #Model 1
  DataW$lag1rWSPD2_T<-c(rep(NA,(0+1)),DataW$rWSPD2_T[(1:(length(DataW$rWSPD2_T)-1))])
  DataW$lag2rWSPD2_T<-c(rep(NA,(0+2)),DataW$rWSPD2_T[(1:(length(DataW$rWSPD2_T)-2))])
  DataW$lag3rWSPD2_T<-c(rep(NA,(0+3)),DataW$rWSPD2_T[(1:(length(DataW$rWSPD2_T)-3))])
  DataW$lag4rWSPD2_T<-c(rep(NA,(0+4)),DataW$rWSPD2_T[(1:(length(DataW$rWSPD2_T)-4))])
  DataW$rWSPDind_T<-as.factor(ifelse(DataW$rWSPD<= 0.25,1,0))
  formula_1<-SWHft ~ lag1SWHft+lag2SWHft+lag3SWHft+lag4SWHft+
    lag1rWSPD2_T + lag2rWSPD2_T + lag3rWSPD2_T +lag4rWSPD2_T+
    z2_T+ z3_T +z4_T+
    I(lag1SWHft*rWSPD2_T)+
    I(lag3SWHft*lag2rWSPD2_T)+
    rWSPDind_T
  formula_2<-~ rWSPD2_T
  runchngptW<-chngptm(formula.1= formula_1, formula.2=formula_2 , family="gaussian", data=DataW, 
                     type="segmented", var.type="bootstrap", weights=NULL,ci.bootstrap.size = 100)
  setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
  load("44039_reg_merged.Rda")
  source("Data_Function_CS.R")
  DataC<-data.in_CS(j,"Calm")
  #Model 1
  DataC$lag1rWSPD2_T<-c(rep(NA,(0+1)),DataC$rWSPD2_T[(1:(length(DataC$rWSPD2_T)-1))])
  DataC$lag2rWSPD2_T<-c(rep(NA,(0+2)),DataC$rWSPD2_T[(1:(length(DataC$rWSPD2_T)-2))])
  DataC$lag3rWSPD2_T<-c(rep(NA,(0+3)),DataC$rWSPD2_T[(1:(length(DataC$rWSPD2_T)-3))])
  DataC$lag4rWSPD2_T<-c(rep(NA,(0+4)),DataC$rWSPD2_T[(1:(length(DataC$rWSPD2_T)-4))])
  DataC$rWSPDind_T<-as.factor(ifelse(DataC$rWSPD<= 0.25,1,0))
  formula_1<-SWHft ~ lag1SWHft+lag2SWHft+lag3SWHft+lag4SWHft+
    lag1rWSPD2_T + lag2rWSPD2_T + lag3rWSPD2_T +lag4rWSPD2_T+
    z2_T+ z3_T +z4_T+
    I(lag1SWHft*rWSPD2_T)+
    I(lag3SWHft*lag2rWSPD2_T)+
    rWSPDind_T
  formula_2<-~ rWSPD2_T
  runchngptC<-chngptm(formula.1= formula_1, formula.2=formula_2 , family="gaussian", data=DataC, 
                     type="segmented", var.type="bootstrap", weights=NULL,ci.bootstrap.size = 100)

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

