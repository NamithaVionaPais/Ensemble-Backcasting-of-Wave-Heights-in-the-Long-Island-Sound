#Train test
#Ensemble_TM-GARCH
setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
load("44039_reg_merged.Rda")
source("Data_Function_CS_BC.R")
library(fGarch)
library(chngpt)

#### WINDY Season #####

Kseq=c(1974:1975,1977:1979,1981:1983,1985:1987,1989:1991,1993:1995,1997:1999,2001:2003)
#Kseq=c(1997)
PredSeqWindy<-list()
PredSeqWindy_p10sd<-list()
PredSeqWindy_m10sd<-list()
for(k in  Kseq)
{
  setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
  load("44039_reg_merged.Rda")
  source("Data_Function_CS_BC.R")
  Data<-data.in_CS_bc(k,"Windy")
  Pred_2005_TM_Windy<-matrix(nrow=nrow(Data),ncol=9)
  Pred_2005_TM_Windy_garch<-matrix(nrow=nrow(Data),ncol=9)
  Pred_2005_TM_Windy_garch_p10sd<-matrix(nrow=nrow(Data),ncol=9)
  Pred_2005_TM_Windy_garch_m10sd<-matrix(nrow=nrow(Data),ncol=9)
  
  
  for(j in c(2005,2006,2007,2009,2010,2011,2013))
  {
    
    setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
    load("44039_reg_merged.Rda")
    source("Data_Function_CS.R")
    Data_train<-data.in_CS(j,"Windy")
    Data_train$lag1rWSPD2_T<-c(rep(NA,(0+1)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-1))])
    Data_train$lag2rWSPD2_T<-c(rep(NA,(0+2)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-2))])
    Data_train$lag3rWSPD2_T<-c(rep(NA,(0+3)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-3))])
    Data_train$lag4rWSPD2_T<-c(rep(NA,(0+4)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-4))])
    Data_train$rWSPDind_T<-as.factor(ifelse(Data_train$rWSPD<= 0.25,1,0))
    formula_1<-SWHft ~ lag1SWHft+lag2SWHft+lag3SWHft+lag4SWHft+
      lag1rWSPD2_T + lag2rWSPD2_T + lag3rWSPD2_T +lag4rWSPD2_T+
      z2_T+ z3_T +z4_T+
      I(lag1SWHft*rWSPD2_T)+
      I(lag3SWHft*lag2rWSPD2_T)+
      rWSPDind_T
    formula_2<-~ rWSPD2_T
    runLag2_Int<-chngptm(formula.1= formula_1, formula.2=formula_2 , family="gaussian", data=Data_train, 
                         type="segmented", var.type="bootstrap", weights=NULL,ci.bootstrap.size = 100)
    
    ResWindy<-residuals(runLag2_Int)
    ModelRes_Windy<-garchFit(~garch(1,2), data = ResWindy,trace=FALSE,include.mean =FALSE)
    Ind_Res_Windy<-ifelse(ResWindy<0,-1,1)
    sigW<-ModelRes_Windy@sigma.t
    RhatWindy<-(sigW*Ind_Res_Windy)
    setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
    load("44039_reg_merged.Rda")
    source("Data_Function_CS_BC.R")
    Data<-data.in_CS_bc(k,"Windy")
    Data$lag1rWSPD2_T<-c(rep(NA,(0+1)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-1))])
    Data$lag2rWSPD2_T<-c(rep(NA,(0+2)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-2))])
    Data$lag3rWSPD2_T<-c(rep(NA,(0+3)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-3))])
    Data$lag4rWSPD2_T<-c(rep(NA,(0+4)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-4))])
    Data$rWSPDind_T<-ifelse(Data$rWSPD_T<= 0.25,1,0)
    vec<-ifelse((Data$rWSPD2_T-runLag2_Int$coefficients[18])>0,(Data$rWSPD2_T-runLag2_Int$coefficients[18]),0)
    Data1<-cbind(rep(1,nrow(Data))
                 ,Data$lag1rWSPD2_T             
                 ,Data$lag2rWSPD2_T      
                 ,Data$lag3rWSPD2_T           
                 ,Data$lag4rWSPD2_T
                 ,Data$z2_T
                 ,Data$z3_T
                 ,Data$z4_T
                 ,(as.numeric(Data$rWSPDind)-1)
                 ,Data$rWSPD2_T
                 ,vec)
    Coef1<-runLag2_Int$coefficients
    Coef11<-Coef1[c(1,6,7,8,9,10,11,12,15,16,17)]
    Coef12<-Coef1[c(2,3,4,5,13,14)]
    X_mat<-matrix(,nrow=nrow(Data),ncol=17)
    X_mat[,c(1,6,7,8,9,10,11,12,15,16,17)]<-Data1
    Pred<-c()
    Pred[1:5]<-0
    Pred_2005_TM_Windy[1:5]<-0
    Pred_2005_TM_Windy_garch[1:5]<-0
    RhatWindy1<-RhatWindy
    for(t in 6:nrow(Data1))
    {
      Pred[t]<-Coef11%*%Data1[t,]+(Coef12[1]*Pred[(t-1)])+(Coef12[2]*Pred[(t-2)])+(Coef12[3]*Pred[(t-3)])+
        (Coef12[4]*Pred[(t-4)])+(Coef12[5]*Pred[(t-1)]*Data$rWSPD2_T[t])+(Coef12[6]*Pred[(t-3)]*Data$lag2rWSPD2_T[t])
      X_mat[t,c(2,3,4,5,13,14)]<-c(Pred[(t-1)],Pred[(t-2)],Pred[(t-3)],Pred[(t-4)],(Pred[(t-1)]*Data$rWSPD2_T[t]),
                                   (Pred[(t-3)]*Data$lag2rWSPD2_T[t]))
    }
    Pred_2005_TM_Windy[6:nrow(Data),(j-2004)]<-Pred[-(1:5)]
    Pred_2005_TM_Windy_garch[6:nrow(Data),(j-2004)]<-Pred[-(1:5)]+RhatWindy1
    V_yhat<-as.vector(diag(X_mat[-(1:5),]%*%vcov(runLag2_Int)[-18,-18]%*%t(X_mat[-(1:5),])))
    Pred_2005_TM_Windy_garch_m10sd[6:nrow(Data),(j-2004)]<-Pred[-(1:5)]+RhatWindy1-(10*sqrt(V_yhat))
    Pred_2005_TM_Windy_garch_p10sd[6:nrow(Data),(j-2004)]<-Pred[-(1:5)]+RhatWindy1+(10*sqrt(V_yhat))
    
    
  }
  
  for(j in c(2008,2012))
  {
    
    setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
    load("44039_reg_merged.Rda")
    source("Data_Function_CS.R")
    Data_train<-data.in_CS(j,"Windy")
    Data_train$lag1rWSPD2_T<-c(rep(NA,(0+1)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-1))])
    Data_train$lag2rWSPD2_T<-c(rep(NA,(0+2)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-2))])
    Data_train$lag3rWSPD2_T<-c(rep(NA,(0+3)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-3))])
    Data_train$lag4rWSPD2_T<-c(rep(NA,(0+4)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-4))])
    Data_train$rWSPDind_T<-as.factor(ifelse(Data_train$rWSPD<= 0.25,1,0))
    formula_1<-SWHft ~ lag1SWHft+lag2SWHft+lag3SWHft+lag4SWHft+
      lag1rWSPD2_T + lag2rWSPD2_T + lag3rWSPD2_T +lag4rWSPD2_T+
      z2_T+ z3_T +z4_T+
      I(lag1SWHft*rWSPD2_T)+
      I(lag3SWHft*lag2rWSPD2_T)+
      rWSPDind_T
    formula_2<-~ rWSPD2_T
    runLag2_Int<-chngptm(formula.1= formula_1, formula.2=formula_2 , family="gaussian", data=Data_train, 
                         type="segmented", var.type="bootstrap", weights=NULL,ci.bootstrap.size = 100)
    
    ResWindy<-residuals(runLag2_Int)
    ModelRes_Windy<-garchFit(~garch(1,2), data = ResWindy,trace=FALSE,include.mean =FALSE)
    Ind_Res_Windy<-ifelse(ResWindy<0,-1,1)
    sigW<-ModelRes_Windy@sigma.t
    RhatWindy<-(sigW*Ind_Res_Windy)
    setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
    load("44039_reg_merged.Rda")
    source("Data_Function_CS_BC.R")
    Data<-data.in_CS_bc(k,"Windy")
    Data$lag1rWSPD2_T<-c(rep(NA,(0+1)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-1))])
    Data$lag2rWSPD2_T<-c(rep(NA,(0+2)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-2))])
    Data$lag3rWSPD2_T<-c(rep(NA,(0+3)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-3))])
    Data$lag4rWSPD2_T<-c(rep(NA,(0+4)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-4))])
    Data$rWSPDind_T<-ifelse(Data$rWSPD<= 0.25,1,0)
    vec<-ifelse((Data$rWSPD2_T-runLag2_Int$coefficients[18])>0,(Data$rWSPD2_T-runLag2_Int$coefficients[18]),0)
    Data1<-cbind(rep(1,nrow(Data))
                 ,Data$lag1rWSPD2_T             
                 ,Data$lag2rWSPD2_T      
                 ,Data$lag3rWSPD2_T           
                 ,Data$lag4rWSPD2_T
                 ,Data$z2_T
                 ,Data$z3_T
                 ,Data$z4_T
                 ,(as.numeric(Data$rWSPDind)-1)
                 ,Data$rWSPD2_T
                 ,vec)
    Coef1<-runLag2_Int$coefficients
    Coef11<-Coef1[c(1,6,7,8,9,10,11,12,15,16,17)]
    Coef12<-Coef1[c(2,3,4,5,13,14)]
    X_mat<-matrix(,nrow=nrow(Data),ncol=17)
    X_mat[,c(1,6,7,8,9,10,11,12,15,16,17)]<-Data1
    
    Pred<-c()
    Pred[1:5]<-0
    Pred_2005_TM_Windy[1:5]<-0
    Pred_2005_TM_Windy_garch[1:5]<-0
    RhatWindy1<-c()
    RhatWindy1[1:1411]<-RhatWindy[1:1411]
    RhatWindy1[1412:3619]<-RhatWindy[1436:3643]
    for(t in 6:nrow(Data))
    {
      Pred[t]<-Coef11%*%Data1[t,]+(Coef12[1]*Pred[(t-1)])+(Coef12[2]*Pred[(t-2)])+(Coef12[3]*Pred[(t-3)])+
        (Coef12[4]*Pred[(t-4)])+(Coef12[5]*Pred[(t-1)]*Data$rWSPD2[t])+(Coef12[6]*Pred[(t-3)]*Data$lag2rWSPD2[t])
      X_mat[t,c(2,3,4,5,13,14)]<-c(Pred[(t-1)],Pred[(t-2)],Pred[(t-3)],Pred[(t-4)],(Pred[(t-1)]*Data$rWSPD2_T[t]),
                                   (Pred[(t-3)]*Data$lag2rWSPD2_T[t]))
    }
    Pred_2005_TM_Windy[6:nrow(Data),(j-2004)]<-Pred[-(1:5)]
    Pred_2005_TM_Windy_garch[6:nrow(Data),(j-2004)]<-Pred[-(1:5)]+RhatWindy1
    V_yhat<-as.vector(diag(X_mat[-(1:5),]%*%vcov(runLag2_Int)[-18,-18]%*%t(X_mat[-(1:5),])))
    Pred_2005_TM_Windy_garch_m10sd[6:nrow(Data),(j-2004)]<-Pred[-(1:5)]+RhatWindy1-(10*sqrt(V_yhat))
    Pred_2005_TM_Windy_garch_p10sd[6:nrow(Data),(j-2004)]<-Pred[-(1:5)]+RhatWindy1+(10*sqrt(V_yhat))
    
  }
  Pred_2005_W<-c()
  Pred_2005_W_p10sd<-c()
  Pred_2005_W_m10sd<-c()
  Pred_2005_W<-rowMeans(Pred_2005_TM_Windy_garch)
  Pred_2005_W_p10sd<-rowMeans(Pred_2005_TM_Windy_garch_p10sd)
  Pred_2005_W_m10sd<-rowMeans(Pred_2005_TM_Windy_garch_m10sd)
  PredSeqWindy[[(k-1973)]]=ifelse(Pred_2005_W<0,0,Pred_2005_W)
  PredSeqWindy_p10sd[[(k-1973)]]<-ifelse(Pred_2005_W_p10sd<0,0,Pred_2005_W_p10sd)
  PredSeqWindy_m10sd[[(k-1973)]]<-ifelse(Pred_2005_W_m10sd<0,0,Pred_2005_W_m10sd)
}

Kseq1=c(1976,1980,1984,1988,1992,1996,2000,2004)
#Kseq1<-c(1976)
for(k in  Kseq1)
{
  setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
  load("44039_reg_merged.Rda")
  source("Data_Function_CS_BC.R")
  Data<-data.in_CS_bc(k,"Windy")
  Pred_2005_TM_Windy<-matrix(nrow=nrow(Data),ncol=9)
  Pred_2005_TM_Windy_garch<-matrix(nrow=nrow(Data),ncol=9)
  Pred_2005_TM_Windy_garch_p10sd<-matrix(nrow=nrow(Data),ncol=9)
  Pred_2005_TM_Windy_garch_m10sd<-matrix(nrow=nrow(Data),ncol=9)
  for(j in c(2005,2006,2007,2009,2010,2011,2013))
  {
    setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
    load("44039_reg_merged.Rda")
    source("Data_Function_CS.R")
    Data_train<-data.in_CS(j,"Windy")
    Data_train$lag1rWSPD2_T<-c(rep(NA,(0+1)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-1))])
    Data_train$lag2rWSPD2_T<-c(rep(NA,(0+2)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-2))])
    Data_train$lag3rWSPD2_T<-c(rep(NA,(0+3)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-3))])
    Data_train$lag4rWSPD2_T<-c(rep(NA,(0+4)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-4))])
    Data_train$rWSPDind_T<-as.factor(ifelse(Data_train$rWSPD<= 0.25,1,0))
    formula_1<-SWHft ~ lag1SWHft+lag2SWHft+lag3SWHft+lag4SWHft+
      lag1rWSPD2_T + lag2rWSPD2_T + lag3rWSPD2_T +lag4rWSPD2_T+
      z2_T+ z3_T +z4_T+
      I(lag1SWHft*rWSPD2_T)+
      I(lag3SWHft*lag2rWSPD2_T)+
      rWSPDind_T
    formula_2<-~ rWSPD2_T
    runLag2_Int<-chngptm(formula.1= formula_1, formula.2=formula_2 , family="gaussian", data=Data_train, 
                         type="segmented", var.type="bootstrap", weights=NULL,ci.bootstrap.size = 100)
    ResWindy<-residuals(runLag2_Int)
    ModelRes_Windy<-garchFit(~garch(1,2), data = ResWindy,trace=FALSE,include.mean = FALSE)
    Ind_Res_Windy<-ifelse(ResWindy<0,-1,1)
    sigW<-ModelRes_Windy@sigma.t
    RhatWindy<-(sigW*Ind_Res_Windy)
    setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
    load("44039_reg_merged.Rda")
    source("Data_Function_CS_BC.R")
    Data<-data.in_CS_bc(k,"Windy")
    Data$lag1rWSPD2_T<-c(rep(NA,(0+1)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-1))])
    Data$lag2rWSPD2_T<-c(rep(NA,(0+2)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-2))])
    Data$lag3rWSPD2_T<-c(rep(NA,(0+3)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-3))])
    Data$lag4rWSPD2_T<-c(rep(NA,(0+4)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-4))])
    Data$rWSPDind_T<-ifelse(Data$rWSPD<= 0.25,1,0)
    vec<-ifelse((Data$rWSPD2_T-runLag2_Int$coefficients[18])>0,(Data$rWSPD2_T-runLag2_Int$coefficients[18]),0)
    Data1<-cbind(rep(1,nrow(Data))
                 ,Data$lag1rWSPD2_T             
                 ,Data$lag2rWSPD2_T      
                 ,Data$lag3rWSPD2_T           
                 ,Data$lag4rWSPD2_T
                 ,Data$z2_T
                 ,Data$z3_T
                 ,Data$z4_T
                 ,(as.numeric(Data$rWSPDind)-1)
                 ,Data$rWSPD2_T
                 ,vec)
    Coef1<-runLag2_Int$coefficients
    Coef11<-Coef1[c(1,6,7,8,9,10,11,12,15,16,17)]
    Coef12<-Coef1[c(2,3,4,5,13,14)]
    X_mat<-matrix(,nrow=nrow(Data),ncol=17)
    X_mat[,c(1,6,7,8,9,10,11,12,15,16,17)]<-Data1
    Pred<-c()
    Pred[1:5]<-0
    Pred_2005_TM_Windy[1:5]<-0
    Pred_2005_TM_Windy_garch[1:5]<-0
    #RhatWindy1<-RhatWindy
    RhatWindy1<-c()
    RhatWindy1[1:1411]<-RhatWindy[1:1411]
    RhatWindy1[1412:1435]<-rep(0,24)
    RhatWindy1[1436:3643]<-RhatWindy[1412:3619]
    for(t in 6:nrow(Data1))
    {
      Pred[t]<-Coef11%*%Data1[t,]+(Coef12[1]*Pred[(t-1)])+(Coef12[2]*Pred[(t-2)])+(Coef12[3]*Pred[(t-3)])+
        (Coef12[4]*Pred[(t-4)])+(Coef12[5]*Pred[(t-1)]*Data$rWSPD2_T[t])+(Coef12[6]*Pred[(t-3)]*Data$lag2rWSPD2_T[t])
      X_mat[t,c(2,3,4,5,13,14)]<-c(Pred[(t-1)],Pred[(t-2)],Pred[(t-3)],Pred[(t-4)],(Pred[(t-1)]*Data$rWSPD2_T[t]),
                                   (Pred[(t-3)]*Data$lag2rWSPD2_T[t]))
    }
    Pred_2005_TM_Windy[6:nrow(Data),(j-2004)]<-Pred[-(1:5)]
    Pred_2005_TM_Windy_garch[6:nrow(Data),(j-2004)]<-Pred[-(1:5)]+RhatWindy1
    V_yhat<-as.vector(diag(X_mat[-(1:5),]%*%vcov(runLag2_Int)[-18,-18]%*%t(X_mat[-(1:5),])))
    Pred_2005_TM_Windy_garch_m10sd[6:nrow(Data),(j-2004)]<-Pred[-(1:5)]+RhatWindy1-(10*sqrt(V_yhat))
    Pred_2005_TM_Windy_garch_p10sd[6:nrow(Data),(j-2004)]<-Pred[-(1:5)]+RhatWindy1+(10*sqrt(V_yhat))
  }
  
  for(j in c(2008,2012))
  {
    setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
    load("44039_reg_merged.Rda")
    source("Data_Function_CS.R")
    Data_train<-data.in_CS(j,"Windy")
    Data_train$lag1rWSPD2_T<-c(rep(NA,(0+1)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-1))])
    Data_train$lag2rWSPD2_T<-c(rep(NA,(0+2)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-2))])
    Data_train$lag3rWSPD2_T<-c(rep(NA,(0+3)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-3))])
    Data_train$lag4rWSPD2_T<-c(rep(NA,(0+4)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-4))])
    Data_train$rWSPDind_T<-as.factor(ifelse(Data_train$rWSPD<= 0.25,1,0))
    formula_1<-SWHft ~ lag1SWHft+lag2SWHft+lag3SWHft+lag4SWHft+
      lag1rWSPD2_T + lag2rWSPD2_T + lag3rWSPD2_T +lag4rWSPD2_T+
      z2_T+ z3_T +z4_T+
      I(lag1SWHft*rWSPD2_T)+
      I(lag3SWHft*lag2rWSPD2_T)+
      rWSPDind_T
    formula_2<-~ rWSPD2_T
    runLag2_Int<-chngptm(formula.1= formula_1, formula.2=formula_2 , family="gaussian", data=Data_train, 
                         type="segmented", var.type="bootstrap", weights=NULL,ci.bootstrap.size = 100)
    ResWindy<-residuals(runLag2_Int)
    ModelRes_Windy<-garchFit(~garch(1,2), data = ResWindy,trace=FALSE,include.mean = FALSE)
    Ind_Res_Windy<-ifelse(ResWindy<0,-1,1)
    sigW<-ModelRes_Windy@sigma.t
    RhatWindy<-(sigW*Ind_Res_Windy)
    setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
    load("44039_reg_merged.Rda")
    source("Data_Function_CS_BC.R")
    Data<-data.in_CS_bc(k,"Windy")
    Data$lag1rWSPD2_T<-c(rep(NA,(0+1)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-1))])
    Data$lag2rWSPD2_T<-c(rep(NA,(0+2)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-2))])
    Data$lag3rWSPD2_T<-c(rep(NA,(0+3)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-3))])
    Data$lag4rWSPD2_T<-c(rep(NA,(0+4)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-4))])
    Data$rWSPDind_T<-ifelse(Data$rWSPD<= 0.25,1,0)
    vec<-ifelse((Data$rWSPD2_T-runLag2_Int$coefficients[18])>0,(Data$rWSPD2_T-runLag2_Int$coefficients[18]),0)
    Data1<-cbind(rep(1,nrow(Data))
                 ,Data$lag1rWSPD2_T             
                 ,Data$lag2rWSPD2_T      
                 ,Data$lag3rWSPD2_T           
                 ,Data$lag4rWSPD2_T
                 ,Data$z2_T
                 ,Data$z3_T
                 ,Data$z4_T
                 ,(as.numeric(Data$rWSPDind)-1)
                 ,Data$rWSPD2_T
                 ,vec)
    Coef1<-runLag2_Int$coefficients
    Coef11<-Coef1[c(1,6,7,8,9,10,11,12,15,16,17)]
    Coef12<-Coef1[c(2,3,4,5,13,14)]
    X_mat<-matrix(,nrow=nrow(Data),ncol=17)
    X_mat[,c(1,6,7,8,9,10,11,12,15,16,17)]<-Data1
    Pred<-c()
    Pred[1:5]<-0
    Pred_2005_TM_Windy[1:5]<-0
    Pred_2005_TM_Windy_garch[1:5]<-0
    RhatWindy1<-c()
    #RhatWindy1[1:1411]<-RhatWindy[1:1411]
    #RhatWindy1[1412:3619]<-RhatWindy[1436:3643]
    RhatWindy1<-RhatWindy
    for(t in 6:nrow(Data))
    {
      Pred[t]<-Coef11%*%Data1[t,]+(Coef12[1]*Pred[(t-1)])+(Coef12[2]*Pred[(t-2)])+(Coef12[3]*Pred[(t-3)])+
        (Coef12[4]*Pred[(t-4)])+(Coef12[5]*Pred[(t-1)]*Data$rWSPD2[t])+(Coef12[6]*Pred[(t-3)]*Data$lag2rWSPD2[t])
      X_mat[t,c(2,3,4,5,13,14)]<-c(Pred[(t-1)],Pred[(t-2)],Pred[(t-3)],Pred[(t-4)],(Pred[(t-1)]*Data$rWSPD2_T[t]),
                                   (Pred[(t-3)]*Data$lag2rWSPD2_T[t]))
    }
    Pred_2005_TM_Windy[6:nrow(Data),(j-2004)]<-Pred[-(1:5)]
    Pred_2005_TM_Windy_garch[6:nrow(Data),(j-2004)]<-Pred[-(1:5)]+RhatWindy1
    V_yhat<-as.vector(diag(X_mat[-(1:5),]%*%vcov(runLag2_Int)[-18,-18]%*%t(X_mat[-(1:5),])))
    Pred_2005_TM_Windy_garch_m10sd[6:nrow(Data),(j-2004)]<-Pred[-(1:5)]+RhatWindy1-(10*sqrt(V_yhat))
    Pred_2005_TM_Windy_garch_p10sd[6:nrow(Data),(j-2004)]<-Pred[-(1:5)]+RhatWindy1+(10*sqrt(V_yhat))
  }
  Pred_2005_W<-c()
  Pred_2005_W_p10sd<-c()
  Pred_2005_W_m10sd<-c()
  Pred_2005_W<-rowMeans(Pred_2005_TM_Windy_garch)
  Pred_2005_W_p10sd<-rowMeans(Pred_2005_TM_Windy_garch_p10sd)
  Pred_2005_W_m10sd<-rowMeans(Pred_2005_TM_Windy_garch_m10sd)
  PredSeqWindy[[(k-1973)]]=ifelse(Pred_2005_W<0,0,Pred_2005_W)
  PredSeqWindy_p10sd[[(k-1973)]]<-ifelse(Pred_2005_W_p10sd<0,0,Pred_2005_W_p10sd)
  PredSeqWindy_m10sd[[(k-1973)]]<-ifelse(Pred_2005_W_m10sd<0,0,Pred_2005_W_m10sd)
}

names(PredSeqWindy)=1974:2004
names(PredSeqWindy_p10sd)=1974:2004
names(PredSeqWindy_m10sd)=1974:2004

Nov_mean=Dec_mean=Jan_mean=Feb_mean=March_mean=c()
for(i in Kseq )
{
  Nov_mean[(i-1973)]<-mean(na.omit(PredSeqWindy[[(i-1973)]][1:720]))
  Dec_mean[(i-1973)]<-mean(na.omit(PredSeqWindy[[(i-1973)]][721:1464]))
  Jan_mean[(i-1973)]<-mean(na.omit(PredSeqWindy[[(i-1973)]][1465:2208]))
  Feb_mean[(i-1973)]<-mean(na.omit(PredSeqWindy[[(i-1973)]][2209:2880]))
  March_mean[(i-1973)]<-mean(na.omit(PredSeqWindy[[(i-1973)]][2881:3624]))
}

for(i in Kseq1)
{
  Nov_mean[(i-1973)]<-mean(na.omit(PredSeqWindy[[(i-1973)]][1:720]))
  Dec_mean[(i-1973)]<-mean(na.omit(PredSeqWindy[[(i-1973)]][721:1464]))
  Jan_mean[(i-1973)]<-mean(na.omit(PredSeqWindy[[(i-1973)]][1465:2208]))
  Feb_mean[(i-1973)]<-mean(na.omit(PredSeqWindy[[(i-1973)]][2209:2904]))
  March_mean[(i-1973)]<-mean(na.omit(PredSeqWindy[[(i-1973)]][2905:3648]))
}


Data_Buoy <- read.csv("~/Documents/Documents - Namitha’s MacBook Pro/MarineScience/NVP-spring22/Sreeram_NR/Data_Buoy.csv")
head(Data_Buoy)
kseq=c(2005,2006,2007,2009,2010,2011,2013)
kseq1=c(2008,2012)
for(i in kseq )
{
  D1<-Data_Buoy[which(Data_Buoy$year==i),]
  D11<-Data_Buoy[which(Data_Buoy$year==(i-1)),]
  D2<-c(D11$SWHft[which(D11$month==11|D11$month==12)],D1$SWHft[which(D1$month==1|D1$month==2|D1$month==3)])
  Nov_mean[(i-1973)]<-mean(na.omit(D2[1:720]))
  Dec_mean[(i-1973)]<-mean(na.omit(D2[721:1464]))
  Jan_mean[(i-1973)]<-mean(na.omit(D2[1465:2208]))
  Feb_mean[(i-1973)]<-mean(na.omit(D2[2209:2880]))
  March_mean[(i-1973)]<-mean(na.omit(D2[2881:3624]))
}

for(i in kseq1)
{
  D1<-Data_Buoy[which(Data_Buoy$year==i),]
  D11<-Data_Buoy[which(Data_Buoy$year==(i-1)),]
  D2<-c(D11$SWHft[which(D11$month==11|D11$month==12)],D1$SWHft[which(D1$month==1|D1$month==2|D1$month==3)])
  Nov_mean[(i-1973)]<-mean(na.omit(D2[1:720]))
  Dec_mean[(i-1973)]<-mean(na.omit(D2[721:1464]))
  Jan_mean[(i-1973)]<-mean(na.omit(D2[1465:2208]))
  Feb_mean[(i-1973)]<-mean(na.omit(D2[2209:2904]))
  March_mean[(i-1973)]<-mean(na.omit(D2[2905:3648]))
}

boxplot(Nov_mean,Dec_mean,Jan_mean,Feb_mean,March_mean,names=c("Nov","Dec","Jan","Feb","March"))



Nov_max=Dec_max=Jan_max=Feb_max=March_max=c()
for(i in Kseq )
{
  Nov_max[(i-1973)]<-max(na.omit(PredSeqWindy[[(i-1973)]][1:720]))
  Dec_max[(i-1973)]<-max(na.omit(PredSeqWindy[[(i-1973)]][721:1464]))
  Jan_max[(i-1973)]<-max(na.omit(PredSeqWindy[[(i-1973)]][1465:2208]))
  Feb_max[(i-1973)]<-max(na.omit(PredSeqWindy[[(i-1973)]][2209:2880]))
  March_max[(i-1973)]<-max(na.omit(PredSeqWindy[[(i-1973)]][2881:3624]))
}

for(i in Kseq1)
{
  Nov_max[(i-1973)]<-max(na.omit(PredSeqWindy[[(i-1973)]][1:720]))
  Dec_max[(i-1973)]<-max(na.omit(PredSeqWindy[[(i-1973)]][721:1464]))
  Jan_max[(i-1973)]<-max(na.omit(PredSeqWindy[[(i-1973)]][1465:2208]))
  Feb_max[(i-1973)]<-max(na.omit(PredSeqWindy[[(i-1973)]][2209:2904]))
  March_max[(i-1973)]<-max(na.omit(PredSeqWindy[[(i-1973)]][2905:3648]))
}
kseq=c(2005,2006,2007,2009,2010,2011,2013)
kseq1=c(2008,2012)
for(i in kseq )
{
  D1<-Data_Buoy[which(Data_Buoy$year==i),]
  D11<-Data_Buoy[which(Data_Buoy$year==(i-1)),]
  D2<-c(D11$SWHft[which(D11$month==11|D11$month==12)],D1$SWHft[which(D1$month==1|D1$month==2|D1$month==3)])
  Nov_max[(i-1973)]<-max(na.omit(D2[1:720]))
  Dec_max[(i-1973)]<-max(na.omit(D2[721:1464]))
  Jan_max[(i-1973)]<-max(na.omit(D2[1465:2208]))
  Feb_max[(i-1973)]<-max(na.omit(D2[2209:2880]))
  March_max[(i-1973)]<-max(na.omit(D2[2881:3624]))
}

for(i in kseq1)
{
  D1<-Data_Buoy[which(Data_Buoy$year==i),]
  D11<-Data_Buoy[which(Data_Buoy$year==(i-1)),]
  D2<-c(D11$SWHft[which(D11$month==11|D11$month==12)],D1$SWHft[which(D1$month==1|D1$month==2|D1$month==3)])
  Nov_max[(i-1973)]<-max(na.omit(D2[1:720]))
  Dec_max[(i-1973)]<-max(na.omit(D2[721:1464]))
  Jan_max[(i-1973)]<-max(na.omit(D2[1465:2208]))
  Feb_max[(i-1973)]<-max(na.omit(D2[2209:2904]))
  March_max[(i-1973)]<-max(na.omit(D2[2905:3648]))
}

boxplot(Nov_max,Dec_max,Jan_max,Feb_max,March_max,names=c("Nov","Dec","Jan","Feb","March"))


#### CALM Season #####

PredSeqCalm<-list()
PredSeqCalm_p10sd<-list()
PredSeqCalm_m10sd<-list()
Kseq=1974:2004
#Kseq=c(1997)
for(k in Kseq)
{
  setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
  load("Sikorsky-Bridgeport_merged.Rda")
  source("Data_Function_CS_BC.R")
  Data<-data.in_CS_bc(k,"Calm")
  Pred_2005_TM_Calm<-matrix(nrow=nrow(Data),ncol=9)
  Pred_2005_TM_Calm_garch<-matrix(nrow=nrow(Data),ncol=9)
  Pred_2005_TM_Calm_garch_p10sd<-matrix(nrow=nrow(Data),ncol=9)
  Pred_2005_TM_Calm_garch_m10sd<-matrix(nrow=nrow(Data),ncol=9)
  for(j in 2005:2013)
  {
    setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
    load("44039_reg_merged.Rda")
    source("Data_Function_CS.R")
    Data_train<-data.in_CS(j,"Calm")
    Data_train$lag1rWSPD2_T<-c(rep(NA,(0+1)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-1))])
    Data_train$lag2rWSPD2_T<-c(rep(NA,(0+2)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-2))])
    Data_train$lag3rWSPD2_T<-c(rep(NA,(0+3)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-3))])
    Data_train$lag4rWSPD2_T<-c(rep(NA,(0+4)),Data_train$rWSPD2_T[(1:(length(Data_train$rWSPD2_T)-4))])
    Data_train$rWSPDind_T<-as.factor(ifelse(Data_train$rWSPD<= 0.25,1,0))
    formula_1<-SWHft ~ lag1SWHft+lag2SWHft+lag3SWHft+lag4SWHft+
      lag1rWSPD2_T + lag2rWSPD2_T + lag3rWSPD2_T +lag4rWSPD2_T+
      z2_T+ z3_T +z4_T+
      I(lag1SWHft*rWSPD2_T)+
      I(lag3SWHft*lag2rWSPD2_T)+
      rWSPDind_T
    formula_2<-~ rWSPD2_T
    runLag2_Int<-chngptm(formula.1= formula_1, formula.2=formula_2 , family="gaussian", data=Data_train, 
                         type="segmented", var.type="bootstrap", weights=NULL,ci.bootstrap.size = 100)
    
    ResCalm<-residuals(runLag2_Int)
    ModelRes_Calm<-garchFit(~garch(1,2), data = ResCalm,trace=FALSE,include.mean =FALSE)
    Ind_Res_Calm<-ifelse(ResCalm<0,-1,1)
    sigW<-ModelRes_Calm@sigma.t
    RhatCalm<-(sigW*Ind_Res_Calm)
    setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
    load("44039_reg_merged.Rda")
    source("Data_Function_CS_BC.R")
    Data<-data.in_CS_bc(k,"Calm")
    Data$lag1rWSPD2_T<-c(rep(NA,(0+1)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-1))])
    Data$lag2rWSPD2_T<-c(rep(NA,(0+2)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-2))])
    Data$lag3rWSPD2_T<-c(rep(NA,(0+3)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-3))])
    Data$lag4rWSPD2_T<-c(rep(NA,(0+4)),Data$rWSPD2_T[(1:(length(Data$rWSPD2_T)-4))])
    Data$rWSPDind_T<-as.factor(ifelse(Data$rWSPD<= 0.25,1,0))
    vec<-ifelse((Data$rWSPD2_T-runLag2_Int$coefficients[18])>0,(Data$rWSPD2_T-runLag2_Int$coefficients[18]),0)
    Data1<-cbind(rep(1,nrow(Data))
                 ,Data$lag1rWSPD2_T             
                 ,Data$lag2rWSPD2_T      
                 ,Data$lag3rWSPD2_T           
                 ,Data$lag4rWSPD2_T
                 ,Data$z2_T
                 ,Data$z3_T
                 ,Data$z4_T
                 ,(as.numeric(Data$rWSPDind)-1)
                 ,Data$rWSPD2_T
                 ,vec)
    Coef1<-runLag2_Int$coefficients
    Coef11<-Coef1[c(1,6,7,8,9,10,11,12,15,16,17)]
    Coef12<-Coef1[c(2,3,4,5,13,14)]
    X_mat<-matrix(,nrow=nrow(Data),ncol=17)
    X_mat[,c(1,6,7,8,9,10,11,12,15,16,17)]<-Data1
    Pred<-c()
    Pred[1:4]<-0
    Pred_2005_TM_Calm[1:4]<-0
    Pred_2005_TM_Calm_garch[1:4]<-0
    RhatCalm1<-c()
    RhatCalm1<-c()
    RhatCalm1<-RhatCalm
    for(t in 5:nrow(Data1))
    {
      Pred[t]<-Coef11%*%Data1[t,]+(Coef12[1]*Pred[(t-1)])+(Coef12[2]*Pred[(t-2)])+(Coef12[3]*Pred[(t-3)])+
        (Coef12[4]*Pred[(t-4)])+(Coef12[5]*Pred[(t-1)]*Data$rWSPD2_T[t])+(Coef12[6]*Pred[(t-3)]*Data$lag2rWSPD2_T[t])
      X_mat[t,c(2,3,4,5,13,14)]<-c(Pred[(t-1)],Pred[(t-2)],Pred[(t-3)],Pred[(t-4)],(Pred[(t-1)]*Data$rWSPD2_T[t]),
                                   (Pred[(t-3)]*Data$lag2rWSPD2_T[t]))
    }
    Pred_2005_TM_Calm[5:nrow(Data),(j-2004)]<-Pred[-(1:4)]
    V_yhat<-as.vector(diag(X_mat[-(1:4),]%*%vcov(runLag2_Int)[-18,-18]%*%t(X_mat[-(1:4),])))
    Pred_2005_TM_Calm_garch[5:nrow(Data),(j-2004)]<-Pred[-(1:4)]+RhatCalm1
    Pred_2005_TM_Calm_garch_m10sd[5:nrow(Data),(j-2004)]<-Pred[-(1:4)]+RhatCalm1-(10*sqrt(V_yhat))
    Pred_2005_TM_Calm_garch_p10sd[5:nrow(Data),(j-2004)]<-Pred[-(1:4)]+RhatCalm1+(10*sqrt(V_yhat))
  }
  
  Pred_2005_C<-c()
  Pred_2005_C_p10sd<-c()
  Pred_2005_C_m10sd<-c()
  
  Pred_2005_C<-rowMeans(Pred_2005_TM_Calm_garch)
  Pred_2005_C_p10sd<-rowMeans(Pred_2005_TM_Calm_garch_p10sd)
  Pred_2005_C_m10sd<-rowMeans(Pred_2005_TM_Calm_garch_m10sd)
  
  PredSeqCalm[[(k-1973)]]<-ifelse(Pred_2005_C<0,0,Pred_2005_C)
  PredSeqCalm_p10sd[[(k-1973)]]<-ifelse(Pred_2005_C_p10sd<0,0,Pred_2005_C_p10sd)
  PredSeqCalm_m10sd[[(k-1973)]]<-ifelse(Pred_2005_C_m10sd<0,0,Pred_2005_C_m10sd)
}


names(PredSeqCalm)=1974:2004
names(PredSeqCalm_p10sd)=1974:2004
names(PredSeqCalm_p10sd)=1974:2004

April_mean=May_mean=June_mean=July_mean=Aug_mean=Sept_mean=Oct_mean=c()
for(i in 1:length(PredSeqCalm) )
{
  April_mean[i]<-mean(na.omit(PredSeqCalm[[i]][1:720]))
  May_mean[i]<-mean(na.omit(PredSeqCalm[[i]][721:1464]))
  June_mean[i]<-mean(na.omit(PredSeqCalm[[i]][1465:2184]))
  July_mean[i]<-mean(na.omit(PredSeqCalm[[i]][2185:2928]))
  Aug_mean[i]<-mean(na.omit(PredSeqCalm[[i]][2929:3672]))
  Sept_mean[i]<-mean(na.omit(PredSeqCalm[[i]][3673:4392]))
  Oct_mean[i]<-mean(na.omit(PredSeqCalm[[i]][4393:5136]))
}

kseq=2005:2013
for(i in kseq )
{
  D1<-Data_Buoy[which(Data_Buoy$year==i),]
  D2<-D1$SWHft[which(D1$month==4|D1$month==5|D1$month==6|D1$month==7|D1$month==8|D1$month==9|D1$month==10)]
  April_mean[(i-1973)]<-mean(na.omit(D2[1:720]))
  May_mean[(i-1973)]<-mean(na.omit(D2[721:1464]))
  June_mean[(i-1973)]<-mean(na.omit(D2[1465:2184]))
  July_mean[(i-1973)]<-mean(na.omit(D2[2185:2928]))
  Aug_mean[(i-1973)]<-mean(na.omit(D2[2929:3672]))
  Sept_mean[(i-1973)]<-mean(na.omit(D2[3673:4392]))
  Oct_mean[(i-1973)]<-mean(na.omit(D2[4393:5136]))
}
boxplot(April_mean,May_mean,June_mean,July_mean,Aug_mean,Sept_mean,Oct_mean,names=c("April","May","June","July","Aug","Sept","Oct"))


April_max=May_max=June_max=July_max=Aug_max=Sept_max=Oct_max=c()
for(i in 1:length(PredSeqCalm) )
{
  April_max[i]<-max(na.omit(PredSeqCalm[[i]][1:720]))
  May_max[i]<-max(na.omit(PredSeqCalm[[i]][721:1464]))
  June_max[i]<-max(na.omit(PredSeqCalm[[i]][1465:2184]))
  July_max[i]<-max(na.omit(PredSeqCalm[[i]][2185:2928]))
  Aug_max[i]<-max(na.omit(PredSeqCalm[[i]][2929:3672]))
  Sept_max[i]<-max(na.omit(PredSeqCalm[[i]][3673:4392]))
  Oct_max[i]<-max(na.omit(PredSeqCalm[[i]][4393:5136]))
}

kseq=2005:2013
for(i in kseq )
{
  D1<-Data_Buoy[which(Data_Buoy$year==i),]
  D2<-D1$SWHft[which(D1$month==4|D1$month==5|D1$month==6|D1$month==7|D1$month==8|D1$month==9|D1$month==10)]
  April_max[(i-1973)]<-max(na.omit(D2[1:720]))
  May_max[(i-1973)]<-max(na.omit(D2[721:1464]))
  June_max[(i-1973)]<-max(na.omit(D2[1465:2184]))
  July_max[(i-1973)]<-max(na.omit(D2[2185:2928]))
  Aug_max[(i-1973)]<-max(na.omit(D2[2929:3672]))
  Sept_max[(i-1973)]<-max(na.omit(D2[3673:4392]))
  Oct_max[(i-1973)]<-max(na.omit(D2[4393:5136]))
}
May_max[5]<-3.8
boxplot(April_max,May_max,June_max,July_max,Aug_max,Sept_max,Oct_max,names=c("April","May","June","July","Aug","Sept","Oct"))

DataBP<-data.frame(matrix(nrow=480,ncol=5))
#provide column names
colnames(DataBP) <- c('MonthG', 'SWHft','SWHft_ind',"Year","SWHft1")
DataBP$MonthG<-factor(rep(c("Jan","Feb","Mar",
                            "Apr","May","Jun",
                            "Jul","Aug","Sep",
                            "Oct","Nov","Dec"),each=40),levels=c("Jan","Feb","Mar",
                                                                 "Apr","May","Jun",
                                                                 "Jul","Aug","Sep",
                                                                 "Oct","Nov","Dec"))
DataBP$Year<-rep(1974:2013,12)
DataBP$SWHft_ind<-NA
DataBP$SWHft[which(DataBP$MonthG=="Nov")]<-Nov_mean
DataBP$SWHft[which(DataBP$MonthG=="Dec")]<-Dec_mean
DataBP$SWHft[which(DataBP$MonthG=="Jan")]<-Jan_mean
DataBP$SWHft[which(DataBP$MonthG=="Feb")]<-Feb_mean
DataBP$SWHft[which(DataBP$MonthG=="Mar")]<-March_mean
DataBP$SWHft[which(DataBP$MonthG=="Apr")]<-April_mean
DataBP$SWHft[which(DataBP$MonthG=="May")]<-May_mean
DataBP$SWHft[which(DataBP$MonthG=="Jun")]<-June_mean
DataBP$SWHft[which(DataBP$MonthG=="Jul")]<-July_mean
DataBP$SWHft[which(DataBP$MonthG=="Aug")]<-Aug_mean
DataBP$SWHft[which(DataBP$MonthG=="Sep")]<-Sept_mean
DataBP$SWHft[which(DataBP$MonthG=="Oct")]<-Oct_mean
ind1<-c(32:40,72:80,112:120,152:160,192:200,232:240,272:280,312:320,352:360,392:400,432:440,472:480)
DataBP$SWHft_ind[ind1]<-DataBP$SWHft[ind1]
DataBP$SWHft1<-DataBP$SWHft
DataBP$SWHft[which(is.na(DataBP$SWHft_ind)==FALSE)]=NA

boxplot(SWHft ~ MonthG, data = DataBP, col = "aliceblue",xlab="Months",pch=16,cex=0.5,ylim=c(0,4))

# Points
stripchart(SWHft_ind ~ MonthG,
           data = DataBP,
           method = "jitter",
           pch=21,
           bg="coral3",
           vertical = TRUE,add=TRUE,cex=0.5)
library(RColorBrewer)
n <- 12
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))


DataBP<-data.frame(matrix(nrow=480,ncol=5))
#provide column names
colnames(DataBP) <- c('MonthG', 'SWHft','SWHft_ind',"Year","SWHft1")
DataBP$MonthG<-factor(rep(c("Jan","Feb","Mar",
                            "Apr","May","Jun",
                            "Jul","Aug","Sep",
                            "Oct","Nov","Dec"),each=40),levels=c("Jan","Feb","Mar",
                                                                 "Apr","May","Jun",
                                                                 "Jul","Aug","Sep",
                                                                 "Oct","Nov","Dec"))
DataBP$Year<-rep(1974:2013,12)
DataBP$SWHft_ind<-NA
DataBP$SWHft[which(DataBP$MonthG=="Nov")]<-Nov_max
DataBP$SWHft[which(DataBP$MonthG=="Dec")]<-Dec_max
DataBP$SWHft[which(DataBP$MonthG=="Jan")]<-Jan_max
DataBP$SWHft[which(DataBP$MonthG=="Feb")]<-Feb_max
DataBP$SWHft[which(DataBP$MonthG=="Mar")]<-March_max
DataBP$SWHft[which(DataBP$MonthG=="Apr")]<-April_max
DataBP$SWHft[which(DataBP$MonthG=="May")]<-May_max
DataBP$SWHft[which(DataBP$MonthG=="Jun")]<-June_max
DataBP$SWHft[which(DataBP$MonthG=="Jul")]<-July_max
DataBP$SWHft[which(DataBP$MonthG=="Aug")]<-Aug_max
DataBP$SWHft[which(DataBP$MonthG=="Sep")]<-Sept_max
DataBP$SWHft[which(DataBP$MonthG=="Oct")]<-Oct_max
ind1<-c(32:40,72:80,112:120,152:160,192:200,232:240,272:280,312:320,352:360,392:400,432:440,472:480)
DataBP$SWHft_ind[ind1]<-DataBP$SWHft[ind1]
DataBP$SWHft1<-DataBP$SWHft
DataBP$SWHft[which(is.na(DataBP$SWHft_ind)==FALSE)]=NA

boxplot(SWHft ~ MonthG, data = DataBP, col = "aliceblue",xlab="Months",pch=16,cex=0.5)

# Points
stripchart(SWHft_ind ~ MonthG,
           data = DataBP,
           method = "jitter",
           pch=21,
           bg="coral3",
           vertical = TRUE,add=TRUE,cex=0.5)
library(RColorBrewer)
n <- 12
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))





#Question: Comparing with buoy WH for 2004-2013. Profile plots.

##################             POT Analysis              ##################

library(lubridate)
Data<-data.frame(seq(ymd_hms("1973-11-1 00:00:00"),ymd_hms(" 2004-10-31 23:00:00"),by="hours"))
colnames(Data)="TimeStamp"
Data$YYYY<-year(Data$TimeStamp)
Data$MM<-month(Data$TimeStamp)
Data$DD<-day(Data$TimeStamp)
Data$hh<-hour(Data$TimeStamp)
Data$EB_SWHft<-NA
summary(Data)
Kseq=c(1974:1975,1977:1979,1981:1983,1985:1987,1989:1991,1993:1995,1997:1999,2001:2003)
Kseq1=c(1976,1980,1984,1988,1992,1996,2000,2004)
for(i in Kseq)
{
  Data$EB_SWHft[which(Data$YYYY==(i-1) & Data$MM==11)]=PredSeqWindy[[(i-1973)]][1:720]
  Data$EB_SWHft[which(Data$YYYY==(i-1) & Data$MM==12)]=PredSeqWindy[[(i-1973)]][721:1464]
  Data$EB_SWHft[which(Data$YYYY==i & Data$MM==1)]=PredSeqWindy[[(i-1973)]][1465:2208]
  Data$EB_SWHft[which(Data$YYYY==i & Data$MM==2)]=PredSeqWindy[[(i-1973)]][2209:2880]
  Data$EB_SWHft[which(Data$YYYY==i & Data$MM==3)]=PredSeqWindy[[(i-1973)]][2881:3624]
}


for(i in Kseq1)
{
  Data$EB_SWHft[which(Data$YYYY==(i-1) & Data$MM==11)]=PredSeqWindy[[(i-1973)]][1:720]
  Data$EB_SWHft[which(Data$YYYY==(i-1) & Data$MM==12)]=PredSeqWindy[[(i-1973)]][721:1464]
  Data$EB_SWHft[which(Data$YYYY==i & Data$MM==1)]=PredSeqWindy[[(i-1973)]][1465:2208]
  Data$EB_SWHft[which(Data$YYYY==i & Data$MM==2)]=PredSeqWindy[[(i-1973)]][2209:2904]
  Data$EB_SWHft[which(Data$YYYY==i & Data$MM==3)]=PredSeqWindy[[(i-1973)]][2905:3648]
}

Kseq=1974:2004

for(i in Kseq)
{
  Data$EB_SWHft[which(Data$YYYY==i & Data$MM==4)]=PredSeqCalm[[(i-1973)]][1:720]
  Data$EB_SWHft[which(Data$YYYY==i & Data$MM==5)]=PredSeqCalm[[(i-1973)]][721:1464]
  Data$EB_SWHft[which(Data$YYYY==i & Data$MM==6)]=PredSeqCalm[[(i-1973)]][1465:2184]
  Data$EB_SWHft[which(Data$YYYY==i & Data$MM==7)]=PredSeqCalm[[(i-1973)]][2185:2928]
  Data$EB_SWHft[which(Data$YYYY==i & Data$MM==8)]=PredSeqCalm[[(i-1973)]][2929:3672]
  Data$EB_SWHft[which(Data$YYYY==i & Data$MM==9)]=PredSeqCalm[[(i-1973)]][3673:4392]
  Data$EB_SWHft[which(Data$YYYY==i & Data$MM==10)]=PredSeqCalm[[(i-1973)]][4393:5136]
}



summary(Data)
getwd()
write.csv(Data,"EnesembleBackcast_1974-2004_1.csv")

Data_Buoy <- read.csv("~/Documents/Documents - Namitha’s MacBook Pro/MarineScience/NVP-spring22/Sreeram_NR/Data_Buoy.csv")
head(Data_Buoy)
library(lubridate)
Data1<-data.frame(seq(ymd_hms("2004-11-1 00:00:00"),ymd_hms(" 2013-10-31 23:00:00"),by="hours"))
colnames(Data1)="TimeStamp"
Data1$YYYY<-year(Data1$TimeStamp)
Data1$MM<-month(Data1$TimeStamp)
Data1$DD<-day(Data1$TimeStamp)
Data1$hh<-hour(Data1$TimeStamp)
Data1$EB_SWHft<-Data_Buoy$SWHft
summary(Data1)
#write.csv(Data1,"Databuoy_2004-2013.csv")
####################              GPD              #################### 
library(ismev)
library(extRemes)
library(POT)
getwd()
setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
Data<-read.csv("EnesembleBackcast_1974-2004_1.csv")
head(Data)
Data1<-read.csv("Databuoy_2004-2013.csv")
head(Data1)
DataComb<-rbind(Data,Data1)
summary(DataComb)
summary(DataComb$EB_SWHft)
#length(which(DataComb$EB_SWHft<0 | DataComb$EB_SWHft>15))
EV_fn<-function(t)
{
  Vec1<-c()
  YYYY=MM=DD<-c()
  n<-length(which(DataComb$hh==0))
  k<-t/24
  for(i in 1:(n/k))
  {
    a<-(t*(i-1))+1
    b<-(t*(i-1))+t
    Vec1[i]<-max(na.omit(DataComb$EB_SWHft[a:b]))
    YYYY[i]<-unique(DataComb$YYYY[a])
    MM[i]<-unique(DataComb$MM[a])
    DD[i]<-unique(DataComb$DD[a])
  }
  data1<-data.frame(matrix(,nrow=n/k,ncol=4))
  data1[,1]<-YYYY
  data1[,2]<-MM
  data1[,3]<-DD
  data1[,4]<-Vec1
  colnames(data1)<-c("YYYY","MM","DD","SWHft_max")
  return(data1)
}

A1<-EV_fn(24)
head(A1)
library(POT)
mle_est <- fitgpd(A1$SWHft_max, 5, "mle")
npy1=length(mle_est$exceed)/length(unique(A1$YYYY));npy1
par(mfrow=c(1,1))
rl.fun<-retlev(mle_est,xlimsup=100,col="red",pch=19,npy=npy1,ci=TRUE)
rl.fun(10)
rl.fun(50)
rl.fun(100)

#getAnywhere(retlev.uvpot)
npy<-length(mle_est$exceed)/length(unique(A1$YYYY))
data <- mle_est$exceed
loc <- mle_est$threshold[1]
scale <- mle_est$param["scale"]
shape <- mle_est$param["shape"]
n <- mle_est$nat
pot.fun <- function(T) {
  p <- rp2prob(T, npy)[, "prob"]
  return(qgpd(p, loc, scale, shape))
}
eps <- 10^(-3)
p_emp <- (1:n - 0.35)/n
xlimsup <- prob2rp((n - 0.35)/n, npy)[, "retper"]
#xlimsup<-100
#117.1429
samp <- rgpd(1000 * n, loc, scale, shape)
samp <- matrix(samp, n, 1000)
samp <- apply(samp, 2, sort)
samp <- apply(samp, 1, sort)
ci_inf <- samp[25, ]
ci_sup <- samp[975, ]
plot(pot.fun, from = 1/npy + eps, to = xlimsup, log = "x",col="red",xlim=c(1/npy + eps,xlimsup),
     main ="Return Level Plot",ylab="Return Level", xlab= "Return Period (Years)",ylim=c(0,32))
points(1/(npy * (1 - p_emp)), sort(data), pch = 20,col="orange",cex=0.5)
lines(1/(npy * (1 - p_emp)), ci_sup, lty = 2)
lines(1/(npy * (1 - p_emp)), ci_inf, lty = 2)
legend("topleft",legend=c("Return Values for SWH(ft)","10-yr: 11.0227 (9.5,13.0)",
                          "50-yr: 14.7095 (10.7,21.6) ","100-yr: 16.6351 (11.2,26.8)"),cex=0.5)
abline(h=rl.fun(10))
abline(v=10)
abline(h=9.5)
abline(h=13.0)

abline(h=rl.fun(50))
abline(v=50)
abline(h=10.7)
abline(h=21.6)

abline(h=rl.fun(100))
abline(v=100)
abline(h=11.0)
abline(h=26.8)



abline(h=11.27071)
abline(v=50)
abline(h=11.2)
abline(h=25.4)
#legend("topleft",legend=c("10-yr"))

abline(h=9.41915)
abline(v=10)
abline(h=8.4)
abline(h=10.7)

rl.fun(10)
rl.fun(50)
rl.fun(100)
#abline(v=118)
################

########Extreme value analysis only with observed WH
setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")

Data1<-read.csv("Databuoy_2004-2013_Observed_Pred.csv")
head(Data1)
#Data2<-cbind(Data1,PredData)
#write.csv(Data2,"Databuoy_2004-2013_Observed_Pred.csv")
DataComb<-Data1
summary(DataComb)
summary(DataComb$EB_SWHft)
#length(which(DataComb$EB_SWHft<0 | DataComb$EB_SWHft>15))

EV_fn1<-function(t)
{
  Vec1<-c()
  YYYY=MM=DD<-c()
  n<-length(which(DataComb$hh==0))
  k<-t/24
  for(i in 1:(n/k))
  {
    a<-(t*(i-1))+1
    b<-(t*(i-1))+t
    Vec1[i]<-max(na.omit(DataComb$EB_SWHft[a:b]))
    YYYY[i]<-unique(DataComb$YYYY[a])
    MM[i]<-unique(DataComb$MM[a])
    DD[i]<-unique(DataComb$DD[a])
  }
  data1<-data.frame(matrix(,nrow=n/k,ncol=4))
  data1[,1]<-YYYY
  data1[,2]<-MM
  data1[,3]<-DD
  data1[,4]<-Vec1
  colnames(data1)<-c("YYYY","MM","DD","SWHft_max")
  return(data1)
}

A1<-EV_fn1(24)
summary(A1)
####
#Largest value
max(A1$SWHft_max)
length(A1$SWHft_max)
#Exceedance Prob per day
1/3287   #1/168
#0.0003
#Exceedance Prob per year
(1/3287)*365
#0.11
1/0.11
#Largest return interval we can report using data is 9.

#### 
library(POT)
mle_est <- fitgpd(A1$SWHft_max, 5, "mle")

npy1=length(mle_est$exceed)/length(unique(A1$YYYY));npy1
par(mfrow=c(1,1))
rl.fun<-retlev(mle_est,col="red",pch=8,npy=npy1,ci=TRUE,ylim=c(5,12))
abline(v=29)
round(rl.fun(10),4)
rl.fun(50)
rl.fun(100)

#getAnywhere(retlev.uvpot)
npy<-length(mle_est$exceed)/length(unique(A1$YYYY))
data <- mle_est$exceed
loc <- mle_est$threshold[1]
scale <- mle_est$param["scale"]
shape <- mle_est$param["shape"]
n <- mle_est$nat
pot.fun <- function(T) {
  p <- rp2prob(T, npy)[, "prob"]
  return(qgpd(p, loc, scale, shape))
}
eps <- 10^(-3)
p_emp <- (1:n - 0.35)/n
xlimsup <- prob2rp((n - 0.35)/n, npy)[, "retper"];xlimsup
##############
1/(npy1*(1-0.9979167))
#28.57189
################
#xlimsup<-100
samp <- rgpd(1000 * n, loc, scale, shape)
samp <- matrix(samp, n, 1000)
samp <- apply(samp, 2, sort)
samp <- apply(samp, 1, sort)
ci_inf <- samp[25, ]
ci_sup <- samp[975, ]

plot(pot.fun, from = 1/npy + eps, to = xlimsup, log = "x",col="red",xlim=c(1/npy + eps,100),
     main ="Return Level Plot",ylab="Return Level", xlab= "Return Period (Years)",ylim=c(0,32))
points(1/(npy * (1 - p_emp)), sort(data), pch = 20,col="orange",cex=0.5)
lines(1/(npy * (1 - p_emp)), ci_sup, lty = 2)
lines(1/(npy * (1 - p_emp)), ci_inf, lty = 2)
legend("topleft",legend=c("Return Values for SWH(ft)","10-yr: 9.3625",
                          "50-yr: 10.0601 ","100-yr: 10.3008"),cex=0.5)
abline(v=10)
abline(h=10.19)



plot(pot.fun, from = 1/npy + eps, to = xlimsup, log = "x",col="red",xlim=c(1/npy + eps,xlimsup),
     main ="Return Level Plot",ylab="Return Level", xlab= "Return Period (Years)",ylim=c(mle_est$threshold[1],max(data)))
points(1/(npy * (1 - p_emp)), sort(data), pch = 8)
lines(1/(npy * (1 - p_emp)), ci_sup, lty = 2)
lines(1/(npy * (1 - p_emp)), ci_inf, lty = 2)
##############################

EV_fn2<-function(t)
{
  Vec1<-c()
  YYYY=MM=DD<-c()
  n<-length(which(DataComb$hh==0))
  k<-t/24
  for(i in 1:(n/k))
  {
    a<-(t*(i-1))+1
    b<-(t*(i-1))+t
    Vec1[i]<-max(na.omit(DataComb$PredData[a:b]))
    YYYY[i]<-unique(DataComb$YYYY[a])
    MM[i]<-unique(DataComb$MM[a])
    DD[i]<-unique(DataComb$DD[a])
  }
  data1<-data.frame(matrix(,nrow=n/k,ncol=4))
  data1[,1]<-YYYY
  data1[,2]<-MM
  data1[,3]<-DD
  data1[,4]<-Vec1
  colnames(data1)<-c("YYYY","MM","DD","SWHft_max_pred")
  return(data1)
}

A2<-EV_fn2(24)
summary(A2)
library(POT)
mle_est <- fitgpd(A2$SWHft_max_pred, 5, "mle")
npy1=length(mle_est$exceed)/length(unique(A2$YYYY));npy1
par(mfrow=c(1,1))
rl.fun<-retlev(mle_est,col="red",npy=npy1,ci=TRUE)
abline(v=29)
rl.fun(10)
rl.fun(50)
rl.fun(100)

#getAnywhere(retlev.uvpot)
npy<-length(mle_est$exceed)/length(unique(A2$YYYY))
data <- mle_est$exceed
loc <- mle_est$threshold[1]
scale <- mle_est$param["scale"]
shape <- mle_est$param["shape"]
n <- mle_est$nat
pot.fun <- function(T) {
  p <- rp2prob(T, npy)[, "prob"]
  return(qgpd(p, loc, scale, shape))
}
eps <- 10^(-3)
p_emp <- (1:n - 0.35)/n
xlimsup <- prob2rp((n - 0.35)/n, npy)[, "retper"]
##############
1/(npy1*(1-0.9984018))
#28.57096
################
#xlimsup<-100
samp <- rgpd(1000 * n, loc, scale, shape)
samp <- matrix(samp, n, 1000)
samp <- apply(samp, 2, sort)
samp <- apply(samp, 1, sort)
ci_inf <- samp[25, ]
ci_sup <- samp[975, ]
plot(pot.fun, from = 1/npy + eps, to = xlimsup, log = "x",col="red",xlim=c(1/npy + eps,xlimsup),
     main ="Return Level Plot",ylab="Return Level", xlab= "Return Period (Years)")
points(1/(npy * (1 - p_emp)), sort(data), pch = 8)
lines(1/(npy * (1 - p_emp)), ci_sup, lty = 2)
lines(1/(npy * (1 - p_emp)), ci_inf, lty = 2)
abline(v=29)

#rl.fun(10)
#npy  #16.8  #21.9
#scale
#9.362469  #9.72181 
#rl.fun(50)
#scale 
#10.06013 #10.71723 
#rl.fun(100)
#scale 
#10.3008  #11.1005 



