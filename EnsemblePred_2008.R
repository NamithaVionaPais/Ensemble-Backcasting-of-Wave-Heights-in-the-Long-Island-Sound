#Ensemble_TM-GARCH
setwd('/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs/Threshold_Model')
load("44039_reg_merged.Rda")
source("Data_Function.R")
library(fGarch)
library(chngpt)
euclidean <- function(a, b) sqrt(sum((a - b)^2))/length(a)
k=2008
Data<-data.in("central",k,"Windy")
Pred_2008_TM_Windy<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Windy_garch<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Windy_garch_p1sd<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Windy_garch_m1sd<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Windy_garch_p2sd<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Windy_garch_m2sd<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Windy_garch_p3sd<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Windy_garch_m3sd<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Windy_garch_p5sd<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Windy_garch_m5sd<-matrix(nrow=nrow(Data),ncol=9)

for(j in c(2005,2006,2007,2009,2010,2011,2013))
{
  Data_train<-data.in("central",j,"Windy")
  formula_1<-SWHft ~ lag1SWHft+lag2SWHft+lag3SWHft+lag4SWHft+
    lag1rWSPD2 + lag2rWSPD2 + lag3rWSPD2 +lag4rWSPD2+
    z2+ z3 +z4+
    I(lag1SWHft*rWSPD2)+
    I(lag3SWHft*lag2rWSPD2)+
    rWSPDind
  formula_2<-~ rWSPD2
  runLag2_Int<-chngptm(formula.1= formula_1, formula.2=formula_2 , family="gaussian", data=Data_train, 
                       type="segmented", var.type="bootstrap", weights=NULL,ci.bootstrap.size = 100)
  #runLag2_Int<-readRDS( paste('TM_',j,'_Windy_slurm','.RDS',sep=""))
  
  ResWindy<-residuals(runLag2_Int)
  ModelRes_Windy<-garchFit(~garch(1,2), data = ResWindy,trace=FALSE,include.mean = FALSE)
  Ind_Res_Windy<-ifelse(ResWindy<0,-1,1)
  sigW<-ModelRes_Windy@sigma.t
  RhatWindy<-(sigW*Ind_Res_Windy)
  vec<-ifelse((Data$rWSPD2-runLag2_Int$coefficients[18])>0,(Data$rWSPD2-runLag2_Int$coefficients[18]),0)
  Data1<-cbind(rep(1,nrow(Data))
               ,Data$lag1rWSPD2             
               ,Data$lag2rWSPD2      
               ,Data$lag3rWSPD2            
               ,Data$lag4rWSPD2
               ,Data$z2
               ,Data$z3
               ,Data$z4
               ,(as.numeric(Data$rWSPDind)-1)
               ,Data$rWSPD2
               ,vec)
  
  
  Coef1<-runLag2_Int$coefficients
  Coef11<-Coef1[c(1,6,7,8,9,10,11,12,15,16,17)]
  Coef12<-Coef1[c(2,3,4,5,13,14)]
  X_mat<-matrix(,nrow=nrow(Data),ncol=17)
  X_mat[,c(1,6,7,8,9,10,11,12,15,16,17)]<-Data1
  
  Pred<-c()
  Pred[1:5]<-0
  Pred_2008_TM_Windy[1:5]<-0
  Pred_2008_TM_Windy_garch[1:5]<-0
  #RhatWindy1<-RhatWindy
  RhatWindy1<-c()
  RhatWindy1[1:1411]<-RhatWindy[1:1411]
  RhatWindy1[1412:1435]<-rep(0,24)
  RhatWindy1[1436:3643]<-RhatWindy[1412:3619]
  for(t in 6:nrow(Data1))
  {
    Pred[t]<-Coef11%*%Data1[t,]+(Coef12[1]*Pred[(t-1)])+(Coef12[2]*Pred[(t-2)])+(Coef12[3]*Pred[(t-3)])+
      (Coef12[4]*Pred[(t-4)])+(Coef12[5]*Pred[(t-1)]*Data$rWSPD2[t])+(Coef12[6]*Pred[(t-3)]*Data$lag2rWSPD2[t])
    X_mat[t,c(2,3,4,5,13,14)]<-c(Pred[(t-1)],Pred[(t-2)],Pred[(t-3)],Pred[(t-4)],(Pred[(t-1)]*Data$rWSPD2[t]),
                                 (Pred[(t-3)]*Data$lag2rWSPD2[t]))
  }
  Pred_2008_TM_Windy[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]
  V_yhat<-as.vector(diag(X_mat[-(1:5),]%*%vcov(runLag2_Int)[-18,-18]%*%t(X_mat[-(1:5),])))
  Pred_2008_TM_Windy_garch[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1
  Pred_2008_TM_Windy_garch_m1sd[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1-sqrt(V_yhat)
  Pred_2008_TM_Windy_garch_p1sd[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1+sqrt(V_yhat)
  Pred_2008_TM_Windy_garch_m2sd[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1-(2*sqrt(V_yhat))
  Pred_2008_TM_Windy_garch_p2sd[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1+(2*sqrt(V_yhat))
  Pred_2008_TM_Windy_garch_m3sd[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1-(3*sqrt(V_yhat))
  Pred_2008_TM_Windy_garch_p3sd[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1+(3*sqrt(V_yhat))
  Pred_2008_TM_Windy_garch_m5sd[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1-(5*sqrt(V_yhat))
  Pred_2008_TM_Windy_garch_p5sd[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1+(5*sqrt(V_yhat))
}


for(j in c(2012))
{
  Data_train<-data.in("central",j,"Windy")
  formula_1<-SWHft ~ lag1SWHft+lag2SWHft+lag3SWHft+lag4SWHft+
    lag1rWSPD2 + lag2rWSPD2 + lag3rWSPD2 +lag4rWSPD2+
    z2+ z3 +z4+
    I(lag1SWHft*rWSPD2)+
    I(lag3SWHft*lag2rWSPD2)+
    rWSPDind
  formula_2<-~ rWSPD2
  runLag2_Int<-chngptm(formula.1= formula_1, formula.2=formula_2 , family="gaussian", data=Data_train, 
                       type="segmented", var.type="bootstrap", weights=NULL,ci.bootstrap.size = 100)
  
  ResWindy<-residuals(runLag2_Int)
  ModelRes_Windy<-garchFit(~garch(1,2), data = ResWindy,trace=FALSE,include.mean = FALSE)
  Ind_Res_Windy<-ifelse(ResWindy<0,-1,1)
  sigW<-ModelRes_Windy@sigma.t
  RhatWindy<-(sigW*Ind_Res_Windy)
  vec<-ifelse((Data$rWSPD2-runLag2_Int$coefficients[18])>0,(Data$rWSPD2-runLag2_Int$coefficients[18]),0)
  Data1<-cbind(rep(1,nrow(Data))
               ,Data$lag1rWSPD2             
               ,Data$lag2rWSPD2      
               ,Data$lag3rWSPD2            
               ,Data$lag4rWSPD2
               ,Data$z2
               ,Data$z3
               ,Data$z4
               ,(as.numeric(Data$rWSPDind)-1)
               ,Data$rWSPD2
               ,vec)
  
  
  Coef1<-runLag2_Int$coefficients
  Coef11<-Coef1[c(1,6,7,8,9,10,11,12,15,16,17)]
  Coef12<-Coef1[c(2,3,4,5,13,14)]
  X_mat<-matrix(,nrow=nrow(Data),ncol=17)
  X_mat[,c(1,6,7,8,9,10,11,12,15,16,17)]<-Data1
  
  Pred<-c()
  Pred[1:5]<-0
  Pred_2008_TM_Windy[1:5]<-0
  Pred_2008_TM_Windy_garch[1:5]<-0
  Pred_2008_TM_Windy_garch_m1sd[1:5]<-0
  Pred_2008_TM_Windy_garch_p1sd[1:5]<-0
  #RhatWindy1<-c()
  RhatWindy1<-RhatWindy
  #RhatWindy1[1:1411]<-RhatWindy[1:1411]
  #RhatWindy1[1412:3619]<-RhatWindy[1436:3643]
  for(t in 6:nrow(Data))
  {
    Pred[t]<-Coef11%*%Data1[t,]+(Coef12[1]*Pred[(t-1)])+(Coef12[2]*Pred[(t-2)])+(Coef12[3]*Pred[(t-3)])+
      (Coef12[4]*Pred[(t-4)])+(Coef12[5]*Pred[(t-1)]*Data$rWSPD2[t])+(Coef12[6]*Pred[(t-3)]*Data$lag2rWSPD2[t])
    X_mat[t,c(2,3,4,5,13,14)]<-c(Pred[(t-1)],Pred[(t-2)],Pred[(t-3)],Pred[(t-4)],(Pred[(t-1)]*Data$rWSPD2[t]),
                                 (Pred[(t-3)]*Data$lag2rWSPD2[t]))
  }
  Pred_2008_TM_Windy[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]
  V_yhat<-as.vector(diag(X_mat[-(1:5),]%*%vcov(runLag2_Int)[-18,-18]%*%t(X_mat[-(1:5),])))
  Pred_2008_TM_Windy_garch[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1
  Pred_2008_TM_Windy_garch_m1sd[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1-sqrt(V_yhat)
  Pred_2008_TM_Windy_garch_p1sd[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1+sqrt(V_yhat)
  Pred_2008_TM_Windy_garch_m2sd[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1-(2*sqrt(V_yhat))
  Pred_2008_TM_Windy_garch_p2sd[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1+(2*sqrt(V_yhat))
  Pred_2008_TM_Windy_garch_m3sd[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1-(3*sqrt(V_yhat))
  Pred_2008_TM_Windy_garch_p3sd[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1+(3*sqrt(V_yhat))
  Pred_2008_TM_Windy_garch_m5sd[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1-(5*sqrt(V_yhat))
  Pred_2008_TM_Windy_garch_p5sd[6:nrow(Data),(j-k+4)]<-Pred[-(1:5)]+RhatWindy1+(5*sqrt(V_yhat))
}

Pred_2008_TM_Windy_garch<-ifelse(Pred_2008_TM_Windy_garch<0,0,Pred_2008_TM_Windy_garch)
Pred_2008_TM_Windy_garch_m1sd<-ifelse(Pred_2008_TM_Windy_garch_m1sd<0,0,Pred_2008_TM_Windy_garch_m1sd)
Pred_2008_TM_Windy_garch_p1sd<-ifelse(Pred_2008_TM_Windy_garch_p1sd<0,0,Pred_2008_TM_Windy_garch_p1sd)
Pred_2008_TM_Windy_garch_m2sd<-ifelse(Pred_2008_TM_Windy_garch_m2sd<0,0,Pred_2008_TM_Windy_garch_m2sd)
Pred_2008_TM_Windy_garch_p2sd<-ifelse(Pred_2008_TM_Windy_garch_p2sd<0,0,Pred_2008_TM_Windy_garch_p2sd)
Pred_2008_TM_Windy_garch_m3sd<-ifelse(Pred_2008_TM_Windy_garch_m3sd<0,0,Pred_2008_TM_Windy_garch_m3sd)
Pred_2008_TM_Windy_garch_p3sd<-ifelse(Pred_2008_TM_Windy_garch_p3sd<0,0,Pred_2008_TM_Windy_garch_p3sd)
Pred_2008_TM_Windy_garch_m5sd<-ifelse(Pred_2008_TM_Windy_garch_m5sd<0,0,Pred_2008_TM_Windy_garch_m5sd)
Pred_2008_TM_Windy_garch_p5sd<-ifelse(Pred_2008_TM_Windy_garch_p5sd<0,0,Pred_2008_TM_Windy_garch_p5sd)

setwd('/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs/Threshold_Model')
load("44039_reg_merged.Rda")
source("Data_Function.R")
library(fGarch)
library(chngpt)
k=2008
Data<-data.in("central",k,"Windy")
ED_Year_W=ED_Year_W_p5sd=ED_Year_W_m5sd=c()
for(i in 1:9)
{
  ED_Year_W[i]<-round(rmse(Pred_2008_TM_Windy_garch[6:nrow(Data),i],Data$SWHft[-(1:5)]),4)
  ED_Year_W_p5sd[i]<-round(rmse(Pred_2008_TM_Windy_garch_p5sd[6:nrow(Data),i],Data$SWHft[-(1:5)]),4)
  ED_Year_W_m5sd[i]<-round(rmse(Pred_2008_TM_Windy_garch_m5sd[6:nrow(Data),i],Data$SWHft[-(1:5)]),4)
}

Pred_2008_W=Pred_2008_W_m1sd=Pred_2008_W_p1sd=c()
Pred_2008_W_m2sd=Pred_2008_W_p2sd=c()
Pred_2008_W_m3sd=Pred_2008_W_p3sd=c()
Pred_2008_W<-rowMeans(Pred_2008_TM_Windy_garch[,-4])
Pred_2008_W_m1sd<-rowMeans(Pred_2008_TM_Windy_garch_m1sd[,-4])
Pred_2008_W_p1sd<-rowMeans(Pred_2008_TM_Windy_garch_p1sd[,-4])
Pred_2008_W_m2sd<-rowMeans(Pred_2008_TM_Windy_garch_m2sd[,-4])
Pred_2008_W_p2sd<-rowMeans(Pred_2008_TM_Windy_garch_p2sd[,-4])
Pred_2008_W_m3sd<-rowMeans(Pred_2008_TM_Windy_garch_m3sd[,-4])
Pred_2008_W_p3sd<-rowMeans(Pred_2008_TM_Windy_garch_p3sd[,-4])
Pred_2008_W_m5sd<-rowMeans(Pred_2008_TM_Windy_garch_m5sd[,-4])
Pred_2008_W_p5sd<-rowMeans(Pred_2008_TM_Windy_garch_p5sd[,-4])


MonthWindy<-c("Nov","Dec","Jan","Feb","Mar")
MonthWindyNum<-c(11,12,1,2,3)
j=2008
par(mar = c(2, 2, 2, 2))
par(mfrow=c(3,2))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)

legend("center",c("Observed-SWHft","Predicted-SWHft","NA-Indicator","5-sd CI"),
       col=c("black","indianred3","yellow2","orange"),pch=19)

i=1
plot(Data$SWHft[which(Data$MM==MonthWindyNum[i])][-(1:5)],type="l",ylim=c(0,10),
     main =paste(MonthWindy[i],(j-1),sep="-"),xaxt='n',lwd=0.5)
lines(Data$SWHft[which(Data$MM==MonthWindyNum[i])][-(1:5)],col=Data$colSWHft[which(Data$MM==MonthWindyNum[i])][-(1:5)],
      pch=20,cex=0.2,type="p",ylim=c(0,10),main =paste(MonthWindy[i],(j-1),sep="-"))
lines(Pred_2008_W[which(Data$MM==MonthWindyNum[i])][-(1:5)],col="indianred3",lwd=0.5)
lines(Pred_2008_W[which(Data$MM==MonthWindyNum[i])][-(1:5)],col="indianred3",type="p",pch=20,cex=0.2)
lines(Pred_2008_W_m5sd[which(Data$MM==MonthWindyNum[i])][-(1:5)],col="orange",type="l",lwd=0.5)
lines(Pred_2008_W_p5sd[which(Data$MM==MonthWindyNum[i])][-(1:5)],col="orange",type="l",lwd=0.5)
#lines(Pred_2008_W_m5sd[which(Data$MM==MonthWindyNum[i])][-(1:5)],col="darkblue",type="l",lwd=0.5)
#lines(Pred_2008_W_p5sd[which(Data$MM==MonthWindyNum[i])][-(1:5)],col="darkblue",type="l",lwd=0.5)


i=2
plot(Data$SWHft[which(Data$MM==MonthWindyNum[i])],type="l",ylim=c(0,10),
     xaxt='n',main =paste(MonthWindy[i],(j-1),sep="-"),lwd=0.5)
lines(Data$SWHft[which(Data$MM==MonthWindyNum[i])],col=Data$colSWHft[which(Data$MM==MonthWindyNum[i])],
      pch=20,cex=0.2,type="p")

lines(Pred_2008_W[which(Data$MM==MonthWindyNum[i])] ,col="indianred3",lwd=0.5)
lines(Pred_2008_W[which(Data$MM==MonthWindyNum[i])] ,col="indianred3",pch=20,cex=0.2,type="p")
lines(Pred_2008_W_m5sd[which(Data$MM==MonthWindyNum[i])],col="orange",type="l",lwd=0.7)
lines(Pred_2008_W_p5sd[which(Data$MM==MonthWindyNum[i])],col="orange",type="l",lwd=0.7)
#lines(Pred_2008_W_m5sd[which(Data$MM==MonthWindyNum[i])],col="darkblue",type="l",lwd=0.5)
#lines(Pred_2008_W_p5sd[which(Data$MM==MonthWindyNum[i])],col="darkblue",type="l",lwd=0.5)

for(i in 3:length(MonthWindy))
{
  plot(Data$SWHft[which(Data$MM==MonthWindyNum[i])],type="l",ylim=c(0,10),main =paste(MonthWindy[i],j,sep="-"),xaxt='n',lwd=0.5)
  lines(Data$SWHft[which(Data$MM==MonthWindyNum[i])],col=Data$colSWHft[which(Data$MM==MonthWindyNum[i])],pch=19,cex=0.2,type="p")
  lines(Pred_2008_W[which(Data$MM==MonthWindyNum[i])] ,col="indianred3",lwd=0.5)
  lines(Pred_2008_W[which(Data$MM==MonthWindyNum[i])] ,col="indianred3",pch=20,cex=0.2,type="p")
  lines(Pred_2008_W_m5sd[which(Data$MM==MonthWindyNum[i])],col="orange",type="l",lwd=0.7)
  lines(Pred_2008_W_p5sd[which(Data$MM==MonthWindyNum[i])],col="orange",type="l",lwd=0.7)
  #lines(Pred_2008_W_m5sd[which(Data$MM==MonthWindyNum[i])],col="darkblue",type="l",lwd=0.5)
  #lines(Pred_2008_W_p5sd[which(Data$MM==MonthWindyNum[i])],col="darkblue",type="l",lwd=0.5)
  
}

library(Metrics)
round(rmse(Pred_2008_W[-(1:5)],Data$SWHft[-(1:5)]),4)

k=2008
Data<-data.in("central",k,"Calm")
Pred_2008_TM_Calm<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Calm_garch<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Calm_garch_p1sd<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Calm_garch_m1sd<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Calm_garch_p2sd<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Calm_garch_m2sd<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Calm_garch_p3sd<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Calm_garch_m3sd<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Calm_garch_p5sd<-matrix(nrow=nrow(Data),ncol=9)
Pred_2008_TM_Calm_garch_m5sd<-matrix(nrow=nrow(Data),ncol=9)

for(j in c(2005:2007,2009:2013))
{
  Data_train<-data.in("central",j,"Calm")
  formula_1<-SWHft ~ lag1SWHft+lag2SWHft+lag3SWHft+lag4SWHft+
    lag1rWSPD2 + lag2rWSPD2 + lag3rWSPD2 +lag4rWSPD2+
    z2+ z3 +z4+
    I(lag1SWHft*rWSPD2)+
    I(lag3SWHft*lag2rWSPD2)+
    rWSPDind
  formula_2<-~ rWSPD2
  runLag2_Int<-chngptm(formula.1= formula_1, formula.2=formula_2 , family="gaussian", data=Data_train, 
                       type="segmented", var.type="bootstrap", weights=NULL,ci.bootstrap.size = 100)
  ResCalm<-residuals(runLag2_Int)
  ModelRes_Calm<-garchFit(~garch(1,2), data = ResCalm,trace=FALSE,include.mean = FALSE)
  Ind_Res_Calm<-ifelse(ResCalm<0,-1,1)
  sigC<-ModelRes_Calm@sigma.t
  RhatCalm<-(sigC*Ind_Res_Calm)
  vec<-ifelse((Data$rWSPD2-runLag2_Int$coefficients[18])>0,(Data$rWSPD2-runLag2_Int$coefficients[18]),0)
  Data1<-cbind(rep(1,nrow(Data))
               ,Data$lag1rWSPD2             
               ,Data$lag2rWSPD2      
               ,Data$lag3rWSPD2            
               ,Data$lag4rWSPD2
               ,Data$z2
               ,Data$z3
               ,Data$z4
               ,(as.numeric(Data$rWSPDind)-1)
               ,Data$rWSPD2
               ,vec)
  Coef1<-runLag2_Int$coefficients
  Coef11<-Coef1[c(1,6,7,8,9,10,11,12,15,16,17)]
  Coef12<-Coef1[c(2,3,4,5,13,14)]
  X_mat<-matrix(,nrow=nrow(Data),ncol=17)
  X_mat[,c(1,6,7,8,9,10,11,12,15,16,17)]<-Data1
  Pred<-c()
  Pred[1:4]<-0
  Pred_2008_TM_Calm[1:4]<-0
  Pred_2008_TM_Calm_garch[1:4]<-0
  RhatCalm1<-c()
  RhatCalm1<-RhatCalm
  for(t in 5:nrow(Data))
  {
    Pred[t]<-Coef11%*%Data1[t,]+(Coef12[1]*Pred[(t-1)])+(Coef12[2]*Pred[(t-2)])+(Coef12[3]*Pred[(t-3)])+
      (Coef12[4]*Pred[(t-4)])+(Coef12[5]*Pred[(t-1)]*Data$rWSPD2[t])+(Coef12[6]*Pred[(t-3)]*Data$lag2rWSPD2[t])
    X_mat[t,c(2,3,4,5,13,14)]<-c(Pred[(t-1)],Pred[(t-2)],Pred[(t-3)],Pred[(t-4)],(Pred[(t-1)]*Data$rWSPD2[t]),
                                 (Pred[(t-3)]*Data$lag2rWSPD2[t]))
  }
  Pred_2008_TM_Calm[5:nrow(Data),(j-k+4)]<-Pred[-(1:4)]
  V_yhat<-as.vector(diag(X_mat[-(1:4),]%*%vcov(runLag2_Int)[-18,-18]%*%t(X_mat[-(1:4),])))
  Pred_2008_TM_Calm_garch[5:nrow(Data),(j-k+4)]<-Pred[-(1:4)]+RhatCalm1
  Pred_2008_TM_Calm_garch_m1sd[5:nrow(Data),(j-k+4)]<-Pred[-(1:4)]+RhatCalm1-sqrt(V_yhat)
  Pred_2008_TM_Calm_garch_p1sd[5:nrow(Data),(j-k+4)]<-Pred[-(1:4)]+RhatCalm1+sqrt(V_yhat)
  Pred_2008_TM_Calm_garch_m2sd[5:nrow(Data),(j-k+4)]<-Pred[-(1:4)]+RhatCalm1-(2*sqrt(V_yhat))
  Pred_2008_TM_Calm_garch_p2sd[5:nrow(Data),(j-k+4)]<-Pred[-(1:4)]+RhatCalm1+(2*sqrt(V_yhat))
  Pred_2008_TM_Calm_garch_m3sd[5:nrow(Data),(j-k+4)]<-Pred[-(1:4)]+RhatCalm1-(3*sqrt(V_yhat))
  Pred_2008_TM_Calm_garch_p3sd[5:nrow(Data),(j-k+4)]<-Pred[-(1:4)]+RhatCalm1+(3*sqrt(V_yhat))
  Pred_2008_TM_Calm_garch_m5sd[5:nrow(Data),(j-k+4)]<-Pred[-(1:4)]+RhatCalm1-(5*sqrt(V_yhat))
  Pred_2008_TM_Calm_garch_p5sd[5:nrow(Data),(j-k+4)]<-Pred[-(1:4)]+RhatCalm1+(5*sqrt(V_yhat))
  
}



Pred_2008_TM_Calm_garch<-ifelse(Pred_2008_TM_Calm_garch<0,0,Pred_2008_TM_Calm_garch)
Pred_2008_TM_Calm_garch_m1sd<-ifelse(Pred_2008_TM_Calm_garch_m1sd<0,0,Pred_2008_TM_Calm_garch_m1sd)
Pred_2008_TM_Calm_garch_p1sd<-ifelse(Pred_2008_TM_Calm_garch_p1sd<0,0,Pred_2008_TM_Calm_garch_p1sd)
Pred_2008_TM_Calm_garch_m2sd<-ifelse(Pred_2008_TM_Calm_garch_m2sd<0,0,Pred_2008_TM_Calm_garch_m2sd)
Pred_2008_TM_Calm_garch_p2sd<-ifelse(Pred_2008_TM_Calm_garch_p2sd<0,0,Pred_2008_TM_Calm_garch_p2sd)
Pred_2008_TM_Calm_garch_m3sd<-ifelse(Pred_2008_TM_Calm_garch_m3sd<0,0,Pred_2008_TM_Calm_garch_m3sd)
Pred_2008_TM_Calm_garch_p3sd<-ifelse(Pred_2008_TM_Calm_garch_p3sd<0,0,Pred_2008_TM_Calm_garch_p3sd)
Pred_2008_TM_Calm_garch_m5sd<-ifelse(Pred_2008_TM_Calm_garch_m5sd<0,0,Pred_2008_TM_Calm_garch_m5sd)
Pred_2008_TM_Calm_garch_p5sd<-ifelse(Pred_2008_TM_Calm_garch_p5sd<0,0,Pred_2008_TM_Calm_garch_p5sd)


k=2008
Data<-data.in("central",k,"Calm")
ED_Year_C=ED_Year_C_p5sd=ED_Year_C_m5sd=c()
for(i in 1:9)
{
  ED_Year_C[i]<-round(rmse(Pred_2008_TM_Calm_garch[6:nrow(Data),i],Data$SWHft[-(1:5)]),4)
  ED_Year_C_p5sd[i]<-round(rmse(Pred_2008_TM_Calm_garch_p5sd[6:nrow(Data),i],Data$SWHft[-(1:5)]),4)
  ED_Year_C_m5sd[i]<-round(rmse(Pred_2008_TM_Calm_garch_m5sd[6:nrow(Data),i],Data$SWHft[-(1:5)]),4)
}

Pred_2008_C=Pred_2008_C_m1sd=Pred_2008_C_p1sd=c()
Pred_2008_C_m2sd=Pred_2008_C_p2sd=c()
Pred_2008_C_m3sd=Pred_2008_C_p3sd=c()
Pred_2008_C<-rowMeans(Pred_2008_TM_Calm_garch[,-4])
Pred_2008_C_m1sd<-rowMeans(Pred_2008_TM_Calm_garch_m1sd[,-4])
Pred_2008_C_p1sd<-rowMeans(Pred_2008_TM_Calm_garch_p1sd[,-4])
Pred_2008_C_m2sd<-rowMeans(Pred_2008_TM_Calm_garch_m2sd[,-4])
Pred_2008_C_p2sd<-rowMeans(Pred_2008_TM_Calm_garch_p2sd[,-4])
Pred_2008_C_m3sd<-rowMeans(Pred_2008_TM_Calm_garch_m3sd[,-4])
Pred_2008_C_p3sd<-rowMeans(Pred_2008_TM_Calm_garch_p3sd[,-4])
Pred_2008_C_m5sd<-rowMeans(Pred_2008_TM_Calm_garch_m5sd[,-4])
Pred_2008_C_p5sd<-rowMeans(Pred_2008_TM_Calm_garch_p5sd[,-4])


j=2008
MonthCalm<-c("Apr","May","Jun","Jul","Aug","Sep","Oct")
MonthCalmNum<-c(4,5,6,7,8,9,10)
par(mar = c(2, 2, 2, 2))
par(mfrow=c(4,2))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
#legend("center",c("Observed-SWHft","Predicted-SWHft","NA-Indicator","3-sd CI"),
       #col=c("black","indianred3","yellow2","orange"),pch=19)
legend("center",c("Observed-SWHft","Predicted-SWHft","NA-Indicator","5-sd CI"),
       col=c("black","indianred3","yellow2","orange"),pch=19)
i=1
plot(Data$SWHft[which(Data$MM==MonthCalmNum[i])][-(1:4)],type="l",
     ylim=c(0,10),main =paste(MonthCalm[i],j,sep="-"),xaxt='n',lwd=0.5)
lines(Data$SWHft[which(Data$MM==MonthCalmNum[i])][-(1:4)],
      col=Data$colSWHft[which(Data$MM==MonthCalmNum[i])][-(1:5)],pch=20,cex=0.2,type="p")
lines(Pred_2008_C[which(Data$MM==MonthCalmNum[i])][-(1:4)],col="indianred3")
lines(Pred_2008_C[which(Data$MM==MonthCalmNum[i])][-(1:4)],col="indianred3",pch=20,cex=0.2,type="p")
lines(Pred_2008_C_m5sd[which(Data$MM==MonthCalmNum[i])][-(1:4)],col="orange",type="l",lwd=0.5)
lines(Pred_2008_C_p5sd[which(Data$MM==MonthCalmNum[i])][-(1:4)],col="orange",type="l",lwd=0.5)

for(i in 2:length(MonthCalm))
{
  plot(Data$SWHft[which(Data$MM==MonthCalmNum[i])],type="l",ylim=c(0,10),
       main =paste(MonthCalm[i],j,sep="-"),xaxt='n',lwd=0.5)
  lines(Data$SWHft[which(Data$MM==MonthCalmNum[i])],col=Data$colSWHft[which(Data$MM==MonthCalmNum[i])],
        pch=20,cex=0.2,type="p",lwd=0.5)
  lines(Pred_2008_C[which(Data$MM==MonthCalmNum[i])] ,col="indianred3")
  lines(Pred_2008_C[which(Data$MM==MonthCalmNum[i])] ,col="indianred3",pch=20,cex=0.2,type="p")
  lines(Pred_2008_C_m5sd[which(Data$MM==MonthCalmNum[i])] ,col="orange",type="l",lwd=0.5)
  lines(Pred_2008_C_p5sd[which(Data$MM==MonthCalmNum[i])] ,col="orange",type="l",lwd=0.5)
}


round(rmse(Pred_2008_C[-(1:5)],Data$SWHft[-(1:5)]),4)
