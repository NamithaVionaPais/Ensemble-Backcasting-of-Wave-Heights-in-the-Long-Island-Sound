library(lubridate)
Datap10<-data.frame(seq(ymd_hms("1973-11-1 00:00:00"),ymd_hms(" 2004-10-31 23:00:00"),by="hours"))
colnames(Datap10)="TimeStamp"
Datap10$YYYY<-year(Datap10$TimeStamp)
Datap10$MM<-month(Datap10$TimeStamp)
Datap10$DD<-day(Datap10$TimeStamp)
Datap10$hh<-hour(Datap10$TimeStamp)
Datap10$EB_SWHft<-NA
summary(Datap10)
Kseq=c(1974:1975,1977:1979,1981:1983,1985:1987,1989:1991,1993:1995,1997:1999,2001:2003)
Kseq1=c(1976,1980,1984,1988,1992,1996,2000,2004)

for(i in Kseq)
{
  Datap10$EB_SWHft[which(Datap10$YYYY==(i-1) & Datap10$MM==11)]=PredSeqWindy_p10sd[[(i-1973)]][1:720]
  Datap10$EB_SWHft[which(Datap10$YYYY==(i-1) & Datap10$MM==12)]=PredSeqWindy_p10sd[[(i-1973)]][721:1464]
  Datap10$EB_SWHft[which(Datap10$YYYY==i & Datap10$MM==1)]=PredSeqWindy_p10sd[[(i-1973)]][1465:2208]
  Datap10$EB_SWHft[which(Datap10$YYYY==i & Datap10$MM==2)]=PredSeqWindy_p10sd[[(i-1973)]][2209:2880]
  Datap10$EB_SWHft[which(Datap10$YYYY==i & Datap10$MM==3)]=PredSeqWindy_p10sd[[(i-1973)]][2881:3624]
}


for(i in Kseq1)
{
  Datap10$EB_SWHft[which(Datap10$YYYY==(i-1) & Datap10$MM==11)]=PredSeqWindy_p10sd[[(i-1973)]][1:720]
  Datap10$EB_SWHft[which(Datap10$YYYY==(i-1) & Datap10$MM==12)]=PredSeqWindy_p10sd[[(i-1973)]][721:1464]
  Datap10$EB_SWHft[which(Datap10$YYYY==i & Datap10$MM==1)]=PredSeqWindy_p10sd[[(i-1973)]][1465:2208]
  Datap10$EB_SWHft[which(Datap10$YYYY==i & Datap10$MM==2)]=PredSeqWindy_p10sd[[(i-1973)]][2209:2904]
  Datap10$EB_SWHft[which(Datap10$YYYY==i & Datap10$MM==3)]=PredSeqWindy_p10sd[[(i-1973)]][2905:3648]
}

Kseq=1974:2004

for(i in Kseq)
{
  Datap10$EB_SWHft[which(Datap10$YYYY==i & Datap10$MM==4)]=PredSeqCalm_p10sd[[(i-1973)]][1:720]
  Datap10$EB_SWHft[which(Datap10$YYYY==i & Datap10$MM==5)]=PredSeqCalm_p10sd[[(i-1973)]][721:1464]
  Datap10$EB_SWHft[which(Datap10$YYYY==i & Datap10$MM==6)]=PredSeqCalm_p10sd[[(i-1973)]][1465:2184]
  Datap10$EB_SWHft[which(Datap10$YYYY==i & Datap10$MM==7)]=PredSeqCalm_p10sd[[(i-1973)]][2185:2928]
  Datap10$EB_SWHft[which(Datap10$YYYY==i & Datap10$MM==8)]=PredSeqCalm_p10sd[[(i-1973)]][2929:3672]
  Datap10$EB_SWHft[which(Datap10$YYYY==i & Datap10$MM==9)]=PredSeqCalm_p10sd[[(i-1973)]][3673:4392]
  Datap10$EB_SWHft[which(Datap10$YYYY==i & Datap10$MM==10)]=PredSeqCalm_p10sd[[(i-1973)]][4393:5136]
}

#a1<-which(Datap10$EB_SWHft %in% sort(Datap10$EB_SWHft,decreasing =TRUE)[1:20])
#summary(Datap10$EB_SWHft[a1])
#Datap10$EB_SWHft[a1]<-25.61  
#summary(Datap10$EB_SWHft)
#getwd()
write.csv(Datap10,"EnesembleBackcast_1974-2004_p10.csv")




library(ismev)
library(extRemes)
library(POT)
getwd()
setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/MarineScience/HPC_Outputs_CS/HPC_Files/ThresholdModel_CS")
Data<-read.csv("EnesembleBackcast_1974-2004_p10.csv")
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
rl.fun(10)
rl.fun(50)
rl.fun(100)
plot(pot.fun, from = 1/npy + eps, to = xlimsup, log = "x",col="red",xlim=c(1/npy + eps,xlimsup),
     main ="Return Level Plot",ylab="Return Level", xlab= "Return Period (Years)",ylim=c(0,32))
points(1/(npy * (1 - p_emp)), sort(data), pch = 20,col="orange",cex=0.5)
lines(1/(npy * (1 - p_emp)), ci_sup, lty = 2)
lines(1/(npy * (1 - p_emp)), ci_inf, lty = 2)
legend("topleft",legend=c("Return Values for SWH(ft)","10-yr: 17.7961 (14.3,23.0)",
                          "50-yr: 25.9501 (17.2,42.5)  ","100-yr: 30.5825 (18.1,54.5)"),cex=0.5)


abline(v=10)
abline(h=23.0)
abline(h=14.3)


plot(pot.fun, from = 1/npy + eps, to = xlimsup, log = "x",col="red",xlim=c(1/npy + eps,xlimsup),
     main ="Return Level Plot",ylab="Return Level", xlab= "Return Period (Years)",ylim=c(0,55))
points(1/(npy * (1 - p_emp)), sort(data), pch = 20,col="orange",cex=0.5)
lines(1/(npy * (1 - p_emp)), ci_sup, lty = 2)
lines(1/(npy * (1 - p_emp)), ci_inf, lty = 2)
legend("topleft",legend=c("Return Values for SWH(ft)","10-yr: 17.79612 (14.3,23.0)",
                          "50-yr: 25.9501 (17.2,42.5)  ","100-yr: 30.58251 (18.1,54.5)"),cex=0.5)
abline(v=50)
abline(h=17.2)
abline(h=42.5)


abline(v=100)
abline(h=18.1)
abline(h=54.5)

#legend("topleft",legend=c("Return Values for SWH(ft)","10-yr: 9.42 (8.4,10.7)",
#"50-yr: 11.27 (9.2,13.7)","100-yr: 12.07 (9.4,14.6)"),cex=0.5)
abline(h=rl.fun(10))
abline(v=10)
abline(h=rl.fun(50))
abline(v=50)
abline(h=rl.fun(100))
abline(v=100)

abline(h=11.27071)
abline(v=50)
abline(h=9.2)
abline(h=13.7)
#legend("topleft",legend=c("10-yr"))

abline(h=9.41915)
abline(v=10)
abline(h=8.4)
abline(h=10.7)

rl.fun(10)
rl.fun(50)
rl.fun(100)
#abline(v=118)
