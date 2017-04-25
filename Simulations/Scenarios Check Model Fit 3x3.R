############################# Check model fit
setwd("Z:/Wetland Birds PostDoc/Occupancy models/Simulations")
colMax <- function(data){
       v<-suppressWarnings( apply(data, 2, max, na.rm = TRUE) )
       v<-ifelse(v=="-Inf",NA,v)
       return(v)
}
sumF<-function(x){ ## Function for discrepancy measures
   wNNA<-which(!is.na(x))
   fit<-sum(as.vector(x[wNNA]))
   return(fit)
}

years<-2005:2014
years<-scale(years)[,1]
yearsplot<-c(2005:2014) #c(1:10) #
nyear<-length(years)
nsite<-5
days<-90:180
days<-scale(days)[,1]
nday<-length(days)
nscn<-9 #H,M,L occupancy; +,- occ trends; +,- no. obs; +,- trends in detection
nit<-10000  #iterations
scn.lbls<-expression("1: High Site Use","2: Medium Site Use","3: Low Site Use", "4: " %up% "Site Use", "5: " %down%  "Site Use",
                     "6: " %up%  "no. Daily Visits", "7: " %down% "no. Daily Visits", "8: " %up% "Detectability", "9: " %down% "Detectability")
exp1<-expression(bar(psi[t])==bar(over(sum(u[d,t,i], n, i == 1), n)) ) # Occupancy across sites
exp2<-expression(bar(z[t,i])==sum(over(u[d,t,i], n), d == 1, n), "Simulated State","Observed Data","Re-Simulated Data")

  
load("PropObs.rdata")
load(file="Simulated State.Rdata")

pdf(paste0("Occupancy year 3x3 Scn.pdf"),14,9) #windows(14,9)
par(mfrow=c(3,3), mar=c(2,2,2,1),oma=c(2,2,0,0))
for(scn in 1:9){
  load(paste("OccParamScn",scn,".RData",sep=""))
  load(paste("Simulated Obs Scn",scn,".RData", sep=""))
  
  yM<-matrix(NA, nyear, nsite)
  y<-Objs$simOcc
  
 for(i in 1:nsite){
    yM[,i]<-colMeans(apply(y[,,,i],3,colMax), na.rm=TRUE)
 }# end nsite
  
# plot psi years
yOcc<-apply(U[,,,scn,999],2,colMeans)
zM<-array(NA, dim=c(3,nyear,nsite))
psifs<-array(NA,dim=c(3,nday,nyear))
psiM<-matrix(NA, 3,nyear)
for(t in 1:nyear){
       for(d in 1:nday){
         psifs[2,d,t]<-median(OccVar[,paste("psi.fs[",d,",",t,"]",sep="")])
         psifs[1,d,t]<-quantile(OccVar[,paste("psi.fs[",d,",",t,"]",sep="")], p=0.975)
         psifs[3,d,t]<-quantile(OccVar[,paste("psi.fs[",d,",",t,"]",sep="")], p=0.025)
       }
  for(i in 1:nsite){
    zM[2,t,i]<-median(OccVar[,paste("zM[",t,",",i,"]",sep="")])
    zM[1,t,i]<-quantile(OccVar[,paste("zM[",t,",",i,"]",sep="")], p=0.975)
    zM[3,t,i]<-quantile(OccVar[,paste("zM[",t,",",i,"]",sep="")], p=0.025)
  }
  psiM[,t]<-rowMeans(zM[,t,])
}

plot(yearsplot,psiM[2,], type="n", xlab="", ylab="", ylim=c(0,1))
polygon( c(yearsplot,rev(yearsplot)),c(apply(psifs[1,,],2,mean),rev(apply(psifs[3,,],2,mean))), col="gray80", border=0)
lines(yearsplot,apply(psifs[2,,],2,mean),lty=1)
points(yearsplot,rowMeans(yM), pch=20)
points(yearsplot,colMeans(yOcc), pch=19, col="blue")
mtext(scn.lbls[scn],3, font=2, outer=FALSE, line=0, adj=0.5 )
if(scn==4) mtext("Mean Site Use",2, font=2, outer=TRUE, line=0.5, adj=0.5, las=0 )
if(scn==8) mtext("Year",1, outer=TRUE, font=2, line= 0.5, adj=0.5)
} # end scn
dev.off()

###################################################### 
###Plot Fit of P(detection)
load("Simulation Parameters.RData") #ParamsSim
load("PropObs.rdata")

par(mfrow=c(3,3), mar=c(2,2,2,1),oma=c(1.5,2,0,0), las=1)
for(scn in 1:9){
  load(paste("OccParamScn",scn,".RData",sep=""))
  load(paste("Simulated Obs Scn",scn,".RData", sep=""))
  SPPL = Objs$SPPLsim[,,,]
  
  n=nrow(OccVar)/2
  ns<- sample(1:(n*2),n)
  dCoef1<-OccVar[ns,grepl("dCoef1",colnames(OccVar))]
  dCoef2<-OccVar[ns,grepl("dCoef2",colnames(OccVar))]
  
  delta<-array(NA,dim=c(n,nyear))
  deltasim<-array(NA,dim=c(nyear))
  
  delta[,1]<- exp(rowMeans(dCoef1[,]) + dCoef2*PropObs[3,1])
  delta[,10]<- exp(rowMeans(dCoef1[,]) + dCoef2*PropObs[3,10])
  deltasim[1]<-exp(mean(ParamsSim$deltaInt[]) + ParamsSim$deltaPrOb * PropObs[3,1]  + ParamsSim$deltaTrd[scn]*(1-1))
  deltasim[10]<-exp(mean(ParamsSim$deltaInt[]) + ParamsSim$deltaPrOb * PropObs[3,10] + ParamsSim$deltaTrd[scn]*(10-1))
  
  curve(1 - median(delta[,1])/(x + median(delta[,1])), from=0, to=70, ylim=c(0,1), ylab="")
  dqH<-quantile(delta[,1], 0.975)
  dqL<-quantile(delta[,1], 0.025)
  polygon(c(0:70,70:0), c(1 - dqH/(c(0:70) + dqH),1 - dqL/(c(70:0) + dqL)), col="gray50", border=NA)
  dqH<-quantile(delta[,10], 0.975)
  dqL<-quantile(delta[,10], 0.025)
  polygon(c(0:70,70:0), c(1 - dqH/(c(0:70) + dqH),1 - dqL/(c(70:0) + dqL)), col=ifelse(scn==3,"#E5E5E580","gray80"), border=NA)
  
  curve(1 - median(delta[,1])/(x + median(delta[,1])), add=TRUE)
  curve(1 - median(delta[,10])/(x + median(delta[,10])), add=TRUE, lty=3)
  
  curve(1 - deltasim[1]/(x + deltasim[1]), col="blue", add=TRUE, lwd=1)
  curve(1 - deltasim[10]/(x + deltasim[10]), col="blue", add=TRUE, lwd=1, lty=3)
  
  x0<-20
  y0<-1-deltasim[1]/(x0 + deltasim[1])
  x1<-20
  y1<-1-deltasim[10]/(x1 + deltasim[10])
  arrows(x0,y0,x1,y1, length=0.1, col=2, lwd=2)
  
  mtext(scn.lbls[scn],3, font=2, outer=FALSE, line=0, adj=0.5 )
  if(scn==8)mtext("Species List Length (no.)",1, outer=TRUE, font=2, line=0.5, adj=0.55, cex=0.8)
  if(scn==4)mtext("Detection Probability",2, outer=TRUE, font=2, adj=0.5, cex=0.8, las=0, line=0.5)
  rug(jitter(as.vector(SPPL[,,1,]), amount=0.25), side=3)
  rug(jitter(as.vector(SPPL[,,10,]),amount=0.25), col="gray")
}
