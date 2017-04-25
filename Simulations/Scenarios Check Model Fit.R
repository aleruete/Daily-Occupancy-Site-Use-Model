############################# Check model fit
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
nyear<-length(years)
nsite<-5
days<-90:180
days<-scale(days)[,1]
nday<-length(days)
nscn<-9 #H,M,L occupancy; +,- occ trends; +,- no. obs; +,- trends in detection
nit<-10000  #iterations

load("PropObs.rdata")
load(file="Simulated State.Rdata")
for(scn in 1:9){
  load(paste("OccParamScn",scn,".RData",sep=""))
  load(paste("Simulated Obs Scn",scn,".RData", sep=""))
  
  ## Bayesian p-values
  ## Assess model fit using a sums-of-squares-type discrepancy
  n=nrow(OccVar)/2
  ns<- sample(1:(n*2),n)
  dCoef1<-OccVar[ns,grepl("dCoef1",colnames(OccVar))]
  dCoef2<-OccVar[ns,grepl("dCoef2",colnames(OccVar))]

  delta<-array(NA,dim=c(n,nyear,nsite))

  for(i in 1:nsite){
    for(t in 1:nyear){
         delta[,t,i]<- exp(dCoef1[,i] + dCoef2*PropObs[3,t])
    }
  }

  yM<-matrix(NA, nyear, nsite)
  y.newM<-array(NA, dim=c(4,nyear, nsite))
  y.newS<-array(NA, dim=c(n,nyear, nsite)) #all samples of seasonal means
  predictedS<-array(NA, dim=c(n,nyear, nsite))
  sq<-array(NA,dim=c(n,dim(yM)))
  sq.new<-array(NA,dim=c(n,dim(yM)))
  sq.newSMP<-array(NA,dim=c(n,dim(yM)))
  y<-Objs$simOcc
  
  system.time(
  for(i in 1:nsite){
    yM[,i]<-colMeans(apply(y[,,,i],3,colMax), na.rm=T)
     for(t in 1:nyear){
       y.new<-array(NA,dim=c(n,50,nday))
       predicted<-array(NA,dim=c(n,50,nday))

       for(d in 1:nday){
         u<-OccVar[1:n,paste("u[",d,",",t,",",i,"]",sep="")]
         if(Objs$nrep[d,t,i]>0){
           for(j in 1:Objs$nrep[d,t,i]){ #nrepMax to avoid NAs
             p<-1 - delta[,t,i]/(Objs$SPPLsim[j,d,t,i] + delta[,t,i])
             predicted[,j,d]<-u*p
             ## Generate new data
             y.new[,j,d] <- rbinom(n,1,predicted[,j,d])
           }# end nrep
         }else{
         y.new[,,d]<-NA
         predicted[,,d]<-NA
         }
       }# end nday
       y.newS[,t,i]<-colMeans(apply(y.new[,,],1,colMax), na.rm=T)
       predictedS[,t,i]<-colMeans(apply(predicted[,,],1,colMax), na.rm=T)
       #residuals
       sq[,t,i]<-(yM[t,i] - predictedS[,t,i])^2
       sq.new[,t,i]<-(y.newS[,t,i] - predictedS[,t,i])^2
     }# end nyear
     print(i)
     y.newM[1:3,,i]<-apply(y.newS[,,i],2,quantile,probs=c(0.025,0.5,0.975))
     y.newM[4,,i]  <-apply(y.newS[,,i],2,sample,1)
  })# end nsite

  site.fit<-matrix(NA,n,nsite)
  site.fit.new<-matrix(NA,n,nsite)
  site.bpv<-numeric()
  for(i in 1:nsite){
    site.fit[,i]<-rowSums(sq[,,i])#,na.rm=T
    site.fit.new[,i]<-rowSums(sq.new[,,i]) #,na.rm=T
    site.bpv[i]<-mean(site.fit.new[,i]>site.fit[,i])
  }
  
  fit<-apply(sq[,,],1,sum)#,na.rm=T
  fit.new<-apply(sq.new[,,],1,sum) #,na.rm=T

  print(bpvalue <- mean(fit.new>fit))			# Bayesian p-value; should be close to 0.5
  print(site.bpv)
  print(paste("Scenario",scn))
  rm(sq, sq.new)
  gc()


################################################################################
#plot psi years
pdf(paste("Occupancy year Scn",scn,".pdf"),14,9) #windows(14,9)
yOcc<-apply(U[,,,scn,999],2,colMeans)
par(mfrow=c(2,4), mar=c(3,4,1,1),oma=c(0,0,3,0))
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

exp1<-expression(bar(psi[t])==bar(sum(over(u[w,t,i], n), i == 1, n)), bar(bar(z[t,i])))
exp2<-expression(bar(z[t,i])==sum(over(u[w,t,i], n), d == 1, n), "Simulated State","Observed Data","Re-Simulated Data")

yearsplot<-c(2005:2014)
plot(yearsplot,psiM[2,], type="n", xlab="Year", ylab="Occupancy", ylim=c(0,1))
mtext("Simulated Species",3, font=3, outer=TRUE, adj=0.5)
polygon( c(yearsplot,rev(yearsplot)),c(psiM[1,],rev(psiM[3,])), col="#BEBEBE50", border=0)
lines(yearsplot, psiM[2,])
lines(yearsplot,apply(psifs[2,,],2,mean),lty=2)
lines(yearsplot,apply(psifs[1,,],2,mean),lty=3)
lines(yearsplot,apply(psifs[3,,],2,mean),lty=3)
points(yearsplot,rowMeans(y.newM[2,,]), pch=20, col=2)
arrows(yearsplot,rowMeans(y.newM[1,,]),yearsplot,rowMeans(y.newM[3,,]), length = 0.05, angle = 90, code = 3, col=2)
points(yearsplot,rowMeans(yM), pch=20)
points(yearsplot,colMeans(yOcc), pch=19, col="blue")
title(main="Region")
legend("topright",legend=exp1, lty=c(2,1), bty="n")

#Plot occurrence per site
for(i in 1:nsite){
     plot(yearsplot,zM[2,,i], type="n", xlab="Year", ylab="Occurrence", ylim=c(0,1), 
          xaxp=c(yearsplot[1],yearsplot[nyear],nyear-1))
     polygon( c(yearsplot,rev(yearsplot)),c(zM[1,,i],rev(zM[3,,i])), col="#BEBEBE50", border=0)
     lines(yearsplot, zM[2,,i])
     points(yearsplot, y.newM[2,,i], pch=20, col=2) #, col=2
     arrows(yearsplot, y.newM[1,,i], yearsplot, y.newM[3,,i], length = 0.05, angle = 90, code = 3, col=2) #, col=2
     points(yearsplot, yM[,i], pch=20) #, col=2
     points(yearsplot,yOcc[i,], pch=19, col="blue") ## Simulated State
     title(main=paste("Site",i))
     text(2006,1,paste("bpv =",round(site.bpv[i],2)))
     legend("topright",legend=exp2, #"Yearly Occupancy (mean z[t,1:nsite])"
                     lty=c(1,NA,NA,NA), pch=c(NA,19,20,20),col=c(1,4,1,2), bty="n", cex=0.8)
}
## Data vs predictions
par(mar=c(4,4,1,1))
plot(yM,y.newM[4,,], ylim=c(0,max(y.newM[4,,])), xlim=c(0,max(yM)), col="red",xlab="Observed annual occupancy", ylab="Predicted annual occupancy")
abline(0,1)

par(mar=c(4,4,1,1))
plot(fit,fit.new, xlab="SSQ discrepancy fo actual data set", ylab="SSQ discrepancy for predicted data set",
     xlim=c(0,max(fit)), ylim=c(0,max(fit.new)))
abline(0,1)
text(min(fit),max(fit.new),paste("bpv =",round(bpvalue,2)), pos=4)
for(i in 1:nsite){
  points(site.fit[,i],site.fit.new[,i], col=i+1)
  text(mean(site.fit[,i]),mean(site.fit.new[,i]),round(site.bpv[i],2))
}
dev.off()

##############################
#Plot Probability of detection
pdf(paste("p(det)Scn",scn,".pdf"),14,9) #windows(14,9)
par(mfrow=c(2,3), mar=c(4,4,1,1), oma=c(0,0,3,0))
n=nrow(OccVar)
dCoef1<-OccVar[,grepl("dCoef1",colnames(OccVar))]
dCoef2<-OccVar[,grepl("dCoef2",colnames(OccVar))]

delta<-array(NA,dim=c(n,nyear,nsite),dimnames=list(c(1:n), 1:10, 1:nsite))

for(i in 1:nsite){
  for(t in 1:nyear){
    delta[,t,i]<- exp(dCoef1[,i] + dCoef2 * PropObs[3,t]) 
  }
}

for(i in 1:nsite){
  plot(1, 1,xlab="Species list length", ylab="p(Detection)",xlim=c(0,100), ylim=c(0,1), type="n")
  title(main=paste("Site",i))
  if(i==1) mtext(paste("Scenario", scn),3, font=2, outer=TRUE, adj=0.5 )
  deltaQ<-quantile(delta[,1,i],probs=c(0.025, 0.975))
  polygon(c(0:100,100:0),c(1 - deltaQ[1]/(c(0:100) + deltaQ[1]),1 - deltaQ[2]/(c(100:0) + deltaQ[2])), col="gray90", border=0)
  deltaQ<-quantile(delta[,10,i],probs=c(0.025, 0.975))
  polygon(c(0:100,100:0),c(1 - deltaQ[1]/(c(0:100) + deltaQ[1]),1 - deltaQ[2]/(c(100:0) + deltaQ[2])), col="gray80", border=0)
  curve(1 - median(delta[,1,i])/(x + median(delta[,1,i])), lwd=2, add=T)
  curve(1 - median(delta[,10,i])/(x + median(delta[,10,i])), lwd=1, add=T)
  legend("topleft", title="No. Long Lists / No. Obs", legend=c("0.06", "0.25"), lwd=c(2,1),
         fill=c("gray90","gray80"), border=0, bty="n")
}
dev.off()

################################################################################
#plot colonization and extinction rates
pdf(paste("ColExt rate Scn",scn,".pdf"),14,9) #windows(14,9)
pCoef1<-OccVar[,grepl("pCoef",colnames(OccVar))][,1]
pCoef2<-OccVar[,grepl("pCoef",colnames(OccVar))][,2]
pCoef3<-OccVar[,grepl("pCoef",colnames(OccVar))][,3]
p.tau <-OccVar[,grepl("p.tau",colnames(OccVar))]
p.sigma<-sqrt(1  / p.tau)
eta.p<-rnorm(nrow(OccVar), 0, p.sigma)
pt.tau <-OccVar[,grepl("pt.tau",colnames(OccVar))]
pt.sigma<-sqrt(1  / pt.tau)
eta.pT<-rnorm(nrow(OccVar), 0, pt.sigma)

gCoef1<-OccVar[,grepl("gCoef",colnames(OccVar))][,1]
gCoef2<-OccVar[,grepl("gCoef",colnames(OccVar))][,2]
gCoef3<-OccVar[,grepl("gCoef",colnames(OccVar))][,3]
g.tau <-OccVar[,grepl("g.tau",colnames(OccVar))]
g.sigma<-sqrt(1  / g.tau)
eta.g<-rnorm(nrow(OccVar), 0, g.sigma)
gt.tau <-OccVar[,grepl("gt.tau",colnames(OccVar))]
gt.sigma<-sqrt(1  / gt.tau)
eta.gT<-rnorm(nrow(OccVar), 0, gt.sigma)


t=1
plot(c(1,nday), c(0,1), type="n", ylab="Colonization/Persistance rate", xlab="Day")
if(t==1) mtext(paste("Scenario", scn),3, font=2, outer=TRUE, adj=0.5 )
phi<-matrix(NA, nrow=n, ncol=nday-1)
gamma<-matrix(NA, nrow=n, ncol=nday-1)
for(d in 1:(nday-1)){
  phi[,d]  <- pnorm(pCoef1 + pCoef2 * days[d] + pCoef3 * (days[d])^2 )
  gamma[,d]<- pnorm(gCoef1 + gCoef2 * days[d] + gCoef3 * (days[d])^2 )

}
boxplot(gamma, col="blue",add = TRUE, at = (1:(nday-1))+0.25, yaxt="n",xaxt="n",outline = FALSE)
boxplot(phi, col="red",add = TRUE, at = (1:(nday-1))-0.25, yaxt="n",xaxt="n",outline = FALSE)

dev.off()
} #end scn