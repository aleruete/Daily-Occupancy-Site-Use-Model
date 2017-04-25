require(foreach)
require(doParallel)
colMax <- function(data){
       v<-suppressWarnings( apply(data, 2, max, na.rm = TRUE) )
       v<-ifelse(v=="-Inf",NA,v)
       return(v)
}

sppList<-read.csv("StudySpecies.csv")
spp10<-c(1,9,12,17,30,33,45,53,59)

load("OccDataUppland.RData")
nrep<-nrepIND[,,]
Nobs<-apply(nrepMAX,2,sum)
NobsSL<-apply(SL,3,sum)
NobsCL<-apply(CL,3,sum)

Obs<-matrix(NA,3,10,dimnames=list(c("Single","Short","Long"), c(2005:2014)))
Obs[1,]<-Nobs-NobsSL-NobsCL
Obs[2,]<-NobsSL
Obs[3,]<-NobsCL

PropObs<-matrix(NA,3,10,dimnames=list(c("Single","Short","Long"), c(2005:2014)))
PropObs[1,]<-Obs[1,]/Nobs
PropObs[2,]<-Obs[2,]/Nobs
PropObs[3,]<-Obs[3,]/Nobs

#pdf(paste("Occupancy year 3x3.pdf"),14,9) #windows(14,9)
par(mfrow=c(3,3), mar=c(2,2,2,1),oma=c(2,2,0,0))
for(spp in spp10){
load("OccDataUppland.RData")
y<-OccData[spp,,,,]

load(paste0("Results/", sppList[spp,4]," OccVar.RData"))

n=nrow(OccVar)/2
ns<- sample(1:(n*2),n)
dCoef1<-OccVar[ns,grepl("dCoef1",colnames(OccVar))]
dCoef2<-OccVar[ns,grepl("dCoef2",colnames(OccVar))]

delta<-array(NA,dim=c(n,nyear,nsite),dimnames=list(c(1:n), years, sites))

for(i in 1:nsite){
  for(t in 1:nyear){
    delta[,t,i]<- exp(dCoef1[,i] + dCoef2*PropObs[3,t]) 
  }
}

yM<-matrix(NA, nyear, nsite)
y.newM<-array(NA, dim=c(4,nyear, nsite))
y.newS<-array(NA, dim=c(n,nyear, nsite)) #all samples of seasonal means
predictedS<-array(NA, dim=c(n,nyear, nsite))

system.time(
for(i in 1:nsite){
  yM[,i]<-colMeans(apply(y[,,,i],3,colMax), na.rm=TRUE)
  for(t in 1:nyear){
    y.new<-array(NA,dim=c(n,40,nday))
    predicted<-array(NA,dim=c(n,40,nday))
    for(d in 1:nday){
      u<-OccVar[1:n,paste("u[",d,",",t,",",i,"]",sep="")]
      if(nrepMAX[d,t,i]>0){
        for(j in 1:nrepMAX[d,t,i]){ #nrepMax to avoid NAs
          p<-1 - delta[,t,i]/(SPPL[j,d,t,i] + delta[,t,i])
          predicted[,j,d]<-u*p
          ## Generate new data
          y.new[,j,d] <- rbinom(n,1,predicted[,j,d])
        }# end nrep
      }
    }# end nday
    y.newS[,t,i]<-colMeans(apply(y.new[,,],1,colMax), na.rm=TRUE)
    predictedS[,t,i]<-colMeans(apply(predicted[,,],1,colMax), na.rm=TRUE)
  }# end nyear
  print(i)
  y.newM[1:3,,i]<-apply(y.newS[,,i],2,quantile,probs=c(0.025,0.5,0.975))
  y.newM[4,,i]  <-apply(y.newS[,,i],2,sample,1)
})# end nsite


################################################################################
#plot psi years
zM<-array(NA, dim=c(3,nyear,nsite))
psifs<-array(NA,dim=c(3,nday,nyear))
psiM<-matrix(NA, 3,nyear)
for(t in 1:nyear){
  for(d in 1:nday){
    psifs[2,d,t]<-median(OccVar[,paste("psi.fs[",d,",",t,"]",sep="")])    # psi.fs is the daily average of all sites.
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

exp1<-expression(bar(psi[t])==bar(over(sum(u[d,t,i], n, i == 1), n)) ) # Occupancy across sites
exp2<-expression(bar(z[t,i])==sum(over(u[w,t,i], n), d == 1, n), "Observed","Simulated") #     lines(years, Occ25days[,i], col=4)

plot(years,psiM[2,], type="n", xlab="", ylab="", ylim=c(0,1), xaxp=c(years[1],years[nyear],nyear-1))
polygon( c(years,rev(years)),c(apply(psifs[1,,],2,mean),rev(apply(psifs[3,,],2,mean))), col="gray80", border=0)
lines(years,apply(psifs[2,,],2,mean),lty=1)

points(years,rowMeans(y.newM[2,,]), pch=20, col=2)
arrows(years,rowMeans(y.newM[1,,]),years,rowMeans(y.newM[3,,]), length = 0.05, angle = 90, code = 3, col=2)
points(years,rowMeans(yM), pch=20)
mtext(sppList[spp,4],3, font=3, outer=F, line=0, adj=0.5, cex=0.8 )

if(spp==spp10[4]) mtext("Mean Site Use", 2, font=2, outer=TRUE, line=0.5, las=0, cex=0.8)
if(spp==spp10[8]) mtext("Year", 1, outer=TRUE, font=2, line=0.5, adj=0.5, cex=0.8)
}# end spp

################################################################################
#pdf(paste("ColExt rate 3x3.pdf"),14,9) #windows(14,9)
par(mfrow=c(3,3), mar=c(2,2,2,1),oma=c(3,3,0,0))
eg.site<-5
eg.year<-10
for(spp in spp10){
  y<-OccData[spp,,,,]
  load(paste0("Results/", sppList[spp,4]," OccVar.RData"))
  #plot colonization and extinction rates
  pCoef1<-OccVar[,grepl("pCoef",colnames(OccVar))][,1]
  pCoef2<-OccVar[,grepl("pCoef",colnames(OccVar))][,2]
  pCoef3<-OccVar[,grepl("pCoef",colnames(OccVar))][,3]
  p.tau <-OccVar[,grepl("p.tau",colnames(OccVar))]
  p.sigma<-sqrt(1  / p.tau)
  eta.p<-OccVar[1:n,paste("eta.p[",eg.site,"]",sep="")]
  pt.tau <-OccVar[,grepl("pt.tau",colnames(OccVar))]
  pt.sigma<-sqrt(1  / pt.tau)
  eta.pT<-eta.p<-OccVar[1:n,paste("eta.pT[",eg.year,"]",sep="")]
  
  gCoef1<-OccVar[,grepl("gCoef",colnames(OccVar))][,1]
  gCoef2<-OccVar[,grepl("gCoef",colnames(OccVar))][,2]
  gCoef3<-OccVar[,grepl("gCoef",colnames(OccVar))][,3]
  g.tau <-OccVar[,grepl("g.tau",colnames(OccVar))]
  g.sigma<-sqrt(1  / g.tau)
  eta.g<-OccVar[1:n,paste("eta.g[",eg.site,"]",sep="")]
  gt.tau <-OccVar[,grepl("gt.tau",colnames(OccVar))]
  gt.sigma<-sqrt(1  / gt.tau)
  eta.gT<-OccVar[1:n,paste("eta.gT[",eg.year,"]",sep="")]
  
  plot(days[c(1,nday)], c(0,1), type="n")
  phi<-matrix(NA, nrow=nrow(OccVar),ncol=nday-1)
  gamma<-matrix(NA, nrow=nrow(OccVar),ncol=nday-1)
  for(d in 1:(nday-1)){
    phi[,d]  <- pnorm(pCoef1 + pCoef2 * scale(days)[,1][d] + pCoef3 * (scale(days)[,1][d])^2 + eta.p + eta.pT )
    gamma[,d]<- pnorm(gCoef1 + gCoef2 * scale(days)[,1][d] + gCoef3 * (scale(days)[,1][d])^2 + eta.g + eta.gT )
  }
  boxplot(gamma, col="lightblue",add = TRUE, at = (2:nday)+89+0.25, yaxt="n",xaxt="n",outline = FALSE, border = "lightblue")
  boxplot(phi, col="pink",add = TRUE, at = (2:nday)+89-0.25, yaxt="n",xaxt="n",outline = FALSE, border = "pink")
  
  gamma.sm<-apply(gamma,2,quantile, c(0.025,0.5,0.975))
  phi.sm<-apply(phi,2,quantile, c(0.025,0.5,0.975))
  lines(days[-1],gamma.sm[2,], col="blue", lwd=2)
  lines(days[-1],phi.sm[2,], col="red", lwd=2)
  
  u<-matrix(NA,2000,90)
  psi<-matrix(NA,2000,90)
  psi[,1]<-round(mean(OccVar[1:n,paste("u[",1,",",eg.year,",",eg.site,"]",sep="")])) #max(y[,1,eg.year,eg.site],na.rm=TRUE)
  u[,1]<-rbinom(2000,1,psi[,1])
  for(d in 2:90){
    psi[,d]<-u[,d-1]*phi[,d-1] + (1-u[,d-1])* gamma[,d-1]
    u[,d]<-rbinom(2000,1,psi[,d])
  }
  lines(days[1:90],colMeans(psi), lwd=2)
  
  points(days,apply(y[,,,eg.site],3,colMax)[,eg.year], ylim=c(0,1),pch=19, col="palegreen")
  sim.occ<-round(colMeans(OccVar[1:n,paste("u[",1:91,",",eg.year,",",eg.site,"]",sep="")]))
  points(days,sim.occ)
  mtext(sppList[spp,4],3, font=3, outer=F, line=0, adj=0.5 )
  if(spp == spp10[4]) mtext("Occupancy / Probability",2, font=1, outer=F, las=0, line=2.5, adj=0.5 )
  if(spp == spp10[8]) mtext("Day",1, font=1, outer=F, las=0, line=2.5, adj=0.5 )
}# end spp
#dev.off()

###################################################### 
###Plot Fit of P(detection)
par(mfrow=c(3,3), mar=c(2,2,2,1),oma=c(2,2,0,0), las=1)
spplim=71 # number of study species
i=0
for(spp in spp10){
  i=i+1
  load(paste("Results/", sppList[spp,4],"OccVar.RData"))
  
  n=nrow(OccVar)
  dCoef1<-OccVar[,grepl("dCoef1",colnames(OccVar))]
  dCoef2<-OccVar[,grepl("dCoef2",colnames(OccVar))]
  
  delta<-array(NA,dim=c(n,nyear))
  
  for(t in 1:nyear){
    delta[,t]<- exp(dCoef1[,1] + dCoef2*PropObs[3,t]) 
  }

  curve(1 - median(delta[,1])/(x + median(delta[,1])), from=0, to=spplim, ylim=c(0,1), ylab="")
  dqH<-quantile(delta[,1], 0.975)
  dqL<-quantile(delta[,1], 0.025)
  polygon(c(0:100,100:0), c(1 - dqH/(c(0:100) + dqH),1 - dqL/(c(100:0) + dqL)), col="gray50", border=NA)
  dqH<-quantile(delta[,10], 0.975)
  dqL<-quantile(delta[,10], 0.025)
  polygon(c(0:100,100:0), c(1 - dqH/(c(0:100) + dqH),1 - dqL/(c(100:0) + dqL)), col=ifelse(spp==3,"#CCCCCC95","gray80"), border=NA)
  
  curve(1 - median(delta[,1])/(x + median(delta[,1])), add=TRUE)
  curve(1 - median(delta[,10])/(x + median(delta[,10])), add=TRUE, lty=3)
  
  x0<-20
  y0<-1-median(delta[,1])/(x0 + median(delta[,1]))
  x1<-20
  y1<-1-median(delta[,10])/(x1 + median(delta[,10]))
  arrows(x0,y0,x1,y1, length=0.1, col=2)
  
  mtext(sppList[spp,4],3, font=3, outer=F, line=0, adj=0.5, cex=0.8 )
  if(i==8)mtext("Species List Length (No.)",1, outer=TRUE, font=2, line=0.5, adj=0.5, cex=0.8)
  if(i==4)mtext("Detection Probability",2, outer=TRUE, font=2, adj=0.5, cex=0.8, las=0, line=0.5)
  rug(jitter(as.vector(SPPL[,,1,]), amount=0.2), side=3)
  rug(jitter(as.vector(SPPL[,,10,]),amount=0.2),col=2)
  
} #end of spp
