##Simulate data according to different assumptions
set.seed(1234)

######## Simulate Real Occupancy data
years<-2005:2014
years<-scale(years)[,1]
nyear<-length(years)
nsite<-5
days<-90:180
days<-scale(days)[,1]
nday<-length(days)
nscn<-9 #H,M,L occupancy; +,- occ trends; +,- no. obs; +,- trends in detection
nit<-10000  #iterations

## Occupancy levels
occ<-matrix(NA,2,3)
occ[,1]<-c(0.7,0.9) # HIGH
occ[,2]<-c(0.3,0.6) # MEDIUM
occ[,3]<-c(0.05,0.2) # LOW

psi<-matrix(NA,nsite,nscn)
for(scn in 1:3){
   psi[,scn]<-round(runif(nsite, occ[1,scn], occ[2,scn]),3) #initial occupancy probability
}   ## end scn
psi[,4:9]<-psi[,2]

### Parameters for linear changes in col/ext coefficients
pCf1<- c(2.13,2.3,-0.4,rep(2.3,6))
pCf2<- c(0,0,0,-1,1,0,0,0,0) #change in dayly persistance over years (for non linear reasons + and - coef have to be inversed)
pCf3<- c(-0.34,1.3,0.67,rep(1.3,6))  #1.12
pCf4<- c(-0.3,-1.26,-0.64,rep(-1.26,6))
eta.p<-rnorm(nsite, 0, 0.7)
gCf1<- c(-0.21,1.4,-2.1,rep(1.4,6))
gCf2<- c(0,0,0,-1,1,0,0,0,0) #change in colonization over years
gCf3<- c(-0.03,2.1,0.4,rep(2.1,6))
gCf4<- c(0.01,-2.55,0.04,rep(-2.55,6))
eta.g<-rnorm(nsite, 0, 0.7) #

U<-array(NA,dim=c(nday,nyear,nsite,nscn,nit))
muU<-array(NA,dim=c(nday,nyear,nsite,nscn,nit))
phi<-array(NA,dim=c(nday,nyear,nsite,nscn))
gamma<-array(NA,dim=c(nday,nyear,nsite,nscn))

for(scn in 1:9){
for(i in 1:nsite){
  for(t in 1:nyear){
    for(d in 1:nday){
        phi[d,t,i,scn]  <-pnorm(pCf1[scn] + pCf2[scn] * years[t] + pCf3[scn] * days[d] + pCf4[scn] * (days[d])^2 + eta.p[i])
        gamma[d,t,i,scn]<-pnorm(gCf1[scn] + gCf2[scn] * years[t] + gCf3[scn] * days[d] + gCf4[scn] * (days[d])^2 + eta.g[i])
    }
     muU[1,t,i,scn,]<-0 #psi[i,scn]
     U[1,t,i,scn,]<-rbinom(nit,1,muU[1,t,i,scn,])
     for(d in 2:nday){
        muU[d,t,i,scn,]<-U[d-1,t,i,scn,]*phi[d-1,t,i,scn] + (1-U[d-1,t,i,scn,])*gamma[d-1,t,i,scn]
        U[d,t,i,scn,]<-rbinom(nit,1,muU[d,t,i,scn,])
     }# end nday
  }# end nyear
}# end nsite
} # end scn
scn.lbls<-expression("1: High Site Use","2: Medium Site Use","3: Low Site Use", "4: " %up% "Site Use", "5: " %down%  "Site Use",
                     "6: " %up%  "no. Daily Visits", "7: " %down% "no. Daily Visits", "8: " %up% "Detectability", "9: " %down% "Detectability")

pdf("Col-Ext.pdf",7,7) #windows(14,9)
par(mfrow=c(3,3), mar=c(4,4,3,0.5), las=1)
for(s in 1:9){
  curve(pnorm(pCf1[s]+pCf2[s]*0   +pCf3[s]*x+pCf4[s]*x^2), from=-1.7,to=1.7, col="red", lwd=2, ylim=c(0,1), ylab="Probability", xlab="Julian day", xaxt="n")
  axis(1,at= days[seq(1,91,10)],labels=c(90:180)[seq(1,91,10)], las=2)
  if(s%in%c(4,5))curve(pnorm(pCf1[s]+pCf2[s]*1.5 +pCf3[s]*x+pCf4[s]*x^2), col="red", lwd=1, add=T)
  if(s%in%c(4,5))curve(pnorm(pCf1[s]+pCf2[s]*-1.5+pCf3[s]*x+pCf4[s]*x^2), col="red", lwd=3, add=T)
  
  curve(pnorm(gCf1[s]+gCf2[s]*0   +gCf3[s]*x+gCf4[s]*x^2), add=T, col="blue", lwd=2)
  if(s%in%c(4,5))curve(pnorm(gCf1[s]+gCf2[s]*1.5 +gCf3[s]*x+gCf4[s]*x^2), add=T, col="blue", lwd=1)
  if(s%in%c(4,5))curve(pnorm(gCf1[s]+gCf2[s]*-1.5+gCf3[s]*x+gCf4[s]*x^2), add=T, col="blue", lwd=3)
  mtext(scn.lbls[s], 3, line=0.5)
  if(s==1)legend("bottomleft", c("Persistance","Colonization"), lwd=2,col=c("red", "blue"), bty="n", cex=0.8)
  if(s==1)legend("bottomright", c("Year 0","Year 5", "Year 10"), lwd=c(1,2,3), bty="n", cex=0.8)
  points(days, U[,1,2,s,999],pch=19, cex=0.8)
}
dev.off()

save(U, file="Simulated State.Rdata")

par(mfrow=c(2,3))
for(i in 1:nsite){
   plot(1:nyear, colMeans(U[,,i,1,999]), type="l", ylim=c(0,1))
   for(scn in 2:nscn)lines(1:nyear,colMeans(U[,,i,scn,999]), col=scn)
}

############ Simulate Observation Process

nrepInt<-c(rep(5,5),0,30,5,5)
nrepSlp<-c(0,0,0,0,0,3,-3,0,0) ## Fixed increment in number of replicates over time
maxrep<-50 #maximum number of replicates to keep it limited
eta.a<-rpois(nsite, 2) #site variability for number of samples

deltaInt<-runif(nsite,3,4) #3
deltaTrd<-c(rep(0,7),-0.15,0.15)  ## Fixed Trend in Detectability
deltaPrOb<- -3 ## increasing detection, LOWERING delta

## It has been observed that number of Comphrehensive lists (Long Lists) increased over time since Svalan was open
load("PropObs.rdata")

for(scn in 1:9){
deltasat<-array(NA,dim=c(nyear,nsite))
nrep<-array(NA,dim=c(nday,nyear,nsite))
simOcc<-array(NA,dim=c(maxrep,nday,nyear,nsite))
SPPLsim<-array(0,dim=c(maxrep,nday,nyear,nsite))
p<-array(NA,dim=c(maxrep,nday,nyear,nsite))
for(i in 1:nsite){
  for(t in 1:nyear){
    if(scn!=8)deltasat[t,i]<-exp(deltaInt[i] + deltaPrOb * PropObs[3,t] + deltaTrd[scn]*(t-1))
    if(scn==8)deltasat[t,i]<-exp(deltaInt[i]+0.2 + deltaPrOb * PropObs[3,t] + deltaTrd[scn]*(t-1))
    for(d in 1:nday){
     ## number of replicates increasing over time
     nrep[d,t,i]<-min(max(1,rpois(1,nrepInt[scn] + nrepSlp[scn]*t + eta.a[i])),maxrep) # at least 1, max 50

     ## Number of samples that are either single records, short lists or comprehensive list
     n10<-trunc(nrep[d,t,i]*PropObs[3,t])
     n5<-trunc(nrep[d,t,i]*PropObs[2,t])
     n1<-nrep[d,t,i]-n10-n5    ## This means that if the number of observations is lower than 3, there will only be single records
     SPPLsim[1:nrep[d,t,i],d,t,i]<-c(sample(c(1:4),n1, replace = T), sample(c(5:9),n5, replace=T), sample(c(10:45),n10, replace=T)) #45 is the longest observed spplist c(10:78),

     # assuming saturation
     p[1:nrep[d,t,i],d,t,i]<- 1 - deltasat[t,i]/ (SPPLsim[1:nrep[d,t,i],d,t,i] + deltasat[t,i])
     simOcc[1:nrep[d,t,i],d,t,i]<-U[d,t,i,scn,999] * rbinom(nrep[d,t,i],1,p[1:nrep[d,t,i],d,t,i])  #rbinom DID observe or not
   }# end nday
 }# end nyear
}# end nsite

# pdf(paste("Nobs Scn",scn,".pdf"),7,7) #windows(14,9)
# par(mfrow=c(3,2), mar=c(3,4,3,1), las=1)
# for(i in 1:nsite){
#   boxplot(nrep[,,i], type="l", ylim=c(0,50), xlab="Year", ylab="No. daily replicates")
#   mtext(paste("Site",i),3,line=0.5)
# }
# dev.off()

# pdf(paste("Detect50 Scn",scn,".pdf"),7,7) #windows(14,9)
# par(mfrow=c(3,2), mar=c(3,4,3,1), las=1)
# for(i in 1:nsite){
#   plot(1:nyear, deltasat[,i], type="l", ylim=c(0,70), xlab="Year", ylab="O50% (List Length for p(det)=0.5)")
#   mtext(paste("Site",i),3,line=0.5)
# }
# dev.off()

Objs<-list(days=days,deltasat=deltasat, nday=nday,nrep=nrep,nscn=nscn,nsite=nsite,
            nyear=nyear,p=p,simOcc=simOcc,SPPLsim=SPPLsim, years=years)
save(Objs, file=paste("Simulated Obs Scn",scn,".RData", sep=""))
} # end scn

ParamsSim<-list(deltaInt=deltaInt,deltaPrOb=deltaPrOb,deltaTrd=deltaTrd,
eta.a=eta.a,eta.g=eta.g,eta.p=eta.p,
gCf1=gCf1,gCf2=gCf2,gCf3=gCf3,gCf4=gCf4,
pCf1=pCf1,pCf2=pCf2,pCf3=pCf3,pCf4=pCf4,
nrepInt=nrepInt,nrepSlp=nrepSlp,psi=psi)
save(ParamsSim, file="Simulation Parameters.RData")
