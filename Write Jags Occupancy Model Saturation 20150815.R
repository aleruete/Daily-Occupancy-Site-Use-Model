### Write occupancy model
require(R2WinBUGS)
path<-"C:/temp/"

#d: day; t: years 2005-2014; i: site/locality; j: observation/replicate
modelstring = "
model{
## Observation Model
for(i in 1:nsite){
 for(t in 1:nyear){
   for(d in 1:nday){
    for(j in 1:nrep[d,t,i]){ #this is nrepIND
      y[j,d,t,i] ~ dbern(Py[j,d,t,i])
      Py[j,d,t,i]<- u[d,t,i]*p[j,d,t,i]
      p[j,d,t,i]<- 1- delta[t,i]/ (SPPL[j,d,t,i] + delta[t,i]) # Species List Length saturation curve Half-ignorance
    }# end nrep
  }# end nday
  zM[t,i]<-mean(u[2:nday,t,i])
  
  ## Detectability coeficients
  log(delta[t,i]) <-dCoef1[i] + dCoef2 * TDVar[t] #Time Dedendant Variable = PLL[t]
 }# end nyear
}# end nsite

## Process model
for(i in 1:nsite){
  for(t in 1:nyear){
    u[1,t,i]~dbern(psiD[i])
    for(d in 2:nday){
      u[d,t,i]~dbern(muU[d,t,i])
      muU[d,t,i]<- u[d-1,t,i]*phi[d-1,t,i] + (1-u[d-1,t,i]) * gamma[d-1,t,i]
      probit(phi[d-1,t,i])   <- pCoef[1]  + pCoef[2] * days[d-1] + pCoef[3] * days[d-1] * days[d-1]
                                          + eta.p[i] + eta.pT[t]
      probit(gamma[d-1,t,i]) <- gCoef[1]  + gCoef[2] * days[d-1] + gCoef[3] * days[d-1] * days[d-1]
                                          + eta.g[i] + eta.gT[t]
   }# end nday
  }# end nyear
}# end nsite


## Summaries
for(t in 1:nyear){
  for(d in 1:nday){
    psi.fs[d,t]<-sum(u[d,t,1:nsite])/nsite
  }
}

#########
## Priors
for(i in 1:nsite){
  psiD[i]~dunif(0,1)
  dCoef1[i]~dnorm(0, 0.001)
  eta.p[i]~dnorm(0, p.tau)
  eta.g[i]~dnorm(0, g.tau)
}
for(t in 1:nyear){
  eta.pT[t]~dnorm(0, pt.tau)
  eta.gT[t]~dnorm(0, gt.tau)
}

dCoef2~dnorm(0, 0.001)

pCoef[1:3]~dmnorm(pC0[], PC0[,])
gCoef[1:3]~dmnorm(gC0[], GC0[,])

p.tau~dgamma(0.001, 0.001)
g.tau~dgamma(0.001, 0.001)
pt.tau~dgamma(0.001, 0.001)
gt.tau~dgamma(0.001, 0.001)


} #end of model
" # close quote for modelstring
writeLines(modelstring,con=paste(path,"WetlandOccSaturation.txt",sep = ""))

