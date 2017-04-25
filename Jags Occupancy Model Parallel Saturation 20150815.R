#########################  Compile Model #######################################
path<-"C:/temp/"
require(rjags)
require(dclone)
set.seed(1234)
require(coda)
#require(mcmcplots)

require(foreach)
require(doParallel)
#source("Z:/R generic scripts/Estimate Mode.R")
invlogit<-function(x){(p<-1/(1+exp(-x)))}
# colMax <- function(data){
#        temp<-suppressWarnings( apply(data, 2, max, na.rm = TRUE) )
#        v<-ifelse(temp=="-Inf",NA,temp)
#        return(v)
# }
# sumF<-function(x){ ## Function for discrepancy measures
#    wNNA<-which(!is.na(x))
#    fit<-sum(as.vector(x[wNNA]))
# #   fit<-sum(as.vector(x),rm.na=T)
#    return(fit)
# }

#################################################################################
#### Occupancy dynamic models
sppList<-read.csv("StudySpecies.csv")
spp10<-c(1:77)

for(spp in spp10){
load("OccDataUppland.RData")
  y<-OccData[spp,,,,]
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

rm(OccData); gc()

uinit<-array(rbinom(nday*nyear*nsite,1,1), dim=c(nday,nyear,nsite))
pC0<-rep(0,3)
gC0<-rep(0,3)
PC0<-diag(0.001,3,3)
GC0<-diag(0.001,3,3)

########## Fit the model in JAGS
mydata<-list(nsite=nsite,nday=nday,nyear=nyear, nrep=nrep, TDVar=PropObs[3,],
             years=scale(years)[,1],days=scale(days)[,1], SPPL=SPPL, y=y,
             pC0=pC0, gC0=gC0, PC0=PC0, GC0=GC0) #
inits<-list(u=uinit)
modelfile<- paste(path,"WetlandOccSaturation.txt",sep = "")


## set cluster
cl <- makeCluster(4, type = "SOCK")
tmp <- clusterEvalQ(cl, library(dclone))
parLoadModule(cl, "glm")
## initiate, update and run JAGS
Occmodel<-jags.parfit(cl, dat = mydata, params=c("dCoef1","dCoef2",
                  "pCoef", "gCoef", "p.tau","g.tau","pt.tau","gt.tau",
                  "eta.g", "eta.p", "eta.gT", "eta.pT",
                  "zM","psiD","psi.fs","u"),
                  model = modelfile, inits=inits, n.chains=4,
                  n.update=5000,n.iter=15000,thin=30)
stopCluster(cl)

## Get posteriors
OccVar<-as.matrix(Occmodel, chains=TRUE)
save(OccVar, file=paste0("Results/", sppList[spp,4]," OccVar.RData"))
rm(Occmodel, OccVar)
gc()
print(paste("I'm done with", as.character(sppList[spp,4])))
}### end spp



