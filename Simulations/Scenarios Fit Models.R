#########################################
path<-"C:/temp/"
require(rjags)
require(dclone)
set.seed(1234)

load("PropObs.rdata")
years<-2005:2014
years<-scale(years)[,1]
nyear<-length(years)
nsite<-5
days<-90:180
days<-scale(days)[,1]
nday<-length(days)
nscn<-9 #H,M,L occupancy; +,- occ trends; +,- no. obs; +,- trends in detection
nit<-10000  #iterations

## Run occupancy model for simulated data
for(scn in 1:9){
  load(paste("Simulated Obs Scn",scn,".RData", sep=""))

  uinit<-array(rbinom(nday*nyear*nsite,1,1), dim=c(nday,nyear,nsite))
  pC0<-rep(0,3)
  gC0<-rep(0,3)
  PC0<-diag(0.001,3,3)
  GC0<-diag(0.001,3,3)

  ########## Fit the model in JAGS
  mydata<-list(nsite=nsite,nday=nday,nyear=nyear, nrep=Objs$nrep, TDVar=PropObs[3,], #scale(years)[,1]
               years=years,days=days, SPPL=Objs$SPPLsim[,,,], y=Objs$simOcc[,,,],
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
                    "psiD","psi.fs","zM","u"),
                    model = modelfile, inits=inits, n.chains=4,
                    n.update=5000,n.iter=15000,thin=30)
  stopCluster(cl)
  OccVar<-as.matrix(Occmodel, chains=TRUE)
  save(OccVar=OccVar, file=paste("OccParamScn",scn,".RData",sep=""))
}