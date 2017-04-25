###Jags complie Data for Species
sppNames<-as.character(read.csv("StudySpecies.csv")[,4])

## the DATA per Wetland
WS<-read.csv("Uppland90180.csv")[,c(2,5,10,11:13)] #only those columns needed
WS$Julian<-as.character(WS$Julian)
WS$Week<-as.character(WS$Week)
WS<-as.matrix(WS)
gc()
####################################
### Compile data in the needed format
####################################
nCL<-10 #complete list parameter
nSL<-5 #short list parameter
nsite<-length(unique(WS[,"WetlandSite"]))
sites<-unique(WS[,"WetlandSite"])

years<-seq(2005,2014)
nyear<-length(years)
days<-seq(90,180)
nday<-length(days)

# Maximum number of visits per day
max.nrep<-40 

#### Create empty arrays
# the real number of visits per day
nrepMAX<-array(0, dim=c(nday,nyear,nsite), dimnames=list(days,years,sites))
# the index for loops; that is at least 1 if that day there ware no visits
nrepIND<-array(1, dim=c(nday,nyear,nsite), dimnames=list(days,years,sites))
SPPL<-array(0, dim=c(max.nrep,nday,nyear,nsite), dimnames=list(c(1:max.nrep),days,years,sites)) 
OccData<-array(NA, dim=c(length(sppNames),max.nrep,nday,nyear,nsite), dimnames=list(sppNames,c(1:max.nrep),days,years,sites))

for(i in sites){
  for(t in years){
    for(d in days){
      wDYS<-which(WS[,"WetlandSite"] == i & WS[,"Year"] == t & WS[,"Julian"] == d)
      D<-as.character(d);T<-as.character(t); I<-as.character(i)
      if(length(wDYS)){
        # temporal indices
        ## all observers names
        obstemp<-unique(WS[,"RecordedBy"][wDYS])
        nrep<-length(obstemp)
        sppUL<-array(NA, dim=c(length(sppNames),nrep), dimnames=list(c(1:length(sppNames)),c(1:nrep)))  ## Unique species observed
        sppLLength<-numeric(nrep)
        for(j in 1:nrep){
          # a visit
          wRDYS<-which(WS[,"WetlandSite"] == i & WS[,"Year"] == t & WS[,"Julian"] == d & WS[,"RecordedBy"] == obstemp[j])
          sppLLength[j]<-length(as.character(unique(WS[,"ScientificName"][wRDYS])))
          sppUL[1:sppLLength[j],j]<-as.character(unique(WS[,"ScientificName"][wRDYS])) # unique species
        }#end j nrep
  
        ## if there are too many replicates I only take the first max.nrep longest species list. This is to avoid acctually one extreme day with 100 rep
        nrepMAX[D,T,I]<-min(nrep,max.nrep)
        nrepIND[D,T,I]<-min(nrep,max.nrep)
        wReps<-order(sppLLength,decreasing=TRUE)[1:nrepMAX[D,T,I]]
        SPPL[1:nrepMAX[D,T,I],D,T,I]<-sppLLength[wReps]
  
        ## for each spp in the list put 1 or zero if it was present on each replicate
        for(spp in sppNames){
          occV<-numeric()
          for(j in 1:nrepMAX[D,T,I]){
            occV[j]<-ifelse(spp %in% sppUL[,wReps[j]],1,0) ## Order observations from longer to shorter Spp List
          }#end j wReps
          OccData[spp,1:nrepMAX[D,T,I],D,T,I]<-occV ##
        } # end spp
      }# end if wDYS
    }# end d
    print(t)
  }# end t
  print(i)
}# end i

SL<-ifelse(SPPL >= nSL & SPPL < nCL,1,0)
CL<-ifelse(SPPL >= nCL,1,0)

rm(WS)
gc()

x<-c("sites","nsite","days","nday","years","nyear","nrepMAX", "nrepIND", "SPPL", "SL", "CL", "OccData")
save(list=x, file="OccDataUppland.RData")

