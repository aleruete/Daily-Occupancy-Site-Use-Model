require(foreach)
require(doParallel)
colMax <- function(data){
       v<-suppressWarnings( apply(data, 2, max, na.rm = TRUE) )
       v[v=="-Inf"]<-NA #ifelse(v=="-Inf",NA,v)
       return(v)
}

sumF<-function(x){ ## Function for discrepancy measures
  wNNA<-which(!is.na(x))
  fit<-sum(as.vector(x[wNNA]))
  return(fit)
}

## Prob at least one
pALO<-function(pdet,n) return(1-(1-pdet)^n)
pALOv<-function(pdet) return(1-prod(1-pdet, na.rm = TRUE))

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

  load(paste0("Results/", sppList[spp,4]," OccVar.RData"))
  SPPL<-ifelse(SPPL[,,,]==0, NA,SPPL[,,,])
  SPPL[1,,,]<-ifelse(is.na(SPPL[1,,,]), 0,SPPL[1,,,])
  
  n=nrow(OccVar)
  ns<- 1:n
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
  ydS<-array(NA, dim=c(nyear, nsite)) # observed site means
  sq<-array(NA,dim=c(n,nyear,nsite))
  sq.new<-array(NA,dim=c(n,nyear,nsite))
  predictedS<-array(NA,dim=c(n,nyear,nsite))
  pti<-array(NA,dim=c(n,nyear,nsite))
  
  ## Bayesian p-values
  ## Assess model fit using a sums-of-squares-type discrepancy
  for(i in 1:nsite){
    yM[,i]<-colMeans(apply(y[,,,i],3,colMax), na.rm=TRUE)
    for(t in 1:nyear){
      y.new<-array(NA,dim=c(n,40,nday))
      y.newMax<-array(NA,dim=c(n,nday))
      ptid<-array(NA,dim=c(n,nday))
      predicted<-array(NA,dim=c(n,nday))
      ydM<-numeric()
      for(d in 1:nday){
        u<-OccVar[1:n,paste("u[",d,",",t,",",i,"]",sep="")]
        ydM[d]<-suppressWarnings(max(y[,d,t,i], na.rm=TRUE)) # max detection per day
        if(ydM[d]=="-Inf")ydM[d]<-NA
        if(exists("p")) rm(p)
        if(nrepMAX[d,t,i]>0){
          p<-matrix(NA, n, nrepMAX[d,t,i])
          for(j in 1:nrepMAX[d,t,i]){ #nrepMax to avoid NAs
            p[,j]<-1 - delta[,t,i]/(SPPL[j,d,t,i] + delta[,t,i])
            
            ## Generate new data
            y.new[,j,d] <- rbinom(n,1,u*p[,j]) #rbinom(n,1,predicted[,j,d])
          }# end nrep
          ptid[,d]<-apply(p,1,pALOv)
          predicted[,d]<-ptid[,d]*u 
        } # end if
        y.newMax[,d] <- suppressWarnings(apply(y.new[,,d],1,max, na.rm=TRUE))
        
      }# end nday
      y.newMax[y.newMax=="-Inf"]<-NA 
      y.newS[,t,i]<-rowMeans(y.newMax, na.rm=TRUE)
      ydS[t,i]<-mean(ydM, na.rm=TRUE)
      predictedS[,t,i]<-rowMeans(predicted[,], na.rm=TRUE)
      pti[,t,i]<-rowMeans(ptid[,], na.rm=TRUE)
      
      #residuals
      sq[,t,i]<-((yM[t,i] - predictedS[,t,i])/sqrt(predictedS[,t,i]*(1-pti[,t,i])))^2 
      sq.new[,t,i]<-((y.newS[,t,i] - predictedS[,t,i])/sqrt(predictedS[,t,i]*(1-pti[,t,i])))^2
    }# end nyear
    print(i)
    y.newM[1:3,,i]<-apply(y.newS[,,i],2,quantile,probs=c(0.025,0.5,0.975))
    y.newM[4,,i]  <-apply(y.newS[,,i],2,sample,1)
  }#)# end nsite
  
  gc()

  site.fit<-matrix(NA,n,nsite)
  site.fit.new<-matrix(NA,n,nsite)
  site.bpv<-numeric()
  for(i in 1:nsite){
    site.fit[,i]<-apply(sq[,,i],1,sum, na.rm=TRUE)
    site.fit.new[,i]<-apply(sq.new[,,i],1,sum, na.rm=TRUE)
    site.bpv[i]<-mean(site.fit.new[,i]>site.fit[,i])
  }

  fit<-apply(sq[,,],1,sum, na.rm=TRUE)
  fit.new<-apply(sq.new[,,],1,sum, na.rm=TRUE)
  bpvalue <- mean(fit.new>fit)			# Bayesian p-value; should be close to 0.5

  gc()
  ################################################################################
  #plot psi years
  pdf(paste0("Results/", sppList[spp,4]," Occupancy year.pdf"),14,9) #windows(14,9)
  par(mfrow=c(3,4), mar=c(3,4,1,1),oma=c(0,0,3,0))
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
  
  exp1<-expression(Z[t]==bar(over(sum(u[d,t,i], n, i == 1), n)) )
  exp2<-expression(z[t,i]==sum(over(u[w,t,i], n), d == 1, n), "Observed","Simulated")
  
  plot(years,psiM[2,], type="n", xlab="Year", ylab="Mean Site Use", ylim=c(0,1), xaxp=c(years[1],years[nyear],nyear-1))
  mtext(sppList[spp,4],3, font=3, outer=TRUE, adj=0.5)
  polygon( c(years,rev(years)),c(apply(psifs[1,,],2,mean),rev(apply(psifs[3,,],2,mean))), col="gray80", border=0)
  lines(years,apply(psifs[2,,],2,mean),lty=1)

  points(years,rowMeans(y.newM[2,,]), pch=20, col=2)
  arrows(years,rowMeans(y.newM[1,,]),years,rowMeans(y.newM[3,,]), length = 0.05, angle = 90, code = 3, col=2)
  points(years,rowMeans(yM), pch=20)
  title(main="Uppland wetlands")
  legend("topright",legend=c("Estimated","Observed","Replicated"), lty=c(1,NA,NA), pch=c(NA,20,20),col=c(1,1,2),bty="n")
  text(2005,1,paste("bpv =",round(bpvalue,2)), pos=4)

  #Plot occurrence per site
  for(i in 1:nsite){
       plot(years,zM[2,,i], type="n", xlab="Year", ylab="Site Use", ylim=c(0,1), xaxp=c(years[1],years[nyear],nyear-1))
       polygon( c(years,rev(years)),c(zM[1,,i],rev(zM[3,,i])), col="gray80", border=0)
       lines(years, zM[2,,i])
       points(years, y.newM[2,,i], pch=20, col=2) #, col=2
       arrows(years, y.newM[1,,i], years, y.newM[3,,i], length = 0.05, angle = 90, code = 3, col=2) #, col=2
       points(years, yM[,i], pch=20) #, col=2
       title(main=sites[i])
       text(2005,1,paste("bpv =",round(site.bpv[i],3)), pos=4)
  }
  ## Data vs predictions
  par(mar=c(4,4,1,1))
  plot(yM,y.newM[4,,], col="red",xlab="Observed mean site use", ylab="Replicated mean site use")
  abline(0,1)

  ## BPV  
  par(mar=c(4,4,1,1))
  plot(fit,fit.new, xlab="SSQ discrepancy fo actual data set", ylab="SSQ discrepancy for replicated data set",
                    #xlim=c(0,max(fit)), ylim=c(0,max(fit.new)))
                    xlim=c(0,max(site.fit)), ylim=c(0,max(site.fit.new)), type="n")
  abline(0,1)
  for(i in 1:nsite){
    points(site.fit[,i],site.fit.new[,i], col=i+1)
    text(mean(site.fit[,i]),mean(site.fit.new[,i]),round(site.bpv[i],3))
  }
  
  dev.off()
rm(OccVar)
print(as.character(sppList[spp,4]))
gc()

################################################################################
  #Plot Psi days
  pdf(paste0("Results/", sppList[spp,4]," Occupancy Season.pdf"),14,9) #windows(14,9)
  par(mfrow=c(3,4), mar=c(3,4,1,1), oma=c(0,0,3,0))
  for(t in 1:nyear){
    plot(days,psifs[2,,t], type="n", xlab="Day", ylab="Occupancy", ylim=c(0,1))
    if(t==1) mtext(sppList[spp,4],3, font=3, outer=TRUE, adj=0.5 )
    polygon( c(days,rev(days)),c(psifs[1,,t],rev(psifs[3,,t])), col="gray80", border=0)
    lines(days, psifs[2,,t])
    ###plot(1:nday, colMeans(u), type="l", ylim=c(0,1))
    OccSite<-matrix(NA,nday,nsite)
    for(i in 1:nsite){
          OccSite[,i]<-colMax(y[,,t,i])
    }
    points(days,rowMeans(OccSite[1:nday,],na.rm=TRUE),pch=20)
    title(main=years[t])
  }
  exp3<-expression(psi[t]==sum(over(u[w,t,i], n), i == 1, n), bar(plain(Obs)))
  legend("bottomleft",legend=exp3, #"Yearly Occupancy (mean z[t,1:nsite])"
                       lty=c(1,NA), pch=c(NA,20), bty="n")
  dev.off()
  
  ################################################################################
  #Plot Probability of detection
  pdf(paste0("Results/", sppList[spp,4]," p(det).pdf"),14,9) #windows(14,9)
  par(mfrow=c(3,4), mar=c(4,4,1,1), oma=c(0,0,3,0))
  n=nrow(OccVar)
  dCoef1<-OccVar[1:n,grepl("dCoef1",colnames(OccVar))]
  dCoef2<-OccVar[1:n,grepl("dCoef2",colnames(OccVar))]
  
  delta<-array(NA,dim=c(n,nyear,nsite),dimnames=list(c(1:n), years, sites))
  
  for(i in 1:nsite){
    for(t in 1:nyear){
         delta[,t,i]<- exp(dCoef1[,i] + dCoef2 * PropObs[3,t]) # scale(years)[t,1]# + eta.a[i]
         #delta[,t,i]<-rpois(n,d.mu[,t,i])
    }
  }
  
  for(i in 1:nsite){
     plot(1, 1,xlab="Species list length", ylab="p(Detection)",xlim=c(0,100), ylim=c(0,1), type="n")
     title(main=sites[i])
     if(i==1) mtext(sppList[spp,4],3, font=3, outer=TRUE, adj=0.5 )
     deltaQ<-quantile(delta[,1,i],probs=c(0.025, 0.975))
     polygon(c(0:100,100:0),c(1 - deltaQ[1]/(c(0:100) + deltaQ[1]),1 - deltaQ[2]/(c(100:0) + deltaQ[2])), col="gray90", border=0)
     deltaQ<-quantile(delta[,10,i],probs=c(0.025, 0.975))
     polygon(c(0:100,100:0),c(1 - deltaQ[1]/(c(0:100) + deltaQ[1]),1 - deltaQ[2]/(c(100:0) + deltaQ[2])), col="gray80", border=0)
     curve(1 - median(delta[,1,i])/(x + median(delta[,1,i])), lwd=2, add=TRUE)
     curve(1 - median(delta[,10,i])/(x + median(delta[,10,i])), lwd=1, add=TRUE)
     legend("topleft", title="No. Long Lists / No. Obs", legend=c("0.06", "0.25"), lwd=c(2,1),
                    fill=c("gray90","gray80"), border=0, bty="n")
  }
  dev.off()
  ################################################################################
  pdf(paste0("Results/", sppList[spp,4]," ColExt rate.pdf"),14,9) #windows(14,9)
  #plot colonization and extinction rates
  pCoef1<-OccVar[,grepl("pCoef",colnames(OccVar))][,1]
  pCoef2<-OccVar[,grepl("pCoef",colnames(OccVar))][,2]
  pCoef3<-OccVar[,grepl("pCoef",colnames(OccVar))][,3]
  #pCoef4<-OccVar[,grepl("pCoef",colnames(OccVar))][,4]
  p.tau <-OccVar[,grepl("p.tau",colnames(OccVar))]
  p.sigma<-sqrt(1  / p.tau)
  eta.p<-rnorm(nrow(OccVar), 0, p.sigma)
  pt.tau <-OccVar[,grepl("pt.tau",colnames(OccVar))]
  pt.sigma<-sqrt(1  / pt.tau)
  eta.pT<-rnorm(nrow(OccVar), 0, pt.sigma)
  
  gCoef1<-OccVar[,grepl("gCoef",colnames(OccVar))][,1]
  gCoef2<-OccVar[,grepl("gCoef",colnames(OccVar))][,2]
  gCoef3<-OccVar[,grepl("gCoef",colnames(OccVar))][,3]
  #gCoef4<-OccVar[,grepl("gCoef",colnames(OccVar))][,4]
  g.tau <-OccVar[,grepl("g.tau",colnames(OccVar))]
  g.sigma<-sqrt(1  / g.tau)
  eta.g<-rnorm(nrow(OccVar), 0, g.sigma)
  gt.tau <-OccVar[,grepl("gt.tau",colnames(OccVar))]
  gt.sigma<-sqrt(1  / gt.tau)
  eta.gT<-rnorm(nrow(OccVar), 0, gt.sigma)
  
  t=1
  plot(days[c(1,nday)], c(0,1), type="n", ylab="Colonization (Red) /Persistance rate (Blue)", xlab="Day")
  if(t==1) mtext(sppList[spp,4],3, font=3, outer=TRUE, adj=0.5 )
  phi<-matrix(NA, nrow=nrow(OccVar),ncol=nday-1)
  gamma<-matrix(NA, nrow=nrow(OccVar),ncol=nday-1)
  for(d in 1:(nday-1)){
    phi[,d]  <- pnorm(pCoef1 + pCoef2 * scale(days)[,1][d] + pCoef3 * (scale(days)[,1][d])^2 )
    gamma[,d]<- pnorm(gCoef1 + gCoef2 * scale(days)[,1][d] + gCoef3 * (scale(days)[,1][d])^2 )
  }
  boxplot(gamma, col="lightblue",add = TRUE, at = (1:(nday-1))+89+0.25, yaxt="n",xaxt="n", outline = FALSE, border = "lightblue")
  boxplot(phi, col="pink",add = TRUE, at = (1:(nday-1))+89-0.25, yaxt="n",xaxt="n", outline = FALSE, border = "pink")
  
  gamma.sm<-apply(gamma,2,quantile, c(0.025,0.5,0.975))
  phi.sm<-apply(phi,2,quantile, c(0.025,0.5,0.975))
  lines(days[-1],gamma.sm[2,], col="blue", lwd=2)
  lines(days[-1],phi.sm[2,], col="red", lwd=2)
  
  dev.off()
  ################################################################################
  rm(OccVar)
  print(paste("I'm done with", as.character(sppList[spp,4])))
  gc()
}### end spp



