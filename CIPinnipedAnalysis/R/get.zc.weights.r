

#' Extracts weights for Zc pups
#' Extracts weights for Zc pups from various tables in the BrandMaster.mdb
#' ACCESS database and returns a dataframe with ancillary fields added.
#' 
#' 
#' @param fdir directory for BrandMaster.mdb
#' @param ENYears values of years that will be flagged as ENSO years
#' @return dataframe of weights for Zc
#' @export
#' @author Jeff Laake
#' @examples
#' 
#' # fdir should be set to location of  BrandMaster.mdb
#' # edir should be set to location of environmental.data.mdb
#' # each defaults below to current working directory
#' # Use "J:/" if it should be Calcurr/Databases
#' # edir="J:/"
#' # fdir="J:/"
#' fdir="." 
#' edir="."
#' plot.weight.series=function (time, predictions, ...) 
#' {
#'	plotCI(time, predictions$fit, predictions$fit - 1.96 * predictions$se, 
#'			predictions$fit + 1.96 * predictions$se, xlab = "Year", 
#'			ylab = "Predicted average weight (kg) 1 Oct", 
#'			...)
#'	lines(time, predictions$fit)
#'	invisible()
#' }
#' 
#' #################################################################################
#' # Cross-sectional analysis
#' #################################################################################
#' # get zc weight values from database
#' zcweights=get.zc.weights(fdir=fdir)
#' #
#' #  exclude brand eval weights and captures in April
#' #
#' zcweights.all=zcweights[zcweights$days>-30&zcweights$days<150&
#'    zcweights$sitecode=="SMI"&zcweights$type=="Initial",]
#' #
#' # compute observed averages - using data within 1 Sept to 15 Nov
#' #
#' zcweights=zcweights.all[zcweights.all$days>=-31&zcweights.all$days<=45,]
#' female.observed=sapply(split(zcweights$weight[zcweights$sex=="F"],
#'  factor(zcweights$cohort[zcweights$sex=="F"],levels=levels(factor(zcweights.all$cohort)))),mean)
#' male.observed=sapply(split(zcweights$weight[zcweights$sex=="M"],
#'  factor(zcweights$cohort[zcweights$sex=="M"],levels=levels(factor(zcweights.all$cohort)))),mean)
#' zcweights=zcweights.all
#' zcweights$batch=factor(paste(zcweights$cohort,zcweights$days))
#' #
#' #  fit growth model
#' #
#' require(nlme)
#' zc.weight.model=lme(weight~sex*days,random=list(~days|cohort,~1|batch),data=zcweights)
#' print(summary(zc.weight.model))
#' # define bootstrap function to compute std error for predicted sex-cohort means
#' bootstrap.se=function(x,nreps)
#' {
#'    pmat=matrix(0,nrow=nreps,ncol=nrow(unique(data.frame(cohort=x$cohort,sex=x$sex))))
#'    i=0
#'    while(i<nreps)
#'    {
#'      xsamp=lapply(split(x,list(x$batch,x$sex)),function(x) if(nrow(x)>0) x[sample(1:nrow(x),replace=TRUE),] else NULL)
#'      xsamp=do.call("rbind",xsamp)
#'      mod=try(lme(formula(zc.weight.model),random=as.formula(zc.weight.model$call$random),data=as.data.frame(xsamp)))
#'      if(class(mod)!="try-error")
#'      {
#'        i=i+1
#'        xsamp=data.frame(days=0,cohort=rep(sort(unique(x$cohort)),2),sex=rep(c("F","M"),each=length(unique(x$cohort))))
#'        pp=predict(mod,newdata=xsamp,level=1)
#'        pmat[i,]=unique(data.frame(predict=as.vector(pp),cohort=xsamp$cohort,sex=xsamp$sex))$predict
#'      }
#'    }
#'    return(sqrt(apply(pmat,2,var)))
#' }
#' # use 100 reps to compute std error
#' stderrors=bootstrap.se(zcweights,100)
#' # Compute predictions and construct dataframes for female and male averages with std errors
#' pp=data.frame(days=0,cohort=rep(sort(unique(zcweights$cohort)),2),sex=rep(c("F","M"),each=length(unique(zcweights$cohort))))
#' pp$predict=predict(zc.weight.model,newdata=pp,level=1)
#' female.averages=data.frame(fit=pp$predict[pp$sex=="F"],se=stderrors[as.numeric(row.names(pp[pp$sex=="F",]))])
#' male.averages=data.frame(fit=pp$predict[pp$sex=="M"],se=stderrors[as.numeric(row.names(pp[pp$sex=="M",]))])
#' # create dataframe with normal conf interval
#' zc.female.averages=female.averages
#' ymin=min(c(female.averages$fit-1.96*female.averages$se, male.averages$fit-1.96*male.averages$se))*.9
#' ymax=max(c(female.averages$fit+1.96*female.averages$se, male.averages$fit+1.96*male.averages$se))*1.1
#' zc.male.averages=male.averages
#' ZCWeight.df=data.frame(female.observed.mean=female.observed,
#' female.adjusted.mean=zc.female.averages$fit,female.adjusted.mean.se=zc.female.averages$se,
#' male.observed.mean=male.observed,
#' male.adjusted.mean=zc.male.averages$fit,male.adjusted.mean.se=zc.male.averages$se)
#' options(width=150)
#' print(ZCWeight.df)
#' #
#' # Create function for plotting and plot female and male averages
#' #
#' pdf("ZCPredictedWeights.pdf",pointsize=10)
#' par(mfrow=c(2,1))
#' maxyear=max(zcweights$cohort)
#' plot.weight.series(1975:maxyear,female.averages,main="California sea lion Female Pups",ylim=c(ymin,ymax))
#' plot.weight.series(1975:maxyear,male.averages,main="California sea lion Male Pups",ylim=c(ymin,ymax))
#' dev.off()
#' 
#' ################################################################################
#' # Construct environmental values for models
#' ################################################################################
#' # Get SST Anomalies at 8 locations using 1994:1996 and 1998:2008 for averages
#' 
#' anomalies=create.SST.anomalies(c(1994:1996,1998:2008),fdir=edir)
#' dir=file.path(edir,"environmental.data.mdb")
#' connection=odbcConnectAccess(dir)
#' # Use Locations 1-5 (ESB,WSB,PtArg,PtSM,PtSL) for pup weight predictions
#' SSTAnomalies=t(apply(anomalies[,,1:5],c(2,1),mean,na.rm=TRUE))
#' maxyear= max(as.numeric(row.names(SSTAnomalies)))
#' minyear= min(as.numeric(row.names(SSTAnomalies)))
#' numyears=maxyear-minyear+1
#' # Compute averages across month grouping so use with period specific growth rate
#' JantoMayAnomalies=rowMeans(SSTAnomalies[,c("Jan","Feb","Mar","Apr","May")])
#' OcttoFebAnomalies=cbind(SSTAnomalies[,c("Oct","Nov","Dec")],rbind(SSTAnomalies[2:nrow(SSTAnomalies),c("Jan","Feb")],data.frame(Jan=NA,Feb=NA)))
#' OcttoFebAnomalies[is.nan(OcttoFebAnomalies)]=NA
#' OcttoFebAnomalies=rowMeans(OcttoFebAnomalies,na.rm=TRUE)
#' JunetoSeptAnomalies=rowMeans(SSTAnomalies[,c("June","July","Aug","Sept")],na.rm=TRUE)
#' JunetoSeptAnomalies[is.nan(JunetoSeptAnomalies)]=NA
#' OcttoDecAnomalies=SSTAnomalies[,c("Oct","Nov","Dec")]
#' OcttoDecAnomalies[is.nan(OcttoDecAnomalies)]=NA
#' OcttoDecAnomalies=rowMeans(OcttoDecAnomalies,na.rm=TRUE)
#' 
#' JunetoFebAnomalies=rowMeans(SSTAnomalies[,c("June","July","Aug","Sept","Oct","Nov","Dec")])
#' JunetoFebAnomalies=cbind(SSTAnomalies[,c("June","July","Aug","Sept","Oct","Nov","Dec")],rbind(SSTAnomalies[2:nrow(SSTAnomalies),c("Jan","Feb")],data.frame(Jan=NA,Feb=NA)))
#' JunetoFebAnomalies[is.nan(JunetoFebAnomalies)]=NA
#' JunetoFebAnomalies=rowMeans(JunetoFebAnomalies,na.rm=TRUE)
#' 
#' x=rbind(data.frame(Year=minyear:maxyear,Season=rep("Spring",numyears),SSTAnomaly=JantoMayAnomalies),
#'   data.frame(Year=minyear:maxyear,Season=rep("Summer",numyears),SSTAnomaly=JunetoSeptAnomalies),
#'   data.frame(Year=minyear:maxyear,Season=rep("Fall",numyears),SSTAnomaly=OcttoDecAnomalies) )
#' x$Season=factor(x$Season,levels=c("Spring","Summer","Fall"))
#' x=x[order(x$Year,x$Season),]
#' MEI=sqlFetch(connection,"MEI")
#' UWI=sqlFetch(connection,"UWIAnomaly")
#' UWI=UWI[order(UWI$Year,UWI$Month),]
#' UWImeansJunetoSept=with(UWI[UWI$Month%in%6:9&UWI$Year>=1975,], tapply(UWIAnomaly,list(Location,Year),mean))
#' 
#' UWIOcttoDec=with(UWI[UWI$Month%in%10:12,], tapply(UWIAnomaly,list(Month,Year,Location),mean))
#' UWIJantoFeb=with(UWI[UWI$Month%in%1:2,], tapply(UWIAnomaly,list(Month,Year,Location),mean))
#' 
#' UWImeansOcttoFeb=NULL
#' for(i in 1:2)
#'    UWImeansOcttoFeb=rbind(UWImeansOcttoFeb,colMeans(rbind(UWIOcttoDec[,,i],UWIJantoFeb[,,i]),na.rm=TRUE))
#'    
#' #
#' # Fit some models to predicted averages
#' #
#' #
#' # Compute correlations between MEI and SST to find the best lag
#' SSTAnomalies.db=data.frame(SSTAnomaly=as.vector(t(SSTAnomalies[-(1:2),])))
#' MEIcor=vector("numeric",8)
#' for(lag in 0:7)
#' {
#'   MEIcor[lag+1]=cor(MEI$MEI[1:(length(MEI$MEI)-lag)],SSTAnomalies.db$SSTAnomaly[(lag+1):length(MEI$MEI)],use="complete.obs")
#'   cat("\nlag = ",lag,"cor = ",MEIcor[lag+1])
#' }
#' lag=which(MEIcor==max(MEIcor))-1
#' #
#' # This MEI is comparable to JunetoSeptSST with defined lag
#' #
#' average.MEI=function(x,months)return(tapply(x$MEI[x$Month%in%months],x$Year[x$Month%in%months],mean))
#' LaggedMEIJunetoSept=average.MEI(MEI,(6:9-lag))
#' LaggedMEIOcttoFeb=average.MEI(MEI,8:12) # assumes 2 month lag to avoid Dec/Jan break
#' 
#' JunetoSeptAnomalies=JunetoSeptAnomalies[4:numyears]
#' OcttoFebAnomalies=OcttoFebAnomalies[4:numyears]
#' JunetoFebAnomalies=JunetoFebAnomalies[4:numyears]
#' zcweights.environ=merge(zcweights,data.frame(cohort=1975:maxyear,SST1=JunetoSeptAnomalies,SST2=OcttoFebAnomalies,SST3=JunetoFebAnomalies,
#'                             MEI=LaggedMEIJunetoSept[-1],MEI1=LaggedMEIOcttoFeb[-1],UWI33=UWImeansJunetoSept[1,],UWI36=UWImeansJunetoSept[2,],
#'                             UWI331=UWImeansOcttoFeb[1,7:41],UWI361=UWImeansOcttoFeb[2,7:41]))
#' zcweights.environ$SST=zcweights.environ$SST1
#' zcweights.environ$SST[zcweights.environ$days>90]=zcweights.environ$SST3[zcweights.environ$days>90]
#' zcweights.environ$MEI[zcweights.environ$days>90]=zcweights.environ$MEI1[zcweights.environ$days>90]
#' zcweights.environ$UWI33[zcweights.environ$days>90]=zcweights.environ$UWI331[zcweights.environ$days>90]
#' zcweights.environ$UWI36[zcweights.environ$days>90]=zcweights.environ$UWI361[zcweights.environ$days>90]
#' 
#' #############################################################################################################
#' # Environmental Cross-sectional analysis of weights
#' #
#' # Evaluate best random effect model with most complex fixed-effect model
#' zc.weight.model=lme(weight~sex*days*SST+UWI33+cohort+MEI,random=list(~days|cohort,~1|batch),data=zcweights.environ)
#' summary(zc.weight.model)$AIC
#' zc.weight.model=lme(weight~sex*days*SST+UWI33+cohort+MEI,random=list(~1|cohort,~1|batch),data=zcweights.environ)
#' summary(zc.weight.model)$AIC
#' zc.weight.model=lme(weight~sex*days*SST+UWI33+cohort+MEI,random=list(~days|cohort),data=zcweights.environ)
#' summary(zc.weight.model)$AIC
#' #
#' # Next evaluate sequence of fixed-effect models with the chosen random effect model
#' #
#' # SST
#' zc.weight.model=lme(weight~sex*days*SST+cohort,random=list(~days|cohort,~1|batch),data=zcweights.environ,method="ML")
#' summary(zc.weight.model)$AIC
#' zc.weight.model=lme(weight~sex*days*SST,random=list(~days|cohort,~1|batch),data=zcweights.environ,method="ML")
#' summary(zc.weight.model)$AIC
#' zc.weight.model=lme(weight~sex*days+days*SST+sex*SST+cohort,random=list(~days|cohort,~1|batch),data=zcweights.environ,method="ML")
#' summary(zc.weight.model)$AIC
#' zc.weight.model=lme(weight~sex*days+days*SST+sex*SST,random=list(~days|cohort,~1|batch),data=zcweights.environ,method="ML")
#' summary(zc.weight.model)$AIC
#' # UWI
#' zc.weight.model=lme(weight~sex*days*UWI33+cohort,random=list(~days|cohort,~1|batch),data=zcweights.environ,method="ML")
#' summary(zc.weight.model)$AIC
#' zc.weight.model=lme(weight~sex*days*UWI33,random=list(~days|cohort,~1|batch),data=zcweights.environ,method="ML")
#' summary(zc.weight.model)$AIC
#' zc.weight.model=lme(weight~sex*days+days*UWI33+sex*UWI33+cohort,random=list(~days|cohort,~1|batch),data=zcweights.environ,method="ML")
#' summary(zc.weight.model)$AIC
#' zc.weight.model=lme(weight~sex*days+days*UWI33+sex*UWI33,random=list(~days|cohort,~1|batch),data=zcweights.environ,method="ML")
#' summary(zc.weight.model)$AIC
#' # MEI
#' zc.weight.model=lme(weight~sex*days*MEI+cohort,random=list(~days|cohort,~1|batch),data=zcweights.environ,method="ML")
#' summary(zc.weight.model)$AIC
#' zc.weight.model=lme(weight~sex*days*MEI,random=list(~days|cohort,~1|batch),data=zcweights.environ,method="ML")
#' summary(zc.weight.model)$AIC
#' zc.weight.model=lme(weight~sex*days+days*MEI+sex*MEI+cohort,random=list(~days|cohort,~1|batch),data=zcweights.environ,method="ML")
#' summary(zc.weight.model)$AIC
#' zc.weight.model=lme(weight~sex*days+days*MEI+sex*MEI,random=list(~days|cohort,~1|batch),data=zcweights.environ,method="ML")
#' summary(zc.weight.model)$AIC
#' # Get predicted means and compute std errors and add to ZCWeight.df table
#' best.zc.weight.model=lme(weight~sex*days*SST+cohort,random=list(~days|cohort,~1|batch),data=zcweights.environ)
#' summary(best.zc.weight.model)
#' zc.weight.model=best.zc.weight.model   
#' bootstrap.se=function(x,nreps)
#' {
#'    pmat=matrix(0,nrow=nreps,ncol=nrow(unique(data.frame(cohort=x$cohort,sex=x$sex))))
#'    i=0
#'    while(i<nreps)
#'    {
#'      xsamp=lapply(split(x,list(x$batch,x$sex)),function(x) if(nrow(x)>0) x[sample(1:nrow(x),replace=TRUE),] else NULL)
#'      xsamp=do.call("rbind",xsamp)
#'      mod=try(lme(formula(zc.weight.model),random=as.formula(zc.weight.model$call$random),data=as.data.frame(xsamp)))
#'      if(class(mod)!="try-error")
#'      {
#'        i=i+1
#'        xsamp$days=0
#'        pp=predict(mod,newdata=xsamp)
#'        pmat[i,]=as.vector(tapply(as.vector(pp),list(x$cohort,x$sex),mean))
#'      }
#'    }
#'    return(sqrt(apply(pmat,2,var)))
#' }
#' # use 100 reps to compute std error
#' stderrors=bootstrap.se(zcweights.environ,100)
#' # Compute predictions and construct dataframes for female and male averages with std errors
#' pp=zcweights.environ
#' pp$days=0
#' pp=predict(zc.weight.model,newdata=pp)
#' pp=tapply(as.vector(pp),list(zcweights.environ$cohort,zcweights.environ$sex),mean)
#' female.averages=data.frame(fit=pp[,1],se=stderrors[1:length(pp[,1])])
#' male.averages=data.frame(fit=pp[,2],se=stderrors[(length(pp[,1])+1):(2*length(pp[,1]))])
#' ZCWeight.df$female.eviron.mean=female.averages$fit
#' ZCWeight.df$female.environ.mean.se=female.averages$se
#' ZCWeight.df$male.environ.mean=male.averages$fit
#' ZCWeight.df$male.environ.mean.se=male.averages$se
#' # Plot predictions and observed
#' win.graph()
#' par(mfrow=c(2,1))
#' with(ZCWeight.df,plot(1975:2009,female.eviron.mean,pch="F",type="b",ylim=c(12,26)))
#' with(ZCWeight.df,points(1975:2009,female.observed.mean,pch="O"))
#' with(ZCWeight.df,lines(1975:2009,female.observed.mean,lty=2))
#' with(ZCWeight.df,plot(1975:2009,male.environ.mean,pch="M",type="b",ylim=c(12,26)))
#' with(ZCWeight.df,points(1975:2009,male.observed.mean,pch="O"))
#' with(ZCWeight.df,lines(1975:2009,male.observed.mean,lty=2))
#' 
#' # Plot predictions and predictions at SST=0
#' pp=zcweights.environ
#' pp$days=0
#' pp$SST=0
#' pp=predict(zc.weight.model,newdata=pp)
#' pp=tapply(as.vector(pp),list(zcweights.environ$cohort,zcweights.environ$sex),mean)
#' female.averages=data.frame(fit=pp[,1],se=stderrors[1:length(pp[,1])])
#' male.averages=data.frame(fit=pp[,2],se=stderrors[(length(pp[,1])+1):(2*length(pp[,1]))])
#' win.graph()
#' par(mfrow=c(2,1))
#' with(ZCWeight.df,plot(1975:2009,female.eviron.mean,pch="F",type="b",ylim=c(12,26)))
#' points(1975:2009,female.averages$fit,pch="S")
#' lines(1975:2009,female.averages$fit,pch="S",lty=2)
#' with(ZCWeight.df,plot(1975:2009,male.environ.mean,pch="M",type="b",ylim=c(12,26)))
#' points(1975:2009,male.averages$fit,pch="S")
#' lines(1975:2009,male.averages$fit,pch="S",lty=2)
#' 
#' 
#' ################################################################################
#' # Longitudinal Environmental Growth Analysis
#' ################################################################################
#' # get zc weight values from database
#' zcweights.lon=get.zc.weights(fdir=fdir)
#' #
#' #  Use SMI from 1 Sept to end of Feb
#' #
#' zcweights.lon=zcweights.lon[zcweights.lon$days>-30&zcweights.lon$days<150&zcweights.lon$sitecode=="SMI",]
#' # Extract the pups with more than one measurement
#' long.id=table(zcweights.lon$AnimalID)
#' long.id=names(long.id[long.id>1])
#' zcweights.lon=zcweights.lon[zcweights.lon$AnimalID%in%long.id,]
#' # Create initial and recap dataframes and exclude any with more than one initial
#' initial=zcweights.lon[zcweights.lon$type=="Initial",]
#' recaps=zcweights.lon[zcweights.lon$type=="Recap",]
#' long.id=table(initial$AnimalID)
#' long.id=names(long.id[long.id==1])
#' initial=initial[initial$AnimalID%in%long.id,]
#' zcweights.lon=rbind(initial,recaps)
#' # Compute growth rates
#' gr=sapply(split(zcweights.lon,list(zcweights.lon$AnimalID)),function(x) diff(x$weight[order(x$days)])/diff(x$days[order(x$days)]))
#' # Merge with data about animal
#' grdata=data.frame(gr=gr,AnimalID=names(gr))
#' grdata=merge(initial,grdata)
#' # Merge growth data with environmental data
#' zcweights.environ.lon=merge(grdata,data.frame(cohort=1975:maxyear,SST=OcttoFebAnomalies,MEI=LaggedMEIOcttoFeb[-1],
#'                                          UWI33=UWImeansOcttoFeb[1,7:41],UWI36=UWImeansOcttoFeb[2,7:41]))
#' # Create model
#' mod=lme(gr~sex*UWI33,random=~1|cohort,data=zcweights.environ.lon,method="ML")
#' summary(mod)$AIC
#' mod=lme(gr~sex*MEI,random=~1|cohort,data=zcweights.environ.lon,method="ML")
#' summary(mod)$AIC
#' mod=lme(gr~sex*SST,random=~1|cohort,data=zcweights.environ.lon,method="ML")
#' summary(mod)$AIC
#' mod=lme(gr~sex*SST,random=~1|cohort,data=zcweights.environ.lon)
#' plot(zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="F"],predict(mod)[zcweights.environ.lon$sex=="F"],pch="F",type="b",ylim=c(0,.12))
#' points(zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="M"],predict(mod)[zcweights.environ.lon$sex=="M"],pch="M",type="b")
#' lines(zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="M"],predict(mod)[zcweights.environ.lon$sex=="M"],lty=2)
#' # Observed and predicted values
#' win.graph()
#' par(mfrow=c(1,2))
#' plot(zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="F"],predict(mod)[zcweights.environ.lon$sex=="F"],pch="F",type="b",ylim=c(0,.14),ylab="Female daily growth rate (kg/day)",xlab="Cohort")
#' points(c(1989,1990,1995,1997:max(zcweights.environ.lon$cohort)),tapply(zcweights.environ.lon$gr[zcweights.environ.lon$sex=="F"],zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="F"],mean),pch="O")
#' plot(zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="M"],predict(mod)[zcweights.environ.lon$sex=="M"],pch="M",type="b",ylim=c(0,.14),ylab="Male daily growth rate (kg/day)",xlab="Cohort")
#' points(c(1989,1990,1995,1997:max(zcweights.environ.lon$cohort)),tapply(zcweights.environ.lon$gr[zcweights.environ.lon$sex=="M"],zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="M"],mean),pch="O")
#' 
#' 
#' ################################################################################
#' #  Brand vs control at brand eval
#' ################################################################################
#' zcweights.be=get.zc.weights(fdir=fdir)
#' #
#' #  Use SMI from 1 Dec to end of Feb; use brand and unmarked only
#' #
#' zcweights.be=zcweights.be[zcweights.be$days>60&zcweights.be$days<150&zcweights.be$sitecode=="SMI"&zcweights.be$tag.type!="Tag",]
#' zcweights.be$tag.type=factor(zcweights.be$tag.type)
#' # lme model of tag.type effect
#' mod=lme(weight~tag.type*sex,random=~1|cohort,data=zcweights.be,method="ML")
#' summary(mod)
#' mod0=lme(weight~sex,random=~1|cohort,data=zcweights.be,method="ML")
#' summary(mod0)
#' mod1=lme(weight~tag.type+sex,random=~1|cohort,data=zcweights.be,method="ML")
#' summary(mod1)
#' anova(mod0,mod1,mod)
#' # kruskal test by sex
#' kruskal.test(weight~tag.type*cohort,data=zcweights.be[zcweights.be$sex=="F",])
#' kruskal.test(weight~tag.type*cohort,data=zcweights.be[zcweights.be$sex=="M",])
#' # paired t-test
#' be.fem=with(zcweights.be,tapply(weight,list(cohort,tag.type,sex),mean))[,,"F"]
#' be.m=with(zcweights.be,tapply(weight,list(cohort,tag.type,sex),mean))[,,"M"]
#' t.test(x=c(be.fem[,1],be.m[,1]),y=c(be.fem[,2],be.m[,2]),paired=TRUE)
#' # Close ACCESS tables
#' odbcCloseAll()
#' 
get.zc.weights=function(fdir="./",ENYears=c(1976,1983,1984,1987,1992,1997,1998))
{
fdir=file.path(fdir,"Brandmaster.mdb")
connection=odbcConnectAccess2007(fdir)
#
# read in weights from ZcBrand table of the Access database
#
zcweights=sqlFetch(connection,"ZCBRAND")
maxyear=max(zcweights$cohort)
#
#  Exclude those with missing weight and unknown sex
#
zcweights=zcweights[!zcweights$sex=="U" & !is.na(zcweights$weight)& !is.na(zcweights$sex),]
zcweights$sex=factor(zcweights$sex)
#
#  read in the weights from the unmarkedpupweights table.  Exclude those with unknown sex or weight and from SCI and SNI.
#
zcweights.unmarked=sqlFetch(connection,"UnmarkedPupWeights")
zcweights.unmarked=zcweights.unmarked[!zcweights.unmarked$sex=="U"& !is.na(zcweights.unmarked$sex),]
#zcweights.unmarked=zcweights.unmarked[!zcweights.unmarked$sitecode=="SNI"& !is.na(zcweights.unmarked$sitecode),]
#zcweights.unmarked=zcweights.unmarked[!zcweights.unmarked$sitecode=="SCI"& !is.na(zcweights.unmarked$sitecode),]
zcweights.unmarked$sex=factor(zcweights.unmarked$sex)
#
#  read in the tag only pup weights but exclude those with missing sex or weight
#
zcTagOnly=sqlFetch(connection,"TagInitial")
zcTagOnly=zcTagOnly[!zcTagOnly$sex=="U" & !is.na(zcTagOnly$weight)& !is.na(zcTagOnly$sex),]
zcTagOnly$sex=factor(zcTagOnly$sex)
#
#  read in the brand recaptures
#
zcRecap=sqlFetch(connection,"Recaptures")
zcRecap=zcRecap[!zcRecap$sex=="U" & !zcRecap$sex=="" &!is.na(zcRecap$weight)& !is.na(zcRecap$sex),]
zcRecap$sex=factor(zcRecap$sex)
zcRecap=merge(zcRecap,subset(zcweights,select=c("cohort","brand")),by="brand")
#
#  read in the tag recaptures
#
zcTagRecap=sqlFetch(connection,"TagRecapture")
zcTagRecap=zcTagRecap[!zcTagRecap$Sex=="U" & !is.na(zcTagRecap$Weight)& !is.na(zcTagRecap$Sex),]
zcTagRecap$sex=factor(zcTagRecap$Sex)
zcTagRecap=merge(zcTagRecap,subset(zcTagOnly,select=c("cohort","AnimalID")),by="AnimalID")
#
#   Add a days field which is the number of days from 1 Oct of the year
#
zcweights$days=floor(as.numeric((zcweights$branddate-as.POSIXct(paste(as.character(zcweights$cohort),"-10-01",sep="")))/(24*3600)))
zcTagOnly$days=floor(as.numeric((zcTagOnly$capturedate-as.POSIXct(paste(as.character(zcTagOnly$cohort),"-10-01",sep="")))))
zcweights.unmarked$days=floor(as.numeric((zcweights.unmarked$capturedate-as.POSIXct(paste(as.character(zcweights.unmarked$cohort),"-10-01",sep="")))/(24*3600)))
zcRecap$days=floor(as.numeric((zcRecap$capturedate-as.POSIXct(paste(as.character(zcRecap$cohort),"-10-01",sep="")))))
zcTagRecap$days=floor(as.numeric((zcTagRecap$capturedate-as.POSIXct(paste(as.character(zcTagRecap$cohort),"-10-01",sep="")))))
zcRecap=zcRecap[zcRecap$days>0,]
#
#  Combine data tables
#
zcweights.unmarked=subset(zcweights.unmarked,select=c("AnimalID","sitecode","cohort","sex","weight","length","girth","days"))
zcweights.unmarked=zcweights.unmarked[!is.na(zcweights.unmarked$weight),]
zcweights$sitecode="SMI"
zcweights=subset(zcweights,select=c("AnimalID","sitecode","cohort","sex","weight","length","girth","days"))
zcTagOnly$sitecode="SMI"
zcTagOnly=subset(zcTagOnly,select=c("AnimalID","sitecode","cohort","sex","weight","length","girth","days"))
zcRecap=subset(zcRecap,select=c("AnimalID","sitecode","cohort","sex","weight","length","girth","days"))
#zcTagRecap$sitecode="SMI"
zcTagRecap=subset(zcTagRecap,select=c("AnimalID","sitecode","cohort","sex","Weight","Length","Girth","days"))
names(zcTagRecap)=c("AnimalID","sitecode","cohort","sex","weight","length","girth","days")
zcweights$AnimalID=as.character(zcweights$AnimalID)
zcTagOnly$AnimalID=as.character(zcTagOnly$AnimalID)
zcTagRecap$AnimalID=as.character(zcTagRecap$AnimalID)
zcRecap$AnimalID=as.character(zcRecap$AnimalID)
zcweights.unmarked$AnimalID=as.character(zcweights.unmarked$AnimalID)
zcweights$tag.type="Brand"
zcweights.unmarked$tag.type="Unmarked"
zcRecap$tag.type="Brand"
zcTagRecap$tag.type="Tag"
zcTagOnly$tag.type="Tag"
zcweights$type="Initial"
zcweights.unmarked$type="Initial"
zcTagOnly$type="Initial"
zcTagRecap$type="Recap"
zcRecap$type="Recap"
zcweights=rbind(zcTagOnly,zcweights,zcweights.unmarked,zcRecap,zcTagRecap)
zcweights$type=factor(zcweights$type)
zcweights$sitecode=factor(zcweights$sitecode)
zcweights$tag.type=factor(zcweights$tag.type)
#
# Define EN field
#
zcweights$EN=0
zcweights$EN[zcweights$cohort%in%ENYears]=1
return(zcweights)
}
