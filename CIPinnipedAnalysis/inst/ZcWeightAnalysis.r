plot.weight.series=function (time, predictions, date,...)
{
	plotCI(time, predictions$fit, 1.96 * predictions$se,
			1.96 * predictions$se, xlab = "Cohort",
			ylab = paste("Predicted average weight (kg)",date),
			...)
	lines(time, predictions$fit)
	invisible()
}
fdir=NULL
dir=NULL
# use "" to use databases in package directory
library(CIPinnipedAnalysis)
fdir=""
dir=""
#################################################################################
# Cross-sectional analysis
#################################################################################
# get zc weight values from database
zcweights=get.zc.weights(fdir=fdir)
#
#  exclude brand eval weights and captures in April
#
zcweights.all=zcweights[zcweights$days>-30&zcweights$days<150&
				zcweights$sitecode=="SMI",]
#
# compute observed averages - using data within 1 Sept to 15 Nov
#
zcweights=zcweights.all[zcweights.all$days>=-31&zcweights.all$days<=45,]
female.observed=sapply(split(zcweights$weight[zcweights$sex=="F"],
				factor(zcweights$cohort[zcweights$sex=="F"],levels=levels(factor(zcweights.all$cohort)))),mean)
male.observed=sapply(split(zcweights$weight[zcweights$sex=="M"],
				factor(zcweights$cohort[zcweights$sex=="M"],levels=levels(factor(zcweights.all$cohort)))),mean)
zcweights=zcweights.all
zcweights$batch=factor(paste(zcweights$cohort,zcweights$days))
#
#  fit growth model
#
require(nlme)
zc.weight.model=lme(weight~sex*days,random=list(~sex*days|cohort,~1|batch),data=zcweights)
zc.weight.model=lme(weight~-1+sex+days+sex:days,random=list(~sex*days|cohort,~1|batch),data=zcweights)
print(summary(zc.weight.model))
# define bootstrap function to compute std error for predicted sex-cohort means
bootstrap.se=function(x,nreps,days)
{
	pmat=matrix(0,nrow=nreps,ncol=nrow(unique(data.frame(cohort=x$cohort,sex=x$sex))))
	i=0
	while(i<nreps)
	{
		cat("Bootstrap ",i,"\n")
		xsamp=lapply(split(x,list(x$batch,x$sex)),function(x) if(nrow(x)>0) x[sample(1:nrow(x),replace=TRUE),] else NULL)
		xsamp=do.call("rbind",xsamp)
		mod=try(lme(formula(zc.weight.model),random=as.formula(zc.weight.model$call$random),data=as.data.frame(xsamp),control=lmeControl(maxIter=100, msMaxIter=100)))
		if(class(mod)!="try-error")
		{
			i=i+1
			xsamp=data.frame(days=days,cohort=rep(sort(unique(x$cohort)),2),sex=rep(c("F","M"),each=length(unique(x$cohort))))
			pp=predict(mod,newdata=xsamp,level=1)
			pmat[i,]=unique(data.frame(predict=as.vector(pp),cohort=xsamp$cohort,sex=xsamp$sex))$predict
		}
	}
	return(sqrt(apply(pmat,2,var)))
}
# use 100 reps to compute std error
stderrors=bootstrap.se(zcweights,100,days=0)
# Compute predictions and construct dataframes for female and male averages with std errors
pp=data.frame(days=0,cohort=rep(sort(unique(zcweights$cohort)),2),sex=rep(c("F","M"),each=length(unique(zcweights$cohort))))
pp$predict=predict(zc.weight.model,newdata=pp,level=1)
female.averages=data.frame(fit=pp$predict[pp$sex=="F"],se=stderrors[as.numeric(row.names(pp[pp$sex=="F",]))])
male.averages=data.frame(fit=pp$predict[pp$sex=="M"],se=stderrors[as.numeric(row.names(pp[pp$sex=="M",]))])
# create dataframe with normal conf interval
zc.female.averages=female.averages
ymin=min(c(female.averages$fit-1.96*female.averages$se, male.averages$fit-1.96*male.averages$se))*.9
ymax=max(c(female.averages$fit+1.96*female.averages$se, male.averages$fit+1.96*male.averages$se))*1.1
zc.male.averages=male.averages
ZCWeight.df=data.frame(female.observed.mean=female.observed,
		female.adjusted.mean=zc.female.averages$fit,female.adjusted.mean.se=zc.female.averages$se,
		male.observed.mean=male.observed,
		male.adjusted.mean=zc.male.averages$fit,male.adjusted.mean.se=zc.male.averages$se)
options(width=150)
print(ZCWeight.df)
#
# Create function for plotting and plot female and male averages
#
pdf("ZCPredictedWeights.pdf",pointsize=10)
par(mfrow=c(2,1))
maxyear=max(zcweights$cohort)
plot.weight.series(1975:maxyear,female.averages,main="California sea lion Female Pups",ylim=c(ymin,ymax),date="1 Oct")
plot.weight.series(1975:maxyear,male.averages,main="California sea lion Male Pups",ylim=c(ymin,ymax),date="1 Oct")
dev.off()

jpeg("ZCPredictedWeightsCalcofi.jpg",,height=600,width=600,quality=100,pointsize=12)
maxyear=max(zcweights$cohort)
par(lty=1)
plot.weight.series(1997:maxyear,female.averages[23:(maxyear-1974),],ylim=c(ymin,ymax),xaxp=c(1998,2012,7),date="1 Oct")
abline(h=17.4)
par(lty=2)
plot.weight.series(1997:maxyear,male.averages[23:(maxyear-1974),],pch=2,add=TRUE,slty=1,date="1 Oct")
abline(h=17.4+2.76)
points(2009,25,pch=2)
lines(x=c(2008.75,2009.25),y=c(25,25),pch=2,lty=2)
points(2009,23,pch=1)
lines(x=c(2008.75,2009.25),y=c(23,23),pch=1,lty=1)
text(2009.4,25,"Males",pos=4)
text(2009.4,23,"Females",pos=4)
par(lty=1)
dev.off()

# 1 Feb predictions
pp=data.frame(days=123,cohort=rep(sort(unique(zcweights$cohort)),2),sex=rep(c("F","M"),each=length(unique(zcweights$cohort))))
pp$predict=predict(zc.weight.model,newdata=pp,level=1)
stderrors=bootstrap.se(zcweights,100,days=123)
female.averages=data.frame(fit=pp$predict[pp$sex=="F"],se=stderrors[as.numeric(row.names(pp[pp$sex=="F",]))])
male.averages=data.frame(fit=pp$predict[pp$sex=="M"],se=stderrors[as.numeric(row.names(pp[pp$sex=="M",]))])
ymin=min(c(female.averages$fit-1.96*female.averages$se, male.averages$fit-1.96*male.averages$se))*.9
ymax=max(c(female.averages$fit+1.96*female.averages$se, male.averages$fit+1.96*male.averages$se))*1.1

jpeg("ZCPredictedWeightsFebCalcofi.jpg",height=600,width=600,quality=100,pointsize=12)
maxyear=max(zcweights$cohort)
par(lty=1)
plot.weight.series(1997:maxyear,female.averages[23:(maxyear-1974),],ylim=c(ymin,ymax),xaxp=c(1998,2012,7),date="1 Feb")
abline(h=25.69)
par(lty=2)
plot.weight.series(1997:maxyear,male.averages[23:(maxyear-1974),],pch=2,add=TRUE,slty=1,date="1 Feb")
abline(h=30.55)
points(2009,40,pch=2)
lines(x=c(2008.75,2009.25),y=c(40,40),pch=2,lty=2)
points(2009,38,pch=1)
lines(x=c(2008.75,2009.25),y=c(38,38),pch=1,lty=1)
text(2009.4,40,"Males",pos=4)
text(2009.4,38,"Females",pos=4)
par(lty=1)
dev.off()


################################################################################
# Construct environmental values for models
################################################################################
# Get SST Anomalies at 8 locations using 1994:1996 and 1998:2008 for averages

anomalies=create.SST.anomalies(c(1994:1996,1998:2008),fdir=fdir)
MEI=getCalcurData("Environ","MEI",dir=dir)
UWI=getCalcurData("Environ","UWIAnomaly",dir=dir)
# Use Locations 1-5 (ESB,WSB,PtArg,PtSM,PtSL) for pup weight predictions
SSTAnomalies=t(apply(anomalies[,,1:5],c(2,1),mean,na.rm=TRUE))
SSTAnomalies[is.nan(SSTAnomalies)]=NA
maxyear= max(as.numeric(row.names(SSTAnomalies)))
minyear= min(as.numeric(row.names(SSTAnomalies)))
#if(is.na(SSTAnomalies[maxyear-minyear+1,11]) &is.na(SSTAnomalies[maxyear-minyear+1,12]))maxyear=maxyear-1
numyears=maxyear-minyear+1
# Compute averages across month grouping so use with period specific growth rate
JantoMayAnomalies=rowMeans(SSTAnomalies[,c("Jan","Feb","Mar","Apr","May")])[1:(maxyear-minyear+1)]
OtoD=SSTAnomalies[,c("Oct","Nov","Dec")][1:(maxyear-minyear+1),]
JtoF=SSTAnomalies[2:nrow(SSTAnomalies),c("Jan","Feb")]
if(nrow(JtoF)<nrow(OtoD))JtoF=rbind(JtoF,c(NA,NA))
OcttoFebAnomalies=as.matrix(cbind(OtoD,JtoF))
OcttoFebAnomalies[is.nan(OcttoFebAnomalies)]=NA
OcttoFebAnomalies=rowMeans(OcttoFebAnomalies,na.rm=TRUE)[1:(maxyear-minyear+1)]
JunetoSeptAnomalies=SSTAnomalies[,c("June","July","Aug","Sept")]
JunetoSeptAnomalies[is.nan(JunetoSeptAnomalies)]=NA
JunetoSeptAnomalies=rowMeans(JunetoSeptAnomalies,na.rm=TRUE)[1:(maxyear-minyear+1)]
OcttoDecAnomalies=SSTAnomalies[,c("Oct","Nov","Dec")]
OcttoDecAnomalies[is.nan(OcttoDecAnomalies)]=NA
OcttoDecAnomalies=rowMeans(OcttoDecAnomalies,na.rm=TRUE)[1:(maxyear-minyear+1)]
JtoD=SSTAnomalies[,c("June","July","Aug","Sept","Oct","Nov","Dec")][1:(maxyear-minyear+1),]
JtoF=SSTAnomalies[2:nrow(SSTAnomalies),c("Jan","Feb")]
if(nrow(JtoF)<nrow(JtoD))JtoF=rbind(JtoF,c(NA,NA))
JunetoFebAnomalies=as.matrix(cbind(JtoD,JtoF))
JunetoFebAnomalies[is.nan(JunetoFebAnomalies)]=NA
JunetoFebAnomalies=rowMeans(JunetoFebAnomalies,na.rm=TRUE)

x=rbind(data.frame(Year=minyear:maxyear,Season=rep("Spring",numyears),SSTAnomaly=as.vector(JantoMayAnomalies)),
		data.frame(Year=minyear:maxyear,Season=rep("Summer",numyears),SSTAnomaly=as.vector(JunetoSeptAnomalies)),
		data.frame(Year=minyear:maxyear,Season=rep("Fall",numyears),SSTAnomaly=as.vector(OcttoDecAnomalies)) )
x$Season=factor(x$Season,levels=c("Spring","Summer","Fall"))
x=x[order(x$Year,x$Season),]

UWI33YrxMonth=with(UWI[UWI$Location=="33N119W",],tapply(UWIAnomaly,list(Year,Month),mean))
UWI33JunetoSept=apply(UWI33YrxMonth[,6:9],1,mean)[1:(maxyear-minyear+4)]
UWI33OcttoFeb=apply(cbind(UWI33YrxMonth[1:(nrow(UWI33YrxMonth)-1),10:12],UWI33YrxMonth[2:nrow(UWI33YrxMonth),1:2]),1,mean)[1:(maxyear-minyear+4)]
UWI33JunetoFeb=apply(cbind(UWI33YrxMonth[1:(nrow(UWI33YrxMonth)-1),6:12],UWI33YrxMonth[2:nrow(UWI33YrxMonth),1:2]),1,mean)[1:(maxyear-minyear+4)]
UWI36YrxMonth=with(UWI[UWI$Location=="36N122W",],tapply(UWIAnomaly,list(Year,Month),mean))
UWI36JunetoSept=apply(UWI36YrxMonth[,6:9],1,mean)[1:(maxyear-minyear+4)]
UWI36OcttoFeb=apply(cbind(UWI36YrxMonth[1:(nrow(UWI36YrxMonth)-1),10:12],UWI36YrxMonth[2:nrow(UWI36YrxMonth),1:2]),1,mean)[1:(maxyear-minyear+4)]
UWI36JunetoFeb=apply(cbind(UWI36YrxMonth[1:(nrow(UWI36YrxMonth)-1),6:12],UWI36YrxMonth[2:nrow(UWI36YrxMonth),1:2]),1,mean)[1:(maxyear-minyear+4)]


#
# Fit some models to predicted averages
#
#
# Compute correlations between MEI and SST to find the best lag
SSTAnomalies.db=data.frame(SSTAnomaly=as.vector(t(SSTAnomalies[-(1:2),])))
MEIcor=vector("numeric",8)
for(lag in 0:7)
{
	MEIcor[lag+1]=cor(MEI$MEI[1:(length(MEI$MEI)-lag)],SSTAnomalies.db$SSTAnomaly[(lag+1):length(MEI$MEI)],use="complete.obs")
	cat("\nlag = ",lag,"cor = ",MEIcor[lag+1])
}
lag=which(MEIcor==max(MEIcor))-1
#
# This MEI is comparable to JunetoSeptSST with defined lag
#
MEIYrxMonth=with(MEI,tapply(MEI,list(Year,Month),mean))

JtoS=6:9-lag
JtoF=6:14-lag
OtoF=10:14-lag
LaggedMEIJunetoSept=apply(MEIYrxMonth[,JtoS],1,mean)[1:(maxyear-minyear-1)]
LaggedMEIOcttoFeb=apply(cbind(MEIYrxMonth[1:(nrow(MEIYrxMonth)-1),OtoF[OtoF<=12]],MEIYrxMonth[2:nrow(MEIYrxMonth),OtoF[OtoF>12]-12]),1,mean)[1:(maxyear-minyear-1)]
LaggedMEIJunetoFeb=apply(cbind(MEIYrxMonth[1:(nrow(MEIYrxMonth)-1),JtoF[JtoF<=12]],MEIYrxMonth[2:nrow(MEIYrxMonth),JtoF[JtoF>12]-12]),1,mean)[1:(maxyear-minyear-1)]

# Create dataframe with weights and environmental variables
JunetoSeptAnomalies=JunetoSeptAnomalies[4:numyears]
OcttoFebAnomalies=OcttoFebAnomalies[4:numyears]
JunetoFebAnomalies=JunetoFebAnomalies[4:numyears]
zcweights.environ=merge(zcweights,data.frame(cohort=1975:maxyear,SST=JunetoSeptAnomalies,SST1=OcttoFebAnomalies,SST2=JunetoFebAnomalies,
				MEI=LaggedMEIJunetoSept[-1],MEI1=LaggedMEIOcttoFeb[-1],MEI2=LaggedMEIJunetoFeb[-1],UWI33=UWI33JunetoSept[-(1:6)],UWI36=UWI36JunetoSept[-(1:6)],
				UWI331=UWI33OcttoFeb[-(1:6)],UWI361=UWI36OcttoFeb[-(1:6)],UWI332=UWI33JunetoFeb[-(1:6)],UWI362=UWI36JunetoFeb[-(1:6)]))
zcweights.environ$cohort.factor=factor(ifelse(zcweights.environ$cohort<1990,0,1),labels=c("<=1989",">=1990"))
zcweights.environ$cohort=zcweights.environ$cohort-min(zcweights.environ$cohort)

#############################################################################################################
# Environmental Cross-sectional analysis of weights
#
# Evaluate best random effect model with most complex fixed-effect model
zc.weight.model=lme(weight~sex*SST+sex:days+SST1:days+sex:SST1:days+UWI33+cohort+MEI,random=list(~days|cohort,~1|batch),data=zcweights.environ)
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*SST+sex:days+SST1:days+sex:SST1:days+UWI33+cohort+MEI,random=list(~1|cohort,~1|batch),data=zcweights.environ)
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*SST+sex:days+SST1:days+sex:SST1:days+UWI33+cohort+MEI,random=list(~days|cohort),data=zcweights.environ)
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*SST+sex:days+SST1:days+sex:SST1:days+UWI33+cohort+MEI,random=list(~sex:days|cohort,~1|batch),data=zcweights.environ)
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*SST+sex:days+SST1:days+sex:SST1:days+UWI33+cohort+MEI,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ)
summary(zc.weight.model)$AIC
#
# Next evaluate sequence of fixed-effect models with the chosen random effect model
#
# SST
zc.weight.model=lme(weight~sex*SST+sex:days+SST1:days+sex:SST1:days+cohort.factor,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*SST+sex:days+SST1:days+sex:SST1:days+cohort,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*SST+sex:days+SST1:days+sex:SST1:days,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*SST+sex:days+SST1:days+cohort.factor,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*SST+sex:days+SST1:days+cohort,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*SST+sex:days+SST1:days,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
# UWI
zc.weight.model=lme(weight~sex*UWI36+sex:days+UWI361:days+sex:UWI361:days+cohort.factor,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*UWI36+sex:days+UWI361:days+sex:UWI361:days+cohort,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*UWI36+sex:days+UWI361:days+sex:UWI361:days,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*UWI36+sex:days+UWI361:days+cohort.factor,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*UWI36+sex:days+UWI361:days+cohort,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*UWI36+sex:days+UWI361:days,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
# MEI
zc.weight.model=lme(weight~sex*MEI+sex:days+MEI1:days+sex:MEI1:days+cohort.factor,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*MEI+sex:days+MEI1:days+sex:MEI1:days+cohort,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*MEI+sex:days+MEI1:days+sex:MEI1:days,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*MEI+sex:days+MEI1:days+cohort.factor,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*MEI+sex:days+MEI1:days+cohort,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
zc.weight.model=lme(weight~sex*MEI+sex:days+MEI1:days,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(zc.weight.model)$AIC
# Get predicted means and compute std errors and add to ZCWeight.df table
best.zc.weight.model=lme(weight~sex*MEI+sex:days+MEI1:days,random=list(~sex*days|cohort,~1|batch),data=zcweights.environ,method="ML")
summary(best.zc.weight.model)
zc.weight.model=best.zc.weight.model

pp=zcweights.environ
pp$days=0
pp1=predict(zc.weight.model,newdata=pp)

pp0=predict(zc.weight.model,newdata=pp,level=0)

pp0=tapply(as.vector(pp0),list(zcweights.environ$cohort,zcweights.environ$sex),mean)
pp1=tapply(as.vector(pp1),list(zcweights.environ$cohort,zcweights.environ$sex),mean)

female.averages=data.frame(fit=pp1[,1],se=stderrors[1:length(pp1[,1])])
male.averages=data.frame(fit=pp1[,2],se=stderrors[(length(pp1[,1])+1):(2*length(pp1[,1]))])

expected.female.averages=data.frame(fit=pp0[,1],se=stderrors[1:length(pp0[,1])])
expected.male.averages=data.frame(fit=pp0[,2],se=stderrors[(length(pp0[,1])+1):(2*length(pp0[,1]))])

#ZCWeight.df=ZCWeight.df[rownames(ZCWeight.df)%in%rownames(pp1),]

ZCWeight.df$female.eviron.mean=female.averages$fit
ZCWeight.df$female.environ.mean.se=female.averages$se
ZCWeight.df$male.environ.mean=male.averages$fit
ZCWeight.df$male.environ.mean.se=male.averages$se
ZCWeight.df$Year=as.numeric(rownames(ZCWeight.df))



# Plot predictions and observed
win.graph()
par(mfrow=c(2,1))
with(ZCWeight.df,plot(Year,female.eviron.mean,pch="F",type="b",ylim=c(12,26)))
with(ZCWeight.df,points(Year,female.observed.mean,pch="O"))
with(ZCWeight.df,lines(Year,female.observed.mean,lty=2))
with(ZCWeight.df,plot(Year,male.environ.mean,pch="M",type="b",ylim=c(12,26)))
with(ZCWeight.df,points(Year,male.observed.mean,pch="O"))
with(ZCWeight.df,lines(Year,male.observed.mean,lty=2))


# Plot residuals of predictions based on fixed effects only versus mixed effects


win.graph()
plot(ZCWeight.df$Year,male.averages$fit-expected.male.averages$fit,pch="M",type="b")
abline(0,0)

win.graph()
plot(ZCWeight.df$Year,female.averages$fit-expected.female.averages$fit,pch="F",type="b",xlab="Year",ylab="Observed- model predicted weight")
abline(0,0)

win.graph()
plot(ZCWeight.df$Year,female.averages$fit,type="b",ylim=c(12,22),xlab="Year",ylab="Female pup weight (kg)")
lines(ZCWeight.df$Year,expected.female.averages$fit,type="b",lty=2,pch=2)
points(2000,14,pch=1)
points(2000,13,pch=2)
text(2000,14,"1 Oct mean weight",pos=4)
text(2000,13,"Model predicted mean weight",pos=4)


bootstrap.se=function(x,nreps)
{
	pmat=matrix(0,nrow=nreps,ncol=nrow(unique(data.frame(cohort=x$cohort,sex=x$sex))))
	i=0
	while(i<nreps)
	{
		cat("Bootstrap ",i,"\n")
		xsamp=lapply(split(x,list(x$batch,x$sex)),function(x) if(nrow(x)>0) x[sample(1:nrow(x),replace=TRUE),] else NULL)
		xsamp=do.call("rbind",xsamp)
		mod=try(lme(formula(zc.weight.model),random=as.formula(zc.weight.model$call$random),data=as.data.frame(xsamp)))
		if(class(mod)!="try-error")
		{
			i=i+1
			xsamp$days=0
			pp=predict(mod,newdata=xsamp)
			pmat[i,]=as.vector(tapply(as.vector(pp),list(x$cohort,x$sex),mean))
		}
	}
	return(sqrt(apply(pmat,2,var)))
}
# use 100 reps to compute std error
stderrors=bootstrap.se(zcweights.environ,100)

# Compute predictions and construct dataframes for female and male averages with std errors
pp=zcweights.environ
pp$days=0
pp1=predict(zc.weight.model,newdata=pp)

pp0=predict(zc.weight.model,newdata=pp,level=0)

pp0=tapply(as.vector(pp0),list(zcweights.environ$cohort,zcweights.environ$sex),mean)
pp1=tapply(as.vector(pp1),list(zcweights.environ$cohort,zcweights.environ$sex),mean)

female.averages=data.frame(fit=pp1[,1],se=stderrors[1:length(pp1[,1])])
male.averages=data.frame(fit=pp1[,2],se=stderrors[(length(pp1[,1])+1):(2*length(pp1[,1]))])

expected.female.averages=data.frame(fit=pp0[,1],se=stderrors[1:length(pp0[,1])])
expected.male.averages=data.frame(fit=pp0[,2],se=stderrors[(length(pp0[,1])+1):(2*length(pp0[,1]))])

ZCWeight.df$female.eviron.mean=female.averages$fit
ZCWeight.df$female.environ.mean.se=female.averages$se
ZCWeight.df$male.environ.mean=male.averages$fit
ZCWeight.df$male.environ.mean.se=male.averages$se
# Plot predictions and observed
win.graph()
par(mfrow=c(2,1))
with(ZCWeight.df,plot(1975:maxyear,female.eviron.mean,pch="F",type="b",ylim=c(12,26)))
with(ZCWeight.df,points(1975:maxyear,female.observed.mean,pch="O"))
with(ZCWeight.df,lines(1975:maxyear,female.observed.mean,lty=2))
with(ZCWeight.df,plot(1975:maxyear,male.environ.mean,pch="M",type="b",ylim=c(12,26)))
with(ZCWeight.df,points(1975:maxyear,male.observed.mean,pch="O"))
with(ZCWeight.df,lines(1975:maxyear,male.observed.mean,lty=2))


# Plot residuals of predictions based on fixed effects only versus mixed effects

win.graph()
plot(1975:maxyear,female.averages$fit-expected.female.averages$fit,pch="F",type="b")
win.graph()
plot(1975:maxyear,male.averages$fit-expected.male.averages$fit,pch="M",type="b")


# Plot predictions and predictions at SST=0
pp=zcweights.environ
pp$days=0
pp$SST=0
pp=predict(zc.weight.model,newdata=pp)
pp=tapply(as.vector(pp),list(zcweights.environ$cohort,zcweights.environ$sex),mean)
female.averages=data.frame(fit=pp[,1],se=stderrors[1:length(pp[,1])])
male.averages=data.frame(fit=pp[,2],se=stderrors[(length(pp[,1])+1):(2*length(pp[,1]))])
win.graph()
par(mfrow=c(2,1))
with(ZCWeight.df,plot(1975:maxyear,female.eviron.mean,pch="F",type="b",ylim=c(12,26)))
points(1975:maxyear,female.averages$fit,pch="S")
lines(1975:maxyear,female.averages$fit,pch="S",lty=2)
with(ZCWeight.df,plot(1975:maxyear,male.environ.mean,pch="M",type="b",ylim=c(12,26)))
points(1975:maxyear,male.averages$fit,pch="S")
lines(1975:maxyear,male.averages$fit,pch="S",lty=2)


################################################################################
# Longitudinal Environmental Growth Analysis
################################################################################
# get zc weight values from database
zcweights.lon=get.zc.weights(fdir="")
#
#  Use SMI from 1 Sept to end of Feb
#
zcweights.lon=zcweights.lon[zcweights.lon$days>-30&zcweights.lon$days<150&zcweights.lon$sitecode=="SMI",]
# Extract the pups with more than one measurement
long.id=table(zcweights.lon$AnimalID)
long.id=names(long.id[long.id>1])
zcweights.lon=zcweights.lon[zcweights.lon$AnimalID%in%long.id,]
# Create initial and recap dataframes and exclude any with more than one initial
initial=zcweights.lon[zcweights.lon$type=="Initial",]
recaps=zcweights.lon[zcweights.lon$type=="Recap",]
long.id=table(initial$AnimalID)
long.id=names(long.id[long.id==1])
initial=initial[initial$AnimalID%in%long.id,]
zcweights.lon=rbind(initial,recaps)
# Compute growth rates
gr=sapply(split(zcweights.lon,list(zcweights.lon$AnimalID)),function(x) (diff(x$weight[order(x$days)])/diff(x$days[order(x$days)]))[1])

# Merge with data about animal
grdata=data.frame(gr=gr,AnimalID=names(gr))
grdata=merge(initial,grdata)

gr.mod=lme(gr~sex,random=list(~sex|cohort,~1|AnimalID),data=grdata)
pp=data.frame(cohort=rep(sort(unique(zcweights.lon$cohort)),2),sex=rep(c("F","M"),each=length(unique(zcweights.lon$cohort))))
pp$gr=predict(gr.mod,pp,level=1)

bootstrap.se=function(x,nreps)
{
	pmat=matrix(0,nrow=nreps,ncol=nrow(unique(data.frame(cohort=x$cohort,sex=x$sex))))
	i=0
	while(i<nreps)
	{
		cat("Bootstrap ",i,"\n")
		xsamp=lapply(split(x,list(x$cohort,x$sex)),function(x) if(nrow(x)>0) x[sample(1:nrow(x),replace=TRUE),] else NULL)
		xsamp=do.call("rbind",xsamp)
		xsamp$AnimalID=1:nrow(xsamp)
		mod=try(lme(formula(gr.mod),random=as.formula(gr.mod$call$random),data=as.data.frame(xsamp),control=lmeControl(maxIter=100, msMaxIter=100)))
		if(class(mod)!="try-error")
		{
			i=i+1
			pp=predict(mod,newdata=xsamp,level=1)
			pmat[i,]=as.vector(tapply(as.vector(pp),list(x$cohort,x$sex),mean))
		}
	}
	return(sqrt(apply(pmat,2,var)))
}
stderrors=bootstrap.se(grdata,100)
female.averages=data.frame(fit=pp$gr[pp$sex=="F"&pp$cohort>=1997],se=stderrors[as.numeric(row.names(pp[pp$sex=="F"&pp$cohort>=1997,]))])
male.averages=data.frame(fit=pp$gr[pp$sex=="M"&pp$cohort>=1997],se=stderrors[as.numeric(row.names(pp[pp$sex=="M"&pp$cohort>=1997,]))])

plot.growth.series=function (time, predictions,...)
{
	plotCI(time, predictions$fit, 1.96 * predictions$se,
			1.96 * predictions$se, xlab = "Cohort",
			ylab = paste("Average daily growth rate(kg/day)"),
			...)
	lines(time, predictions$fit)
	invisible()
}

jpeg("growth.jpg",pointsize=10,res=600)
par(lty=1)
plot.growth.series(sort(unique(pp$cohort[pp$cohort>=1997])),female.averages,ylim=c(0,.12),xaxp=c(1998,2014,7))
par(lty=2)
plot.weight.series(sort(unique(pp$cohort[pp$cohort>=1997]))+.2,male.averages,pch=2,add=TRUE,slty=1,date="1 Oct")
points(2005,0.04,pch=2)
lines(x=c(2004.75,2005.25),y=c(.040,.040),pch=2,lty=2)
points(2005,0.03,pch=1)
lines(x=c(2004.75,2005.25),y=c(0.03,0.03),pch=1,lty=1)
text(2005.4,0.04,"Males",pos=4)
text(2005.4,0.03,"Females",pos=4)
dev.off()


# Merge growth data with environmental data
zcweights.environ.lon=merge(grdata,data.frame(cohort=1975:maxyear,SST=OcttoFebAnomalies,MEI=LaggedMEIOcttoFeb[-1],
				UWI33=UWImeansOcttoFeb[1,-(1:6)],UWI36=UWImeansOcttoFeb[2,-(1:6)]))
# Create model
mod=lme(gr~sex*UWI33,random=~1|cohort,data=zcweights.environ.lon,method="ML")
summary(mod)$AIC
mod=lme(gr~sex*MEI,random=~1|cohort,data=zcweights.environ.lon,method="ML")
summary(mod)$AIC
mod=lme(gr~sex*SST,random=~1|cohort,data=zcweights.environ.lon,method="ML")
summary(mod)$AIC
mod=lme(gr~sex*SST,random=~1|cohort,data=zcweights.environ.lon)
plot(zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="F"],predict(mod)[zcweights.environ.lon$sex=="F"],pch="F",type="b",ylim=c(0,.12))
points(zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="M"],predict(mod)[zcweights.environ.lon$sex=="M"],pch="M",type="b")
lines(zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="M"],predict(mod)[zcweights.environ.lon$sex=="M"],lty=2)
# Observed and predicted values
win.graph()
par(mfrow=c(1,2))
plot(zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="F"],predict(mod)[zcweights.environ.lon$sex=="F"],pch="F",type="b",ylim=c(0,.14),ylab="Female daily growth rate (kg/day)",xlab="Cohort")
points(c(1989,1990,1995,1997:2010,2012:max(zcweights.environ.lon$cohort)),tapply(zcweights.environ.lon$gr[zcweights.environ.lon$sex=="F"],zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="F"],mean),pch="O")
plot(zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="M"],predict(mod)[zcweights.environ.lon$sex=="M"],pch="M",type="b",ylim=c(0,.14),ylab="Male daily growth rate (kg/day)",xlab="Cohort")
points(c(1989,1990,1995,1997:2010,2012:max(zcweights.environ.lon$cohort)),tapply(zcweights.environ.lon$gr[zcweights.environ.lon$sex=="M"],zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="M"],mean),pch="O")


################################################################################
#  Brand vs control at brand eval
################################################################################
zcweights.be=get.zc.weights(fdir="")
#
#  Use SMI from 1 Dec to end of Feb; use brand and unmarked only
#
zcweights.be=zcweights.be[zcweights.be$days>60&zcweights.be$days<150&zcweights.be$sitecode=="SMI"&zcweights.be$tag.type!="Tag",]
zcweights.be$tag.type=factor(zcweights.be$tag.type)
# lme model of tag.type effect
mod=lme(weight~tag.type*sex,random=~1|cohort,data=zcweights.be,method="ML")
summary(mod)
mod0=lme(weight~sex,random=~1|cohort,data=zcweights.be,method="ML")
summary(mod0)
mod1=lme(weight~tag.type+sex,random=~1|cohort,data=zcweights.be,method="ML")
summary(mod1)
anova(mod0,mod1,mod)
# kruskal test by sex
kruskal.test(weight~tag.type*cohort,data=zcweights.be[zcweights.be$sex=="F",])
kruskal.test(weight~tag.type*cohort,data=zcweights.be[zcweights.be$sex=="M",])
# paired t-test
be.fem=with(zcweights.be,tapply(weight,list(cohort,tag.type,sex),mean))[,,"F"]
be.m=with(zcweights.be,tapply(weight,list(cohort,tag.type,sex),mean))[,,"M"]
t.test(x=c(be.fem[,1],be.m[,1]),y=c(be.fem[,2],be.m[,2]),paired=TRUE)


zcweights.all=zcweights[zcweights$days>-30&zcweights$days<80&zcweights$cohort>=1997&
				 zcweights$sitecode=="SMI"&zcweights$type=="Initial",]
# Distribution of timing of weight samples
batches=apply(table(zcweights.all$cohort,zcweights.all$days),2,function(x) length(x[x>0]))
barplot(batches,xlab="Days from 1 Oct",ylab="Frequency",space=0)
abline(v=which(as.numeric(names(batches))==0),lwd=2)


par(mfrow=c(2,1))
batches=with(zcweights.all[zcweights.all$sex=="M",],aggregate(weight,list(cohort,days),mean))
names(batches)=c("Year","Days","Weight")
plot(batches$Year,batches$Weight,xlab="Year",ylab="Average male weight in batch")
batches=with(zcweights.all[zcweights.all$sex=="F",],aggregate(weight,list(cohort,days),mean))
names(batches)=c("Year","Days","Weight")
plot(batches$Year,batches$Weight,xlab="Year",ylab="Average female weight in batch")

