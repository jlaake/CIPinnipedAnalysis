\name{get.cu.weights}
\alias{get.cu.weights}
\title{Extracts weights for Cu pups
: from various tables in the cutagnew.mdb ACCESS
database and returns a dataframe with ancillary fields added.}
\usage{
  get.cu.weights(fdir = "", ENYears = c(1976, 1983, 1984,
  1987, 1992, 1997, 1998))
}
\arguments{
  \item{fdir}{directory for cutagnew.mdb}

  \item{ENYears}{values of years that will be flagged as
  ENSO years}
}
\value{
  dataframe of weights for Cu
}
\description{
  Extracts weights for Cu pups : from various tables in the
  cutagnew.mdb ACCESS database and returns a dataframe with
  ancillary fields added.
}
\examples{
# fdir should be set to location of  BrandMaster.mdb
# edir should be set to location of environmental.data.mdb
# each defaults below to current working directory
# Use "J:/" if it should be Calcurr/Databases
# edir="J:/"
# fdir="J:/"
fdir=""
edir=""
plot.weight.series=function (time, predictions, ...)
{
plotCI(time, predictions$fit, predictions$fit - 1.96 * predictions$se,
predictions$fit + 1.96 * predictions$se, xlab = "Year",
ylab = "Predicted average weight (kg) 1 Oct", ...)
lines(time, predictions$fit)
invisible()
}

#################################################################################
# Cross-sectional analysis
#################################################################################
# get cu weight values from database
cuweights=get.cu.weights(fdir=fdir)
#
#  exclude brand eval weights and captures in April
#
cuweights.all=cuweights[cuweights$days>-30&cuweights$days<150,]
#
# compute observed averages - using data within 1 Sept to 15 Nov
#
cuweights=cuweights.all[cuweights.all$days>=-31&cuweights.all$days<=45,]
female.observed=sapply(split(cuweights$weight[cuweights$sex=="F"],
factor(cuweights$cohort[cuweights$sex=="F"],levels=levels(factor(cuweights.all$cohort)))),mean)
male.observed=sapply(split(cuweights$weight[cuweights$sex=="M"],
factor(cuweights$cohort[cuweights$sex=="M"],levels=levels(factor(cuweights.all$cohort)))),mean)
cuweights=cuweights.all
#
#  fit growth model
#
require(nlme)
cu.weight.model=lme(weight~sex*days,random=~days|cohort,data=cuweights)
print(summary(cu.weight.model))
# define bootstrap function to compute std error for predicted sex-cohort means
bootstrap.se=function(x,nreps)
{
pmat=matrix(0,nrow=nreps,ncol=nrow(unique(data.frame(cohort=x$cohort,sex=x$sex))))
i=0
while(i<nreps)
{
xsamp=lapply(split(x,list(x$sex)),function(x) if(nrow(x)>0) x[sample(1:nrow(x),replace=TRUE),] else NULL)
xsamp=do.call("rbind",xsamp)
mod=try(lme(formula(cu.weight.model),random=as.formula(cu.weight.model$call$random),data=as.data.frame(xsamp)))
if(class(mod)!="try-error")
{
i=i+1
xsamp=data.frame(days=0,cohort=rep(sort(unique(x$cohort)),2),sex=rep(c("F","M"),each=length(unique(x$cohort))))
pp=predict(mod,newdata=xsamp,level=1)
pmat[i,]=unique(data.frame(predict=as.vector(pp),cohort=xsamp$cohort,sex=xsamp$sex))$predict
}
}
return(sqrt(apply(pmat,2,var)))
}
# use 100 reps to compute std error
stderrors=bootstrap.se(cuweights,100)
# Compute predictions and construct dataframes for female and male averages with std errors
pp=data.frame(days=0,cohort=rep(sort(unique(cuweights$cohort)),2),sex=rep(c("F","M"),each=length(unique(cuweights$cohort))))
pp$predict=predict(cu.weight.model,newdata=pp,level=1)
female.averages=data.frame(fit=pp$predict[pp$sex=="F"],se=stderrors[as.numeric(row.names(pp[pp$sex=="F",]))])
male.averages=data.frame(fit=pp$predict[pp$sex=="M"],se=stderrors[as.numeric(row.names(pp[pp$sex=="M",]))])
# create dataframe with normal conf interval
cu.female.averages=female.averages
ymin=min(c(female.averages$fit-1.96*female.averages$se, male.averages$fit-1.96*male.averages$se))*.9
ymax=max(c(female.averages$fit+1.96*female.averages$se, male.averages$fit+1.96*male.averages$se))*1.1
cu.male.averages=male.averages
CUWeight.df=data.frame(female.observed.mean=female.observed,
female.adjusted.mean=cu.female.averages$fit,female.adjusted.mean.se=cu.female.averages$se,
male.observed.mean=male.observed,
male.adjusted.mean=cu.male.averages$fit,male.adjusted.mean.se=cu.male.averages$se)
options(width=150)
print(CUWeight.df)
pdf("CUPredictedWeights.pdf",pointsize=10)
par(mfrow=c(2,1))
maxyear=max(cuweights$cohort)
plot.weight.series(1975:maxyear,female.averages,main="Fur Seal Female Pups",ylim=c(ymin,ymax))
plot.weight.series(1975:maxyear,male.averages,main="Fur Seal Male Pups",ylim=c(ymin,ymax))
dev.off()
################################################################################
# Construct environmental values for models
################################################################################
# Get SST Anomalies at 8 locations using 1994:1996 and 1998:2008 for averages
anomalies=create.SST.anomalies(c(1994:1996,1998:2008),fdir=edir)
if(edir=="")edir=system.file(package="CIPinnipedAnalysis")
dir=file.path(edir,"environmental.data.mdb")
require(RODBC)
connection=odbcConnectAccess(dir)
# Use Locations 1-5 (ESB,WSB,PtArg,PtSM,PtSL) for pup weight predictions
SSTAnomalies=t(apply(anomalies[,,1:5],c(2,1),mean,na.rm=TRUE))
maxyear= max(as.numeric(row.names(SSTAnomalies)))
minyear= min(as.numeric(row.names(SSTAnomalies)))
numyears=maxyear-minyear+1
# Compute averages across month grouping so use with period specific growth rate
JantoMayAnomalies=rowMeans(SSTAnomalies[,c("Jan","Feb","Mar","Apr","May")])
OcttoFebAnomalies=cbind(SSTAnomalies[,c("Oct","Nov","Dec")],rbind(SSTAnomalies[2:nrow(SSTAnomalies),c("Jan","Feb")],data.frame(Jan=NA,Feb=NA)))
OcttoFebAnomalies[is.nan(OcttoFebAnomalies)]=NA
OcttoFebAnomalies=rowMeans(OcttoFebAnomalies,na.rm=TRUE)
JunetoSeptAnomalies=rowMeans(SSTAnomalies[,c("June","July","Aug","Sept")],na.rm=TRUE)
JunetoSeptAnomalies[is.nan(JunetoSeptAnomalies)]=NA
OcttoDecAnomalies=SSTAnomalies[,c("Oct","Nov","Dec")]
OcttoDecAnomalies[is.nan(OcttoDecAnomalies)]=NA
OcttoDecAnomalies=rowMeans(OcttoDecAnomalies,na.rm=TRUE)

JunetoFebAnomalies=rowMeans(SSTAnomalies[,c("June","July","Aug","Sept","Oct","Nov","Dec")])
JunetoFebAnomalies=cbind(SSTAnomalies[,c("June","July","Aug","Sept","Oct","Nov","Dec")],rbind(SSTAnomalies[2:nrow(SSTAnomalies),c("Jan","Feb")],data.frame(Jan=NA,Feb=NA)))
JunetoFebAnomalies[is.nan(JunetoFebAnomalies)]=NA
JunetoFebAnomalies=rowMeans(JunetoFebAnomalies,na.rm=TRUE)

x=rbind(data.frame(Year=minyear:maxyear,Season=rep("Spring",numyears),SSTAnomaly=JantoMayAnomalies),
data.frame(Year=minyear:maxyear,Season=rep("Summer",numyears),SSTAnomaly=JunetoSeptAnomalies),
data.frame(Year=minyear:maxyear,Season=rep("Fall",numyears),SSTAnomaly=OcttoDecAnomalies) )
x$Season=factor(x$Season,levels=c("Spring","Summer","Fall"))
x=x[order(x$Year,x$Season),]
MEI=sqlFetch(connection,"MEI")
UWI=sqlFetch(connection,"UWIAnomaly")
UWI=UWI[order(UWI$Year,UWI$Month),]
UWImeansJunetoSept=with(UWI[UWI$Month\%in\%6:9&UWI$Year>=1975,], tapply(UWIAnomaly,list(Location,Year),mean))

UWIOcttoDec=with(UWI[UWI$Month\%in\%10:12,], tapply(UWIAnomaly,list(Month,Year,Location),mean))
UWIJantoFeb=with(UWI[UWI$Month\%in\%1:2,], tapply(UWIAnomaly,list(Month,Year,Location),mean))

UWImeansOcttoFeb=NULL
for(i in 1:2)
UWImeansOcttoFeb=rbind(UWImeansOcttoFeb,colMeans(rbind(UWIOcttoDec[,,i],UWIJantoFeb[,,i]),na.rm=TRUE))

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
cat("\\nlag = ",lag,"cor = ",MEIcor[lag+1])
}
lag=which(MEIcor==max(MEIcor))-1
#
# This MEI is comparable to JunetoSeptSST with defined lag
#
average.MEI=function(x,months)return(tapply(x$MEI[x$Month\%in\%months],x$Year[x$Month\%in\%months],mean))
LaggedMEIJunetoSept=average.MEI(MEI,(6:9-lag))
LaggedMEIOcttoFeb=average.MEI(MEI,8:12) # assumes 2 month lag to avoid Dec/Jan break

JunetoSeptAnomalies=JunetoSeptAnomalies[4:numyears]
OcttoFebAnomalies=OcttoFebAnomalies[4:numyears]
JunetoFebAnomalies=JunetoFebAnomalies[4:numyears]
cuweights.environ=merge(cuweights,data.frame(cohort=1975:maxyear,SST1=JunetoSeptAnomalies,SST2=OcttoFebAnomalies,SST3=JunetoFebAnomalies,
MEI=LaggedMEIJunetoSept[-1],MEI1=LaggedMEIOcttoFeb[-1],UWI33=UWImeansJunetoSept[1,],UWI36=UWImeansJunetoSept[2,],
UWI331=UWImeansOcttoFeb[1,-(1:6)],UWI361=UWImeansOcttoFeb[2,-(1:6)]))
cuweights.environ$SST=cuweights.environ$SST1
cuweights.environ$SST[cuweights.environ$days>90]=cuweights.environ$SST3[cuweights.environ$days>90]
cuweights.environ$MEI[cuweights.environ$days>90]=cuweights.environ$MEI1[cuweights.environ$days>90]
cuweights.environ$UWI33[cuweights.environ$days>90]=cuweights.environ$UWI331[cuweights.environ$days>90]
cuweights.environ$UWI36[cuweights.environ$days>90]=cuweights.environ$UWI361[cuweights.environ$days>90]

#############################################################################################################
# Environmental Cross-sectional analysis of weights
#
# Evaluate best random effect model with most complex fixed-effect model
cu.weight.model=lme(weight~sex*days*SST+UWI33+cohort+MEI,random=~days|cohort,data=cuweights.environ)
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex*days*SST+UWI33+cohort+MEI,random=list(~1|cohort),data=cuweights.environ)
summary(cu.weight.model)$AIC
#
# Next evaluate sequence of fixed-effect models with the chosen random effect model
#
# SST
cu.weight.model=lme(weight~sex*days*SST+cohort,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex*days*SST,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex*days+days*SST+sex*SST+cohort,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex*days+days*SST+sex*SST,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex*SST+days+cohort,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex*SST+days,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex+days+SST+cohort,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex+days+SST,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
# UWI
cu.weight.model=lme(weight~sex*days*UWI33+cohort,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex*days*UWI33,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex*days+days*UWI33+sex*UWI33+cohort,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex*days+days*UWI33+sex*UWI33,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex+days+UWI33,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex+days+UWI33,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
# MEI
cu.weight.model=lme(weight~sex*days*MEI+cohort,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex*days*MEI,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex*days+days*MEI+sex*MEI+cohort,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex*days+days*MEI+sex*MEI,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex+days+MEI,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
cu.weight.model=lme(weight~sex+days+MEI,random=~days|cohort,data=cuweights.environ,method="ML")
summary(cu.weight.model)$AIC
# Get predicted means and compute std errors and add to CUWeight.df table
best.cu.weight.model=lme(weight~sex*SST+days+cohort,random=~days|cohort,data=cuweights.environ)
summary(best.cu.weight.model)
cu.weight.model=best.cu.weight.model
bootstrap.se=function(x,nreps)
{
pmat=matrix(0,nrow=nreps,ncol=nrow(unique(data.frame(cohort=x$cohort,sex=x$sex))))
i=0
while(i<nreps)
{
xsamp=lapply(split(x,list(x$sex)),function(x) if(nrow(x)>0) x[sample(1:nrow(x),replace=TRUE),] else NULL)
xsamp=do.call("rbind",xsamp)
mod=try(lme(formula(cu.weight.model),random=as.formula(cu.weight.model$call$random),data=as.data.frame(xsamp)))
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
stderrors=bootstrap.se(cuweights.environ,100)
# Compute predictions and construct dataframes for female and male averages with std errors
pp=cuweights.environ
pp$days=0
pp=predict(cu.weight.model,newdata=pp)
pp=tapply(as.vector(pp),list(cuweights.environ$cohort,cuweights.environ$sex),mean)
female.averages=data.frame(fit=pp[,1],se=stderrors[1:length(pp[,1])])
male.averages=data.frame(fit=pp[,2],se=stderrors[(length(pp[,1])+1):(2*length(pp[,1]))])
CUWeight.df$female.eviron.mean=female.averages$fit
CUWeight.df$female.environ.mean.se=female.averages$se
CUWeight.df$male.environ.mean=male.averages$fit
CUWeight.df$male.environ.mean.se=male.averages$se
# Plot predictions and observed
win.graph()
par(mfrow=c(2,1))
with(CUWeight.df,plot(1975:maxyear,female.eviron.mean,pch="F",type="b",ylim=c(4,16)))
with(CUWeight.df,points(1975:maxyear,female.observed.mean,pch="O"))
with(CUWeight.df,lines(1975:maxyear,female.observed.mean,lty=2))
with(CUWeight.df,plot(1975:maxyear,male.environ.mean,pch="M",type="b",ylim=c(4,16)))
with(CUWeight.df,points(1975:maxyear,male.observed.mean,pch="O"))
with(CUWeight.df,lines(1975:maxyear,male.observed.mean,lty=2))

# Plot predictions and predictions at SST=0
pp=cuweights.environ
pp$days=0
pp$SST=0
pp=predict(cu.weight.model,newdata=pp)
pp=tapply(as.vector(pp),list(cuweights.environ$cohort,cuweights.environ$sex),mean)
female.averages=data.frame(fit=pp[,1],se=stderrors[1:length(pp[,1])])
male.averages=data.frame(fit=pp[,2],se=stderrors[(length(pp[,1])+1):(2*length(pp[,1]))])
win.graph()
par(mfrow=c(2,1))
year.seq=1975:max(as.numeric(row.names(CUWeight.df)))
with(CUWeight.df,plot(year.seq,female.eviron.mean,pch="F",type="b",ylim=c(4,16),xlab="Year"))
points(year.seq,female.averages$fit,pch="S")
lines(year.seq,female.averages$fit,pch="S",lty=2)
with(CUWeight.df,plot(year.seq,male.environ.mean,pch="M",type="b",ylim=c(4,16),xlab="Year"))
points(year.seq,male.averages$fit,pch="S")
lines(year.seq,male.averages$fit,pch="S",lty=2)
}
\author{
  Jeff Laake
}

