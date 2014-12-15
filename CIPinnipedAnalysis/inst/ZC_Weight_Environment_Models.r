#if fdir="" it looks for data files in the CalcurData package directory of your R library.
#if fdir=NULL it looks in databases.txt in CalcurData package directory to get the database location
#if fdir is anything else it uses the value of fdir as the directory for database.  
#The scripts check for the value of fdir and if it exists the script will not change the value; otherwise it sets it to NULL
if(!exists("fdir"))fdir=NULL
require(CIPinnipedAnalysis)
if(!exists("lastyear"))lastyear=2013
if(!exists("nboot"))nboot=100
if(!exists("anomalies"))
{
	sdir=system.file(package="CIPinnipedAnalysis")
	source(file.path(sdir,"CreateAnomalies.r"))
}
if(!exists("use.calcofi"))use.calcofi=FALSE
#################################################################################
# Cross-sectional analysis
#################################################################################
# get zc weight values from database
zcweights=get.zc.weights(fdir=fdir)
zcweights=zcweights[zcweights$cohort<=lastyear,]
if(use.calcofi)zcweights=zcweights[zcweights$cohort>=1984,]
#
#  exclude brand eval weights and captures in April and only use SMI
#
zcweights=zcweights[zcweights$days>-30&zcweights$days<150&
				zcweights$sitecode=="SMI",]
zcweights$batch=factor(paste(zcweights$cohort,zcweights$days))

if(use.calcofi)
{
	data(calcofi)
#	calcofi=calcofi[as.numeric(calcofi$Station)<=3,]
	calcofi=sapply(calcofi[,-(1:2)],function(x) tapply(x,calcofi$Year,mean))
	calcofi=t(t(calcofi)-colMeans(calcofi))
	zcweights.environ=merge(zcweights,cbind(calcofi,data.frame(cohort=1984:lastyear,SST=JunetoSeptAnomalies[13:numyears],SST1=OcttoFebAnomalies[13:numyears],SST2=JunetoFebAnomalies[13:numyears],
							MEI=LaggedMEIJunetoSept[-(1:10)],MEI1=LaggedMEIOcttoFeb[-(1:10)],MEI2=LaggedMEIJunetoFeb[-(1:10)],UWI33=UWImeansJunetoSept[1,-(1:15)],UWI36=UWImeansJunetoSept[2,-(1:15)],
							UWI331=UWImeansOcttoFeb[1,-(1:15)],UWI361=UWImeansOcttoFeb[2,-(1:15)],UWI332=UWImeansJunetoFeb[1,-(1:15)],UWI362=UWImeansJunetoFeb[2,-(1:15)])))
	fixed.f=list(weight~sex*SST+sex:days+SST1:days+sex:SST1:days+cohort.factor,
			weight~sex*SST+sex:days+SST1:days+sex:SST1:days+cohort,
			weight~sex*SST+sex:days+SST1:days+sex:SST1:days,
			weight~sex*SST+sex:days+SST1:days+cohort.factor,
			weight~sex*SST+sex:days+SST1:days+cohort,
			weight~sex*SST+sex:days+SST1:days,
			weight~sex*UWI36+sex:days+UWI361:days+sex:UWI361:days+cohort.factor,
			weight~sex*UWI36+sex:days+UWI361:days+sex:UWI361:days+cohort,
			weight~sex*UWI36+sex:days+UWI361:days+sex:UWI361:days,
			weight~sex*UWI36+sex:days+UWI361:days+cohort.factor,
			weight~sex*UWI36+sex:days+UWI361:days+cohort,
			weight~sex*UWI36+sex:days+UWI361:days,
			weight~sex*MEI+sex:days+MEI1:days+sex:MEI1:days+cohort.factor,
			weight~sex*MEI+sex:days+MEI1:days+sex:MEI1:days+cohort,
			weight~sex*MEI+sex:days+MEI1:days+sex:MEI1:days,
			weight~sex*MEI+sex:days+MEI1:days+cohort.factor,
			weight~sex*MEI+sex:days+MEI1:days+cohort,
			weight~sex*MEI+sex:days,
			weight~sex*dynamic_height_0_500m+sex:days+SST1:days+sex:SST1:days+cohort.factor,
			weight~sex*dynamic_height_0_500m+sex:days+SST1:days+sex:SST1:days+cohort,
			weight~sex*dynamic_height_0_500m+sex:days+SST1:days+sex:SST1:days,
			weight~sex*dynamic_height_0_500m+sex:days+SST1:days+cohort.factor,
			weight~sex*dynamic_height_0_500m+sex:days+SST1:days+cohort,
			weight~sex*dynamic_height_0_500m+sex:days,
			weight~sex*R_SIGMA_75m+sex:days+SST1:days+sex:SST1:days+cohort.factor,
			weight~sex*R_SIGMA_75m+sex:days+SST1:days+sex:SST1:days+cohort,
			weight~sex*R_SIGMA_75m+sex:days+SST1:days+sex:SST1:days,
			weight~sex*R_SIGMA_75m+sex:days+SST1:days+cohort.factor,
			weight~sex*R_SIGMA_75m+sex:days+SST1:days+cohort,
			weight~sex*R_SIGMA_75m+sex:days,
			weight~sex*R_POTEMP_25m+sex:days+R_POTEMP_25m:days+sex:R_POTEMP_25m:days+cohort.factor,
			weight~sex*R_POTEMP_25m+sex:days+R_POTEMP_25m:days+sex:R_POTEMP_25m:days+cohort,
			weight~sex*R_POTEMP_25m+sex:days+R_POTEMP_25m:days+sex:R_POTEMP_25m:days,
			weight~sex*R_POTEMP_25m+sex:days+R_POTEMP_25m:days+cohort.factor,
			weight~sex*R_POTEMP_25m+sex:days+R_POTEMP_25m:days+cohort,
			weight~sex*R_POTEMP_25m+sex:days+R_POTEMP_25m:days,
			weight~sex*R_SIGMA_75m+sex:days+SST1:days+cohort+stratification,
			weight~sex*R_SIGMA_75m+sex:days+SST1:days+cohort+pycnocline_depth,
			weight~sex*R_SIGMA_75m+sex:days+SST1:days+cohort+R_POTEMP_25m,
			weight~sex*R_SIGMA_75m+sex:days+SST1:days+cohort+R_POTEMP_75m,
			weight~sex*R_SIGMA_75m+sex:days+SST1:days+cohort+R_SIGMA_25m,
			weight~sex*R_SIGMA_75m+sex:days+SST1:days+cohort+dynamic_height_0_500m,
			weight~sex*R_SIGMA_75m+sex:days+SST1:days+cohort+R_O2_25m,
			weight~sex*R_SIGMA_75m+sex:days+SST1:days+cohort+R_O2_75m,
			weight~sex*R_SIGMA_75m+sex:days+SST1:days+cohort+R_NO3_25m,
			weight~sex*R_SIGMA_75m+sex:days+SST1:days+cohort+R_NO3_75m)
} else {
	zcweights.environ=merge(zcweights,data.frame(cohort=1975:lastyear,SST=JunetoSeptAnomalies[4:numyears],SST1=OcttoFebAnomalies[4:numyears],SST2=JunetoFebAnomalies[4:numyears],
					MEI=LaggedMEIJunetoSept[-1],MEI1=LaggedMEIOcttoFeb[-1],MEI2=LaggedMEIJunetoFeb[-1],UWI33=UWImeansJunetoSept[1,-(1:6)],UWI36=UWImeansJunetoSept[2,-(1:6)],
					UWI331=UWImeansOcttoFeb[1,-(1:6)],UWI361=UWImeansOcttoFeb[2,-(1:6)],UWI332=UWImeansJunetoFeb[1,-(1:6)],UWI362=UWImeansJunetoFeb[2,-(1:6)]))
	fixed.f=list(weight~sex*SST+sex:days+SST1:days+sex:SST1:days+cohort.factor,
			weight~sex*SST+sex:days+SST1:days+sex:SST1:days+cohort,
			weight~sex*SST+sex:days+SST1:days+sex:SST1:days,
			weight~sex*SST+sex:days+SST1:days+cohort.factor,
			weight~sex*SST+sex:days+SST1:days+cohort,
			weight~sex*SST+sex:days+SST1:days,
			weight~sex*UWI36+sex:days+UWI361:days+sex:UWI361:days+cohort.factor,
			weight~sex*UWI36+sex:days+UWI361:days+sex:UWI361:days+cohort,
			weight~sex*UWI36+sex:days+UWI361:days+sex:UWI361:days,
			weight~sex*UWI36+sex:days+UWI361:days+cohort.factor,
			weight~sex*UWI36+sex:days+UWI361:days+cohort,
			weight~sex*UWI36+sex:days+UWI361:days,
			weight~sex*MEI+sex:days+MEI1:days+sex:MEI1:days+cohort.factor,
			weight~sex*MEI+sex:days+MEI1:days+sex:MEI1:days+cohort,
			weight~sex*MEI+sex:days+MEI1:days+sex:MEI1:days,
			weight~sex*MEI+sex:days+MEI1:days+cohort.factor,
			weight~sex*MEI+sex:days+MEI1:days+cohort,
			weight~sex*MEI+sex:days)
}

zcweights.environ$cohort.factor=factor(ifelse(zcweights.environ$cohort<1990,0,1),labels=c("<=1989",">=1990"))
zcweights.environ$cohort=zcweights.environ$cohort-min(zcweights.environ$cohort)

# First fit a sequence of random effect models with REML to assess best random model with same fixed model
random.f=list(list(~1|cohort,~1|batch),list(~days|cohort,~1|batch),list(~sex:days|cohort,~1|batch),list(~sex*days|cohort,~1|batch),list(~-1+sex+days|cohort,~1|batch))
fixed.f1=list(weight~sex*SST+sex:days+SST1:days+sex:SST1:days+UWI33+cohort+MEI)
res.environ=fitmixed(fixed.f1,random.f,data=zcweights.environ) 

# Using that random model, fit a sequence of fixed effect models with ML and use AIC to assess best fixed model
random.f=list(res.environ$best.r)
res.environ=fitmixed(fixed.f,random.f,data=zcweights.environ) 

# Finally fit best fixed/random model with REML
zc.weight.model=lme(fixed=res.environ$best.f,random=res.environ$best.r,data=zcweights.environ,method="REML",control=lmeControl(opt="optim"))
print(summary(zc.weight.model))


bootstrap.se=function(x,nreps,days=0)
{
	pmat=matrix(0,nrow=nreps,ncol=nrow(unique(data.frame(cohort=x$cohort,sex=x$sex))))
	i=0
	while(i<nreps)
	{
		cat("Bootstrap ",i,"\n")
		xsamp=lapply(split(x,list(x$batch,x$sex)),function(x) if(nrow(x)>0) x[sample(1:nrow(x),replace=TRUE),] else NULL)
		xsamp=do.call("rbind",xsamp)
		mod=try(lme(fixed=res.environ$best.f,random=res.environ$best.r,data=as.data.frame(xsamp),method="REML",control=lmeControl(opt="optim")))
		if(class(mod)!="try-error")
		{
			i=i+1
			xsamp$days=days
			pp=predict(mod,newdata=xsamp)
			pmat[i,]=as.vector(tapply(as.vector(pp),list(x$cohort,x$sex),mean))
		}
	}
	return(sqrt(apply(pmat,2,var)))
}

#################  1 Oct Predictions ####################
# use 100 reps to compute std error
stderrors=bootstrap.se(zcweights.environ,nboot,days=0)
# Compute fall 1 Oct predictions and construct dataframes for female and male averages with std errors
pp=zcweights.environ
pp$days=0
# predictions at 1 Oct with random effects
pp1=predict(zc.weight.model,newdata=pp)
# predictions at 1 Oct with fixed effects only
pp0=predict(zc.weight.model,newdata=pp,level=0)
# compute mean values which essentially acts as unique
pp0=tapply(as.vector(pp0),list(zcweights.environ$cohort,zcweights.environ$sex),mean)
pp1=tapply(as.vector(pp1),list(zcweights.environ$cohort,zcweights.environ$sex),mean)
# create list with estimates and std errors for averages
female.averages=data.frame(fit=pp1[,1],se=stderrors[1:length(pp1[,1])])
male.averages=data.frame(fit=pp1[,2],se=stderrors[(length(pp1[,1])+1):(2*length(pp1[,1]))])
# create list with estimates and std errors for expected averages
expected.female.averages=data.frame(fit=pp0[,1],se=stderrors[1:length(pp0[,1])])
expected.male.averages=data.frame(fit=pp0[,2],se=stderrors[(length(pp0[,1])+1):(2*length(pp0[,1]))])

if(exists("ZCWeight.df"))
{
	ZCWeight.df$female.environ.mean.fall=NA
	ZCWeight.df$female.environ.mean.fall.se=NA
	ZCWeight.df$male.environ.mean.fall=NA
	ZCWeight.df$male.environ.mean.fall.se=NA
	add=0
	if(use.calcofi)add=9
	ZCWeight.df$female.environ.mean.fall[add+1:length(female.averages$fit)]=female.averages$fit
	ZCWeight.df$female.environ.mean.fall.se[add+1:length(female.averages$fit)]=female.averages$se
	ZCWeight.df$male.environ.mean.fall[add+1:length(female.averages$fit)]=male.averages$fit
	ZCWeight.df$male.environ.mean.fall.se[add+1:length(female.averages$fit)]=male.averages$se
	
	ZCWeight.df$female.environ.mean.fall.fixed[add+1:length(female.averages$fit)]=expected.female.averages$fit
	ZCWeight.df$male.environ.mean.fall.fixed[add+1:length(female.averages$fit)]=expected.male.averages$fit
	
    # Plot predictions and observed
	jpeg("ZCEnvironObserved&Predicted.jpg",height=600,width=600,quality=100,pointsize=12)
	par(mfrow=c(2,1))
	with(ZCWeight.df,plot(ZCWeight.df$Year,female.environ.mean.fall,pch="F",type="b",ylim=c(12,26)))
	with(ZCWeight.df,points(ZCWeight.df$Year,female.observed.mean.fall,pch="O"))
	with(ZCWeight.df,lines(ZCWeight.df$Year,female.observed.mean.fall,lty=2))
	with(ZCWeight.df,plot(ZCWeight.df$Year,male.environ.mean.fall,pch="M",type="b",ylim=c(12,26)))
	with(ZCWeight.df,points(ZCWeight.df$Year,male.observed.mean.fall,pch="O"))
	with(ZCWeight.df,lines(ZCWeight.df$Year,male.observed.mean.fall,lty=2))
	dev.off()
	
	# Plot residuals of predictions based on fixed effects only versus mixed effects
	jpeg("ZCEnvironFixedEffectResiduals.jpg",height=600,width=600,quality=100,pointsize=12)
	par(mfrow=c(2,1))
	plot(ZCWeight.df$Year[add+1:length(male.averages$fit)],male.averages$fit-expected.male.averages$fit,pch="M",type="b")
	abline(0,0)
	plot(ZCWeight.df$Year[add+1:length(male.averages$fit)],female.averages$fit-expected.female.averages$fit,pch="F",type="b",xlab="Year",ylab="Observed- model predicted weight")
	abline(0,0)
	dev.off()
	
	jpeg("ZCEnvironFixedEffectFemalePredictions&Observed.jpg",height=600,width=600,quality=100,pointsize=12)
	plot(ZCWeight.df$Year[add+1:length(male.averages$fit)],female.averages$fit,type="b",ylim=c(12,22),xlab="Year",ylab="Female pup weight (kg)")
	lines(ZCWeight.df$Year[add+1:length(male.averages$fit)],expected.female.averages$fit,type="b",lty=2,pch=2) 
	points(2000,14,pch=1)
	points(2000,13,pch=2)
	text(2000,14,"1 Oct mean weight",pos=4)
	text(2000,13,"Model predicted mean weight",pos=4)
	dev.off()
}

#################  1 Feb Predictions ####################
# use 100 reps to compute std error
stderrors=bootstrap.se(zcweights.environ,nboot,days=123)
# Compute fall 1 Feb predictions and construct dataframes for female and male averages with std errors
pp=zcweights.environ
pp$days=123
# predictions at 1 Oct with random effects
pp1=predict(zc.weight.model,newdata=pp)
# predictions at 1 Oct with fixed effects only
pp0=predict(zc.weight.model,newdata=pp,level=0)
# compute mean values which essentially acts as unique
pp0=tapply(as.vector(pp0),list(zcweights.environ$cohort,zcweights.environ$sex),mean)
pp1=tapply(as.vector(pp1),list(zcweights.environ$cohort,zcweights.environ$sex),mean)
# create list with estimates and std errors for averages
female.averages=data.frame(fit=pp1[,1],se=stderrors[1:length(pp1[,1])])
male.averages=data.frame(fit=pp1[,2],se=stderrors[(length(pp1[,1])+1):(2*length(pp1[,1]))])
# create list with estimates and std errors for expected averages
expected.female.averages=data.frame(fit=pp0[,1],se=stderrors[1:length(pp0[,1])])
expected.male.averages=data.frame(fit=pp0[,2],se=stderrors[(length(pp0[,1])+1):(2*length(pp0[,1]))])
if(exists("ZCWeight.df"))
{
	ZCWeight.df$female.environ.mean.winter=NA
	ZCWeight.df$female.environ.mean.winter.se=NA
	ZCWeight.df$male.environ.mean.winter=NA
	ZCWeight.df$male.environ.mean.winter.se=NA
	
	ZCWeight.df$female.environ.mean.winter[add+1:length(female.averages$fit)]=female.averages$fit
	ZCWeight.df$female.environ.mean.winter.se[add+1:length(female.averages$fit)]=female.averages$se
	ZCWeight.df$male.environ.mean.winter[add+1:length(female.averages$fit)]=male.averages$fit
	ZCWeight.df$male.environ.mean.winter.se[add+1:length(female.averages$fit)]=male.averages$se
}
	

	
		
		
		
		
