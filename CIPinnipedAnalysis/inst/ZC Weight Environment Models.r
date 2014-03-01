require(CIPinnipedAnalysis)
if(!exists("fdir"))fdir=""
if(!exists("nboot"))nboot=100
# use "" to use databases in Calcur installed package directory; use NULL to use default Databases directory J:/Master  or specify directory
#fdir=NULL
if(!exists("anomalies"))
{
	sdir=system.file(package="CIPinnipedAnalysis")
	source(file.path(sdir,"CreateAnomalies.r"))
}
#################################################################################
# Cross-sectional analysis
#################################################################################
# get zc weight values from database
zcweights=get.zc.weights(fdir=fdir)
#
#  exclude brand eval weights and captures in April and only use SMI
#
zcweights=zcweights[zcweights$days>-30&zcweights$days<150&
				zcweights$sitecode=="SMI",]
zcweights$batch=factor(paste(zcweights$cohort,zcweights$days))


# Create dataframe with weights and environmental variables
# use 1975 and on for weight data; remove 1972-1974
JunetoSeptAnomalies=JunetoSeptAnomalies[4:numyears]
OcttoFebAnomalies=OcttoFebAnomalies[4:numyears]
JunetoFebAnomalies=JunetoFebAnomalies[4:numyears]
zcweights.environ=merge(zcweights,data.frame(cohort=1975:lastyear,SST=JunetoSeptAnomalies,SST1=OcttoFebAnomalies,SST2=JunetoFebAnomalies,
				MEI=LaggedMEIJunetoSept[-1],MEI1=LaggedMEIOcttoFeb[-1],MEI2=LaggedMEIJunetoFeb[-1],UWI33=UWImeansJunetoSept[1,-(1:6)],UWI36=UWImeansJunetoSept[2,-(1:6)],
				UWI331=UWImeansOcttoFeb[1,-(1:6)],UWI361=UWImeansOcttoFeb[2,-(1:6)],UWI332=UWImeansJunetoFeb[1,-(1:6)],UWI362=UWImeansJunetoFeb[2,-(1:6)]))
zcweights.environ$cohort.factor=factor(ifelse(zcweights.environ$cohort<1990,0,1),labels=c("<=1989",">=1990"))
zcweights.environ$cohort=zcweights.environ$cohort-min(zcweights.environ$cohort)

# First fit a sequence of random effect models with REML to assess best random model with same fixed model
random.f=list(list(~1|cohort,~1|batch),list(~days|cohort,~1|batch),list(~sex:days|cohort,~1|batch),list(~sex*days|cohort,~1|batch),list(~-1+sex+days|cohort,~1|batch))
fixed.f=list(weight~sex*SST+sex:days+SST1:days+sex:SST1:days+UWI33+cohort+MEI)
res.environ=fitmixed(fixed.f,random.f,data=zcweights.environ) 

# Using that random model, fit a sequence of fixed effect models with ML and use AIC to assess best fixed model
random.f=list(res.environ$best.r)
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
	ZCWeight.df$female.eviron.mean.fall=female.averages$fit
	ZCWeight.df$female.environ.mean.fall..se=female.averages$se
	ZCWeight.df$male.environ.mean.fall=male.averages$fit
	ZCWeight.df$male.environ.mean.fall.se=male.averages$se


    # Plot predictions and observed
    win.graph()
    par(mfrow=c(2,1))
    with(ZCWeight.df,plot(ZCWeight.df$Year,female.eviron.mean.fall,pch="F",type="b",ylim=c(12,26)))
    with(ZCWeight.df,points(ZCWeight.df$Year,female.observed.mean.fall,pch="O"))
    with(ZCWeight.df,lines(ZCWeight.df$Year,female.observed.mean.fall,lty=2))
    with(ZCWeight.df,plot(ZCWeight.df$Year,male.environ.mean.fall,pch="M",type="b",ylim=c(12,26)))
    with(ZCWeight.df,points(ZCWeight.df$Year,male.observed.mean.fall,pch="O"))
    with(ZCWeight.df,lines(ZCWeight.df$Year,male.observed.mean.fall,lty=2))


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
	ZCWeight.df$female.eviron.winter.mean=female.averages$fit
	ZCWeight.df$female.environ.winter.mean.se=female.averages$se
	ZCWeight.df$male.environ.winter.mean=male.averages$fit
	ZCWeight.df$male.environ.winter.mean.se=male.averages$se
}
	

# Plot predictions and predictions at SST=0
#pp=zcweights.environ
#pp$days=0
#pp$SST=0
#pp=predict(zc.weight.model,newdata=pp)
#pp=tapply(as.vector(pp),list(zcweights.environ$cohort,zcweights.environ$sex),mean)
#female.averages=data.frame(fit=pp[,1],se=stderrors[1:length(pp[,1])])
#male.averages=data.frame(fit=pp[,2],se=stderrors[(length(pp[,1])+1):(2*length(pp[,1]))])
#win.graph()
#par(mfrow=c(2,1))
#with(ZCWeight.df,plot(ZCWeight.df$Year,female.eviron.mean,pch="F",type="b",ylim=c(12,26)))
#points(ZCWeight.df$Year,female.averages$fit,pch="S")
#lines(ZCWeight.df$Year,female.averages$fit,pch="S",lty=2)
#with(ZCWeight.df,plot(ZCWeight.df$Year,male.environ.mean,pch="M",type="b",ylim=c(12,26)))
#points(ZCWeight.df$Year,male.averages$fit,pch="S")
#lines(ZCWeight.df$Year,male.averages$fit,pch="S",lty=2)




