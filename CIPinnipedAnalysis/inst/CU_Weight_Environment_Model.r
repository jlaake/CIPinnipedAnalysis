#if fdir="" it looks for data files in the CalcurData package directory of your R library.
#if fdir=NULL it looks in databases.txt in CalcurData package directory to get the database location
#if fdir is anything else it uses the value of fdir as the directory for database.  
#The scripts check for the value of fdir and if it exists the script will not change the value; otherwise it sets it to NULL
if(!exists("fdir"))fdir=NULL
if(!exists("nboot"))nboot=100
####################################
# Set this value; be aware that all of the environmental data has to be entered through Feb of lastyear+1 
# for the growth script to work properly
if(!exists("lastyear"))lastyear=2014
################################################################################
# Construct environmental values for models
################################################################################
if(!exists("anomalies"))
{
	sdir=system.file(package="CIPinnipedAnalysis")
	source(file.path(sdir,"CreateAnomalies.r"))
}


# get cu weight values from database
cuweights=get.cu.weights(fdir=fdir)
#
#  exclude brand eval weights and captures in April
#
cuweights=cuweights[cuweights$days>-30&cuweights$days<150,]

cuweights.environ=merge(cuweights,data.frame(cohort=1975:lastyear,SST=JunetoSeptAnomalies[4:numyears],SST1=OcttoFebAnomalies[4:numyears],SST2=JunetoFebAnomalies[4:numyears],
				MEI=LaggedMEIJunetoSept[-1],MEI1=LaggedMEIOcttoFeb[-1],MEI2=LaggedMEIJunetoFeb[-1],UWI33=UWImeansJunetoSept[1,-(1:6)],UWI36=UWImeansJunetoSept[2,-(1:6)],
				UWI331=UWImeansOcttoFeb[1,-(1:6)],UWI361=UWImeansOcttoFeb[2,-(1:6)],UWI332=UWImeansJunetoFeb[1,-(1:6)],UWI362=UWImeansJunetoFeb[2,-(1:6)]))
cuweights.environ$cohort.factor=factor(ifelse(cuweights.environ$cohort<1990,0,1),labels=c("<=1989",">=1990"))
cuweights.environ$SST=cuweights.environ$SST1
cuweights.environ$SST[cuweights.environ$days>90]=cuweights.environ$SST3[cuweights.environ$days>90]
cuweights.environ$MEI[cuweights.environ$days>90]=cuweights.environ$MEI1[cuweights.environ$days>90]
cuweights.environ$UWI33[cuweights.environ$days>90]=cuweights.environ$UWI331[cuweights.environ$days>90]
cuweights.environ$UWI36[cuweights.environ$days>90]=cuweights.environ$UWI361[cuweights.environ$days>90]

#############################################################################################################
# Environmental Cross-sectional analysis of weights
#
# First fit a sequence of random effect models with REML to assess best random model with same fixed model
random.f=list(list(~1|cohort),list(~days|cohort),list(~sex:days|cohort),list(~sex*days|cohort),list(~-1+sex+days|cohort))
fixed.f=list(weight~sex*SST+sex:days+SST1:days+sex:SST1:days+UWI33+cohort+MEI)
res.environ=fitmixed(fixed.f,random.f,data=cuweights.environ) 

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
res.environ=fitmixed(fixed.f,random.f,data=cuweights.environ) 

# Finally fit best fixed/random model with REML
cu.weight.model=lme(fixed=res.environ$best.f,random=res.environ$best.r,data=cuweights.environ,method="REML",control=lmeControl(opt="optim"))
print(summary(cu.weight.model))

bootstrap.se=function(x,nreps)
{
	pmat=matrix(0,nrow=nreps,ncol=nrow(unique(data.frame(cohort=x$cohort,sex=x$sex))))
	i=0
	while(i<nreps)
	{
		xsamp=lapply(split(x,list(x$sex)),function(x) if(nrow(x)>0) x[sample(1:nrow(x),replace=TRUE),] else NULL)
		xsamp=do.call("rbind",xsamp)
		mod=try(lme(formula(cu.weight.model),random=as.formula(cu.weight.model$call$random),data=as.data.frame(xsamp),control=lmeControl(opt="optim")))
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
stderrors=bootstrap.se(cuweights.environ,nboot)
# Compute predictions and construct dataframes for female and male averages with std errors
pp=cuweights.environ
pp$days=0
# predictions with random effects
pp1=predict(cu.weight.model,newdata=pp)
pp1=tapply(as.vector(pp1),list(cuweights.environ$cohort,cuweights.environ$sex),mean)
female.averages=data.frame(fit=pp1[,1],se=stderrors[1:length(pp1[,1])])
male.averages=data.frame(fit=pp1[,2],se=stderrors[(length(pp1[,1])+1):(2*length(pp1[,1]))])

if(exists("CUWeight.df"))
{
   CUWeight.df$female.environ.mean=NA
   CUWeight.df$female.environ.mean.se=NA
   CUWeight.df$male.environ.mean=NA
   CUWeight.df$male.environ.mean.se=NA
   numrecs=1:length(female.averages$fit)
   CUWeight.df$female.environ.mean[numrecs]=female.averages$fit
   CUWeight.df$female.environ.mean.se[numrecs]=female.averages$se
   CUWeight.df$male.environ.mean[numrecs]=male.averages$fit
   CUWeight.df$male.environ.mean.se[numrecs]=male.averages$se
  
   # Plot predictions and observed
   jpeg("CUEnvironObserved&Predicted.jpg",height=600,width=600,quality=100,pointsize=12)
   par(mfrow=c(2,1))
   year.seq=1975:max(as.numeric(row.names(CUWeight.df)))
   with(CUWeight.df,plot(year.seq,female.environ.mean,pch="F",type="b",ylim=c(4,16),xlab="Year"))
   with(CUWeight.df,points(year.seq,female.observed.mean,pch="O"))
   with(CUWeight.df,lines(year.seq,female.observed.mean,lty=2))
   with(CUWeight.df,plot(year.seq,male.environ.mean,pch="M",type="b",ylim=c(4,16),xlab="Year"))
   with(CUWeight.df,points(year.seq,male.observed.mean,pch="O"))
   with(CUWeight.df,lines(year.seq,male.observed.mean,lty=2))
   dev.off()

   # predictions with fixed effects only
   pp0=predict(cu.weight.model,newdata=pp,level=0)
   pp0=tapply(as.vector(pp0),list(cuweights.environ$cohort,cuweights.environ$sex),mean)
   expected.female.averages=data.frame(fit=pp0[,1],se=stderrors[1:length(pp0[,1])])
   expected.male.averages=data.frame(fit=pp0[,2],se=stderrors[(length(pp0[,1])+1):(2*length(pp0[,1]))])

   jpeg("CUEnvironFixedEffectFemalePredictions&Observed.jpg",height=600,width=600,quality=100,pointsize=12)
   par(mfrow=c(2,1))
   with(CUWeight.df,plot(year.seq[numrecs],expected.female.averages$fit,pch="F",type="b",ylim=c(4,16),xlab="Year"))
   points(year.seq,CUWeight.df$female.observed.mean,pch="O")
   lines(year.seq[numrecs],expected.female.averages$fit,pch="S",lty=2)
   with(CUWeight.df,plot(year.seq[numrecs],expected.male.averages$fit,pch="F",type="b",ylim=c(4,16),xlab="Year"))
   points(year.seq,CUWeight.df$male.observed.mean,pch="O")
   lines(year.seq[numrecs],expected.male.averages$fit,pch="S",lty=2)
   dev.off()

   # Plot residuals of predictions based on fixed effects only versus mixed effects
   jpeg("CUEnvironFixedEffectResiduals.jpg",height=600,width=600,quality=100,pointsize=12)
   par(mfrow=c(2,1))
   plot(year.seq[numrecs],female.averages$fit-expected.female.averages$fit,pch="F",type="b")
   plot(year.seq[numrecs],male.averages$fit-expected.male.averages$fit,pch="M",type="b")
   dev.off()
}
