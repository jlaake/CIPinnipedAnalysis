#if fdir="" it looks for data files in the CalcurData package directory of your R library.
#if fdir=NULL it looks in databases.txt in CalcurData package directory to get the database location
#if fdir is anything else it uses the value of fdir as the directory for database.  
#The scripts check for the value of fdir and if it exists the script will not change the value; otherwise it sets it to NULL
if(!exists("fdir"))fdir=NULL
if(!exists("nboot"))nboot=100
####################################
# Set this value; be aware that all of the environmental data has to be entered through Feb of lastyear+1 
# for the growth script to work properly
lastyear=2013
####################################
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
stderrors=bootstrap.se(cuweights,nboot)
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
plot_weight.series(1975:lastyear,female.averages,main="Fur Seal Female Pups",ylim=c(ymin,ymax),date="1 Oct")
plot_weight.series(1975:lastyear,male.averages,main="Fur Seal Male Pups",ylim=c(ymin,ymax),date="1 Oct")
dev.off()
################################################################################
# Construct environmental values for models
################################################################################
if(!exists("anomalies"))
{
	sdir=system.file(package="CIPinnipedAnalysis")
	source(file.path(sdir,"CreateAnomalies.r"))
}


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
stderrors=bootstrap.se(cuweights.environ,nboot)
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
pp1=predict(cu.weight.model,newdata=pp)

pp0=predict(cu.weight.model,newdata=pp,level=0)

pp0=tapply(as.vector(pp0),list(cuweights.environ$cohort,cuweights.environ$sex),mean)
pp1=tapply(as.vector(pp1),list(cuweights.environ$cohort,cuweights.environ$sex),mean)

female.averages=data.frame(fit=pp1[,1],se=stderrors[1:length(pp1[,1])])
male.averages=data.frame(fit=pp1[,2],se=stderrors[(length(pp1[,1])+1):(2*length(pp1[,1]))])

expected.female.averages=data.frame(fit=pp0[,1],se=stderrors[1:length(pp0[,1])])
expected.male.averages=data.frame(fit=pp0[,2],se=stderrors[(length(pp0[,1])+1):(2*length(pp0[,1]))])

win.graph()
par(mfrow=c(2,1))
year.seq=1975:max(as.numeric(row.names(CUWeight.df)))
with(CUWeight.df,plot(year.seq,female.eviron.mean,pch="F",type="b",ylim=c(4,16),xlab="Year"))
points(year.seq,female.averages$fit,pch="S")
lines(year.seq,female.averages$fit,pch="S",lty=2)
with(CUWeight.df,plot(year.seq,male.environ.mean,pch="M",type="b",ylim=c(4,16),xlab="Year"))
points(year.seq,male.averages$fit,pch="S")
lines(year.seq,male.averages$fit,pch="S",lty=2)

# Plot residuals of predictions based on fixed effects only versus mixed effects

win.graph()
plot(1975:maxyear,female.averages$fit-expected.female.averages$fit,pch="F",type="b")
win.graph()
plot(1975:maxyear,male.averages$fit-expected.male.averages$fit,pch="M",type="b")

