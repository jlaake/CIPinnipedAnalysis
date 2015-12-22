#if fdir="" it looks for data files in the CalcurData package directory of your R library.
#if fdir=NULL it looks in databases.txt in CalcurData package directory to get the database location
#if fdir is anything else it uses the value of fdir as the directory for database.  
#The scripts check for the value of fdir and if it exists the script will not change the value; otherwise it sets it to NULL
if(!exists("fdir"))fdir=NULL
require(CIPinnipedAnalysis)
if(!exists("nboot"))nboot=100
if(!exists("lastyear"))lastyear=2014
#################################################################################
# Cross-sectional analysis
#################################################################################
# get zc weight values from database
zcweights=get.zc.weights(fdir=fdir)
zcweights=zcweights[zcweights$cohort<=lastyear,]
#
#  exclude brand eval weights and captures in April and only use SMI
#
zcweights.all=zcweights[zcweights$days>-30&zcweights$days<150&zcweights$sitecode=="SMI",]
#
# compute observed fall averages - using data within 1 Sept to 15 Nov
#
zcweights=zcweights.all[zcweights.all$days>=-31&zcweights.all$days<=45,]
female.observed=sapply(split(zcweights$weight[zcweights$sex=="F"],
				factor(zcweights$cohort[zcweights$sex=="F"],levels=levels(factor(zcweights.all$cohort)))),mean)
male.observed=sapply(split(zcweights$weight[zcweights$sex=="M"],
				factor(zcweights$cohort[zcweights$sex=="M"],levels=levels(factor(zcweights.all$cohort)))),mean)
female.observed.n=sapply(split(zcweights$weight[zcweights$sex=="F"],
				factor(zcweights$cohort[zcweights$sex=="F"],levels=levels(factor(zcweights.all$cohort)))),length)
male.observed.n=sapply(split(zcweights$weight[zcweights$sex=="M"],
				factor(zcweights$cohort[zcweights$sex=="M"],levels=levels(factor(zcweights.all$cohort)))),length)

fall.avg.day=tapply(zcweights$days,zcweights$cohort,mean)
all.cohorts=sort(unique(zcweights.all$cohort))
#
# compute observed winter averages - using data > 15 Nov
#
zcweights=zcweights.all[zcweights.all$days>45,]
female.winter.observed=sapply(split(zcweights$weight[zcweights$sex=="F"],
				factor(zcweights$cohort[zcweights$sex=="F"],levels=levels(factor(zcweights.all$cohort)))),mean)
male.winter.observed=sapply(split(zcweights$weight[zcweights$sex=="M"],
				factor(zcweights$cohort[zcweights$sex=="M"],levels=levels(factor(zcweights.all$cohort)))),mean)
male.winter.observed[is.nan(male.winter.observed)]=NA
female.winter.observed[is.nan(female.winter.observed)]=NA
female.winter.observed.n=sapply(split(zcweights$weight[zcweights$sex=="F"],
				factor(zcweights$cohort[zcweights$sex=="F"],levels=levels(factor(zcweights.all$cohort)))),length)
male.winter.observed.n=sapply(split(zcweights$weight[zcweights$sex=="M"],
				factor(zcweights$cohort[zcweights$sex=="M"],levels=levels(factor(zcweights.all$cohort)))),length)

winter.avg.day=tapply(zcweights$days,factor(zcweights$cohort,levels=sort(unique(zcweights.all$cohort))),mean)
# For modelling use all of the data from 1 Sept to 1 Mar
zcweights=zcweights.all
zcweights$batch=factor(paste(zcweights$cohort,zcweights$days))
#
#  fit growth model- only adjusting for timing of collection date in fall
#
# First fit a sequence of random effect models with REML to assess best random model with same fixed model
random.f=list(list(~1|cohort,~1|batch),list(~days|cohort,~1|batch),list(~sex:days|cohort,~1|batch),list(~sex*days|cohort,~1|batch),list(~-1+sex+days|cohort,~1|batch))
fixed.f=list(weight~sex*days)
res.adjust=fitmixed(fixed.f,random.f,data=zcweights) 

# Using that random model, fit a sequence of fixed effect models with ML and use AIC to assess best fixed model
random.f=list(res.adjust$best.r)
fixed.f=list(
		weight~sex*days,
		weight~sex+days,
		weight~sex,
		weight~days,
		weight~1)
res.adjust=fitmixed(fixed.f,random.f,data=zcweights) 

# Finally fit best fixed/random model with REML
zc.weight.model=lme(fixed=res.adjust$best.f,random=res.adjust$best.r,data=zcweights,method="REML",control=lmeControl(opt="optim"))
print(summary(zc.weight.model))
zc.weight.adjust.model=zc.weight.model

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
		mod=try(lme(res.adjust$best.f,random=res.adjust$best.r,data=as.data.frame(xsamp),control=lmeControl(opt="optim")))
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


################# 1 Oct Predictions #####################
stderrors=bootstrap.se(zcweights,nboot,days=0)
# Compute predictions and construct dataframes for female and male averages with bootstrap std errors
pp=data.frame(days=0,cohort=rep(sort(unique(zcweights$cohort)),2),sex=rep(c("F","M"),each=length(unique(zcweights$cohort))))
pp$predict=predict(zc.weight.model,newdata=pp,level=1)
female.averages=data.frame(fit=pp$predict[pp$sex=="F"],se=stderrors[as.numeric(row.names(pp[pp$sex=="F",]))])
male.averages=data.frame(fit=pp$predict[pp$sex=="M"],se=stderrors[as.numeric(row.names(pp[pp$sex=="M",]))])
# create dataframe 
ZCWeight.df=data.frame(female.observed.mean.fall=female.observed,female.fall.n=female.observed.n,
		female.adjusted.mean.fall=female.averages$fit,female.adjusted.mean.fall.se=female.averages$se,
		male.observed.mean.fall=male.observed,male.fall.n=male.observed.n,
		male.adjusted.mean.fall=male.averages$fit,male.adjusted.mean.fall.se=male.averages$se,fall.avg.day=fall.avg.day)
ZCWeight.df=cbind(Year=as.numeric(rownames(ZCWeight.df)),ZCWeight.df)


################# 1 Feb Predictions #####################
pp=data.frame(days=123,cohort=rep(sort(unique(zcweights$cohort)),2),sex=rep(c("F","M"),each=length(unique(zcweights$cohort))))
pp$predict=predict(zc.weight.model,newdata=pp,level=1)
stderrors=bootstrap.se(zcweights,nboot,days=123)
winter.female.averages=data.frame(fit=pp$predict[pp$sex=="F"],se=stderrors[as.numeric(row.names(pp[pp$sex=="F",]))])
winter.male.averages=data.frame(fit=pp$predict[pp$sex=="M"],se=stderrors[as.numeric(row.names(pp[pp$sex=="M",]))])

# add to dataframe
ZCWeight.df$female.observed.mean.winter=female.winter.observed
ZCWeight.df$female.winter.n=female.winter.observed.n
ZCWeight.df$female.adjusted.mean.winter=winter.female.averages$fit
ZCWeight.df$female.adjusted.mean.se.winter=winter.female.averages$se
ZCWeight.df$male.observed.mean.winter=male.winter.observed
ZCWeight.df$male.winter.n=male.winter.observed.n
ZCWeight.df$male.adjusted.mean.winter=winter.male.averages$fit
ZCWeight.df$male.adjusted.mean.se.winter=winter.male.averages$se
ZCWeight.df$winter.avg.day=winter.avg.day

