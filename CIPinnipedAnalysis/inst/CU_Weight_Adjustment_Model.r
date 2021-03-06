#if fdir="" it looks for data files in the CalcurData package directory of your R library.
#if fdir=NULL it looks in databases.txt in CalcurData package directory to get the database location
#if fdir is anything else it uses the value of fdir as the directory for database.  
#The scripts check for the value of fdir and if it exists the script will not change the value; otherwise it sets it to NULL
if(!exists("fdir"))fdir=NULL
if(!exists("nboot"))nboot=100
####################################

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
if(!exists("lastyear"))lastyear=max(cuweights.all$cohort)
cuweights.all=cuweights.all[cuweights.all$cohort<=lastyear,]
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
#
#  fit growth model- only adjusting for timing of collection date in fall
#
# First fit a sequence of random effect models with REML to assess best random model with same fixed model
random.f=list(list(~1|cohort),list(~days|cohort),list(~sex:days|cohort),list(~sex*days|cohort),list(~-1+sex+days|cohort))
fixed.f=list(weight~sex*days)
res.adjust=fitmixed(fixed.f,random.f,data=cuweights) 

# Using that random model, fit a sequence of fixed effect models with ML and use AIC to assess best fixed model
random.f=list(res.adjust$best.r)
fixed.f=list(
		weight~sex*days,
		weight~sex+days,
		weight~sex,
		weight~days,
		weight~1)
res.adjust=fitmixed(fixed.f,random.f,data=cuweights) 

# Finally fit best fixed/random model with REML
cu.weight.model=lme(fixed=res.adjust$best.f,random=res.adjust$best.r,data=res.adjust$data,method="REML",control=lmeControl(opt="optim"))
print(summary(cu.weight.model))


# define bootstrap function to compute std error for predicted sex-cohort means
bootstrap.se=function(x,nreps)
{
	pmat=matrix(0,nrow=nreps,ncol=nrow(unique(data.frame(cohort=x$cohort,sex=x$sex))))
	i=0
	while(i<nreps)
	{
		xsamp=lapply(split(x,list(x$sex,x$cohort)),function(x) if(nrow(x)>0) x[sample(1:nrow(x),replace=TRUE),] else NULL)
		xsamp=do.call("rbind",xsamp)
		mod=try(lme(formula(cu.weight.model),random=as.formula(cu.weight.model$call$random),data=as.data.frame(xsamp),control=lmeControl(opt="optim")))
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
stderrors=bootstrap.se(res.adjust$data,nboot)
# Compute predictions and construct dataframes for female and male averages with std errors
pp=data.frame(days=0,cohort=rep(sort(unique(cuweights$cohort)),2),sex=rep(c("F","M"),each=length(unique(cuweights$cohort))))
pp$predict=predict(cu.weight.model,newdata=pp,level=1)
female.averages=data.frame(fit=pp$predict[pp$sex=="F"],se=stderrors[as.numeric(row.names(pp[pp$sex=="F",]))])
male.averages=data.frame(fit=pp$predict[pp$sex=="M"],se=stderrors[as.numeric(row.names(pp[pp$sex=="M",]))])
# create dataframe with normal conf interval
cu.female.averages=female.averages
cu.male.averages=male.averages
CUWeight.df=data.frame(female.observed.mean.fall=female.observed,
		female.adjusted.mean.fall=cu.female.averages$fit,female.adjusted.mean.fall.se=cu.female.averages$se,
		male.observed.mean.fall=male.observed,
		male.adjusted.mean.fall=cu.male.averages$fit,male.adjusted.meanfall..se=cu.male.averages$se)
CUWeight.df=cbind(Year=as.numeric(rownames(CUWeight.df)),CUWeight.df)
