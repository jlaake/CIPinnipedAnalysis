#if fdir="" it looks for data files in the CalcurData package directory of your R library.
#if fdir=NULL it looks in databases.txt in CalcurData package directory to get the database location
#if fdir is anything else it uses the value of fdir as the directory for database.  
#The scripts check for the value of fdir and if it exists the script will not change the value; otherwise it sets it to NULL
if(!exists("fdir"))fdir=NULL
if(!exists("nboot"))nboot=100
####################################
# Set this value; be aware that all of the environmental data has to be entered through Feb of lastyear+1 
# for the growth script to work properly
if(!exists("lastyear"))lastyear=2013
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
cu.male.averages=male.averages
CUWeight.df=data.frame(female.observed.mean=female.observed,
		female.adjusted.mean=cu.female.averages$fit,female.adjusted.mean.se=cu.female.averages$se,
		male.observed.mean=male.observed,
		male.adjusted.mean=cu.male.averages$fit,male.adjusted.mean.se=cu.male.averages$se)

