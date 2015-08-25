# uses zcweights.diet which can either be assigned to zcweights.environ or zcweights.environ.abun
# the latter limits to use of diet data where fish abundance data is also available.
# also fixed.f can be set to the environment set or to the set with fish abundance as well
#
#if fdir="" it looks for data files in the CalcurData package directory of your R library.
#if fdir=NULL it looks in databases.txt in CalcurData package directory to get the database location
#if fdir is anything else it uses the value of fdir as the directory for database.  
#The scripts check for the value of fdir and if it exists the script will not change the value; otherwise it sets it to NULL
if(!exists("fdir"))fdir=NULL
require(CIPinnipedAnalysis)
if(!exists("lastyear"))lastyear=2014
if(!exists("nboot"))nboot=100
if(!exists("anomalies"))
	source(file.path(system.file(package="CIPinnipedAnalysis"),"CreateAnomalies.r"))
if(!exists("use.calcofi"))use.calcofi=FALSE
# if ZC_Weight_Environment_Model has not been run, run it now
if(!exists("ZCWeight.df"))
{
	source(file.path(system.file(package="CIPinnipedAnalysis"),"ZC_Weight_Environment_Model.r"))
} else 
{
	if(is.null(ZCWeight.df$female.environ.mean.fall))
		source(file.path(system.file(package="CIPinnipedAnalysis"),"ZC_Weight_Environment_Model.r"))
}

# select data with fish abundance
zcweights.environ.diet=zcweights.diet[!is.na(zcweights.diet$sardine),]
zcweights.environ.diet=droplevels(zcweights.environ.diet)

#
# Next evaluate sequence of fixed-effect models with the chosen random effect model
#
random.f=list(res.environ$best.r)

fixed.f.diet=c(fixed.f,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ PC1"))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ squid "))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ hake "))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ rockfish "))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ sardine "))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ sa"))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ anchovy"))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ lo.cal"))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ hi.cal"))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ ratio"))))

res.environ.diet=fitmixed(fixed.f.diet,random.f,data=zcweights.environ.diet) 
res.environ.diet=compute_AICc(res.environ.diet,length(table(zcweights.environ.diet$Year))*2,5)

#  fit best environment-diet growth model with REML
zc.weight.model.environ.diet=lme(res.environ.diet$best.f,random=res.environ.diet$best.r,data=zcweights.environ.diet,control=lmeControl(opt="optim"),method="REML")
summary(zc.weight.model.environ.diet)


bootstrap.se=function(x,nreps,days=0)
{
	pmat=matrix(0,nrow=nreps,ncol=nrow(unique(data.frame(cohort=x$cohort,sex=x$sex))))
	i=0
	while(i<nreps)
	{
		cat("Bootstrap ",i,"\n")
		xsamp=lapply(split(x,list(x$batch,x$sex)),function(x) if(nrow(x)>0) x[sample(1:nrow(x),replace=TRUE),] else NULL)
		xsamp=do.call("rbind",xsamp)
		mod=try(lme(fixed=res.environ.diet$best.f,random=res.environ.diet$best.r,data=as.data.frame(xsamp),method="REML",control=lmeControl(opt="optim")))
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
stderrors=bootstrap.se(zcweights.environ.diet,nboot,days=0)

# Compute fall 1 Oct predictions and construct dataframes for female and male averages with std errors
pp=zcweights.environ.diet
pp$days=0
pp1=predict(zc.weight.model.environ.diet,newdata=pp)
pp0=predict(zc.weight.model.environ.diet,newdata=pp,level=0)

pp0=tapply(as.vector(pp0),list(zcweights.environ.diet$cohort,zcweights.environ.diet$sex),mean)
pp1=tapply(as.vector(pp1),list(zcweights.environ.diet$cohort,zcweights.environ.diet$sex),mean)
# create list with estimates and std errors for averages
female.averages=data.frame(fit=pp1[,1],se=stderrors[1:length(pp1[,1])])
male.averages=data.frame(fit=pp1[,2],se=stderrors[(length(pp1[,1])+1):(2*length(pp1[,1]))])
# create list with estimates and std errors for expected averages
expected.female.averages=data.frame(fit=pp0[,1],se=stderrors[1:length(pp0[,1])])
expected.male.averages=data.frame(fit=pp0[,2],se=stderrors[(length(pp0[,1])+1):(2*length(pp0[,1]))])

ZCWeight.df1=data.frame(Year=sort(unique(zcweights.environ.diet$Year)),female.environ.diet.mean.fall=female.averages$fit,
female.environ.diet.mean.fall.se=female.averages$se,
male.environ.diet.mean.fall=male.averages$fit,
male.environ.diet.mean.fall.se=male.averages$se,
female.environ.diet.mean.fall.fixed=expected.female.averages$fit,
male.environ.diet.mean.fall.fixed=expected.male.averages$fit)

ZCWeight.df=merge(ZCWeight.df,ZCWeight.df1,all.x=TRUE)
